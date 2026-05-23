get_vardata = function(abundance_se, which = "cols") {
  if (which %in% "cols") {
    colData(abundance_se) |> as.data.frame()
  } else {
    rowData(abundance_se) |> as.data.frame()
  }
}

sets_run_pca_experimental_design = function(lc_sets) {
  # lc_sets = tar_read(lc_g01); tar_source("R")
  # pca_sets = purrr::map(lc_sets[1:10], run_pca_experimental_design, .progress = TRUE)
  pca_sets = furrr::future_map(
    lc_sets,
    run_pca_experimental_design
  )
  pca_sets
}

run_pca_experimental_design = function(in_lc) {
  # in_lc = lc_sets[[1]]
  # in_lc = lc_sets$AN002911
  tmp_norm = assays(in_lc$abundance)$normalized
  if (
    any(is.na(in_lc$icikt_cor$cor)) ||
      any(is.infinite(tmp_norm)) ||
      min(tmp_norm, na.rm = TRUE) < 0
  ) {
    return(in_lc)
  }

  all_abundance = in_lc$abundance
  possible_outliers = determine_outliers(
    in_lc$icikt_cor,
    get_vardata(all_abundance)
  )
  just_outlier_info = possible_outliers |>
    dplyr::transmute(
      sample_id = sample_id,
      median_ici_correlation = med_cor,
      outlier_fraction = frac,
      score = score,
      outlier = outlier
    )
  just_outlier_info = dplyr::left_join(
    data.frame(sample_id = all_abundance$sample_id),
    just_outlier_info,
    by = "sample_id"
  )

  colData(all_abundance) = cbind(
    colData(all_abundance),
    just_outlier_info |> dplyr::select(-sample_id)
  )
  all_res = list(
    pca = do_pca(all_abundance),
    design = all_abundance |>
      get_vardata() |>
      determine_experimental_design(),
    order = determine_neigbors_ordering(
      in_lc$icikt_cor,
      get_vardata(all_abundance)
    )
  )

  no_qc = all_abundance[, all_abundance$experiment_or_qc %in% "experiment"]
  not_outlier = possible_outliers |>
    dplyr::filter(!outlier) |>
    dplyr::pull(sample_id)
  no_outlier = all_abundance[,
    all_abundance$sample_id %in% not_outlier
  ]
  not_outlier_qc = base::intersect(
    get_vardata(all_abundance) |>
      dplyr::filter(experiment_or_qc %in% "experiment") |>
      dplyr::pull(sample_id),
    not_outlier
  )

  not_outlier_res = list(
    pca = do_pca(no_outlier),
    design = no_outlier |> get_vardata() |> determine_experimental_design(),
    order = determine_neigbors_ordering(
      in_lc$icikt_cor,
      get_vardata(no_outlier)
    )
  )
  noqc_res = list(
    pca = do_pca(no_qc),
    design = no_qc |> get_vardata() |> determine_experimental_design(),
    order = determine_neigbors_ordering(
      in_lc$icikt_cor,
      get_vardata(no_qc)
    )
  )

  noqc_nooutlier = all_abundance[,
    all_abundance$sample_id %in% not_outlier_qc
  ]

  noqc_nooutlier_res = list(
    pca = do_pca(noqc_nooutlier),
    design = noqc_nooutlier |> get_vardata() |> determine_experimental_design(),
    order = determine_neigbors_ordering(
      in_lc$icikt_cor,
      get_vardata(noqc_nooutlier)
    )
  )

  in_lc$abundance = all_abundance
  in_lc$pca_groups = list(
    all = all_res,
    not_outlier = not_outlier_res,
    not_qc = noqc_res,
    not_outlier_qc = noqc_nooutlier_res
  )
  in_lc
}


determine_neigbors_ordering = function(ici_values, group_data) {
  # ici_values = in_lc$icikt_cor; group_data = get_vardata(in_lc$abundance)
  ici_square = ICIKendallTau::long_df_2_cor_matrix(ici_values)
  ici_square = ici_square[group_data$sample_id, group_data$sample_id]
  med_ici = visualizationQualityControl::median_correlations(
    ici_square,
    sample_classes = group_data$factors,
    between_classes = TRUE
  )

  ici_within = med_ici |>
    dplyr::filter(group %in% "within") |>
    dplyr::mutate(within = med_cor)

  ici_between = med_ici |>
    dplyr::filter(group %in% "between") |>
    dplyr::summarise(between = max(med_cor), .by = sample_id)

  ici_compare = dplyr::left_join(
    ici_within[, c("sample_id", "within", "sample_class")],
    ici_between,
    by = "sample_id"
  ) |>
    dplyr::mutate(w_b_better = (within >= between), factors = sample_class)

  ici_dendrogram = visualizationQualityControl::similarity_reorderbyclass(
    ici_square,
    sample_classes = NULL,
    transform = "sub_1"
  )

  sample_order = ici_dendrogram$names

  actual_rles = rle(sort(group_data$factors)) |> rle_2_df()

  sample_factors = group_data$factors
  names(sample_factors) = group_data$sample_id

  sample_factors_ordered = sample_factors[sample_order]
  names(sample_factors_ordered) = NULL
  distance_rles = rle(sample_factors_ordered) |> rle_2_df()

  join_rles = dplyr::left_join(
    actual_rles,
    distance_rles,
    by = "value",
    suffix = c(".ideal", ".ici")
  )
  join_rles = join_rles |> dplyr::mutate(diff = length.ideal - length.ici)

  list(
    ici = ici_square,
    dendrogram = ici_dendrogram,
    diff_rle = join_rles,
    diff_median_ici = ici_compare
  )
  #return(NULL)
}

rle_2_df = function(rle_data) {
  # rle_data = actual_rles
  tibble::tibble(value = rle_data$values, length = rle_data$lengths)
}

determine_outliers = function(icikt_cor, ssf) {
  # icikt_cor = in_lc$icikt_cor; ssf = get_vardata(in_lc$abundance)

  ici_matrix = ICIKendallTau::long_df_2_cor_matrix(icikt_cor)
  ici_matrix = ici_matrix[ssf$sample_id, ssf$sample_id]

  med_cor = visualizationQualityControl::median_correlations(
    ici_matrix,
    ssf$factors
  )
  ici_outliers = visualizationQualityControl::determine_outliers(
    median_correlations = med_cor
  )
  ici_outliers
}

do_pca = function(abundance_se) {
  # norm_data = in_lc$normalized |>
  #   get_just_values()
  # ssf_wide = in_lc$ssf
  # abundance_se = in_lc$abundance
  norm_matrix = assays(abundance_se)$normalized
  ssf = get_vardata(abundance_se, "cols")
  norm_matrix[
    is.na(norm_matrix)
  ] = 0
  log_matrix = log1p(norm_matrix)
  pca_res = prcomp(t(log_matrix), center = TRUE, scale. = FALSE)

  extra_vars = c(
    "experiment_or_qc",
    "median_ici_correlation",
    "outlier_fraction",
    "score",
    "outlier"
  )
  pca_var = visualizationQualityControl::visqc_score_contributions(pca_res$x)
  pca_anova = visualizationQualityControl::visqc_test_pca_scores(
    pca_res$x,
    ssf |>
      dplyr::select(!(tidyselect::all_of(extra_vars)))
  )
  pca_x = pca_res$x
  pca_df = cbind(as.data.frame(pca_x), ssf)

  n_pc = ncol(pca_x)
  use_pc = colnames(pca_x)
  label_vals = purrr::map(seq(1, n_pc), \(x) {
    pca_var$labels[x]
  })
  names(label_vals) = use_pc

  labelled::var_label(pca_df) = label_vals

  list(pca = pca_df, variance = pca_var, anova = pca_anova)
}

sets_create_qcqa_plots = function(pca_sets) {
  out_plots = purrr::map(pca_sets, create_qcqa_plots)
  out_plots
}

create_qcqa_plots = function(in_pca_set) {
  # in_pca_set = pca_sets[[1]]

  all_factors = unique(in_pca_set$abundance$factors)
  factor_colors = scales::hue_pal()(length(all_factors))
  names(factor_colors) = all_factors

  n_value = 40
  cor_colors = circlize::colorRamp2(
    seq(-1, 1, length.out = n_value),
    scico::scico(n_value, palette = "berlin")
  )
  raw_abundance = log2(assays(in_pca_set$abundance)$raw)
  norm_abundance = log2(assays(in_pca_set$abundance)$normalized)

  sample_metadata = get_vardata(in_pca_set$abundance)

  raw_boxplot = create_boxplot(raw_abundance, sample_metadata, factor_colors)
  norm_boxplot = create_boxplot(norm_abundance, sample_metadata, factor_colors)

  pca_heatmap_median_ici_each = purrr::map(
    in_pca_set$pca_groups,
    create_pca_median_ici_heatmap,
    factor_colors,
    cor_colors
  )

  ordering_scores_each = purrr::map(
    in_pca_set$pca_groups,
    \(in_group) {
      dplyr::bind_rows(
        score_rle(in_group$order$diff_rle),
        score_median_ici(in_group$order$diff_median_ici)
      )
    }
  )

  list(
    raw = raw_boxplot,
    normalized = norm_boxplot,
    iterative_plots = pca_heatmap_median_ici_each,
    scores = ordering_scores_each
  )
}

score_rle = function(rle_df) {
  # rle_df = pca_group$order$diff_rle
  rle_summary_all = rle_df |>
    dplyr::summarise(
      grp_diff = sum(diff),
      grp_ideal = sum(length.ideal),
      fraction_diff = grp_diff / grp_ideal,
      factors = "all"
    )
  rle_summary_group = rle_df |>
    dplyr::summarise(
      grp_diff = sum(diff),
      grp_ideal = sum(length.ideal),
      fraction_diff = grp_diff / grp_ideal,
      .by = value
    ) |>
    dplyr::mutate(factors = value) |>
    dplyr::select(-value)
  dplyr::bind_rows(rle_summary_all, rle_summary_group) |>
    dplyr::mutate(method = "rle")
}

score_median_ici = function(ici_df) {
  # ici_df = pca_group$order$diff_median_ici
  ici_summary_all = ici_df |>
    dplyr::summarise(
      grp_diff = sum(w_b_better) - dplyr::n(),
      grp_ideal = dplyr::n(),
      fraction_diff = grp_diff / grp_ideal,
      factors = "all"
    )
  ici_summary_group = ici_df |>
    dplyr::summarise(
      grp_diff = sum(w_b_better) - dplyr::n(),
      grp_ideal = dplyr::n(),
      fraction_diff = grp_diff / grp_ideal,
      .by = factors
    )
  dplyr::bind_rows(ici_summary_all, ici_summary_group) |>
    dplyr::mutate(method = "median-ici")
}

create_pca_median_ici_heatmap = function(pca_group, factor_colors, cor_colors) {
  # pca_group = in_pca_set$pca_groups$all
  pca_plot = pca_group$pca$pca |>
    ggplot(aes(x = PC1, y = PC2, color = factors, shape = outlier)) +
    geom_point(size = 3) +
    scale_color_manual(values = factor_colors)

  heatmap_data = pca_group$order

  sample_info = pca_group$pca$pca

  sample_info = dplyr::left_join(
    data.frame(sample_id = rownames(heatmap_data$ici)),
    sample_info,
    by = "sample_id"
  )

  heatmap_colors = list(
    factors = factor_colors,
    outlier = c("FALSE" = "#000000", "TRUE" = "#FFFFFF")
  )

  left_var = ComplexHeatmap::HeatmapAnnotation(
    df = sample_info[, c("factors", "outlier")],
    col = heatmap_colors,
    which = "row",
    show_annotation_name = FALSE
  )
  bottom_var = ComplexHeatmap::HeatmapAnnotation(
    df = sample_info[, c("factors", "outlier")],
    col = heatmap_colors,
    which = "column",
    show_legend = FALSE
  )

  out_heatmap = ComplexHeatmap::Heatmap(
    matrix = heatmap_data$ici,
    col = cor_colors,
    name = "ICI-Kt",
    right_annotation = left_var,
    bottom_annotation = bottom_var,
    show_column_names = FALSE,
    cluster_columns = pca_group$order$dendrogram$dendrogram,
    cluster_rows = pca_group$order$dendrogram$dendrogram
  )

  out_sina = pca_group$pca$pca |>
    ggplot(aes(
      x = score,
      y = factors,
      color = factors,
      shape = outlier,
      group = factors
    )) +
    geom_sina(size = 3) +
    scale_color_manual(values = factor_colors) +
    theme(legend.position = "top")
  return(list(pca = pca_plot, heatmap = out_heatmap, median_ici = out_sina))
}

create_boxplot = function(
  abundance_matrix,
  sample_metadata,
  point_colors_fill
) {
  # abundance_matrix = raw_abundance; point_colors_fill = factor_colors
  sample_metadata = sample_metadata |>
    dplyr::arrange(factors, sample_id)
  abundance_matrix = abundance_matrix[, sample_metadata$sample_id]

  sample_metadata$factors = factor(sample_metadata$factors)

  abundance_df = abundance_matrix |>
    as.data.frame() |>
    tidyr::pivot_longer(cols = everything(), names_to = "sample_id")
  abundance_df = dplyr::left_join(
    abundance_df,
    sample_metadata[, c("sample_id", "factors")],
    by = "sample_id"
  )

  sample_metadata$sample_id = factor(
    sample_metadata$sample_id,
    levels = sample_metadata$sample_id
  )
  abundance_df$sample_id = factor(
    abundance_df$sample_id,
    levels = levels(sample_metadata$sample_id)
  )
  out_boxplot = abundance_df |>
    ggplot(aes(x = value, y = sample_id, group = sample_id)) +
    geom_sina(aes(color = factors)) +
    geom_boxplot(alpha = 0) +
    scale_color_manual(values = point_colors_fill) +
    labs(x = "log2(Abundance)") +
    theme(legend.position = "top")
  out_boxplot
}
