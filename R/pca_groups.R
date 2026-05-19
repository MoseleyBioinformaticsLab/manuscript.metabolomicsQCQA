get_vardata = function(abundance_se, which = "cols") {
  if (which %in% "cols") {
    colData(abundance_se) |> as.data.frame()
  } else {
    rowData(abundance_se) |> as.data.frame()
  }
}

sets_run_pca_experimental_design = function(lc_sets) {
  # lc_sets = tar_read(lc_g07)
  pca_sets = furrr::future_map(
    lc_sets,
    run_pca_experimental_design
  )
  pca_sets
}

run_pca_experimental_design = function(in_lc) {
  # in_lc = lc_sets[[1]]
  all_abundance = in_lc$abundance
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

  possible_outliers = determine_outliers(
    in_lc$icikt_cor,
    get_vardata(all_abundance)
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

  rownames(possible_outliers) = possible_outliers$sample_id
  possible_outliers = possible_outliers[all_abundance$sample_id, ]

  all_abundance$med_cor = possible_outliers$med_cor
  all_abundance$score = possible_outliers$score
  all_abundance$outlier = possible_outliers$outlier
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

  ici_dendrogram = visualizationQualityControl::similarity_reorder(
    ici_square,
    transform = "sub_1"
  )
  ici_dendrogram$samples = rownames(ici_square)[ici_dendrogram$indices]
  sample_order = rownames(ici_square)[ici_dendrogram$indices]

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
  norm_matrix[is.na(norm_matrix)] = 0
  log_matrix = log1p(norm_matrix)
  pca_res = prcomp(t(log_matrix), center = TRUE, scale. = FALSE)

  pca_var = visualizationQualityControl::visqc_score_contributions(pca_res$x)
  pca_anova = visualizationQualityControl::visqc_test_pca_scores(
    pca_res$x,
    ssf
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

create_qcqa_plots = function(in_pca_set) {
  # in_pca_set = pca_sets[[3]]

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
}

create_pca_heatmap = function(pca_group, factor_colors) {
  # pca_group = in_pca_set$pca_groups$all
  pca_plot = pca_group$pca$pca |>
    ggplot(aes(x = PC1, y = PC2, color = factors)) +
    geom_point(size = 3) +
    scale_color_manual(values = factor_colors)
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
    ggplot(aes(x = value, y = sample_id, color = factors, group = sample_id)) +
    geom_sina() +
    geom_boxplot(alpha = 0.2) +
    scale_color_manual(values = point_colors_fill) +
    labs(x = "log2(Abundance)")
  out_boxplot
}
