get_vardata = function(se_data, which = "cols") {
  if (which %in% "cols") {
    colData(se_data) |> as.data.frame()
  } else {
    rowData(se_data) |> as.data.frame()
  }
}

add_coldata = function(se_data, additional_data) {
  if (nrow(additional_data) != nrow(colData(se_data))) {
    stop("Number of rows doesn't match!")
  }
  if ("sample_id" %in% colnames(additional_data)) {
    tmp_df = tibble::tibble(sample_id = se_data$sample_id)
    additional_data = dplyr::left_join(
      tmp_df,
      additional_data,
      by = "sample_id"
    )
    additional_data = additional_data |>
      dplyr::select(-sample_id)
  }
  colData(se_data) = cbind(colData(se_data), additional_data)
  se_data
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

get_pca_run_times = function(pca_sets) {
  out_time = purrr::imap(pca_sets, \(x, id) {
    if (!is.null(x[["pca_groups"]])) {
      tibble::tibble(id = id, time = x$pca_groups$wall_time)
    } else {
      NULL
    }
  }) |>
    purrr::list_rbind()
  out_time
}

check_na_normalized = function(cor_sets) {
  has_problem = purrr::map_lgl(cor_sets, \(x) {
    tmp_norm = x$normalized |> get_just_values()
    if (
      any(is.na(x$icikt_cor$cor)) ||
        any(is.infinite(tmp_norm)) ||
        min(tmp_norm, na.rm = TRUE) < 0
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  has_problem
}

get_specific_metadata = function(se_data, item) {
  metadata(se_data)[[item]]
}

run_pca_experimental_design = function(se_data) {
  # se_data = lc_sets[[1]]
  # se_data = lc_sets$AN002911
  # se_data = lc_sets$AN001749
  # se_data = lc_sets$AN004304
  # check lc_g04; AN000777; AN001296
  # lc_sets = tar_read(lc_g10); se_data = lc_sets$AN003604; tar_source("R")
  # se_data = cor_sets[[1]]
  # tar_load(se_pca_g02); se_data = se_pca_g02[["AN004656"]]; tar_source("R")
  start_time = Sys.time()

  if (metadata(se_data) |> check_isstop()) {
    return(se_data)
  }

  se_kept = se_data[, is.na(se_data$removed)]
  icikt_data = get_specific_metadata(se_kept, "icikt")$cor
  possible_outliers = determine_outliers(
    icikt_data,
    se_kept |> get_vardata()
  )
  just_outlier_info = possible_outliers |>
    dplyr::transmute(
      sample_id = sample_id,
      median_ici_correlation = med_cor,
      outlier_fraction = frac,
      score = score,
      outlier = outlier
    )

  # this is here so that things removed during ICI-Kt, and therefore
  # outlier detection, get set to FALSE, as they are NOT actually outliers,
  # but still need to be handled.
  just_outlier_info = dplyr::left_join(
    data.frame(sample_id = se_data$sample_id),
    just_outlier_info,
    by = "sample_id"
  )

  just_outlier_info = just_outlier_info |>
    dplyr::mutate(outlier = dplyr::replace_values(outlier, NA ~ FALSE))

  factor_types = get_specific_metadata(se_data, "factor_types")
  factor_types = dplyr::bind_rows(
    factor_types,
    tibble::tribble(
      ~factor                  , ~type      ,
      "median_ici_correlation" , "computed" ,
      "outlier_fraction"       , "computed" ,
      "score"                  , "computed" ,
      "outlier"                , "computed"
    )
  )

  se_data = add_coldata(se_data, just_outlier_info)
  metadata(se_data)$factor_types = factor_types

  se_data = create_additional_factors_ssf(se_data)
  factor_types = get_specific_metadata(se_data, "factor_types")

  all_res = NULL
  not_outlier_res = NULL
  noqc_res = NULL
  noqc_nooutlier_res = NULL

  # first remove anything that needs to be
  se_kept = se_data[, is.na(se_data$removed), drop = FALSE]

  all_res = list(
    pca = do_pca(se_kept),
    design = se_kept |>
      get_vardata() |>
      determine_experimental_design(factor_types),
    order = determine_neigbors_ordering(
      icikt_data,
      get_vardata(se_kept)
    ),
    rsd = calculate_variance_stuff(se_kept)
  )
  all_res$scores = dplyr::bind_rows(
    score_rle(all_res$order$diff_rle),
    score_median_ici(all_res$order$diff_median_ici)
  )

  has_qc = length(which(se_kept$experiment_or_qc %in% "qc")) > 0
  has_outliers = which(se_kept$outlier)

  if (has_qc) {
    no_qc = se_kept[, se_kept$experiment_or_qc %in% "experiment", drop = FALSE]
    if (ncol(no_qc) >= 2) {
      noqc_res = list(
        pca = do_pca(no_qc),
        design = no_qc |>
          get_vardata() |>
          determine_experimental_design(factor_types),
        order = determine_neigbors_ordering(
          icikt_data,
          get_vardata(no_qc)
        ),
        rsd = calculate_variance_stuff(no_qc)
      )
      noqc_res$scores = dplyr::bind_rows(
        score_rle(noqc_res$order$diff_rle),
        score_median_ici(noqc_res$order$diff_median_ici)
      )
    }
  }

  if (length(has_outliers) > 0) {
    no_outlier = se_kept[, !se_kept$outlier, drop = FALSE]
    if (ncol(no_outlier) >= 2) {
      not_outlier_res = list(
        pca = do_pca(no_outlier),
        design = no_outlier |>
          get_vardata() |>
          determine_experimental_design(factor_types),
        order = determine_neigbors_ordering(
          icikt_data,
          get_vardata(no_outlier)
        ),
        rsd = calculate_variance_stuff(no_outlier)
      )
      not_outlier_res$scores = dplyr::bind_rows(
        score_rle(not_outlier_res$order$diff_rle),
        score_median_ici(not_outlier_res$order$diff_median_ici)
      )
    }
  }

  if (has_qc && (length(has_outliers) > 0)) {
    not_outlier_qc = base::intersect(no_outlier$sample_id, no_qc$sample_id)
    noqc_nooutlier = se_kept[,
      se_kept$sample_id %in% not_outlier_qc,
      drop = FALSE
    ]
    if (ncol(noqc_nooutlier) >= 2) {
      noqc_nooutlier_res = list(
        pca = do_pca(noqc_nooutlier),
        design = noqc_nooutlier |>
          get_vardata() |>
          determine_experimental_design(factor_types),
        order = determine_neigbors_ordering(
          icikt_data,
          get_vardata(noqc_nooutlier)
        ),
        rsd = calculate_variance_stuff(noqc_nooutlier)
      )
      noqc_nooutlier_res$scores = dplyr::bind_rows(
        score_rle(noqc_nooutlier_res$order$diff_rle),
        score_median_ici(noqc_nooutlier_res$order$diff_median_ici)
      )
    }
  }

  end_time = Sys.time()
  wall_time = difftime(end_time, start_time, units = "secs") |> as.numeric()
  metadata(se_data)$pca_groups = list(
    all = all_res,
    not_outlier = not_outlier_res,
    not_qc = noqc_res,
    not_outlier_qc = noqc_nooutlier_res,
    wall_time = wall_time
  )
  metadata(se_data)$qc = add_step(
    metadata(se_data)$qc,
    "succeeded",
    "pca"
  )
  metadata(se_data)$factor_types = factor_types

  se_data
}


determine_neigbors_ordering = function(ici_values, group_data) {
  # ici_values = get_specific_metadata(se_kept, "icikt")$cor; group_data = get_vardata(se_kept)
  ici_square = ICIKendallTau::long_df_2_cor_matrix(ici_values)
  ici_square = ici_square[group_data$sample_id, group_data$sample_id]
  med_ici = visualizationQualityControl::median_correlations(
    ici_square,
    sample_classes = group_data$ssf,
    between_classes = TRUE
  )

  ici_within = med_ici |>
    dplyr::filter(group %in% "within") |>
    dplyr::mutate(within_icikt = med_cor)

  ici_between = med_ici |>
    dplyr::filter(group %in% "between") |>
    dplyr::summarise(between_icikt = max(med_cor), .by = sample_id)

  if (nrow(ici_between) == 0) {
    ici_compare = ici_within |>
      dplyr::mutate(between_icikt = 0) |>
      dplyr::mutate(
        w_b_better = (within_icikt >= between_icikt),
        ssf = sample_class
      )
  } else {
    ici_compare = dplyr::left_join(
      ici_within[, c("sample_id", "within_icikt", "sample_class")],
      ici_between,
      by = "sample_id"
    ) |>
      dplyr::mutate(
        w_b_better = (within_icikt >= between_icikt),
        ssf = sample_class
      )
  }

  ici_dendrogram = visualizationQualityControl::similarity_reorderbyclass(
    ici_square,
    sample_classes = NULL,
    transform = "sub_1"
  )

  sample_order = ici_dendrogram$names

  actual_rles = rle(sort(group_data$ssf)) |> rle_2_df()

  sample_factors = group_data$ssf
  names(sample_factors) = group_data$sample_id

  sample_factors_ordered = sample_factors[sample_order]
  names(sample_factors_ordered) = NULL
  distance_rles = rle(sample_factors_ordered) |>
    rle_2_df() |>
    dplyr::summarise(length = max(length), .by = value)

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
  # icikt_cor = se_data$icikt_cor; ssf = get_vardata(se_data$abundance)

  ici_matrix = ICIKendallTau::long_df_2_cor_matrix(icikt_cor)
  ici_matrix = ici_matrix[ssf$sample_id, ssf$sample_id]

  med_cor = visualizationQualityControl::median_correlations(
    ici_matrix,
    ssf$ssf
  )
  ici_outliers = visualizationQualityControl::determine_outliers(
    median_correlations = med_cor
  )
  ici_outliers
}

check_matches_ssf = function(ssf_df) {
  # ssf_df = ssf[, c("cognitive_status", "ssf")]
  ssf_rle = rle_2_df(rle(ssf_df$ssf))

  other_cols = names(ssf_df)
  other_cols = other_cols[!(other_cols %in% "ssf")]
  not_match_ssf = purrr::map_lgl(ssf_df[, other_cols, drop = FALSE], \(in_col) {
    other_rle = rle_2_df(rle(in_col))
    !identical(other_rle, ssf_rle)
  })
  names(not_match_ssf[not_match_ssf])
}

do_pca = function(se_data) {
  # norm_data = se_data$normalized |>
  #   get_just_values()
  # ssf_wide = se_data$ssf
  # se_data = se_data$abundance
  if (metadata(se_data) |> check_isstop()) {
    return(se_data)
  }
  norm_matrix = assays(se_data)$normalized
  ssf = get_vardata(se_data, "cols")
  norm_matrix[
    is.na(norm_matrix)
  ] = 0
  log_matrix = log1p(norm_matrix)
  pca_res = prcomp(t(log_matrix), center = TRUE, scale. = FALSE)

  factor_types = metadata(se_data)$factor_types

  pca_var = visualizationQualityControl::visqc_score_contributions(pca_res$x)
  pca_90 = which(pca_var$cumulative > 0.9) |> min()
  if ((length(pca_90) == 0) || is.infinite(pca_90) || is.na(pca_90)) {
    use_test = seq_len(ncol(pca_res$x))
  } else {
    use_test = seq_len(pca_90)
  }

  # we test the SSF, the various underlying factors that make up the SSF, as
  # well as the additional-sample-data, and combinations of additional-sample-data and
  # SSF
  test_vars = c(
    factor_types |>
      dplyr::filter(
        type %in%
          c("subject_sample_factor", "additional_sample_data", "additional_ssf")
      ) |>
      dplyr::pull(factor),
    "ssf"
  )

  different_to_ssf = check_matches_ssf(
    ssf |> dplyr::select(tidyselect::all_of(test_vars))
  )
  test_vars2 = c(different_to_ssf, "ssf")
  pca_anova = visualizationQualityControl::visqc_test_pca_scores(
    pca_res$x[, use_test, drop = FALSE],
    ssf |>
      dplyr::select(tidyselect::all_of(test_vars2))
  )
  pca_anova = pca_anova |>
    dplyr::mutate(
      padjust = p.adjust(p.value, method = "bonferroni"),
      .by = variable
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


create_additional_factors_ssf = function(se_data) {
  factor_types = get_specific_metadata(se_data, "factor_types")
  ssf_df = get_vardata(se_data, "cols")

  additional_factors = factor_types |>
    dplyr::filter(type %in% "additional_sample_data") |>
    dplyr::pull(factor)

  if (length(additional_factors) > 0) {
    ssf_factor = ssf_df[["ssf"]]
    additional_df = purrr::map(additional_factors, \(in_factor) {
      use_factor = ssf_df[[in_factor]]
      if (length(unique(use_factor)) == nrow(ssf_df)) {
        return(NULL)
      }
      new_factor = paste0(ssf_factor, "-", use_factor)
      out_tbl = tibble::tibble(new_f = new_factor)
      names(out_tbl) = paste0("ssf_", in_factor)
      out_tbl
    }) |>
      purrr::list_cbind()

    if (ncol(additional_df) > 0) {
      factor_types_add = tibble::tibble(
        factor = colnames(additional_df),
        type = "additional_ssf"
      )
      factor_types = dplyr::bind_rows(factor_types, factor_types_add)
      metadata(se_data)$factor_types = factor_types
      se_data = add_coldata(se_data, additional_df)
    }
  }
  return(se_data)
}

sets_create_qcqa_plots = function(pca_sets) {
  out_plots = purrr::map(pca_sets, create_qcqa_plots)
  out_plots
}


create_qcqa_plots = function(se_data, report_list) {
  # se_data = pca_sets[[1]]
  # report_list = check_qc_outputs(se_data)

  use_samples = is.na(se_data$removed)
  se_data = se_data[, use_samples, drop = FALSE]

  sample_metadata = get_vardata(se_data)
  all_ssf = unique(sample_metadata$ssf)
  ssf_colors = scales::hue_pal()(length(all_ssf))
  names(ssf_colors) = all_ssf

  n_ssf = length(all_ssf)
  if (n_ssf <= 3) {
    guide_nrow = 1
    guide_ncol = 3
  } else {
    guide_ncol = 3
    guide_nrow = ceiling(n_ssf / guide_ncol)
  }

  n_value = 40
  cor_colors = circlize::colorRamp2(
    seq(-1, 1, length.out = n_value),
    scico::scico(n_value, palette = "berlin")
  )
  raw_abundance = log2(assays(se_data)$raw)
  norm_abundance = log2(assays(se_data)$normalized)

  raw_boxplot = create_boxplot(
    raw_abundance,
    sample_metadata,
    ssf_colors,
    guide_nrow,
    guide_ncol
  )
  norm_boxplot = create_boxplot(
    norm_abundance,
    sample_metadata,
    ssf_colors,
    guide_nrow,
    guide_ncol
  )

  dynamic_ranges = as.matrix(log10(sample_metadata$dynamic_range)) |> t()
  rownames(dynamic_ranges) = "dynamic_range"
  colnames(dynamic_ranges) = sample_metadata$sample_id
  dyn_range_sina = create_sina_boxplot(
    dynamic_ranges,
    sample_metadata,
    ssf_colors,
    guide_nrow,
    guide_ncol
  ) +
    labs(x = "Log10(Dynamic Range)")

  if ("pca" %in% report_list$steps) {
    pca_data = metadata(se_data)$pca_groups
    pca_heatmap_median_ici_each = purrr::map(
      pca_data[c(
        "all",
        "not_outlier",
        "not_qc",
        "not_outlier_qc"
      )],
      create_pca_median_ici_heatmap,
      ssf_colors,
      cor_colors,
      guide_nrow,
      guide_ncol
    )
  } else {
    pca_heatmap_median_ici_each = NULL
  }

  list(
    raw = raw_boxplot,
    normalized = norm_boxplot,
    dynamic_range = dyn_range_sina,
    iterative_plots = pca_heatmap_median_ici_each
  )
}

score_rle = function(rle_df) {
  # rle_df = pca_group$order$diff_rle
  rle_summary_all = rle_df |>
    dplyr::summarise(
      n_differences = sum(diff),
      group_size = sum(length.ideal),
      differences_fraction = n_differences / group_size,
      ssf = "all"
    )
  rle_summary_group = rle_df |>
    dplyr::summarise(
      n_differences = sum(diff),
      group_size = sum(length.ideal),
      differences_fraction = n_differences / group_size,
      .by = value
    ) |>
    dplyr::mutate(ssf = value) |>
    dplyr::select(-value)
  dplyr::bind_rows(rle_summary_all, rle_summary_group) |>
    dplyr::mutate(method = "rle")
}

score_median_ici = function(ici_df) {
  # ici_df = pca_group$order$diff_median_ici
  ici_summary_all = ici_df |>
    dplyr::summarise(
      n_differences = sum(!w_b_better),
      group_size = dplyr::n(),
      differences_fraction = n_differences / group_size,
      ssf = "all"
    )
  ici_summary_group = ici_df |>
    dplyr::summarise(
      n_differences = sum(!w_b_better),
      group_size = dplyr::n(),
      differences_fraction = n_differences / group_size,
      .by = ssf
    )
  dplyr::bind_rows(ici_summary_all, ici_summary_group) |>
    dplyr::mutate(method = "median-ici")
}

create_pca_median_ici_heatmap = function(
  pca_group,
  ssf_colors,
  cor_colors,
  guide_nrow,
  guide_ncol
) {
  if (is.null(pca_group)) {
    return(list(
      eval = FALSE,
      pca = NULL,
      heatmap = NULL,
      median_ici = NULL
    ))
  }
  # pca_group = in_pca_set$pca_groups$all
  pca_plot = pca_group$pca$pca |>
    ggplot(aes(x = PC1, y = PC2, color = ssf, shape = outlier)) +
    geom_point(size = 3) +
    scale_color_manual(values = ssf_colors) +
    theme(
      legend.position = "top",
      legend.location = "plot",
      legend.justification.top = "left",
      plot.title.position = "plot"
    ) +
    guides(
      color = guide_legend(ncol = guide_ncol, nrow = guide_nrow),
      shape = guide_legend(ncol = 1)
    )

  heatmap_data = pca_group$order

  sample_info = pca_group$pca$pca

  sample_info = dplyr::left_join(
    data.frame(sample_id = rownames(heatmap_data$ici)),
    sample_info,
    by = "sample_id"
  )

  heatmap_colors = list(
    ssf = ssf_colors,
    outlier = c("FALSE" = "#000000", "TRUE" = "#FFFFFF")
  )

  left_var = ComplexHeatmap::HeatmapAnnotation(
    df = sample_info[, c("ssf", "outlier")],
    col = heatmap_colors,
    which = "row",
    show_annotation_name = FALSE
  )
  bottom_var = ComplexHeatmap::HeatmapAnnotation(
    df = sample_info[, c("ssf", "outlier")],
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
    show_row_names = FALSE,
    cluster_columns = pca_group$order$dendrogram$dendrogram,
    cluster_rows = pca_group$order$dendrogram$dendrogram
  )

  out_sina = pca_group$pca$pca |>
    ggplot(aes(
      x = median_ici_correlation,
      y = ssf,
      color = ssf,
      shape = outlier,
      group = ssf
    )) +
    geom_sina(size = 3) +
    scale_color_manual(values = ssf_colors) +
    theme(
      legend.position = "top",
      legend.location = "plot",
      legend.justification.top = "left",
      plot.title.position = "plot"
    ) +
    guides(
      color = guide_legend(nrow = guide_nrow, ncol = guide_ncol),
      shape = guide_legend(ncol = 1)
    )
  return(list(
    eval = TRUE,
    pca = pca_plot,
    heatmap = out_heatmap,
    median_ici = out_sina
  ))
}

create_boxplot = function(
  abundance_matrix,
  sample_metadata,
  point_colors_fill,
  guide_nrow,
  guide_ncol
) {
  # abundance_matrix = raw_abundance; point_colors_fill = factor_colors
  sample_metadata = sample_metadata |>
    dplyr::arrange(ssf, sample_id)
  abundance_matrix = abundance_matrix[, sample_metadata$sample_id, drop = FALSE]

  sample_metadata$factors = factor(sample_metadata$ssf)

  abundance_df = abundance_matrix |>
    as.data.frame() |>
    tidyr::pivot_longer(cols = everything(), names_to = "sample_id")
  abundance_df = dplyr::left_join(
    abundance_df,
    sample_metadata[, c("sample_id", "ssf")],
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
    geom_sina(aes(color = ssf)) +
    geom_boxplot(alpha = 0) +
    scale_color_manual(values = point_colors_fill) +
    labs(x = "log2(Abundance)") +
    theme(
      legend.position = "top",
      legend.location = "plot",
      legend.justification.top = "left",
      plot.title.position = "plot"
    ) +
    guides(color = guide_legend(nrow = guide_nrow, ncol = guide_ncol))
  out_boxplot
}


create_sina_boxplot = function(
  abundance_matrix,
  sample_metadata,
  point_colors_fill,
  guide_nrow,
  guide_ncol
) {
  # abundance_matrix = raw_abundance; point_colors_fill = factor_colors
  # abundance_matrix = dynamic_ranges; point_colors_fill = ssf_colors
  sample_metadata = sample_metadata |>
    dplyr::arrange(ssf, sample_id)
  abundance_matrix = abundance_matrix[, sample_metadata$sample_id, drop = FALSE]

  sample_metadata$ssf = factor(sample_metadata$ssf)

  abundance_df = abundance_matrix |>
    as.data.frame() |>
    tidyr::pivot_longer(cols = everything(), names_to = "sample_id")
  abundance_df = dplyr::left_join(
    abundance_df,
    sample_metadata[, c("sample_id", "ssf")],
    by = "sample_id"
  )

  out_boxplot = abundance_df |>
    ggplot(aes(x = value, y = ssf, group = ssf)) +
    geom_boxplot(alpha = 0) +
    geom_sina(aes(color = ssf)) +
    scale_color_manual(values = point_colors_fill) +
    labs(x = "Dynamic Range") +
    theme(
      legend.position = "top",
      legend.location = "plot",
      legend.justification.top = "left",
      plot.title.position = "plot"
    ) +
    guides(color = guide_legend(nrow = guide_nrow, ncol = guide_ncol))
  out_boxplot
}


sets_extract_grouping_scores = function(pca_sets) {
  # pca_sets = tar_read(pca_g03); tar_source("R")
  out_scores = purrr::map(pca_sets, extract_grouping_scores) |>
    purrr::list_rbind()
  out_scores
}

extract_grouping_scores = function(in_pca) {
  analysis_id = in_pca$id
  if (is.null(in_pca[["pca_groups"]])) {
    return(NULL)
  }
  out_df = purrr::imap(
    in_pca$pca_groups[c("all", "not_outlier", "not_qc", "not_outlier_qc")],
    \(pca_group, id) {
      # id = "all"
      # pca_group = in_pca$pca_groups[[id]]
      if (!is.null(pca_group)) {
        tmp_score = pca_group$scores
        tmp_score$sample_set = id
        tmp_score$analysis_id = analysis_id
      } else {
        tmp_score = tibble::tibble(
          n_differences = NA,
          group_size = NA,
          differences_fraction = NA,
          factors = NA,
          method = NA,
          sample_set = id,
          analysis_id = analysis_id
        )
      }
      return(tmp_score)
    }
  ) |>
    purrr::list_rbind()
  out_df
}
