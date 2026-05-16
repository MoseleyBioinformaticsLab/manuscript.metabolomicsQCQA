get_vardata = function(abundance_se, which = "cols") {
  if (which %in% "cols") {
    colData(abundance_se) |> as.data.frame()
  } else {
    rowData(abundance_se) |> as.data.frame()
  }
}

sets_run_pca_experimental_design = function(lc_sets) {
  # lc_sets = tar_read(lc_g07)
  pca_sets = purrr::map(lc_sets, run_pca_experimental_design)
}

run_pca_experimental_design = function(in_lc) {
  # in_lc = lc_sets[[1]]
  pca_res = do_pca(in_lc$abundance)
  design = in_lc$abundance |> get_vardata() |> determine_experimental_design()
  in_lc$pca_all = pca_res
  in_lc$design_all = design

  possible_outliers = determine_outliers(
    in_lc$icikt_cor,
    get_vardata(in_lc$abundance)
  )

  no_qc = in_lc$abundance[, in_lc$abundance$experiment_or_qc %in% "experiment"]
  not_outlier = possible_outliers |>
    dplyr::filter(!outlier) |>
    dplyr::pull(sample_id)
  no_outlier = in_lc$abundance[,
    in_lc$abundance$sample_id %in% not_outlier
  ]
  not_outlier_qc = base::intersect(
    get_vardata(in_lc$abundance) |>
      dplyr::filter(experiment_or_qc %in% "experiment") |>
      dplyr::pull(sample_id),
    not_outlier
  )
  pca_noqc = do_pca(no_qc)
  in_lc
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
