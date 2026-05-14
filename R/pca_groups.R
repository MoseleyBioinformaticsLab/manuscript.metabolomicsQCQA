sets_run_pca_experimental_design = function(lc_sets) {
  # lc_sets = tar_read(lc_g07)
  pca_sets = purrr::map(lc_sets, run_pca_experimental_design)
}

run_pca_experimental_design = function(in_lc) {
  # in_lc = lc_sets[[1]]
  pca_res = do_pca(in_lc$abundance)
  design = determine_experimental_design(in_lc$ssf)
  in_lc$pca_all = pca_res
  in_lc$design_all = design
  in_lc
}

do_pca = function(norm_data, ssf_wide) {
  # norm_data = in_lc$normalized |>
  #   get_just_values()
  # ssf_wide = in_lc$ssf
  norm_matrix = as.matrix(norm_data)
  norm_matrix[is.na(norm_matrix)] = 0
  log_matrix = log1p(norm_matrix)
  pca_res = prcomp(t(log_matrix), center = TRUE, scale. = FALSE)

  pca_var = visualizationQualityControl::visqc_score_contributions(pca_res$x)
  pca_anova = visualizationQualityControl::visqc_test_pca_scores(
    pca_res$x,
    ssf_wide
  )
  pca_x = pca_res$x
  pca_df = cbind(as.data.frame(pca_x), ssf_wide)

  n_pc = ncol(pca_x)
  use_pc = colnames(pca_x)
  label_vals = purrr::map(seq(1, n_pc), \(x) {
    pca_var$labels[x]
  })
  names(label_vals) = use_pc

  labelled::var_label(pca_df) = label_vals

  list(pca = pca_df, variance = pca_var, anova = pca_anova)
}
