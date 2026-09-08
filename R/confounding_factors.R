detect_confounding_factors = function(se_data) {
  # sets_pca = tar_read(se_pca_g34); se_data = sets_pca[["AN000923"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g05); se_data = sets_pca[["AN000262"]]; tar_source("R")
  report_list = check_qc_outputs(se_data)

  pca_groups = get_specific_metadata(se_data, "pca_groups")
  if (!is.null(report_list$pca)) {
    use_group = pca_groups[report_list$pca[length(report_list$pca)]]
  } else {
    return(NULL)
  }

  pca_scores = use_group[[1]]$pca$pca
  pca_variance = use_group[[1]]$pca$variance

  # also pull the variances for PC1 and PC2, as I think they
  # will help us determine what is actually *real* and what
  # is not.

  clustering_pcs = purrr::map(c("PC1", "PC2"), \(in_pc) {
    # in_pc = "PC1"
    pc_scores = pca_scores[[in_pc]]
    names(pc_scores) = pca_scores$sample_id

    pc_dist = dist(pc_scores)

    pc_hclust = hclust(pc_dist)
    pc_members = cutree(pc_hclust, k = 2)

    pc_aov = stats::aov(pc_scores ~ pc_members)
    tidy_aov = broom::tidy(pc_aov)[1, ]
    tidy_aov$PC = in_pc
    tidy_aov
  }) |>
    purrr::list_rbind()
  clustering_pcs$analysis_id = metadata(se_data)$mwtab$ID
  clustering_pcs = dplyr::left_join(
    clustering_pcs,
    pca_variance,
    by = c("PC" = "pc")
  )
  clustering_pcs
}

sets_detect_confounding_factors = function(sets_pca) {
  out_factors = purrr::map(sets_pca, detect_confounding_factors) |>
    purrr::list_rbind()
  out_factors
}


add_confounding_factor_status = function(
  sets_qc_color,
  sets_confounding_factors
) {
  sets_qc_confounds = sets_qc_color |>
    dplyr::mutate(confounding_factor_detected = "no")

  sets_yellow_red = sets_qc_color |>
    dplyr::filter(pca_anova_pvalue_status %in% c("yellow", "red"))
  possible_confound = sets_confounding_factors |>
    dplyr::filter(
      analysis_id %in% sets_yellow_red$analysis_id,
      p.value <= 1e-10,
      statistic >= 300,
      percent >= 0.3
    )
}
