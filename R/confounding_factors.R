detect_confounding_factors = function(se_data) {
  # sets_pca = tar_read(se_pca_g34); se_data = sets_pca[["AN000923"]]; tar_source("R")
  report_list = check_qc_outputs(se_data)

  pca_groups = get_specific_metadata(se_data, "pca_groups")
  if (!is.null(report_list$pca)) {
    use_group = pca_groups[report_list$pca[length(report_list$pca)]]
  } else {
    return(NULL)
  }

  pca_scores = use_group[[1]]$pca$pca

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
  clustering_pcs
}

sets_detect_confounding_factors = function(sets_pca) {
  out_factors = purrr::map(sets_pca, detect_confounding_factors) |>
    purrr::list_rbind()
  out_factors
}
