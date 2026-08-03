find_pca_remove_gray = function(sets_qc_evaluation) {
  overall_status = find_green_analyses(sets_qc_evaluation)
  keep_status = overall_status |>
    dplyr::filter(
      completed_steps %in% "se; dynamic_range; correlation; pca",
      !overall_status %in% "gray"
    )
  keep_status
}

parse_statuses = function(keep_status) {}
