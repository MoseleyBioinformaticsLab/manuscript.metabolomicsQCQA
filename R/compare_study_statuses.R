compare_study_statuses = function(not_gray_pca) {
  # tar_load(not_gray_pca)

  not_gray_pca = not_gray_pca |>
    dplyr::filter_out(is.na(study_id))
  by_study = split(
    not_gray_pca[c("internal_id", "overall_status")],
    not_gray_pca$study_id
  )

  n_in_study = purrr::map_int(by_study, nrow)
  keep_studies = by_study[n_in_study > 1]

  comparison_datasets = matrix(0, 9, nrow = 3, ncol = 3)
  rownames(comparison_datasets) = colnames(comparison_datasets) = c(
    "green",
    "yellow",
    "red"
  )

  for (istudy in seq_along(keep_studies)) {
    tmp_study = keep_studies[[istudy]]

    for (irow in seq_len(nrow(tmp_study) - 1)) {
      start_status = tmp_study$overall_status[irow]
      for (jrow in seq((irow + 1), nrow(tmp_study))) {
        end_status = tmp_study$overall_status[jrow]
        comparison_datasets[start_status, end_status] = comparison_datasets[
          start_status,
          end_status
        ] +
          1
      }
    }
  }
  comparison_df = tibble::as_tibble(
    comparison_datasets,
    rownames = "initial_status"
  )
  list(n_studies = length(keep_studies), status_comparisons = comparison_df)
}
