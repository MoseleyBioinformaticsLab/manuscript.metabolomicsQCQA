find_pca_remove_gray = function(sets_qc_color) {
  keep_status = overall_status |>
    dplyr::filter(
      completed_steps %in% "se; dynamic_range; correlation; pca",
      !overall_status %in% "gray"
    )
  keep_status
}

parse_statuses = function(not_gray_pca) {
  # keep_status = tar_read(not_gray_pca); tar_source("R");
  status_cols = grep("status$", names(not_gray_pca), value = TRUE)
  status_cols = grep(
    "^summarized|^has_qc|^score_icikt|^score_rle|^ssf_metadata|^overall",
    status_cols,
    invert = TRUE,
    value = TRUE
  )
  count_statuses = not_gray_pca |>
    dplyr::summarise(N = dplyr::n(), .by = tidyselect::all_of(status_cols)) |>
    dplyr::arrange(dplyr::desc(N))

  count_statuses_with_match = purrr::map(
    seq_len(nrow(count_statuses)),
    \(irow) {
      find_matching_statuses(not_gray_pca, count_statuses[irow, ])
    }
  ) |>
    purrr::list_rbind()
  count_statuses_with_match
}

find_matching_statuses = function(keep_status, status_definition) {
  # keep_status = tar_read(not_gray_pca); status_definition = count_statuses[1, ]; tar_source("R")
  match_status = rep(TRUE, nrow(keep_status))

  use_cols = grep(
    "N",
    names(status_definition),
    ignore.case = FALSE,
    invert = TRUE,
    value = TRUE
  )
  for (icol in use_cols) {
    match_status = match_status &
      (keep_status[[icol]] %in% status_definition[[icol]])
  }

  testthat::expect_equal(sum(match_status), status_definition[["N"]])
  out_status = keep_status[match_status, ]
  status_definition$match = list(out_status)
  status_definition
}
