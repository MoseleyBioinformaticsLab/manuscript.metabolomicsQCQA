find_pca_remove_gray = function(sets_qc_color) {
  keep_status = sets_qc_color |>
    dplyr::filter(
      completed_steps %in% "se; dynamic_range; correlation; pca",
      !overall_status %in% "gray"
    )
  keep_status
}

parse_statuses = function(not_gray_pca) {
  # tar_load(not_gray_pca); tar_source("R");
  status_cols = setup_statuses(names(not_gray_pca))

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
  testthat::expect_equal(length(unique(out_status$overall_status)), 1)
  status_definition$match = list(out_status)
  status_definition$overall = unique(out_status$overall_status)
  status_definition
}

calculate_percentage_overall = function(parsed_grouped_statuses) {
  # tar_load(parsed_grouped_statuses); tar_source("R")
  split_parsed = split(parsed_grouped_statuses, parsed_grouped_statuses$overall)

  overall_percentage = purrr::map(
    split_parsed[c("green", "yellow", "red")],
    \(in_parsed) {
      overall_total = sum(in_parsed$N)
      in_parsed |>
        dplyr::mutate(perc_overall = N / overall_total)
    }
  ) |>
    purrr::list_rbind()
  overall_percentage
}
