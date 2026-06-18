sets_check_was_transformed = function(sets_neg) {
  # tar_load(sets_neg); tar_source("R");
  transform_status = purrr::map(sets_neg, check_transformation) |>
    purrr::list_rbind()
  transform_status
}

check_transformation = function(in_se) {
  # in_se = sets_neg[[1]]
  raw_data = assays(in_se)$raw

  raw_mean = mean(raw_data, na.rm = TRUE)
  raw_median = median(raw_data, na.rm = TRUE)

  med_mean = tibble::tibble(
    median_mean = abs(raw_mean / raw_median),
    id = metadata(in_se)$qc$id
  ) |>
    dplyr::mutate(
      status = dplyr::case_when(
        (median_mean >= 0.95) && (median_mean <= 1.05) ~ "log-space",
        TRUE ~ "raw-space"
      )
    )
  med_mean
}
