determine_experimental_design = function(in_ssf) {
  # in_dataset = tar_read(cor_g06)[[1]]
  extra_vars = c(
    "experiment_or_qc",
    "median_ici_correlation",
    "outlier_fraction",
    "score",
    "outlier"
  )
  in_ssf = in_ssf |>
    dplyr::select(!(tidyselect::all_of(extra_vars)))
  factor_info = count_possible_factors(in_ssf)

  n_rows = nrow(factor_info)

  use_factors = factor_info$factor[seq_len(n_rows - 1)]

  n_by_factor = in_ssf |>
    dplyr::summarise(
      `N-Samples` = dplyr::n(),
      .by = tidyselect::all_of(use_factors)
    )
  all_factors = c(use_factors, "factors")
  n_by_factor_factors = in_ssf |>
    dplyr::select(tidyselect::all_of(all_factors)) |>
    dplyr::summarise(
      `N-Samples` = dplyr::n(),
      .by = tidyselect::all_of(all_factors)
    )

  summary_df = tibble::tibble(
    n_factors = n_rows - 1,
    n_groups = factor_info$n_group[n_rows]
  )
  list(
    summary = summary_df,
    factor_counts = factor_info,
    n_by_factor = n_by_factor,
    n_by_factor_factors = n_by_factor_factors
  )
}

count_possible_factors = function(ssf) {
  left_over = ssf |>
    dplyr::select(c(-subject, -sample_id, -factors))
  n_left_over = ncol(left_over)
  n_each = purrr::imap(left_over, \(in_col, id) {
    tibble::tibble(factor = id, n_group = length(unique(in_col)))
  }) |>
    purrr::list_rbind()
  factor_info = tibble::tibble(
    factor = "factors",
    n_group = length(unique(ssf$factors))
  )
  dplyr::bind_rows(n_each, factor_info)
}
