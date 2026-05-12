determine_experimental_design = function(in_dataset) {
  # in_dataset = tar_read(cor_g01)[[1]]
  factor_info = count_possible_factors(in_dataset$ssf)
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
