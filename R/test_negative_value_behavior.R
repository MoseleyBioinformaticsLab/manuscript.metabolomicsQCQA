test_negative_value_behavior = function() {
  # tar_source("R")
  # qc_run1 = tar_read(run_g01)
  # qc_run2 = tar_read(run_g02)
  # qc_run3 = tar_read(run_g03)

  all_runs = c(qc_run1, qc_run2, qc_run3)
  qc_info = get_qc_outputs(all_runs)

  neg_qc = qc_info |>
    dplyr::filter(why %in% "negative values present") |>
    dplyr::pull(id) |>
    unique()

  neg_list = all_runs[neg_qc]
  out_res = purrr::map(neg_list, get_negative_test) |>
    purrr::list_rbind()
}

get_negative_test = function(in_se) {
  # in_se = neg_list[["AN006242"]]
  raw_data = assays(in_se)$raw
  raw_data[is.na(raw_data)] = 0
  raw_data[raw_data < 0] = NA

  all_test = ICIKendallTau::test_left_censorship(raw_data, global_na = NA)

  ssf_test = try(ICIKendallTau::test_left_censorship(
    raw_data,
    global_na = NA,
    in_se$ssf
  ))
  if (!inherits(ssf_test, "try-error")) {
    binom_res = dplyr::bind_rows(
      broom::tidy(all_test$binomial_test),
      broom::tidy(ssf_test$binomial_test)
    ) |>
      dplyr::mutate(id = metadata(in_se)$id, subset = c("all", "ssf"))
  } else {
    binom_res = broom::tidy(all_test$binomial_test) |>
      dplyr::mutate(id = metadata(in_se)$id, subset = "all")
  }
  binom_res
}
