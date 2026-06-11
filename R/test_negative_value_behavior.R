test_negative_value_behavior = function(sets_neg) {
  # tar_source("R")
  # tar_load(c(neg_g01, neg_g02, neg_g03))
  # sets_neg = c(neg_g01, neg_g02, neg_g03)

  qc_info = get_qc_outputs(sets_neg)

  neg_qc = qc_info |>
    dplyr::filter(why %in% "negative values present") |>
    dplyr::pull(id) |>
    unique()

  neg_list = sets_neg[neg_qc]
  out_res = purrr::map(neg_list, get_negative_test) |>
    purrr::list_rbind()
  out_res
}

get_negative_test = function(in_se) {
  # in_se = neg_list[[1]]
  raw_data = assays(in_se)$raw
  raw_data[raw_data == 0] = NA

  org_tests = do_lc_tests(raw_data, in_se$ssf, ".org")

  neg_data = raw_data
  neg_data[raw_data < 0] = NA

  neg_tests = do_lc_tests(neg_data, in_se$ssf, ".neg")

  out_tests = dplyr::bind_cols(org_tests, neg_tests) |>
    dplyr::mutate(id = metadata(in_se)$qc$id)
  out_tests
}

do_lc_tests = function(raw_data, ssf, suffix = ".org") {
  # ssf = in_se$ssf; suffix = ".org";
  na_test = tibble::tibble(
    estimate = 0,
    statistic = as.numeric(NA),
    p.value = 1,
    parameter = as.numeric(NA),
    conf.low = as.numeric(NA),
    conf.high = as.numeric(NA),
    method = "Exact binomial test",
    alternative = "greater"
  )

  all_test = try(ICIKendallTau::test_left_censorship(
    raw_data,
    global_na = c(0, NA)
  ))
  if (inherits(all_test, "try-error") || is.null(all_test)) {
    all_out = na_test
  } else {
    all_out = broom::tidy(all_test$binomial_test)
  }
  all_out = all_out |>
    dplyr::mutate(subset = "all")

  ssf_test = try(ICIKendallTau::test_left_censorship(
    raw_data,
    global_na = c(0, NA),
    ssf
  ))
  if (inherits(ssf_test, "try-error") || is.null(ssf_test)) {
    ssf_out = na_test
  } else {
    ssf_out = broom::tidy(ssf_test$binomial_test)
  }
  ssf_out = ssf_out |>
    dplyr::mutate(subset = "ssf")
  test_out = dplyr::bind_rows(all_out, ssf_out)
  names(test_out) = paste0(names(test_out), suffix)
  test_out
}


summarize_neg_tests = function(sets_neg_test) {
  sets_neg_test |>
    dplyr::summarise(
      n = dplyr::n(),
      sig = sum(padjust <= 0.05),
      perc = sig / n,
      .by = subset
    )
}
