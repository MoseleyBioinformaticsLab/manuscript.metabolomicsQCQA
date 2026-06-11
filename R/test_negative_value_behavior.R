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

get_contains_negatives = function(sets_neg) {
  qc_info = get_qc_outputs(sets_neg)

  neg_qc = qc_info |>
    dplyr::filter(why %in% "negative values present") |>
    dplyr::pull(id) |>
    unique()

  neg_list = sets_neg[neg_qc]
  neg_list
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


plot_neg_tests = function(sets_neg_test) {
  # tar_source("R"); tar_load(sets_neg_test)
  # we identify and remove those that had an NA statistic so we aren't
  # overinflating anything.
  na_org = sets_neg_test |>
    dplyr::filter(is.na(statistic.org) | is.na(statistic.neg)) |>
    dplyr::pull(id) |>
    unique()

  sets_neg_test_present = sets_neg_test |>
    dplyr::filter(!(id %in% na_org)) |>
    dplyr::mutate(
      p.value.org = dplyr::case_when(
        p.value.org == 0 ~ 1e-300,
        TRUE ~ p.value.org
      ),
      p.value.neg = dplyr::case_when(
        p.value.neg == 0 ~ 1e-300,
        TRUE ~ p.value.neg
      )
    ) |>
    dplyr::mutate(
      p.value_neg_vs_org = -1 * log10(p.value.neg / p.value.org),
      estimate_neg_vs_org = log2(estimate.neg / estimate.org)
    )
  p_value_plot = sets_neg_test_present |>
    ggplot(aes(x = p.value_neg_vs_org)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, color = "red") +
    facet_wrap(~subset.neg, ncol = 1) +
    labs(x = "-1 x Log10(Neg-p / Missing-p)")
  estimate_plot = sets_neg_test_present |>
    ggplot(aes(x = estimate_neg_vs_org)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, color = "red") +
    facet_wrap(~subset.neg, ncol = 1) +
    labs(x = "Log2(Neg-estimate / Missing-estimate)")

  sets_neg_direction = sets_neg_test_present |>
    dplyr::mutate(
      sign_p = sign(p.value_neg_vs_org),
      sign_p2 = dplyr::case_when(
        sign_p == 0 ~ 1,
        TRUE ~ sign_p
      ),
      sign_estimate = sign(estimate_neg_vs_org),
      sign_estimate2 = dplyr::case_when(
        sign_estimate == 0 ~ 1,
        TRUE ~ sign_estimate
      )
    )
  sets_p_test_nozero = sets_neg_direction |>
    dplyr::filter(sign_p != 0) |>
    dplyr::summarise(n = dplyr::n(), .by = c(subset.neg, sign_p)) |>
    dplyr::summarise(
      binom_test = binom.test(n) |> broom::tidy(),
      .by = subset.neg
    )

  sets_p_test_zeropos = sets_neg_direction |>
    dplyr::summarise(n = dplyr::n(), .by = c(subset.neg, sign_p2)) |>
    dplyr::summarise(
      binom_test = binom.test(n) |> broom::tidy(),
      .by = subset.neg
    )

  sets_estimate_test_nozero = sets_neg_direction |>
    dplyr::filter(sign_estimate != 0) |>
    dplyr::summarise(n = dplyr::n(), .by = c(subset.neg, sign_estimate)) |>
    dplyr::summarise(
      binom_test = binom.test(n) |> broom::tidy(),
      .by = subset.neg
    )

  sets_estimate_test_zeropos = sets_neg_direction |>
    dplyr::summarise(n = dplyr::n(), .by = c(subset.neg, sign_estimate2)) |>
    dplyr::summarise(
      binom_test = binom.test(n) |> broom::tidy(),
      .by = subset.neg
    )

  list(
    p_histogram = p_value_plot,
    estimate_histogram = estimate_plot,
    p_tests = list(
      remove_zero = sets_p_test_nozero,
      zero_to_positive = sets_p_test_zeropos
    ),
    estimate_tests = list(
      remove_zero = sets_estimate_test_nozero,
      zero_to_positive = sets_estimate_test_zeropos
    )
  )
}
