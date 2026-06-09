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
  out_res_all = out_res |>
    dplyr::filter(subset %in% "all") |>
    dplyr::mutate(padjust = p.adjust(p.value, method = "BH"))
  out_res_subset = out_res |>
    dplyr::filter(subset %in% "ssf") |>
    dplyr::mutate(padjust = p.adjust(p.value, method = "BH"))
  dplyr::bind_rows(out_res_all, out_res_subset)
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


summarize_neg_tests = function(sets_neg_test) {
  sets_neg_test |>
    dplyr::summarise(
      n = dplyr::n(),
      sig = sum(padjust <= 0.05),
      perc = sig / n,
      .by = subset
    )
}
