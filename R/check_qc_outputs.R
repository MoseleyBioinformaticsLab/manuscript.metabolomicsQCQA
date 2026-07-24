check_qc_outputs = function(se_data) {
  # se_data = readRDS("tests_problems/pca_tempfile_g01_1.rds")
  last_step = extract_qc_last_step(metadata(se_data)$qc)

  if (nchar(last_step$succeeded) > 0) {
    last_succeeded = last_step$succeeded
  } else {
    last_succeeded = "none"
  }
  if (nchar(last_step$failed) > 0) {
    last_failed = last_step$failed
    why_failed = last_step$why
  } else if (last_step$do_next %in% "stop") {
    last_failed = "none"
    why_failed = last_step$why
  } else {
    last_failed = "none"
    why_failed = "none"
  }

  if (last_succeeded %in% "pca") {
    report_steps = c("se", "dynamic_range", "correlation", "pca")
    pca_outputs = metadata(se_data)$pca

    pca_has_results = purrr::map_lgl(
      pca_outputs[c("all", "not_outlier", "not_qc", "not_outlier_qc")],
      \(x) {
        !is.null(x)
      }
    )
    pca_report_results = names(pca_has_results)[pca_has_results]
    report_list = list(
      steps = report_steps,
      pca = pca_report_results,
      failed = c(last = last_failed, why = why_failed)
    )
  }
  if (last_succeeded %in% "calculate_dynamic_range / normalized") {
    report_list = list(
      steps = c("se", "dynamic_range"),
      pca = NULL,
      failed = c(last = last_failed, why = why_failed)
    )
  }
  if (last_succeeded %in% "se") {
    report_list = list(
      steps = "se",
      pca = NULL,
      failed = c(last = last_failed, why = why_failed)
    )
  }
  if (last_succeeded %in% "none") {
    report_list = list(
      steps = "none",
      pca = NULL,
      failed = c(last = last_failed, why = why_failed)
    )
  }
  return(report_list)
}

find_weird_red_status = function(sets_qc_evaluation) {
  red_lc = sets_qc_evaluation |>
    dplyr::filter(
      completed_steps %in% "se; dynamic_range; correlation; pca",
      !pca_anova_pvalue_status %in% "gray"
    ) |>
    dplyr::filter(
      left_censorship_test_status %in% "red",
      pca_anova_pvalue_status %in% "green"
    )
  red_lc
}
