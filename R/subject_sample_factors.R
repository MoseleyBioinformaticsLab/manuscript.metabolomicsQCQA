widen_ssf_json = function(ssf_data2, id_cols = NULL) {
  if (is.null(id_cols)) {
    id_cols = c("subject", "sample_id")
  }
  ssf_data2 |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      id_cols = tidyselect::all_of(id_cols),
      names_from = factor,
      values_from = item
    ) |>
    janitor::clean_names(unique_sep = "...")
}

create_ssf_groups_json = function(ssf_wide, id_cols = NULL) {
  if (is.null(id_cols)) {
    id_cols = c("sample_id", "subject")
  }

  concat_ssf = ssf_wide |>
    tidyr::unite("ssf", !tidyselect::all_of(id_cols), sep = "|")

  concat_output = dplyr::left_join(
    ssf_wide,
    concat_ssf,
    by = id_cols
  )

  concat_output
}

check_ssf_json = function(ssf_data2, min_ssf = 2, min_n = 5) {
  ssf_classes = purrr::map_chr(ssf_data2, class)
  if (any(ssf_classes %in% "list")) {
    return("list factor")
  }
  ssf_n = ssf_data2 |>
    dplyr::summarise(n = dplyr::n(), .by = factors)

  if (any(ssf_n$n >= min_n) && (nrow(ssf_n) >= min_ssf)) {
    return("good")
  } else {
    return("not enough samples")
  }
}
