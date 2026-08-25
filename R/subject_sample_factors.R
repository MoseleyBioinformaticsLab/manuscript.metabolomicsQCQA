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

determine_extra_ssf = function(ssf_wide, id_cols = NULL) {
  if (is.null(id_cols)) {
    id_cols = c("sample_id", "subject")
  }

  other_factors = base::setdiff(colnames(ssf_wide), id_cols)
  if (length(other_factors) == 1) {
    ssf_concat = create_ssf_groups_json(ssf_wide, id_cols = id_cols)
    return(list(ssf_concat = ssf_concat, removed = NULL))
  }

  ssf_concat = create_ssf_groups_json(ssf_wide, id_cols = id_cols)

  ssf_ordered = ssf_concat |>
    dplyr::arrange(ssf)

  n_each = rle_2_df(rle(ssf_ordered$ssf))
  n_1 = sum(n_each$length == 1)
  n_sum = sum(n_each$length)
  n1_sum = n_1 / n_sum
  if (n1_sum > 0.5) {
    alt_groups = purrr::map(other_factors, \(remove_factor) {
      keep_factors = base::setdiff(colnames(ssf_wide), remove_factor)
      ssf_alt = create_ssf_groups_json(
        ssf_wide |> dplyr::select(all_of(c(id_cols, keep_factors))),
        id_cols
      ) |>
        dplyr::arrange(ssf)
      alt_each = rle_2_df(rle(ssf_alt$ssf))
      alt_1 = sum(alt_each$length == 1)
      alt_sum = sum(alt_each$length)
      alt_group = nrow(alt_each)
      alt_g10 = sum(alt_each$length > (alt_sum * 0.1))
      alt_all = sum(alt_each$length == alt_sum)
      tibble::tibble(
        removed = remove_factor,
        n1_sum = alt_1 / alt_sum,
        n_all = alt_all / alt_sum,
        n10 = alt_g10 / alt_group,
        ssf_alt = list(ssf_alt)
      )
    }) |>
      purrr::list_rbind()
  } else {
    return(list(ssf_concat = ssf_concat, remove_ssf = NULL))
  }
  # if nothing fixes the problem, or creates a single group, then just return what we started with
  if (all(alt_groups$n1_sum == 1) || (all(alt_groups$n_all == 1))) {
    return(list(ssf_concat = ssf_concat, remove_ssf = NULL))
  }
  alt_groups_low_n1 = alt_groups |>
    dplyr::filter(n1_sum <= 0.5)

  if (nrow(alt_groups_low_n1) == 1) {
    return(list(
      ssf_concat = alt_groups_low_n1$ssf_alt[[1]],
      remove_ssf = alt_groups_low_n1$removed
    ))
  }
  alt_groups_low_n10 = alt_groups_low_n1 |>
    dplyr::filter(n10 <= 0.5)

  if (nrow(alt_groups_low_n10) == 1) {
    return(list(
      ssf_concat = alt_groups_low_n10$ssf_alt[[1]],
      remove_ssf = alt_groups_low_n10$removed
    ))
  } else if (nrow(alt_groups_low_n10) == 0) {
    alt_groups_return = alt_groups_low_n1 |>
      dplyr::arrange(n10, n1_sum) |>
      dplyr::slice_head(n = 1)
    return(list(
      ssf_concat = alt_groups_return$ssf_alt[[1]],
      remove_ssf = alt_groups_return$removed
    ))
  }
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
