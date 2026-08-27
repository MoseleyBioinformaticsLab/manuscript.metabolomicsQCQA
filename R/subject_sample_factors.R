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

check_if_each = function(in_factor) {
  n_uniq = length(unique(in_factor))
  n_uniq == length(in_factor)
}

check_if_single = function(in_factor) {
  length(unique(in_factor)) == 1
}

rename_no_trailing = function(in_df, id) {
  new_names = c(id, paste0(id, "_rep"))
  names(in_df) = new_names
  in_df
}

# this function checks if a factor has all unique values equal to the
# rows of the ssf data itself, and alternatively is just a single group.
# if equal to number of rows, tries removing the numbers at the end and
# modifying the data
# this should actually be run right at the beginning of things
process_each_and_all_factors = function(
  ssf_data,
  factor_type
) {
  factors_to_check = get_factors_of_type(factor_type)
  is_each = purrr::map_lgl(
    ssf_data[, factors_to_check, drop = FALSE],
    check_if_each
  )
  is_single = purrr::map_lgl(
    ssf_data[, factors_to_check, drop = FALSE],
    check_if_single
  )

  if (sum(is_single) == length(factors_to_check)) {
    return(list(
      ssf_data = ssf_data,
      factor_type = factor_type
    ))
  }
  if (sum(is_each) > 0) {
    factor_is_each = factors_to_check[is_each]
    no_trailing_each = purrr::map(
      factor_is_each,
      \(use_factor) {
        remove_trailing_numbers(ssf_data[[use_factor]], use_factor)
      }
    ) |>
      purrr::list_cbind()

    is_each_2 = purrr::map_lgl(
      no_trailing_each[, factor_is_each, drop = FALSE],
      check_if_each
    )

    if (sum(!is_each_2) > 0) {
      keep_each = paste0(factor_is_each[!is_each_2], c("", "_replicate"))
      no_trailing_each = no_trailing_each[, keep_each, drop = FALSE]
    }

    if (ncol(no_trailing_each) > 0) {
      return(replace_and_add_data_factors(
        ssf_data,
        factor_type,
        no_trailing_each
      ))
    }
  } else {
    return(list(ssf_data = ssf_data, factor_type = factor_type))
  }
}

replace_and_add_data_factors = function(
  ssf_data,
  factor_type,
  new_data
) {
  has_match = base::intersect(colnames(ssf_data), colnames(new_data))
  not_match = base::setdiff(colnames(new_data), colnames(ssf_data))

  match_source = "modified by group detection"
  notmatch_source = "created by group detection"
  for (imatch in has_match) {
    ssf_data[[imatch]] = new_data[[imatch]]
    fac_loc = factor_type$factor %in% imatch
    factor_type$source[fac_loc] = match_source
  }

  if (length(not_match) > 0) {
    extra_data = new_data[, not_match, drop = FALSE]
    extra_factors = tibble::tibble(
      factor = colnames(extra_data),
      type = "additional_sample_data",
      source = notmatch_source
    )
    ssf_data = dplyr::bind_cols(ssf_data, extra_data)
    factor_type = dplyr::bind_rows(factor_type, extra_factors)
  }
  return(list(ssf_data = ssf_data, factor_type = factor_type))
}

remove_trailing_numbers = function(in_factor, id) {
  if (any(stringr::str_detect(in_factor, "\\_[[:digit:]]+$"))) {
    use_split = "\\_"
    together = "_"
  } else if (any(stringr::str_detect(in_factor, "\\-[[:digit:]]+$"))) {
    use_split = "\\-"
    together = "-"
  } else if (any(stringr::str_detect(in_factor, "\\ [[:digit:]]+$"))) {
    use_split = "\\ "
    together = " "
  } else {
    # if there is nothing here, we just return it
    # essentially as is
    new_factor = tibble::tibble(id = in_factor)
    names(new_factor) = id
    return(new_factor)
  }
  split_factor = stringr::str_split(in_factor, use_split)

  new_factor = purrr::map(split_factor, \(in_split) {
    n_in = length(in_split)
    tibble::tibble(
      id = paste0(in_split[1:(n_in - 1)], collapse = together),
      id_rep = in_split[n_in]
    )
  }) |>
    purrr::list_rbind()
  new_names = c(id, paste0(id, "_replicate"))
  names(new_factor) = new_names
  new_factor
}
# this function checks the factors for trailing numbers that may be causing
# the split, and removes them and puts them into their own "replicate" field
# the numbers in this case must be at the end of the string, preceded by a -, _, or " ".
#
# Maybe we should be running this first, then doing the concatenation to SSF
# and checking how many runs of one sample are still present, and then
# try dropping if necessary.

create_ssf_modify_factor_type = function(
  ssf_data,
  id_cols,
  factor_type,
  remove_factors = NULL
) {
  if (length(remove_factors) > 0) {
    for (iremove in remove_factors) {
      match_loc = factor_type$factor %in% iremove
      factor_type[match_loc, "type"] = "additional_sample_data"
      factor_type[match_loc, "source"] = "modified by group detection"
    }
  }
  ssf_factors = get_factors_of_type(factor_type)
  use_concat = c(id_cols, ssf_factors)
  ssf_concat = create_ssf_groups_json(
    ssf_data[, use_concat, drop = FALSE],
    id_cols = id_cols
  )

  return(list(ssf_concat = ssf_concat, factor_type = factor_type))
}

determine_extra_ssf = function(ssf_data, id_cols, factor_type) {
  # ssf_data = ssf_wide
  if (is.null(id_cols)) {
    id_cols = c("sample_id", "subject")
  }
  to_remove_factors = NULL

  ssf_factors = get_factors_of_type(factor_type)
  use_concat = c(id_cols, ssf_factors)
  if (nrow(ssf_data) == 1) {
    ssf_concat = create_ssf_groups_json(
      ssf_data[, use_concat, drop = FALSE],
      id_cols = id_cols
    )
    return(list(
      ssf_concat = ssf_concat,
      factor_type = factor_type,
      data = ssf_data
    ))
  }

  raw_file_match = grep(
    ".*raw_file_name.*|rawfilename",
    ssf_factors,
    value = TRUE
  )
  if (length(raw_file_match) > 0) {
    raw_loc = factor_type$factor %in% raw_file_match
    factor_type[raw_loc, "source"] = "modified by group detection"
    factor_type[raw_loc, "type"] = "additional_sample_data"
    ssf_factors = get_factors_of_type(factor_type)
  }

  if (length(ssf_factors) == 0) {
    tmp_res = create_ssf_modify_factor_type(ssf_data, id_cols, factor_type)
    tmp_res$data = ssf_data
    return(
      tmp_res
    )
  }

  remove_each_all = process_each_and_all_factors(ssf_data, factor_type)
  ssf_data = remove_each_all$ssf_data
  factor_type = remove_each_all$factor_type

  ssf_factors = get_factors_of_type(factor_type)

  use_concat = c(id_cols, ssf_factors)
  ssf_concat = create_ssf_groups_json(
    ssf_data[, use_concat, drop = FALSE],
    id_cols = id_cols
  )

  ssf_ordered = ssf_concat |>
    dplyr::arrange(ssf)

  n_each = rle_2_df(rle(ssf_ordered$ssf))
  n_1 = sum(n_each$length == 1)
  n_sum = sum(n_each$length)
  n1_sum = n_1 / n_sum
  if (n1_sum > 0.5 && length(ssf_factors) > 0) {
    alt_groups = purrr::map(ssf_factors, \(remove_factor) {
      ssf_alt_list = create_ssf_modify_factor_type(
        ssf_data,
        id_cols,
        factor_type,
        remove_factor
      )
      ssf_alt_list$data = ssf_data

      tmp_concat = ssf_alt_list$ssf_concat |> dplyr::arrange(ssf)
      alt_each = rle_2_df(rle(tmp_concat$ssf))
      alt_1 = sum(alt_each$length == 1)
      alt_sum = sum(alt_each$length)
      alt_group = nrow(alt_each)
      alt_g10 = sum(alt_each$length > (alt_sum * 0.1))
      alt_all = sum(alt_each$length == alt_sum)
      alt_med_size = median(alt_each$length)
      alt_n_groups = nrow(alt_each)
      tibble::tibble(
        removed = remove_factor,
        n1_sum = alt_1 / alt_sum,
        n_all = alt_all / alt_sum,
        n10 = alt_g10 / alt_group,
        med_size = alt_med_size,
        n_groups = alt_n_groups,
        ssf_alt = list(ssf_alt_list)
      )
    }) |>
      purrr::list_rbind()

    dist_groups = alt_groups |>
      dplyr::select(n1_sum, n_all, n10, med_size, n_groups) |>
      dplyr::distinct()

    if (
      all(alt_groups$n1_sum == 1) ||
        (all(alt_groups$n_all == 1)) ||
        (nrow(dist_groups) == 1)
    ) {
      return(list(
        ssf_concat = ssf_concat,
        factor_type = factor_type,
        data = ssf_data
      ))
    }
    alt_groups_low_n1 = alt_groups |>
      dplyr::filter(n1_sum <= 0.5)

    if (nrow(alt_groups_low_n1) == 1) {
      return(alt_groups_low_n1$ssf_alt[[1]])
    } else if (nrow(alt_groups_low_n1) == 0) {
      return(list(
        ssf_concat = ssf_concat,
        factor_type = factor_type,
        data = ssf_data
      ))
    }
    alt_groups_low_n10 = alt_groups_low_n1 |>
      dplyr::filter(n10 <= 0.5)

    if (nrow(alt_groups_low_n10) == 1) {
      return(alt_groups_low_n10$ssf_alt[[1]])
    } else if (nrow(alt_groups_low_n10) == 0) {
      alt_groups_return = alt_groups_low_n1 |>
        dplyr::arrange(n10, n1_sum) |>
        dplyr::slice_head(n = 1)
      return(alt_groups_return$ssf_alt[[1]])
    } else {
      alt_groups_return = alt_groups_low_n10 |>
        dplyr::arrange(n10, n1_sum, dplyr::desc(med_size), n_groups) |>
        dplyr::slice_head(n = 1)
      return(alt_groups_return$ssf_alt[[1]])
    }
  } else {
    # message("last")
    return(list(
      ssf_concat = ssf_concat,
      factor_type = factor_type,
      data = ssf_data
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
