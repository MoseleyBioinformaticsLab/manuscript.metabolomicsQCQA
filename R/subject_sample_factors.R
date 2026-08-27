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

# this function checks if a factor has all unique values equal to the
# rows of the ssf data itself, and alternatively is just a single group.
# if equal to number of rows, tries removing the numbers at the end and
# modifying the data
check_single_factors_to_remove = function(
  ssf_ordered,
  other_factors,
  to_remove_factors,
  factor_type
) {
  other_is_single = purrr::map_lgl(other_factors, \(in_factor) {
    n_uniq = length(unique(ssf_ordered[[in_factor]]))
    n_uniq == nrow(ssf_ordered)
  })
  other_is_all = purrr::map_lgl(other_factors, \(in_factor) {
    n_uniq = length(unique(ssf_ordered[[in_factor]]))
    n_uniq == 1
  })

  if (sum(other_is_all) == length(other_factors)) {
    return(list(
      ssf_concat = ssf_ordered,
      other_factors = NULL,
      factor_type = factor_type
    ))
  }
  if (sum(other_is_single) > 0) {
    other_factors = other_factors[!other_factors[other_is_all]]
    out_list = try_removing_trailing_numbers(
      ssf_ordered,
      other_factors,
      factor_type
    )
  } else {
    to_remove_factors = c(to_remove_factors, other_factors[other_is_single])
    other_factors = base::setdiff(other_factors, to_remove_factors)
    return(list(
      other_factors = other_factors,
      to_remove_factors = to_remove_factors
    ))
  }
}


# this function checks the factors for trailing numbers that may be causing
# the split, and removes them and puts them into their own "replicate" field
# the numbers in this case must be at the end of the string, preceded by a -, _, or " ".
#
# Maybe we should be running this first, then doing the concatenation to SSF
# and checking how many runs of one sample are still present, and then
# try dropping if necessary.

# this should be combined with the singleton test, and we only try this on
# things that are necessary.
# We do them both at the beginning, and abort even ealier if we need to
# while keeping at least a single factor to put into the SSF
# finally, we should just keep track of them from the factor_type
# variable, instead of adding and taking away from other factors.
try_removing_trailing_numbers = function(ssf_data, other_factors, factor_type) {
  has_mostly_end_number = purrr::map_lgl(other_factors, \(in_factor) {
    n_has = stringr::str_detect(
      ssf_data[[in_factor]],
      "\\_[[:digit:]]+$|\\-[[:digit:]]+$|\\ [[:digit:]]+$"
    )
    sum(n_has) / nrow(ssf_data) > 0.8
  })
  # we only want to try this process if there are some factors with the condition present
  # otherwise, we just return the data we were given
  if (sum(has_mostly_end_number) == 0) {
    return(list(
      ssf = ssf_data,
      other_factors = other_factors,
      factor_type = factor_type
    ))
  }

  try_to_split = other_factors[has_mostly_end_number]
  out_new_factors = purrr::map(try_to_split, \(in_factor) {
    if (any(stringr::str_detect(ssf_data[[in_factor]], "\\_[[:digit:]]+$"))) {
      use_split = "\\_"
      together = "_"
    } else if (
      any(stringr::str_detect(ssf_data[[in_factor]], "\\-[[:digit:]]+$"))
    ) {
      use_split = "\\-"
      together = "-"
    } else if (
      any(stringr::str_detect(ssf_data[[in_factor]], "\\ [[:digit:]]+$"))
    ) {
      use_split = "\\ "
      together = " "
    } else {
      return(NULL)
    }
    split_factor = stringr::str_split(ssf_data[[in_factor]], use_split)

    new_factor = purrr::map(split_factor, \(in_split) {
      n_in = length(in_split)
      tibble::tibble(
        id = paste0(in_split[1:(n_in - 1)], collapse = together),
        id_rep = in_split[n_in]
      )
    }) |>
      purrr::list_rbind()
    nchar_new_factor = nchar(new_factor)
    if (all(nchar_new_factor == 0)) {
      return(NULL)
    }
    new_id = c(in_factor, paste0(in_factor, "_replicate"))
    names(new_factor) = new_id
    new_factor
  }) |>
    purrr::list_cbind()

  if (ncol(out_new_factors) == 0) {
    return(list(
      ssf = ssf_data,
      other_factors = other_factors,
      factor_type = factor_type
    ))
  }
  ssf_cols = colnames(ssf_data)
  new_cols = colnames(out_new_factors)

  match_cols = base::intersect(ssf_cols, new_cols)
  not_match_cols = base::setdiff(new_cols, ssf_cols)

  for (imatch in match_cols) {
    ssf_data[[imatch]] = new_cols[[imatch]]
    type_loc = factor_type$factor %in% imatch
    factor_type$source[imatch] = "modified by group detection"
  }
  not_match_data = new_cols[, not_match_cols, drop = FALSE]
  ssf_data = dplyr::bind_cols(ssf_data, not_match_data)
  factor_type = dplyr::bind_rows(
    factor_type,
    tibble::tibble(
      factor = colnames(not_match_data),
      type = "additional_sample_data",
      source = "created by group detection"
    )
  )

  return(list(
    ssf = ssf_data,
    other_factors = other_factors,
    factor_type = factor_type
  ))
  # now we need to fix the factor types
}

create_ssf_modify_factor_type = function(
  ssf_data,
  id_cols,
  factor_type,
  remove_factors
) {
  ssf_concat = create_ssf_groups_json(ssf_data, id_cols = id_cols)
  if (length(remove_factors) > 0) {
    for (iremove in remove_factors) {
      match_loc = factor_type$factor %in% iremove
      factor_type[match_loc, "type"] = "additional_sample_data"
      factor_type[match_loc, "source"] = "modified by group detection"
    }
  }
  return(list(ssf_concat = ssf_concat, factor_type = factor_type))
}

determine_extra_ssf = function(ssf_data, id_cols = NULL, factor_type) {
  # ssf_data = ssf_wide |> dplyr::select(tidyselect::all_of(use_for_factors))
  if (is.null(id_cols)) {
    id_cols = c("sample_id", "subject")
  }
  to_remove_factors = NULL

  if (nrow(ssf_data) == 1) {
    ssf_concat = create_ssf_groups_json(ssf_data, id_cols = id_cols)
    return(list(ssf_concat = ssf_concat, removed = NULL))
  }

  other_factors = base::setdiff(colnames(ssf_data), id_cols)

  raw_file_match = grep(
    ".*raw_file_name.*|rawfilename",
    other_factors,
    value = TRUE
  )
  if (length(raw_file_match) > 0) {
    to_remove_factors = raw_file_match
    other_factors = base::setdiff(other_factors, to_remove_factors)
  } else {
    to_remove_factors = NULL
  }

  use_concat = c(id_cols, other_factors)
  if (length(other_factors) == 0) {
    return(
      create_ssf_modify_factor_type(ssf_data |> tidyselect::all_of(use_concat)),
      id_cols,
      factor_type,
      to_remove_factors
    )
  }

  ssf_concat = create_ssf_groups_json(
    ssf_data |> dplyr::select(tidyselect::all_of(use_concat)),
    id_cols = id_cols
  )

  ssf_ordered = ssf_concat |>
    dplyr::arrange(ssf)

  n_each = rle_2_df(rle(ssf_ordered$ssf))
  n_1 = sum(n_each$length == 1)
  n_sum = sum(n_each$length)
  n1_sum = n_1 / n_sum
  if (n1_sum > 0.5) {
    other_factor_list = check_single_factors_to_remove(
      ssf_ordered,
      other_factors,
      to_remove_factors
    )
    other_factors = other_factor_list$other_factors
    to_remove_factors = other_factor_list$to_remove_factors
    if (length(other_factors) > 0) {
      alt_groups = purrr::map(other_factors, \(remove_factor) {
        keep_factors = base::setdiff(colnames(ssf_data), remove_factor)
        ssf_alt = create_ssf_groups_json(
          ssf_data |> dplyr::select(all_of(c(id_cols, keep_factors))),
          id_cols
        ) |>
          dplyr::arrange(ssf)
        alt_each = rle_2_df(rle(ssf_alt$ssf))
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
          ssf_alt = list(ssf_alt)
        )
      }) |>
        purrr::list_rbind()
    } else {
      return(list(ssf_concat = ssf_concat, remove_ssf = to_remove_factors))
    }
  } else {
    return(list(ssf_concat = ssf_concat, remove_ssf = to_remove_factors))
  }
  # if nothing fixes the problem, or creates a single group, then just return what we started with
  dist_groups = alt_groups |>
    dplyr::select(n1_sum, n_all, n10, med_size, n_groups) |>
    dplyr::distinct()

  if (
    all(alt_groups$n1_sum == 1) ||
      (all(alt_groups$n_all == 1)) ||
      (nrow(dist_groups) == 1)
  ) {
    return(list(ssf_concat = ssf_concat, remove_ssf = to_remove_factors))
  }
  alt_groups_low_n1 = alt_groups |>
    dplyr::filter(n1_sum <= 0.5)

  if (nrow(alt_groups_low_n1) == 1) {
    return(list(
      ssf_concat = alt_groups_low_n1$ssf_alt[[1]],
      remove_ssf = c(to_remove_factors, alt_groups_low_n1$removed)
    ))
  }
  alt_groups_low_n10 = alt_groups_low_n1 |>
    dplyr::filter(n10 <= 0.5)

  if (nrow(alt_groups_low_n10) == 1) {
    return(list(
      ssf_concat = alt_groups_low_n10$ssf_alt[[1]],
      remove_ssf = c(to_remove_factors, alt_groups_low_n10$removed)
    ))
  } else if (nrow(alt_groups_low_n10) == 0) {
    alt_groups_return = alt_groups_low_n1 |>
      dplyr::arrange(n10, n1_sum) |>
      dplyr::slice_head(n = 1)
    return(list(
      ssf_concat = alt_groups_return$ssf_alt[[1]],
      remove_ssf = c(to_remove_factors, alt_groups_return$removed)
    ))
  } else {
    alt_groups_return = alt_groups_low_n10 |>
      dplyr::arrange(n10, n1_sum, dplyr::desc(med_size), n_groups) |>
      dplyr::slice_head(n = 1)
    return(list(
      ssf_concat = alt_groups_return$ssf_alt[[1]],
      remove_ssf = c(to_remove_factors, alt_groups_return$removed)
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
