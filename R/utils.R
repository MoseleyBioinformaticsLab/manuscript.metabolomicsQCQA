grab_last_pca_data = function(se_data) {
  pca_data = get_specific_metadata(se_data, "pca_groups")
  last_run = check_qc_outputs(se_data)
  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]
  pca_data_use
}

count_ancillary_missing = function() {
  tar_load(sets_qc_color_confounds)

  ancillary_and_no_data = sets_qc_color_confounds |>
    dplyr::filter(
      data_source %in% c("ancillary", "none")
    )
  n_data = nrow(ancillary_and_no_data)
  n_miss = sum(ancillary_and_no_data$data_source %in% "none")
}

count_code = function() {
  all_r_files = fs::dir_ls(
    "~/Documents/manuscripts/in_progress/rmflight_metabolomicsDatabaseQCQA/R"
  )

  n_r_files = length(all_r_files)

  file_contents = purrr::map(all_r_files, \(in_file) {
    # in_file = all_r_files[25]
    tmp_contents = readLines(in_file)
    no_comment = grep("\\# ", tmp_contents, value = TRUE, invert = TRUE)
    n_line = length(no_comment)
    n_function = length(grep("function\\(", no_comment))
    tibble::tibble(lines = n_line, functions = n_function)
  }) |>
    purrr::list_rbind()
}

find_analysis_in_group = function(analysis_id) {
  group_sets = readRDS("resources/group_sets.rds")
  all_ids = purrr::map(group_sets$set, \(in_set) {
    in_set
  }) |>
    purrr::list_rbind()

  all_ids |>
    dplyr::filter(id %in% analysis_id)
}

fix_sample_ids = function(raw_ids) {
  # tar_source("R")
  # raw_ids = c(
  #   "Bur-0_W-N-_P3ID33_a",
  #   "Bur-0_W-N-_P3ID33_b",
  #   "Bur-0_W-N-_P3ID33_c",
  #   "Bur-0_W-N-_P3ID33_d",
  #   "Bur-0_W-N-_P3ID33_e",
  #   "Bur-0_W-N-_P3ID33_f",
  #   "Bur-0_W-N-_P3ID34_a",
  #   "Bur-0_W-N-_P3ID34_b",
  #   "Bur-0_W-N-_P3ID34_c",
  #   "Bur-0_W-N-_P3ID34_d",
  #   "Bur-0_W-N-_P3ID34_e",
  #   "Bur-0_W-N+_P3ID33_a",
  #   "Bur-0_W-N+_P3ID33_b",
  #   "Bur-0_W-N+_P3ID33_c",
  #   "Bur-0_W-N+_P3ID33_d",
  #   "Bur-0_W-N+_P3ID34_a",
  #   "Bur-0_W-N+_P3ID34_b",
  #   "Bur-0_W-N+_P3ID34_c",
  #   "Bur-0_W-N+_P3ID34_d",
  #   "Bur-0_W+N-_P3ID33_a",
  #   "Bur-0_W+N-_P3ID33_b",
  #   "Bur-0_W+N-_P3ID33_c",
  #   "Bur-0_W+N-_P3ID33_d",
  #   "Bur-0_W+N-_P3ID34_a",
  #   "Bur-0_W+N-_P3ID34_b",
  #   "Bur-0_W+N-_P3ID34_c",
  #   "Bur-0_W+N-_P3ID34_d",
  #   "Bur-0_W+N+_P3ID33_a",
  #   "Bur-0_W+N+_P3ID33_b",
  #   "Bur-0_W+N+_P3ID33_c",
  #   "Bur-0_W+N+_P3ID33_d",
  #   "Bur-0_W+N+_P3ID34_a",
  #   "Bur-0_W+N+_P3ID34_c",
  #   "Bur-0_W+N+_P3ID34_d",
  #   "Col-0_W-N-_P3ID33_a",
  #   "Col-0_W-N-_P3ID33_b",
  #   "Col-0_W-N-_P3ID33_c",
  #   "Col-0_W-N-_P3ID33_d",
  #   "Col-0_W-N-_P3ID33_e",
  #   "Col-0_W-N-_P3ID33_f",
  #   "Col-0_W-N-_P3ID34_a",
  #   "Col-0_W-N-_P3ID34_b",
  #   "Col-0_W-N-_P3ID34_c",
  #   "Col-0_W-N-_P3ID34_d",
  #   "Col-0_W-N-_P3ID34_e",
  #   "Col-0_W-N+_P3ID33_a",
  #   "Col-0_W-N+_P3ID33_b",
  #   "Col-0_W-N+_P3ID33_c",
  #   "Col-0_W-N+_P3ID33_d",
  #   "Col-0_W-N+_P3ID34_a",
  #   "Col-0_W-N+_P3ID34_b",
  #   "Col-0_W-N+_P3ID34_d",
  #   "Col-0_W+N-_P3ID33_a",
  #   "Col-0_W+N-_P3ID33_b",
  #   "Col-0_W+N-_P3ID33_c",
  #   "Col-0_W+N-_P3ID34_a",
  #   "Col-0_W+N-_P3ID34_b",
  #   "Col-0_W+N-_P3ID34_c",
  #   "Col-0_W+N+_P3ID33_a",
  #   "Col-0_W+N+_P3ID33_b",
  #   "Col-0_W+N+_P3ID33_c",
  #   "Col-0_W+N+_P3ID33_d",
  #   "Col-0_W+N+_P3ID34_a",
  #   "Col-0_W+N+_P3ID34_b",
  #   "Col-0_W+N+_P3ID34_c",
  #   "Col-0_W+N+_P3ID34_d",
  #   "Cvi-0_W-N-_P3ID33_a",
  #   "Cvi-0_W-N-_P3ID33_b",
  #   "Cvi-0_W-N-_P3ID33_c",
  #   "Cvi-0_W-N-_P3ID33_d",
  #   "Cvi-0_W-N-_P3ID33_e",
  #   "Cvi-0_W-N-_P3ID33_f",
  #   "Cvi-0_W-N-_P3ID34_a",
  #   "Cvi-0_W-N-_P3ID34_b",
  #   "Cvi-0_W-N-_P3ID34_c",
  #   "Cvi-0_W-N-_P3ID34_d",
  #   "Cvi-0_W-N-_P3ID34_e",
  #   "Cvi-0_W-N+_P3ID33_a",
  #   "Cvi-0_W-N+_P3ID33_b",
  #   "Cvi-0_W-N+_P3ID33_c",
  #   "Cvi-0_W-N+_P3ID33_d",
  #   "Cvi-0_W-N+_P3ID34_a",
  #   "Cvi-0_W-N+_P3ID34_b",
  #   "Cvi-0_W-N+_P3ID34_c",
  #   "Cvi-0_W-N+_P3ID34_d",
  #   "Cvi-0_W+N-_P3ID33_a",
  #   "Cvi-0_W+N-_P3ID33_b",
  #   "Cvi-0_W+N-_P3ID33_c",
  #   "Cvi-0_W+N-_P3ID34_a",
  #   "Cvi-0_W+N-_P3ID34_b",
  #   "Cvi-0_W+N-_P3ID34_c",
  #   "Cvi-0_W+N-_P3ID34_d",
  #   "Cvi-0_W+N+_P3ID33_a",
  #   "Cvi-0_W+N+_P3ID33_b",
  #   "Cvi-0_W+N+_P3ID33_c",
  #   "Cvi-0_W+N+_P3ID33_d",
  #   "Cvi-0_W+N+_P3ID34_a",
  #   "Cvi-0_W+N+_P3ID34_b",
  #   "Cvi-0_W+N+_P3ID34_c",
  #   "Cvi-0_W+N+_P3ID34_d",
  #   "Sha_W-N-_P3ID33_a",
  #   "Sha_W-N-_P3ID33_b",
  #   "Sha_W-N-_P3ID33_c",
  #   "Sha_W-N-_P3ID33_d",
  #   "Sha_W-N-_P3ID33_e",
  #   "Sha_W-N-_P3ID34_a",
  #   "Sha_W-N-_P3ID34_b",
  #   "Sha_W-N-_P3ID34_c",
  #   "Sha_W-N-_P3ID34_e",
  #   "Sha_W-N+_P3ID33_a",
  #   "Sha_W-N+_P3ID33_b",
  #   "Sha_W-N+_P3ID33_c",
  #   "Sha_W-N+_P3ID33_d",
  #   "Sha_W-N+_P3ID34_a",
  #   "Sha_W-N+_P3ID34_b",
  #   "Sha_W-N+_P3ID34_c",
  #   "Sha_W+N-_P3ID33_a",
  #   "Sha_W+N-_P3ID33_b",
  #   "Sha_W+N-_P3ID33_c",
  #   "Sha_W+N-_P3ID34_a",
  #   "Sha_W+N-_P3ID34_b",
  #   "Sha_W+N-_P3ID34_c",
  #   "Sha_W+N-_P3ID34_d",
  #   "Sha_W+N+_P3ID33_a",
  #   "Sha_W+N+_P3ID33_b",
  #   "Sha_W+N+_P3ID33_c",
  #   "Sha_W+N+_P3ID33_d",
  #   "Sha_W+N+_P3ID34_a",
  #   "Sha_W+N+_P3ID34_b",
  #   "Sha_W+N+_P3ID34_d",
  #   "Tsu-0_W-N-_P3ID33_a",
  #   "Tsu-0_W-N-_P3ID33_b",
  #   "Tsu-0_W-N-_P3ID33_c",
  #   "Tsu-0_W-N-_P3ID33_d",
  #   "Tsu-0_W-N-_P3ID33_e",
  #   "Tsu-0_W-N-_P3ID33_f",
  #   "Tsu-0_W-N-_P3ID34_a",
  #   "Tsu-0_W-N-_P3ID34_b",
  #   "Tsu-0_W-N-_P3ID34_c",
  #   "Tsu-0_W-N-_P3ID34_d",
  #   "Tsu-0_W-N-_P3ID34_e",
  #   "Tsu-0_W-N+_P3ID33_a",
  #   "Tsu-0_W-N+_P3ID33_b",
  #   "Tsu-0_W-N+_P3ID33_c",
  #   "Tsu-0_W-N+_P3ID33_d",
  #   "Tsu-0_W-N+_P3ID34_a",
  #   "Tsu-0_W-N+_P3ID34_b",
  #   "Tsu-0_W-N+_P3ID34_c",
  #   "Tsu-0_W+N-_P3ID33_a",
  #   "Tsu-0_W+N-_P3ID33_b",
  #   "Tsu-0_W+N-_P3ID33_c",
  #   "Tsu-0_W+N-_P3ID33_d",
  #   "Tsu-0_W+N-_P3ID34_a",
  #   "Tsu-0_W+N-_P3ID34_b",
  #   "Tsu-0_W+N-_P3ID34_c",
  #   "Tsu-0_W+N-_P3ID34_d",
  #   "Tsu-0_W+N+_P3ID33_a",
  #   "Tsu-0_W+N+_P3ID33_b",
  #   "Tsu-0_W+N+_P3ID33_c",
  #   "Tsu-0_W+N+_P3ID33_d",
  #   "Tsu-0_W+N+_P3ID34_a",
  #   "Tsu-0_W+N+_P3ID34_b",
  #   "Tsu-0_W+N+_P3ID34_c"
  # )
  # tar_source("R")
  # raw_ids = c(
  #   "Metabolite",
  #   "4T1-1_R7+",
  #   "4T1-15_R11+",
  #   "4T1-19_R15+",
  #   "4T1-27_R18+",
  #   "CTL-1_R23+",
  #   "CTL-5_R2+",
  #   "CTL-9_R1+",
  #   "CTL-3_R6+",
  #   "CTL-7_R3+",
  #   "CTL-11_R21+",
  #   "4T1-17_R4+",
  #   "4T1-21_R22+",
  #   "4T1-25_R16+",
  #   "CTL-4_R20+",
  #   "CTL-8_R13+",
  #   "CTL-12_R14+"
  # )

  # raw_ids = c(
  #   "Metabolite",
  #   "4T1-1_R7-",
  #   "4T1-15_R11-",
  #   "4T1-19_R15-",
  #   "4T1-27_R18-",
  #   "CTL-1_R23-",
  #   "CTL-5_R2-",
  #   "CTL-9_R1-",
  #   "CTL-3_R6-",
  #   "CTL-7_R3-",
  #   "CTL-11_R21-",
  #   "4T1-17_R4-",
  #   "4T1-21_R22-",
  #   "4T1-25_R16-",
  #   "CTL-4_R20-",
  #   "CTL-8_R13-",
  #   "CTL-12_R14-"
  # )

  plus_locs = stringr::str_locate_all(raw_ids, "\\+")

  has_plus = purrr::map_lgl(plus_locs, \(x) {
    nrow(x) > 0
  })
  if (any(has_plus)) {
    plus_all = do.call(rbind, plus_locs)
    plus_diff = plus_all[, 2] - plus_all[, 1]
    if (any(plus_diff != 0)) {
      stop("There are multiple plus signs in a row! Stopping.")
    }
    plus_all = plus_all[, 1]
  }

  nchar_all = nchar(raw_ids)
  minus_locs = stringr::str_locate_all(raw_ids, "\\-")
  has_minus = purrr::map_lgl(minus_locs, \(x) {
    nrow(x) > 0
  })
  if (any(has_minus)) {
    minus_all = do.call(rbind, minus_locs)
    minus_diff = minus_all[, 2] - minus_all[, 1]
    if (any(minus_diff != 0)) {
      stop("There are multiple minus signs in a row! Stopping.")
    }
    minus_all = minus_all[, 1]
  }

  if (any(has_plus) && any(has_minus)) {
    plus_rle = rle(sort(plus_all)) |> rle_2_df() |> dplyr::filter(length >= 2)
    minus_rle = rle(sort(minus_all)) |> rle_2_df() |> dplyr::filter(length >= 2)
    match_locs = dplyr::inner_join(plus_rle, minus_rle, by = "value")
  } else {
    match_locs = tibble::tibble(
      value = vector("integer"),
      length = vector("integer")
    )
  }

  # this is the most pain in the butt problem, there are + and - that match
  # up with each other, and you only want to replace them with something when
  # they actually match up.
  if (nrow(match_locs) > 0) {
    seq_index = seq_len(length(raw_ids))

    new_ids = purrr::map_chr(seq_index, \(in_index) {
      replace_match_pm(
        raw_ids[in_index],
        plus_locs[[in_index]],
        minus_locs[[in_index]],
        match_locs$value
      )
    })
  } else if (any(has_plus)) {
    new_ids = gsub("\\+$", "_p_", raw_ids)
  } else if (any(has_minus)) {
    new_ids = gsub("\\-$", "_m_", raw_ids)
  } else {
    new_ids = raw_ids
  }
  n_raw = length(unique(raw_ids))
  n_new = length(unique(new_ids))
  if (n_raw != n_new) {
    stop("Length of new ids doesn't match raw ids! Stopping.")
  }
  return(new_ids)
}

replace_match_pm = function(in_string, in_plus, in_minus, loc_values) {
  # in_string = raw_ids[1]; in_plus = plus_locs[[1]]; in_minus = minus_locs[[1]]; loc_values = c(6, 8, 10)
  # in_string = raw_ids[12]; in_plus = plus_locs[[12]]; in_minus = minus_locs[[12]]; loc_values = c(6, 8, 10)
  in_plus = in_plus[, 1]
  in_plus = in_plus[in_plus %in% loc_values]
  in_minus = in_minus[, 1]
  in_minus = in_minus[in_minus %in% loc_values]

  if (length(in_plus) > 0) {
    plus_df = data.frame(loc = in_plus, replace = "_p_")
  } else {
    plus_df = data.frame(loc = vector("double"), replace = vector("character"))
  }
  if (length(in_minus) > 0) {
    minus_df = data.frame(loc = in_minus, replace = "_m_")
  } else {
    minus_df = data.frame(loc = vector("double"), replace = vector("character"))
  }
  replace_df = rbind(plus_df, minus_df)
  replace_df = replace_df |>
    dplyr::arrange(loc)

  out_string = in_string
  for (irow in seq_len(nrow(replace_df))) {
    n_out = nchar(out_string)
    tmp1 = substring(out_string, 1, replace_df$loc[irow] - 1)
    tmp2 = substring(out_string, replace_df$loc[irow] + 1, n_out)
    out_string = paste0(tmp1, replace_df$replace[irow], tmp2)
    if (irow < nrow(replace_df)) {
      replace_df$loc = replace_df$loc + 2
    }
  }
  out_string
}

create_new_sample_ids = function(in_sample_ids) {
  # in_sample_ids = factor_df$sample_id

  mapping_df = tibble::tibble(original = in_sample_ids)

  unique_df = unique(mapping_df)
  unique_df = unique_df |>
    dplyr::mutate(
      fix1 = fix_sample_ids(original),
      new_id = janitor::make_clean_names(fix1)
    )
  if (length(unique(unique_df$new_id)) != nrow(unique_df)) {
    stop("Something went wrong trying to fix the sample IDs!")
  }
  mapping_df = dplyr::left_join(
    mapping_df,
    unique_df[, c("original", "new_id")],
    by = "original",
    relationship = "many-to-one"
  )
  return(mapping_df$new_id)
}
