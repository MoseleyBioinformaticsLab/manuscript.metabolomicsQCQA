extract_mwtab_ids = function(mwtab_files) {
  split_files = stringr::str_split(fs::path_file(mwtab_files), "%2F")
  mwtab_df = tibble::tibble(
    id = purrr::map_chr(split_files, \(in_file) {
      in_file[7]
    }),
    file = mwtab_files
  )
  mwtab_df
}

check_mwtab = function(processed) {
  check_results = furrr::future_map(
    processed,
    run_mwtab_checks_json,
    min_n = 3,
    min_metabolites = 50,
    max_min_value = 20,
    min_ssf = 1,
    use_ssf_only = "yes",
    remove_pooled = "no"
  )
  check_results
}

process_mwtab = function(mwtab_datasets) {
  processed_results = furrr::future_map2(
    mwtab_datasets$file,
    mwtab_datasets$id,
    \(x, y) {
      parse_json_save(x, y, ancillary_path, save_path)
    }
  )
  processed_results
}

parse_json_save = function(mwtab_file, id, ancillary_path, save_path) {
  #mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000023%2Fmwtab%2Ftxt.txt"  # MS file
  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000693%2Fmwtab%2Ftxt.txt"  # NMR file

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000001%2Fmwtab%2Ftxt.txt"

  #mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000004%2Fmwtab%2Ftxt.txt"

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000038%2Fmwtab%2Ftxt.txt"

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000555%2Fmwtab%2Ftxt.txt"

  # mwtab_file = tar_read(dataset_AN000172)
  # mwtab_file = tar_read(mwtab_AN003177)
  # mwtab_file = tar_read(mwtab_AN002445)
  # mwtab_file = tar_read(mwtab_AN003894)
  # mwtab_file = tar_read(mwtab_AN001238)
  # mwtab_file = tar_read(dataset_AN000023)
  # mwtab_file = tar_read(dataset_AN000555)
  # mwtab_file = tar_read(dataset_AN001154)
  # mwtab_file = tar_read(dataset_AN007151)
  # id = "AN007151"
  # mwtab_file = tar_read(dataset_AN004150)
  # id = "AN004150"
  # mwtab_file = tar_read(dataset_AN003938)
  # id = "AN003938"
  # mwtab_file = tar_read(dataset_AN000428)
  # id = "AN000428"
  # mwtab_file = tar_read(dataset_AN000111)
  # id = "AN000111"
  # mwtab_file = tar_read(dataset_AN000352) # check on this one to see what is going on
  # id = "AN000352"
  # mwtab_file = tar_read(dataset_AN000593)
  # id = "AN000593"
  # mwtab_file = tar_read(dataset_AN000409)
  # id = "AN000409"
  # mwtab_file = tar_read(dataset_AN005567)
  # id = "AN005567"
  # mwtab_file = tar_read(dataset_AN003914)
  # id = "AN003914"
  # mwtab_file = tar_read(dataset_AN007063)
  # id = "AN007063"
  # mwtab_file = tar_read(dataset_AN002549)
  # id = "AN002549"
  # mwtab_file = tar_read(dataset_AN000033)
  # id = "AN000033"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000152%2Fmwtab%2Ftxt.txt.json"
  # id = "AN000152"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000023%2Fmwtab%2Ftxt.txt.json" # MS file
  # id = "AN000023"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000111%2Fmwtab%2Ftxt.txt.json" # MS file
  # id = "AN000111"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000255%2Fmwtab%2Ftxt.txt.json"
  # id = "AN000255"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000344%2Fmwtab%2Ftxt.txt.json"
  # id = "AN000344"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003150%2Fmwtab%2Ftxt.txt.json"
  # id = "AN003150"

  # ancillary_path = "mwtab/ancillary/"
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN005285%2Fmwtab%2Ftxt.txt.json"
  # id = "AN005285"

  mwtab_list = jsonlite::fromJSON(mwtab_file, simplifyVector = FALSE)

  parsed_data = purrr::imap(mwtab_list, \(block, block_id) {
    switch(
      block_id,
      `METABOLOMICS WORKBENCH` = parse_list(block),
      PROJECT = parse_list(block),
      STUDY = parse_list(block),
      SUBJECT = parse_list(block),
      SUBJECT_SAMPLE_FACTORS = parse_factors_json(block),
      COLLECTION = parse_list(block),
      TREATMENT = parse_list(block),
      SAMPLEPREP = parse_list(block),
      CHROMATOGRAPHY = parse_list(block),
      ANALYSIS = parse_list(block),
      MS = parse_list(block),
      MS_METABOLITE_DATA = suppressMessages(parse_ms_data_json(block)),
      NM = parse_list(block),
      NMR_BINNED_DATA = suppressMessages(parse_nmr_data_json(block)),
      NMR_METABOLITE_DATA = suppressMessages(parse_nmr_data_json(block)),
      block
    )
  })

  parsed_data$ID = id
  # if (is.null(parsed_data$SUBJECT_SAMPLE_FACTORS)) {
  #   # browser()
  #   return(parsed_data)
  # }

  if (!is.null(parsed_data$MS_METABOLITE_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$MS_METABOLITE_DATA$DATA
    parsed_data$METABOLITES = parsed_data$MS_METABOLITE_DATA$METABOLITES
    parsed_data$MEASUREMENT_TYPE = "MS"
    parsed_data$DATASOURCE = "mwtab"
  } else if (!is.null(parsed_data$NMR_BINNED_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$NMR_BINNED_DATA$DATA
    parsed_data$METABOLITES = parsed_data$NMR_BINNED_DATA$METABOLITES
    parsed_data$MEASUREMENT_TYPE = "NMR"
    parsed_data$DATASOURCE = "mwtab"
  } else if (!is.null(parsed_data$NMR_METABOLITE_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$NMR_METABOLITE_DATA$DATA
    parsed_data$METABOLITES = parsed_data$NMR_METABOLITE_DATA$METABOLITES
    parsed_data$MEASUREMENT_TYPE = "NMR"
    parsed_data$DATASOURCE = "mwtab"
  } else {
    parsed_data = check_and_parse_ancillary_json(parsed_data, ancillary_path)
  }

  parsed_digest = digest::digest(parsed_data, algo = "md5")
  out_file = fs::path(save_path, paste0(id, ".rds"))
  saveRDS(parsed_data, file = out_file)
  return(list(file = out_file, md5 = parsed_digest))
}

parse_list = function(in_list) {
  tibble::tibble(field = names(in_list), value = unlist(in_list))
}

create_results_txt = function(mwtab_df) {
  study_id = mwtab_df |>
    dplyr::filter(field %in% "STUDY_ID") |>
    dplyr::pull(value)

  analysis_id = mwtab_df |>
    dplyr::filter(field %in% "ANALYSIS_ID") |>
    dplyr::pull(value)

  if ((length(study_id) == 0) || (length(analysis_id) == 0)) {
    return(character(0))
  }
  results_txt = paste0(study_id, "_", analysis_id, "_Results.txt")
  return(results_txt)
}

convert_lgl_to_dbl = function(data_df) {
  col_types = readr::spec(data_df)

  is_lgl = purrr::map_lgl(col_types$cols, \(x) {
    any(grepl("logical", class(x)))
  })
}


check_and_parse_ancillary_json = function(
  parsed_data,
  ancillary_path = "/big_data/data/MetabolomicsWorkbench/ancillary/"
) {
  # parsed_data = tar_read(processed_AN000111)
  # ancillary_path = "/big_data/data/rmflight_icikt_poster_iecm/data/ancillary/"
  if (!is.null(parsed_data$MS)) {
    parsed_data$MEASUREMENT_TYPE = "MS"
  }
  if (!is.null(parsed_data$NM)) {
    parsed_data$MEASUREMENT_TYPE = "NMR"
  }
  results_file = character(0)

  results_file = create_results_txt(parsed_data[["METABOLOMICS WORKBENCH"]])

  if (length(results_file) == 0) {
    return(parsed_data)
  }

  ancillary_files = fs::dir_ls(ancillary_path, recurse = TRUE)
  ancillary_just_file = fs::path_file(ancillary_files)
  names(ancillary_files) = ancillary_just_file

  match_file = ancillary_files[results_file]

  # can't find a match, just return early
  if (is.na(match_file)) {
    return(parsed_data)
  }

  ancillary_data = readr::read_delim(
    match_file,
    delim = "\t",
    show_col_types = FALSE,
    name_repair = "universal_quiet"
  ) |>
    janitor::clean_names(unique_sep = "...")

  check_character = function(data_df) {
    col_types = readr::spec(data_df)
    is_character = purrr::map_lgl(col_types$cols, \(x) {
      any(grepl("character", class(x)))
    })
    names(is_character) = colnames(data_df)
    is_character
  }

  check_character_2 = function(data_df) {
    likely_character = purrr::map_lgl(data_df, \(x) {
      rlang::is_character(x)
    })
    names(likely_character) = colnames(data_df)
    likely_character
  }

  convert_all_numeric = function(data_df) {
    for (icol in seq_len(ncol(data_df))) {
      data_df[[icol]] = suppressWarnings(as.numeric(data_df[[icol]]))
    }
    data_df
  }

  check_row_to_numeric = function(data_df, in_row = 1) {
    num_1 = suppressWarnings(as.numeric(data_df[in_row, ]))
    !is.na(num_1)
  }
  check_character_retype = function(ancillary_data) {
    character_cols = check_character_2(ancillary_data)

    if (sum(character_cols) == ncol(ancillary_data)) {
      ancillary_tmp = suppressMessages(readr::type_convert(ancillary_data[
        -1,
      ]))
      character_cols = check_character_2(ancillary_tmp)
      if (sum(character_cols) == ncol(ancillary_tmp)) {
        ancillary_tmp = suppressMessages(readr::type_convert(ancillary_data[
          c(-1, -2),
        ]))
        character_cols = check_character_2(ancillary_tmp)
        if (sum(character_cols) == ncol(ancillary_tmp)) {
          ancillary_num = convert_all_numeric(ancillary_tmp)
          return(ancillary_num)
        } else {
          return(ancillary_tmp)
        }
      } else {
        return(ancillary_tmp)
      }
    } else {
      return(ancillary_data)
    }
  }

  common_met_fields = c(
    "quant_m_z",
    "m_z",
    "mz",
    "other_id",
    "compound",
    "ret_index",
    "spectra",
    "kegg",
    "pub_chem",
    "reactome",
    "chebi",
    "che_bi",
    "metabolite"
  )

  is_metabolite_name = purrr::map(common_met_fields, \(x) {
    grep(x, names(ancillary_data), value = TRUE)
  }) |>
    purrr::list_c() |>
    unique()

  if (length(is_metabolite_name) > 0) {
    metabolite_data = ancillary_data[, is_metabolite_name]
    sample_ids = base::setdiff(names(ancillary_data), is_metabolite_name)
    ancillary_checked = ancillary_data[, sample_ids]

    character_cols = which(check_character_2(ancillary_checked))
    for (icol in character_cols) {
      ancillary_checked[[icol]] = as.numeric(ancillary_checked[[icol]])
    }

    feature_id = janitor::make_clean_names(
      metabolite_data[[1]],
      unique_sep = "..."
    )
    ancillary_checked$feature_id = feature_id
    metabolite_data$feature_id = feature_id
  } else {
    ancillary_checked = check_character_retype(ancillary_data)
    if (is.null(ancillary_checked)) {
      return(parsed_data)
    }
    likely_character = check_character_2(ancillary_checked)
    if (sum(likely_character) == 0) {
      metabolite_data = tibble::tibble(
        feature_id = janitor::make_clean_names(
          seq_len(nrow(ancillary_checked)),
          unique_sep = "..."
        )
      )
      ancillary_checked$feature_id = metabolite_data$feature_id
    } else {
      metabolite_data = ancillary_checked[, likely_character]
      metabolite_data$feature_id = janitor::make_clean_names(
        metabolite_data[[1]],
        unique_sep = "..."
      )
      ancillary_checked$feature_id = metabolite_data$feature_id
    }
  }

  parsed_data$MEASUREMENTS = ancillary_checked
  parsed_data$METABOLITES = metabolite_data
  parsed_data$DATASOURCE = "ancillary"
  return(parsed_data)
}


get_check = function(checked_mwtab) {
  tibble::as_tibble(checked_mwtab$CHECK)
}


run_mwtab_checks_json = function(
  processed_rds,
  min_n = 3,
  min_metabolites = 50,
  max_min_value = 20,
  min_ssf = 1,
  use_ssf_only = "yes",
  remove_pooled = "no"
) {
  # processed_mwtab = tar_read(processed_AN002428)
  # processed_mwtab = tar_read(processed_AN003177)
  # processed_mwtab = tar_read(processed_AN003224)
  # processed_mwtab = tar_read(processed_AN003894)
  # processed_mwtab = tar_read(processed_AN000176)
  # processed_mwtab = tar_read(processed_AN001959)
  # processed_mwtab = tar_read(processed_AN001936)
  # processed_mwtab = tar_read(processed_AN000409)
  # processed_mwtab = tar_read(processed_AN003719)
  # processed_mwtab = tar_read(processed_AN004603)
  # processed_mwtab = tar_read(processed_AN000033)
  # processed_mwtab = tar_read(processed_AN005159)
  # processed_mwtab = tar_read(processed_AN003251)
  # processed_mwtab = tar_read(processed)[[28]]
  # processed_mwtab = tar_read(processed)[[60]]
  # processed_rds = tar_read(processed)[[575]]
  # min_n = 3
  # min_metabolites = 50
  # max_min_value = 20
  # min_ssf = 1
  # use_ssf_only = "yes"
  # remove_pools = "no"

  processed_mwtab = readRDS(processed_rds$file)
  processed_md5 = digest::digest(processed_mwtab, algo = "md5")

  check_res = list(
    md5_match = FALSE,
    features = "",
    n_features = NA,
    ssf = "",
    range = "",
    datasource = "",
    id = processed_mwtab$ID,
    file = processed_rds$file,
    md5 = processed_rds$md5
  )
  if (is.null(processed_mwtab$ID)) {
    check_res$id = ""
  }

  if (processed_rds$md5 != processed_md5) {
    return(check_res)
  } else {
    check_res$md5_match = TRUE
  }

  if (!is.null(processed_mwtab$DATASOURCE)) {
    check_res$datasource = processed_mwtab$DATASOURCE
  }

  tmp_measures = processed_mwtab$MEASUREMENTS
  if (is.null(tmp_measures)) {
    check_res$features = "NA"
    return(check_res)
  }
  match_features = base::intersect(
    tmp_measures$feature_id,
    processed_mwtab$METABOLITES$feature_id
  )
  tmp_measures = tmp_measures |>
    dplyr::filter(feature_id %in% match_features)
  n_features = nrow(tmp_measures)
  if (is.na(n_features) || is.null(n_features)) {
    check_res$features = "NA"
    return(check_res)
  }

  check_res$n_features = n_features
  if (n_features < min_metabolites) {
    check_res$features = "few"
    return(check_res)
  } else {
    check_res$features = "good"
  }

  if (is.null(processed_mwtab$SUBJECT_SAMPLE_FACTORS)) {
    check_res$ssf = "no ssf data"
    return(check_res)
  }

  ssf_data = processed_mwtab$SUBJECT_SAMPLE_FACTORS

  unique_samples = unique(ssf_data$sample_id)
  match_samples = base::intersect(
    unique_samples,
    names(processed_mwtab$MEASUREMENTS)
  )
  # match_samples = purrr::map(unique_samples, \(x) {
  #   tmp = grep(x, colnames(processed_mwtab$MEASUREMENTS), value = TRUE)
  #   if (length(tmp) == 1) {
  #     return(tibble::tibble(ssf = x, measurements = tmp))
  #   } else {
  #     return(NULL)
  #   }
  # }) |>
  #   purrr::list_rbind()

  if ((is.null(match_samples)) || (length(match_samples) == 0)) {
    check_res$ssf = "no common samples"
    return(check_res)
  }
  ssf_data2 = ssf_data |>
    dplyr::filter(sample_id %in% match_samples)

  if (remove_pooled %in% "yes") {
    is_pooled_qc = ssf_data2 |>
      dplyr::filter(
        grepl(
          ".*pool.*|.*qc.*|.*blank.*",
          sample_id,
          ignore.case = TRUE
        ) |
          grepl(".*pool.*|.*qc.*|.*blank.*", item, ignore.case = TRUE)
      ) |>
      dplyr::pull(sample_id)

    ssf_data2 = ssf_data2 |>
      dplyr::filter(!(sample_id %in% is_pooled_qc))
  }

  if (use_ssf_only %in% "yes") {
    ssf_data_use = ssf_data2 |>
      dplyr::filter(type %in% "subject_sample_factor")
  } else {
    ssf_data_use = ssf_data2
  }

  if (nrow(ssf_data_use) == 0) {
    check_res$ssf = "no samples left"
    return(check_res)
  }

  ssf_wide = widen_ssf_json(ssf_data_use)

  ssf_wide_grouped = create_factor_groups_json(ssf_wide)

  ssf_check_result = check_ssf_json(ssf_wide_grouped, min_ssf, min_n)

  check_res$ssf = ssf_check_result

  match_samples2 = match_samples[match_samples %in% ssf_wide_grouped$sample_id]

  measurements = processed_mwtab$MEASUREMENTS |>
    dplyr::filter(feature_id %in% match_features) |>
    dplyr::select(tidyselect::all_of(match_samples2)) |>
    as.matrix()
  factors = ssf_wide_grouped$factors

  # we check the number of features kept after checking for non-missing
  # in at least a single sample
  range_check_result = check_dynamic_range_single_sample(
    measurements,
    value_cutoff = max_min_value
  )
  check_res$range = range_check_result$range
  check_res$n_features = range_check_result$n
  if (range_check_result$n < min_metabolites) {
    check_res$features = "few"
  } else {
    check_res$features = "good"
  }

  check_res
}

check_dynamic_range_single_sample = function(
  in_values,
  value_cutoff = 20,
  min_range = 10
) {
  n_val = nrow(in_values) * ncol(in_values)
  n_na = sum(is.na(in_values))
  frac_na = n_na / n_val
  if (frac_na > 0.95) {
    return(list(range = "too many na", n = NA))
  }

  value_range = range(in_values, na.rm = TRUE)
  if (max(value_range) < value_cutoff) {
    in_values = 2^in_values
  }

  keep_values = try(visualizationQualityControl::keep_non_missing_percentage(
    in_values,
    keep_num = 1,
    missing_value = c(0, NA)
  ))
  if (inherits(keep_values, "try-error")) {
    return(list(range = "keep-error", n = 0))
  }
  in_values = in_values[keep_values, ]
  value_range = range(in_values, na.rm = TRUE)

  dynamic_range = try(value_range[2] / value_range[1])
  if (inherits(dynamic_range, "try-error")) {
    return(list(range = "error", n = 0))
  }
  if (dynamic_range < min_range) {
    return(list(range = "low", n = nrow(in_values)))
  } else if (dynamic_range >= min_range) {
    return(list(range = "good", n = nrow(in_values)))
  } else {
    return(list(range = "unknown", n = nrow(in_values)))
  }
}

widen_ssf_json = function(ssf_data2) {
  ssf_data2 |>
    tidyr::pivot_wider(
      id_cols = c(subject, sample_id),
      names_from = factor,
      values_from = item
    ) |>
    janitor::clean_names(unique_sep = "...")
}

create_factor_groups_json = function(ssf_wide) {
  other_cols = names(ssf_wide)[
    !(names(ssf_wide) %in% c("sample_id", "subject"))
  ]

  concat_factors = ssf_wide |>
    tidyr::unite("factors", tidyselect::all_of(other_cols), sep = "|")

  if (!is.null(ssf_wide[["factors"]])) {
    ssf_nofact = ssf_wide |> dplyr::select(-factors)
  } else {
    ssf_nofact = ssf_wide
  }
  concat_output = dplyr::left_join(
    ssf_nofact,
    concat_factors,
    by = c("subject", "sample_id")
  )

  concat_output
}

check_ssf_json = function(ssf_data2, min_ssf = 2, min_n = 5) {
  ssf_n = ssf_data2 |>
    dplyr::summarise(n = dplyr::n(), .by = factors)

  if (any(ssf_n$n >= min_n) && (nrow(ssf_n) >= min_ssf)) {
    return("good")
  } else {
    return("not enough samples")
  }
}

check_missingness_ranks_json = function(measurements, factors) {
  if ((sum(is.na(measurements)) == 0) && (sum(measurements == 0) == 0)) {
    return("NO MISSING VALUES")
  }

  if (ncol(measurements) < 2) {
    return("TOO FEW SAMPLES")
  }

  ranked_data = ICIKendallTau::rank_order_data(
    measurements,
    sample_classes = factors
  )

  ranked_null = purrr::map_lgl(ranked_data, is.null)
  if (any(ranked_null)) {
    return("ALL NA IN A GROUP")
  }

  if (length(unique(factors)) == 1) {
    ranked_use = vector("list", 1)
    ranked_use[[1]] = ranked_data
    names(ranked_use) = unique(factors)
    ranked_data = ranked_use
  }
  ranked_cor = purrr::map(ranked_data, \(in_data) {
    min_rank = in_data$n_na_rank |>
      dplyr::summarise(min_rank = min(median_rank), .by = n_na)
    full_cor = cor(
      in_data$n_na_rank$n_na,
      in_data$n_na_rank$median_rank,
      method = "kendall"
    )
    min_cor = cor(
      min_rank$n_na,
      min_rank$min_rank,
      method = "kendall"
    )
    tibble::tibble(median = full_cor, min = min_cor)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(sign_diff = sign(median) * sign(min))

  if (any(is.na(ranked_cor$sign_diff))) {
    return("NA CORRELATION")
  } else if (any(ranked_cor$sign_diff == -1)) {
    return("SIGN DIFFERENCE")
  } else if (any(sign(ranked_cor$median) == 1)) {
    return("POSITIVE MEDIAN")
  } else if (any(sign(ranked_cor$min) == 1)) {
    return("POSITIVE MINIMUM")
  } else {
    return("GOOD")
  }
}

count_factors_replicates_long = function(subject_sample_factors) {
  # subject_sample_factors = parsed_data$SUBJECT_SAMPLE_FACTORS

  other_cols = names(subject_sample_factors)[
    !(names(subject_sample_factors) %in% "sample_id")
  ]

  concat_factors = subject_sample_factors |>
    tidyr::unite("factors", tidyselect::all_of(other_cols), sep = "|")
  subject_sample_factors = dplyr::left_join(
    subject_sample_factors,
    concat_factors,
    by = "sample_id"
  )

  subject_sample_factors$sample_id = janitor::make_clean_names(
    subject_sample_factors$sample_id,
    unique_sep = "..."
  )

  n_samples = subject_sample_factors |>
    dplyr::summarise(n = dplyr::n(), .by = factors)
  list(sample_factors = subject_sample_factors, n = n_samples)
}


parse_nmr_data_json = function(nmr_data) {
  # nmr_data = block_data[["NMR_BINNED_DATA"]]
  # nmr_data = mwtab_list$NMR_BINNED_DATA

  nmr_df_1 = try(tibble::as_tibble(
    nmr_data$Data[[1]],
    .name_repair = "minimal"
  ))
  if (inherits(nmr_df_1, "try-error")) {
    return(NULL)
  }

  nmr_data_df = tryCatch(
    {
      ms_create_tibble(nmr_data$Data)
    },
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(nmr_data_df)) {
    return(NULL)
  }

  nmr_data_df$feature_id = janitor::make_clean_names(
    nmr_data_df[[1]],
    unique_sep = "..."
  )

  if (!is.null(nmr_data$Metabolites)) {
    metabolite_df = tryCatch(
      {
        ms_create_tibble(nmr_data$Metabolites)
      },
      error = function(e) {
        return(NULL)
      }
    )
    metabolite_df$feature_id = janitor::make_clean_names(
      metabolite_df[[1]],
      unique_sep = "..."
    )
    if (is.null(metabolite_df)) {
      metabolite_df = tibble::tibble(feature_id = nmr_data_df$feature_id)
    }
  } else {
    metabolite_df = tibble::tibble(feature_id = nmr_data_df$feature_id)
  }

  list(DATA = nmr_data_df, METABOLITES = metabolite_df)
}

ms_create_tibble = function(ms_data_json) {
  tmp_list = purrr::map(
    ms_data_json,
    tibble::as_tibble,
    .name_repair = "minimal"
  )
  tmp_tibble = do.call(rbind, tmp_list)

  for (icol in seq_len(ncol(tmp_tibble))) {
    tmp_guess = readr::parse_guess(tmp_tibble[[icol]])
    if (inherits(tmp_guess, "logical")) {
      tmp_tibble[[icol]] = readr::parse_number(tmp_tibble[[icol]])
    } else {
      tmp_tibble[[icol]] = tmp_guess
    }
  }

  tmp_tibble |>
    custom_convert_numeric() |>
    janitor::clean_names(unique_sep = "...")
}

custom_convert_numeric = function(in_df) {
  numeric_df = suppressWarnings(purrr::map(in_df, as.numeric)) |>
    dplyr::bind_cols()
  is_character = purrr::map_lgl(numeric_df, \(x) {
    all(is.na(x))
  }) |>
    which()

  if (length(is_character) > 0) {
    numeric_df[, is_character] = in_df[, is_character]
  }
  numeric_df
}

parse_ms_data_json = function(ms_data) {
  # ms_data = block_data[["MS_METABOLITE_DATA"]]
  # message("ms data!")
  # ms_data = mwtab_list$MS_METABOLITE_DATA

  ms_df_1 = try(
    tibble::as_tibble(ms_data$Data[[1]], .name_repair = "minimal")
  )
  if (inherits(ms_df_1, "try-error")) {
    return(NULL)
  }
  ms_data_df = tryCatch(
    {
      ms_create_tibble(ms_data$Data)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(ms_data_df)) {
    return(NULL)
  }
  ms_data_df$feature_id = janitor::make_clean_names(
    ms_data_df[[1]],
    unique_sep = "..."
  )

  if (!is.null(ms_data$Metabolites)) {
    metabolite_df = tryCatch(
      {
        ms_create_tibble(ms_data$Metabolites)
      },
      error = function(e) {
        return(NULL)
      }
    )
    metabolite_df$feature_id = janitor::make_clean_names(
      metabolite_df[[1]],
      unique_sep = "..."
    )
    if (is.null(metabolite_df)) {
      metabolite_df = tibble::tibble(feature_id = ms_data_df$feature_id)
    }
  } else {
    metabolite_df = tibble::tibble(feature_id = ms_data_df$feature_id)
  }

  return(list(DATA = ms_data_df, METABOLITES = metabolite_df))
}


parse_factors_json = function(factor_data) {
  # message("factor data!")
  # factor_data = mwtab_list[["SUBJECT_SAMPLE_FACTORS"]]

  factor_df = purrr::map(factor_data, \(in_sample) {
    # in_sample = factor_data[[1]]
    sub_sample = tibble::tibble(
      subject = in_sample[["Subject ID"]],
      sample_id = janitor::make_clean_names(
        in_sample[["Sample ID"]],
        unique_sep = "..."
      )
    )
    factors = tibble::tibble(
      factor = names(in_sample$Factors),
      item = unlist(in_sample$Factors),
      type = "subject_sample_factor"
    )

    out_df = dplyr::bind_cols(sub_sample, factors)

    if (!is.null(in_sample$`Additional sample data`)) {
      additional = tibble::tibble(
        factor = names(in_sample$`Additional sample data`),
        item = unlist(in_sample$`Additional sample data`),
        type = "additional_sample_data"
      )
      add_df = dplyr::bind_cols(sub_sample, additional)
      out_df = dplyr::bind_rows(out_df, add_df)
    }
    out_df
  }) |>
    purrr::list_rbind()

  if (nrow(factor_df) == 0) {
    return(NULL)
  }
  null_subject = factor_df[["subject"]] == "-"
  if (any(factor_df$factor %in% "sample_id")) {
    return(NULL)
  }

  # these basically can't be empty, if they are, bad things happen.
  empty_factors = factor_df |>
    dplyr::filter(nchar(factor) == 0)
  unique_samples = factor_df$sample_id |> unique()

  if (nrow(empty_factors) > 0) {
    mod_rows = nrow(empty_factors) %% length(unique_samples)
    if (!all(unique_samples %in% empty_factors$sample_id) || !(mod_rows == 0)) {
      return(NULL)
    } else {
      factor_df = factor_df |>
        dplyr::filter(!(nchar(factor) == 0))
    }
  }

  if (any(nchar(factor_df$factor) == 0)) {
    return(NULL)
  }

  if (any(nchar(factor_df$sample_id) == 0)) {
    return(NULL)
  }

  return(factor_df)
}

split_and_return_factors = function(factor_string) {
  # factor_string = "AFTER_MEAL_TIME:30_MIN | gender:female | POST_SURGERY_TIME:12MO | SURGERY:GASTRIC_BAND"
  # factor_string = "Treatment:Media alone"

  pipe_split = stringr::str_split(factor_string, "\\|")[[1]] |> trimws()
  colon_split = stringr::str_split(pipe_split, "\\:")
  factor_df = tibble::tibble(
    factor = purrr::map_chr(colon_split, \(x) {
      x[1]
    }),
    item = purrr::map_chr(colon_split, \(x) {
      x[2]
    })
  )
  factor_df
}


parse_metabolites_block = function(metabolites_block) {
  # metabolites_block = block_data$METABOLITES

  end_loc = which(grepl("METABOLITES_END", metabolites_block)) - 1
  metabolites_range = seq(2, end_loc)
  metabolites_block = metabolites_block[metabolites_range]

  if (length(metabolites_block) < 2) {
    return(NULL)
  }
  has_id_row = grepl(
    "metabolite_name|metabolite_id|id|kegg",
    metabolites_block[1],
    ignore.case = TRUE
  )

  all_split = stringr::str_split(metabolites_block, "\t")
  if (has_id_row) {
    col_id = janitor::make_clean_names(all_split[[1]], unique_sep = "...")
    all_split = all_split[-1]
  } else {
    col_id = janitor::make_clean_names(
      seq(1, length(all_split[[1]])),
      unique_sep = "..."
    )
  }

  split_bind = do.call(rbind, all_split)

  if (ncol(split_bind) < length(col_id)) {
    col_id = col_id[seq_len(ncol(split_bind))]
  } else if (ncol(split_bind) > length(col_id)) {
    col_id = c(
      col_id,
      janitor::make_clean_names(seq_len(ncol(split_bind) - length(col_id)))
    )
  }
  colnames(split_bind) = col_id
  tibble::as_tibble(split_bind)
}


find_possible_mwtab = function() {
  mwtab_results = readRDS("data/processed/mwtab_results.rds")
  n_analyses = nrow(mwtab_results)
  n_na = sum(is.na(mwtab_results$n_features))
  n_less_200 = sum(
    !is.na(mwtab_results$n_features) & (mwtab_results$n_features < 200),
    na.rm = TRUE
  )
  mwtab_nonna = mwtab_results |>
    dplyr::filter(!is.na(n_features)) |>
    dplyr::filter(n_features >= 200)
  n_nonna = nrow(mwtab_nonna)
  keep_reps = purrr::map(mwtab_nonna$reps, \(in_reps) {
    # in_reps = mwtab_nonna$reps[[1]]
    keep_reps = in_reps |>
      dplyr::filter(n >= 3)
    keep_reps
  })
  has_multi = purrr::map_lgl(keep_reps, \(x) {
    (nrow(x) >= 2) && (sum(x$n >= 5) >= 1)
  })
  n_has_multi = sum(has_multi)

  mwtab_nonna$reps = keep_reps

  check_results = mwtab_nonna[has_multi, ]
  n_data = tibble::tibble(
    `N Analyses` = n_analyses,
    `N NA Features` = n_na,
    `N Less 200` = n_less_200,
    `N Few Factors` = sum(!has_multi),
    `N Multiple Factors` = n_has_multi
  )
  out_results = list(N = n_data, check = check_results)
  saveRDS(out_results, "data/processed/mwtab_check.rds")
}


check_ranking_mwtab = function() {
  out_results = readRDS("data/processed/mwtab_check.rds")

  check_results = out_results$check |>
    dplyr::mutate(rds_file = fs::path("data", "processed", id, ext = "rds"))

  check_each = purrr::pmap_chr(
    check_results,
    load_and_check_mwtab,
    .progress = TRUE
  )

  check_results$value_check = check_each
  saveRDS(check_results, "data/processed/mwtab_rank.rds")
  return(invisible(NULL))
}

find_mwtab_keep = function() {
  check_results = readRDS("data/processed/mwtab_rank.rds")

  keep_results = check_results |>
    dplyr::filter(value_check %in% c("good", "sign difference"))

  saveRDS(keep_results, "data/processed/mwtab_keep.rds")
}

load_and_check_mwtab = function(rds_file, reps, ...) {
  # rds_file = check_results$rds_file[70]
  # reps = check_results$reps[[70]]

  # rds_file = check_results$rds_file[7]
  # reps = check_results$reps[[7]]

  keep_factors = reps
  mwtab_data = readRDS(rds_file)

  sample_factors = mwtab_data$SUBJECT_SAMPLE_FACTORS |>
    dplyr::filter(factors %in% keep_factors$factors)

  sample_subjects = mwtab_data$SUBJECT
  if (any(grepl("\\;", sample_subjects$value))) {
    return("multiple species")
  }

  if (
    length(
      base::setdiff(
        sample_factors$sample_id,
        colnames(mwtab_data$MEASUREMENTS)
      ) >
        0
    )
  ) {
    match_samples = base::intersect(
      sample_factors$sample_id,
      colnames(mwtab_data$MEAUREMENTS)
    )
    sample_factors = sample_factors |>
      dplyr::filter(sample_id %in% match_samples)

    n_rep = sample_factors |>
      dplyr::summarise(n = dplyr::n(), .by = factors)
    keep_it = (nrow(n_rep) >= 2) && (sum(n_rep$n >= 5) >= 1)

    if (!keep_it) {
      return("missing or not enough samples")
    }
  }

  if (
    inherits(mwtab_data$MEASUREMENTS[[sample_factors$sample_id[1]]], "list")
  ) {
    return("list data")
  }

  metabolite_values = mwtab_data$MEASUREMENTS[, sample_factors$sample_id] |>
    as.matrix()

  ranked_data = ICIKendallTau::rank_order_data(
    metabolite_values,
    sample_classes = sample_factors$factors
  )
  ranked_null = purrr::map_lgl(ranked_data, is.null)
  if (any(ranked_null)) {
    return("all NA in a group")
  }

  if (
    (sum(is.na(metabolite_values)) == 0) && (sum(metabolite_values == 0) == 0)
  ) {
    return("no missing")
  }
  rank_cor = purrr::map(ranked_data, \(in_data) {
    # in_data = ranked_data[[1]]
    min_rank = in_data$n_na_rank |>
      dplyr::summarise(min_rank = min(median_rank), .by = n_na)
    full_cor = cor(
      in_data$n_na_rank$n_na,
      in_data$n_na_rank$median_rank,
      method = "kendall"
    )
    min_cor = cor(
      min_rank$n_na,
      min_rank$min_rank,
      method = "kendall"
    )
    tibble::tibble(median = full_cor, min = min_cor)
  }) |>
    purrr::list_rbind()

  rank_cor = rank_cor |>
    dplyr::mutate(sign_diff = sign(median) * sign(min))

  if (any(is.na(rank_cor$sign_diff))) {
    return("na corr")
  } else if (any(rank_cor$sign_diff == -1)) {
    return("sign difference")
  } else if (any(sign(rank_cor$median) == 1)) {
    return("positive median")
  } else if (any(sign(rank_cor$min) == 1)) {
    return("positive minimum")
  } else {
    return("good")
  }
}


convert_mwtab_json_smd = function(json_checked) {
  # json_checked = tar_read(checked_AN000033)
  # json_checked = tar_read(checked_AN007099)

  check_res = json_checked$CHECK
  if (
    !(check_res$FEATURE_CHECK %in% "GOOD") ||
      !(check_res$SSF_CHECK %in% "GOOD SSF") ||
      !(check_res$RANK_CHECK %in% c("GOOD", "SIGN DIFFERENCE")) ||
      !(check_res$RANGE_CHECK %in% "GOOD")
  ) {
    return(NULL)
  }

  sample_metadata = json_checked$SUBJECT_METADATA

  if (json_checked$MEASUREMENT_TYPE %in% "MS") {
    json_checked$MEASUREMENT_INFO = json_checked$MS
  } else if (json_checked$MEASUREMENT_TYPE %in% "NMR") {
    json_checked$MEASUREMENT_INFO = json_checked$NMR
  }
  project_data = json_checked[c(
    "METABOLOMICS WORKBENCH",
    "PROJECT",
    "STUDY",
    "SUBJECT",
    "COLLECTION",
    "TREATMENT",
    "SAMPLEPREP",
    "CHROMATOGRAPHY",
    "ANALYSIS",
    "MEASUREMENT_TYPE",
    "MEASUREMENT_INFO",
    "CHECK"
  )]

  match_samples = base::intersect(
    sample_metadata$sample_id,
    colnames(json_checked$MEASUREMENTS)
  )

  sample_metadata = sample_metadata |>
    dplyr::filter(sample_id %in% match_samples)
  measurements = json_checked$MEASUREMENTS[, c(
    "feature_id",
    sample_metadata$sample_id
  )]

  feature_metadata = json_checked$METABOLITES
  feature_DF = DataFrame(feature_metadata)
  rownames(feature_DF) = feature_DF$feature_id

  sample_DF = DataFrame(sample_metadata)
  rownames(sample_DF) = sample_DF$sample_id
  measurement_array = measurements |>
    dplyr::select(-feature_id) |>
    as.matrix()

  rownames(measurement_array) = measurements$feature_id
  keep_metabolites = base::intersect(
    rownames(feature_DF),
    rownames(measurement_array)
  )
  feature_DF = feature_DF[keep_metabolites, , drop = FALSE]
  measurement_array = measurement_array[
    feature_DF$feature_id,
    sample_DF$sample_id
  ]

  if (typeof(measurement_array) %in% "character") {
    mode(measurement_array) = "numeric"
  }
  measurement_array[measurement_array <= 0] = NA

  if (nrow(measurement_array) == 0) {
    return(NULL)
  }

  normalized_array = median_normalize(measurement_array)

  # making sure next steps use metabolites that are present in some subset
  # of the samples, in this case 50% of at least one class
  present_metabolites = visualizationQualityControl::keep_non_missing_percentage(
    measurement_array,
    sample_DF$factors,
    keep_num = 0.5
  )

  if (sum(present_metabolites) <= 20) {
    return(NULL)
  }

  measurement_array = measurement_array[present_metabolites, , drop = FALSE]
  normalized_array = normalized_array[present_metabolites, , drop = FALSE]
  feature_DF = feature_DF[present_metabolites, , drop = FALSE]
  out_smd = SummarizedExperiment(
    assays = SimpleList(
      counts = measurement_array,
      normalized = normalized_array
    ),
    rowData = feature_DF,
    colData = sample_DF,
    metadata = project_data
  )

  out_smd
}


convert_list_to_df = function(measurements) {
  id = "AN000555"
  rds_file = fs::path("data", "processed", id, ext = "rds")
  mwtab_data = readRDS(rds_file)
  measurements = mwtab_data$MEASUREMENTS

  out_vals = purrr::map(measurements, \(in_meas) {
    unlist(in_meas)
  })
}


create_ancillary_file_list = function(repaired_json_dir) {
  repaired_json_dir = "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12/"
  json_files = fs::dir_ls(repaired_json_dir)
  tar_source("R")

  check_external_results = function(in_file) {
    # in_file = "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003030%2Fmwtab%2Ftxt.txt.json"
    # in_file = "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000001%2Fmwtab%2Ftxt.txt.json"
    json_data = jsonlite::fromJSON(in_file, simplifyVector = FALSE)
    results_file = paste0(
      json_data$`METABOLOMICS WORKBENCH`$STUDY_ID,
      "_",
      json_data$`METABOLOMICS WORKBENCH`$ANALYSIS_ID,
      "_Results.txt"
    )
    tmp_id = extract_mwtab_ids(in_file)
    tmp_id$results = results_file
    tmp_id
  }

  results_files = purrr::map(
    json_files,
    check_external_results,
    .progress = TRUE
  ) |>
    purrr::list_rbind()
  cat(
    results_files,
    sep = "\n",
    file = "/big_data/data/MetabolomicsWorkbench/ancillary/ancillary_files.txt"
  )
  return(invisible(NULL))
}

download_ancillary = function(
  ancillary_files = "/big_data/data/MetabolomicsWorkbench/ancillary/ancillary_files.txt",
  download_dir = "/big_data/data/MetabolomicsWorkbench/ancillary/"
) {
  # the list of files can be created by doing:
  # grep -R -h -o "ST.*[[:digit:]].*\\_AN.*[[:digit:]].*.txt" > ancillary_files.txt
  url_loc = "https://www.metabolomicsworkbench.org/studydownload/"
  exist_files = fs::dir_ls(download_dir) |> fs::path_file()
  all_files = readLines(ancillary_files)
  all_files = all_files[!(all_files %in% exist_files)]
  url_full = paste0(url_loc, all_files)

  dest_files = fs::path(download_dir, all_files)
  for (iloc in seq_len(length(all_files))) {
    try(download.file(url = url_full[iloc], destfile = dest_files[iloc]))
    Sys.sleep(2)
  }
  return(invisible(NULL))
}


test_json_parsing = function() {
  tar_source("R")
  json_files = fs::dir_ls(
    "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12"
  )
  ancillary_path = "/big_data/data/MetabolomicsWorkbench/ancillary/"
  mwtab_datasets = extract_mwtab_ids(json_files)

  set.seed(1234)
  test_datasets = mwtab_datasets |>
    dplyr::slice_sample(prop = 0.025) |>
    dplyr::transmute(mwtab_file = file, id = id) |>
    dplyr::mutate(ancillary_path = ancillary_path)

  json_parsed = purrr::pmap(test_datasets, parse_json, .progress = TRUE) |>
    flighttools::ft_notify_success_error()
  check_parsed = purrr::map(
    json_parsed,
    run_mwtab_checks_json,
    min_n = 5,
    min_metabolites = 100,
    min_ssf = 2,
    use_ssf_only = "yes",
    .progress = TRUE
  ) |>
    flighttools::ft_notify_success_error()
}

handle_pooled_samples = function() {
  # this would be code that is useful for "other" types of tasks, or
  # for processing for inclusion in an atlas and then excluding.
  # for this project, we want to get rid of them right now, because
  # they don't add any information, and can cloud our analysis.
  just_subject_sample = ssf_data2 |>
    dplyr::select(subject, sample_id) |>
    dplyr::distinct()
  pooled_qc_status = just_subject_sample |>
    dplyr::mutate(
      factor = "pooled_qc_status",
      item = dplyr::case_when(
        sample_id %in% is_pooled_qc ~ "pooled_or_qc_sample",
        TRUE ~ "experimental"
      ),
      type = "subject_sample_factor"
    )

  ssf_data2 = dplyr::bind_rows(ssf_data2, pooled_qc_status)

  other_cols = names(ssf_wide)[
    !(names(ssf_wide) %in% c("sample_id", "subject", "pooled_qc_status"))
  ]

  concat_factors = ssf_wide |>
    tidyr::unite("factors", tidyselect::all_of(other_cols), sep = "|")
  concat_output = dplyr::left_join(
    ssf_wide,
    concat_factors |> dplyr::select(-pooled_qc_status),
    by = c("subject", "sample_id")
  )

  ssf_wide_grouped_noqc = ssf_wide_grouped |>
    dplyr::filter(pooled_qc_status %in% "experimental")
}

generate_output_lists = function(
  check_data,
  no_data_file = "output_data/no_data.txt",
  no_common_file = "output_data/non_matching_samples.txt"
) {
  no_data = check_data |>
    dplyr::filter(FEATURE_CHECK %in% "NA FEATURES") |>
    dplyr::transmute(CHECK = FEATURE_CHECK, DATASET = ID)
  non_matching = check_data |>
    dplyr::filter(SSF_CHECK %in% "NO COMMON SAMPLES") |>
    dplyr::transmute(CHECK = SSF_CHECK, DATASET = ID)

  write.table(
    no_data,
    file = no_data_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  write.table(
    non_matching,
    file = no_common_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  return(invisible(NULL))
}

filter_good_experiments_splitup = function(check_result, save_file) {
  if (fs::file_exists(save_file)) {
    fs::file_delete(save_file)
  }

  good_result = check_result |>
    dplyr::filter(
      (md5_match & features %in% "good") &
        (ssf %in% "good") &
        (range %in% "good")
    ) |>
    dplyr::mutate(
      group = paste0(
        "g",
        stringr::str_pad(
          sample(seq(1, 20), dplyr::n(), replace = TRUE),
          width = 2,
          side = "left",
          pad = "0"
        )
      )
    )
  split_sets = split(good_result, good_result$group)
  nested_sets = purrr::map(split_sets, \(in_set) {
    tibble::tibble(group = in_set$group[1], set = list(in_set))
  }) |>
    purrr::list_rbind()
  saveRDS(nested_sets, save_file)
  nested_sets
}

bind_check_results = function(checked) {
  purrr::map(checked, \(x) {
    tibble::as_tibble(x)
  }) |>
    purrr::list_rbind()
}
