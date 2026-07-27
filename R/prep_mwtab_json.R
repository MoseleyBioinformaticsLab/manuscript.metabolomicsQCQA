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


process_mwtab = function(mwtab_datasets, ancillary_datasets, save_path) {
  processed_results = purrr::map2(
    mwtab_datasets$file,
    mwtab_datasets$id,
    \(x, y) {
      parse_json_save(
        x,
        y,
        ancillary_datasets = ancillary_datasets,
        save_path = save_path
      )
    },
    .progress = TRUE
  )
  processed_results
}


bind_group_processed = function(
  processed,
  size_group = 200,
  save_file = "resources/group_sets.rds"
) {
  # processed = tar_read(processed); n_group = 30; save_file = "resources/group_sets.rds"

  if (fs::file_exists(save_file)) {
    fs::file_delete(save_file)
  }

  analysis_ids = extract_mwtab_ids(names(processed))
  names(processed) = analysis_ids$id
  processed_df = purrr::imap(processed, \(in_p, id) {
    tibble::tibble(id = id, file = in_p$file, md5 = in_p$md5)
  }) |>
    purrr::list_rbind()
  n_group = ceiling(nrow(processed_df) / size_group)
  processed_group = processed_df |>
    dplyr::mutate(
      group = paste0(
        "g",
        stringr::str_pad(
          sample(seq_len(n_group), dplyr::n(), replace = TRUE),
          width = 2,
          side = "left",
          pad = "0"
        )
      )
    )

  split_sets = split(processed_group, processed_group$group)
  nested_sets = purrr::map(split_sets, \(in_set) {
    tibble::tibble(group = in_set$group[1], set = list(in_set))
  }) |>
    purrr::list_rbind()
  saveRDS(nested_sets, save_file)
  nested_sets
}

parse_json_save = function(mwtab_file, id, ancillary_datasets, save_path) {
  # message(id)
  #mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000032%2Fmwtab%2Ftxt.txt.json"  # MS file
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

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN004699%2Fmwtab%2Ftxt.txt.json"; id = "AN004699"; tar_source("R")

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003604%2Fmwtab%2Ftxt.txt.json"

  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN004188%2Fmwtab%2Ftxt.txt.json"
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000112%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary/"
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003915%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN002713%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN002792%2Fmwtab%2Ftxt.txt.json";ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003159%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003350%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN005090%2Fmwtab%2Ftxt.txt.json" ; ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003606%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary"; tar_source("R")
  # mwtab_file = "mwtab/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN005021%2Fmwtab%2Ftxt.txt.json"; ancillary_path = "mwtab/ancillary"; tar_source("R")
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
    parsed_data = check_and_parse_ancillary_json(
      parsed_data,
      ancillary_datasets
    )
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
  ancillary_datasets
) {
  # parsed_data = tar_read(processed_AN000111)
  # ancillary_path = "/big_data/data/rmflight_icikt_poster_iecm/data/ancillary/"
  if (!is.null(parsed_data$MS)) {
    parsed_data$MEASUREMENT_TYPE = "MS"
  }
  if (!is.null(parsed_data$NM)) {
    parsed_data$MEASUREMENT_TYPE = "NMR"
  }

  ancillary_just_file = ancillary_datasets |>
    dplyr::filter(id %in% parsed_data$ID) |>
    dplyr::pull(file)

  # can't find a match, just return early
  if (length(ancillary_just_file) == 0) {
    return(parsed_data)
  }

  ancillary_data = readr::read_delim(
    ancillary_just_file,
    delim = "\t",
    show_col_types = FALSE,
    name_repair = create_new_sample_ids
  )

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
    "metabolite",
    "mass_rt"
  )

  common_not_met = c(
    "mz_xml",
    "mzxml",
    "xml",
    "mz_ml",
    "mzml",
    "mzdata",
    "mz_data"
  )

  is_metabolite_name = purrr::map(common_met_fields, \(x) {
    grep(x, names(ancillary_data), value = TRUE)
  }) |>
    purrr::list_c() |>
    unique()

  is_not_metabolite = purrr::map(common_not_met, \(x) {
    grep(x, is_metabolite_name, value = TRUE)
  }) |>
    purrr::list_c() |>
    unique()

  is_metabolite_name = base::setdiff(is_metabolite_name, is_not_metabolite)
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


parse_nmr_data_json = function(nmr_data) {
  # nmr_data = block_data[["NMR_BINNED_DATA"]]
  # nmr_data = mwtab_list$NMR_BINNED_DATA

  nmr_df_1 = try(tibble::as_tibble(
    nmr_data$Data[[1]],
    .name_repair = create_new_sample_ids
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
    if (!is.null(metabolite_df)) {
      metabolite_df$feature_id = janitor::make_clean_names(
        metabolite_df[[1]],
        unique_sep = "..."
      )
    } else {
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

  tmp_tibble = tmp_tibble |>
    custom_convert_numeric()
  names(tmp_tibble) = create_new_sample_ids(names(tmp_tibble))
  tmp_tibble
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
  # ms_data = mwtab_list$MS_METABOLITE_DATA

  ms_df_1 = try(
    tibble::as_tibble(ms_data$Data[[1]], .name_repair = create_new_sample_ids)
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
    if (!is.null(metabolite_df)) {
      metabolite_df$feature_id = janitor::make_clean_names(
        metabolite_df[[1]],
        unique_sep = "..."
      )
    } else {
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
      sample_id = in_sample[["Sample ID"]]
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
  # need code here to handle using fix_sample_ids first, and then make_clean_names
  if (nrow(factor_df) == 0) {
    return(NULL)
  }
  null_subject = factor_df[["subject"]] == "-"
  if (any(factor_df$factor %in% c("sample_id", "SAMPLE_ID"))) {
    is_sample_id = factor_df$factor %in% c("sample_id", "SAMPLE_ID")
    factor_df$factor[is_sample_id] = "sample_id_1"
  }

  factor_df = factor_df |>
    dplyr::mutate(sample_id = create_new_sample_ids(sample_id))

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
