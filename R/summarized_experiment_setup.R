sets_readrds_convert_to_se = function(
  set_df,
  processed_group,
  mwtab_raw_files
) {
  if (inherits(set_df, "character")) {
    stop(
      "Run to processed_group first, and then the rest!  tar_make(processed_group)"
    )
  }
  all_files = set_df$file
  all_md5 = set_df$md5
  all_id = set_df$id
  all_se = furrr::future_map2(
    all_files,
    all_md5,
    readrds_convert_to_se,
    mwtab_raw_files
  )
  names(all_se) = all_id
  all_se
}

check_md5 = function(in_data, md5_sum) {
  new_md5 = digest::digest(in_data, algo = "md5")
  stopifnot(new_md5 == md5_sum)
  return(invisible(NULL))
}

# takes the FTP listing, and actually finds the files
# that likely have the raw data in them, to match against
# the study ID.
get_mwtab_ftp_files = function(
  file = "resources/mwbench_ftp_studies_files_2026-06-25.txt"
) {
  all_files = readr::read_lines(file = file)
  all_files = gsub("\\s+", " ", all_files)
  all_files = gsub("->.*", "", all_files)
  file_list = strsplit(all_files, " ")
  file_matrix = do.call(rbind, file_list)
  just_files = file_matrix[, 9]

  likely_raw = just_files[
    !(grepl("supp|_Results|docx$|xlsx$|txt$|csv$", just_files))
  ]
  likely_raw
}

check_convert_character_matrix = function(raw_matrix) {
  if (is.character(raw_matrix[1, 1])) {
    safe_matrix = apply(
      raw_matrix,
      2,
      \(x) {
        y = stringr::str_replace_all(x, "\\<|\\>", "-")
        z = stringr::str_replace_all(y, "\\.E", "E")
        z
      },
      simplify = TRUE
    )
    if (is.null(dim(safe_matrix))) {
      safe_matrix = matrix(
        safe_matrix,
        nrow = nrow(raw_matrix),
        ncol = ncol(raw_matrix),
        byrow = FALSE
      )
      rownames(safe_matrix) = rownames(raw_matrix)
      colnames(safe_matrix) = colnames(raw_matrix)
    }
    mode(safe_matrix) = "numeric"
  } else {
    safe_matrix = raw_matrix
  }
  return(safe_matrix)
}

check_all_na = function(raw_matrix) {
  na_vals = is.na(raw_matrix)
  all(na_vals)
}

readrds_convert_to_se = function(rds_file, md5_sum, mwtab_raw_files) {
  # setup_set = tar_read(setup_g01); tar_source("R")
  # using setup_g01
  # setup_data = setup_set[[4]]
  # setup_data = setup_set[[55]]
  # setup_data = setup_set[["AN000731"]]
  # setup_data = setup_set[["AN003244"]]

  # using setup_g02
  # setup_data = setup_set[["AN004188"]]
  # setup_data = setup_set[["AN001560"]]
  # setup_data = setup_set[["AN002381"]]
  # setup_data = setup_set[["AN002935"]]
  # setup_data = setup_set[["AN003604"]]
  # setup_data = setup_set[["AN003915"]]

  # id = "AN000072"; rds_file = "parsed_rds/AN000072.rds"; md5_sum = "781a19811ab6cf3ca6524224ab2e3507"; tar_source("R")
  # id = "AN003915"; rds_file = "parsed_rds/AN003915.rds"; tar_source("R")
  # id = "AN004851"; rds_file = "parsed_rds/AN004851.rds"; tar_source("R")
  # id = "AN002880"; rds_file = "parsed_rds/AN002880.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN000627.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN001180.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN002196.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN002298.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN001558.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN005958.rds"; tar_load(mwtab_raw_files); tar_source("R")
  # rds_file = "parsed_rds/AN002713.rds"; tar_load(mwtab_raw_files); tar_source("R")
  processed_data = readRDS(rds_file)
  check_md5(processed_data, md5_sum)

  ssf_data = processed_data$SUBJECT_SAMPLE_FACTORS
  abundance_data = processed_data$MEASUREMENTS
  metabolite_data = processed_data$METABOLITES

  study_id = processed_data[["METABOLOMICS WORKBENCH"]] |>
    dplyr::filter(field %in% "STUDY_ID") |>
    dplyr::pull(value)

  if (!is.na(nchar(study_id[1]))) {
    processed_data$RAW_DATA = grep(study_id, mwtab_raw_files, value = TRUE)
  } else {
    processed_data$RAW_DATA = NULL
  }

  if (is.null(abundance_data)) {
    se_data = create_empty_se(processed_data, why = "no feature data")
    return(se_data)
  }
  if (is.null(ssf_data)) {
    se_data = create_empty_se(processed_data, why = "no subject sample factors")
    return(se_data)
  }

  original_abundance_cols = colnames(abundance_data)
  original_ssf_samples = unique(ssf_data$sample_id)
  original_abundance_metabolites = abundance_data$feature_id
  original_metabolites = metabolite_data$feature_id

  abundance_data = intersect_samples(
    abundance_data,
    original_ssf_samples
  )
  ssf_data = ssf_data |>
    dplyr::filter(sample_id %in% colnames(abundance_data))

  match_metabolites = base::intersect(
    original_metabolites,
    abundance_data$feature_id
  )
  metabolite_data = metabolite_data |>
    dplyr::filter(feature_id %in% match_metabolites)
  raw_matrix = get_just_values(abundance_data)
  raw_matrix = raw_matrix[match_metabolites, , drop = FALSE]

  n_feature = nrow(raw_matrix)
  n_sample = ncol(raw_matrix)

  n_filtered_samples = tibble::tribble(
    ~source        , ~original                       , ~kept            , ~filtered                                                                              ,
    "ssf"          , length(original_ssf_samples)    , ncol(raw_matrix) , paste0(base::setdiff(original_ssf_samples, colnames(abundance_data)), collapse = "; ") ,
    "measurements" , length(original_abundance_cols) , ncol(raw_matrix) , paste0(base::setdiff(colnames(raw_matrix), original_ssf_samples), collapse = ";")
  )
  n_filtered_features = tibble::tribble(
    ~source       , ~original                    , ~kept                 ,
    "abundance"   , nrow(abundance_data)         , n_feature             ,
    "metabolites" , length(original_metabolites) , nrow(metabolite_data)
  )

  processed_data$FILTERED_FEATURES = n_filtered_features
  processed_data$FILTERED_SAMPLES = n_filtered_samples
  if (n_feature == 0) {
    if (n_filtered_features$kept[1] == 0) {
      why = "no features; didn't match metabolites"
    } else {
      why = "no features; unknown"
    }

    se_data = create_empty_se(processed_data, why = why)
    return(se_data)
  }

  if (n_sample == 0) {
    if (n_filtered_samples$kept[2] == 0) {
      why = "no samples; didn't match ssf"
    } else {
      why = "no samples; unknown"
    }
    se_data = create_empty_se(processed_data, why = why)
    return(se_data)
  }

  safe_matrix = check_convert_character_matrix(raw_matrix)
  safe_na = is.na(safe_matrix)

  if (all(safe_na)) {
    se_data = create_empty_se(processed_data, why = "no numeric abundances")
    return(se_data)
  }

  use_features = rownames(safe_matrix)
  keep_features = try(visualizationQualityControl::keep_non_missing_percentage(
    safe_matrix,
    keep_num = 1
  ))

  if (!inherits(keep_features, "try-error")) {
    use_features = rownames(safe_matrix)[keep_features]
  }
  safe_matrix = safe_matrix[use_features, , drop = FALSE]

  if (nrow(safe_matrix) == 0) {
    se_data = create_empty_se(
      processed_data,
      why = "no features present in at least 1 sample"
    )
    return(se_data)
  }

  feature_data = metabolite_data |>
    dplyr::filter(feature_id %in% rownames(safe_matrix)) |>
    DataFrame()
  rownames(feature_data) = feature_data$feature_id
  feature_data = feature_data[rownames(safe_matrix), ]

  ssf = clean_ssf_factors(ssf_data)

  if (length(unique(ssf$subject)) == 1) {
    ssf$subject = NULL
    id_cols = "sample_id"
  } else {
    id_cols = c("subject", "sample_id")
  }
  factor_type = ssf |>
    dplyr::select(factor, type) |>
    dplyr::distinct()

  ssf_wide = widen_ssf_json(ssf, id_cols)

  if (nrow(ssf_wide) > length(unique(ssf$sample_id))) {
    se_data = create_empty_se(
      processed_data,
      why = "duplicate samples somewhere"
    )
    return(se_data)
  }
  is_list = purrr::map_lgl(ssf_wide, is.list)

  if (any(is_list)) {
    se_data = create_empty_se(
      processed_data,
      why = "duplicate sample attributes"
    )
    return(se_data)
  }
  ssf_wide = ssf_wide |>
    dplyr::filter(sample_id %in% colnames(safe_matrix))
  safe_matrix = safe_matrix[, ssf_wide$sample_id, drop = FALSE]

  use_for_factors = c(id_cols, grab_factor_types(factor_type))
  ssf_concat = create_ssf_groups_json(
    ssf_wide |> dplyr::select(tidyselect::all_of(use_for_factors)),
    id_cols
  )

  other_factors = c(
    "sample_id",
    grab_factor_types(factor_type, "additional_sample_data")
  )
  if (length(other_factors) > 1) {
    ssf_use = dplyr::left_join(
      ssf_wide[, other_factors],
      ssf_concat,
      by = "sample_id"
    )
  } else {
    ssf_use = ssf_concat
  }
  factor_type = dplyr::bind_rows(
    factor_type,
    data.frame(factor = c("ssf", "experiment_or_qc"), type = "derived")
  )

  pooled_or_qc = determine_qc_samples(
    ssf |> dplyr::filter(sample_id %in% ssf_use$sample_id)
  )

  # check if we accidentally grabbed all the samples as QC
  use_setdiff = base::setdiff(ssf_use$sample_id, pooled_or_qc)

  # try to mitigate by removing the additional data thing from the long SSF.
  if (length(use_setdiff) == 0) {
    ssf_trim = ssf |>
      dplyr::filter(
        type %in% "subject_sample_factor",
        sample_id %in% ssf_use$sample_id
      )
    pooled_or_qc = determine_qc_samples(ssf_trim)
    use_setdiff = base::setdiff(ssf_use$sample_id, pooled_or_qc)
    # and if nothing seems right, then we keep everything as experiment
    if (length(use_setdiff) == 0) {
      pooled_or_qc = vector("character")
    }
  }

  ssf_use = ssf_use |>
    dplyr::mutate(
      experiment_or_qc = dplyr::case_when(
        sample_id %in% pooled_or_qc ~ "qc",
        TRUE ~ "experiment"
      )
    )

  se_out = SummarizedExperiment(
    assays = SimpleList(raw = safe_matrix),
    rowData = feature_data,
    colData = ssf_use,
    metadata = list(
      mwtab = clean_processed_data(processed_data),
      factor_types = factor_type,
      qc = list(
        id = processed_data$ID,
        failed = "",
        succeeded = "create_se",
        do_next = "go"
      )
    )
  )
  se_out
}

clean_processed_data = function(processed_data) {
  if (!is.null(processed_data$MS_METABOLITE_DATA)) {
    processed_data$MS_METABOLITE_DATA = NULL
  }
  if (!is.null(processed_data$NMR_METABOLITE_DATA)) {
    processed_data$NMR_METABOLITE_DATA = NULL
  }
  if (!is.null(processed_data$NMR_BINNED_DATA)) {
    processed_data$NMR_BINNED_DATA = NULL
  }
  processed_data$MEASUREMENTS = NULL
  processed_data$METABOLITES = NULL
  processed_data
}

create_empty_se = function(processed_data, why = "") {
  processed_data = clean_processed_data(processed_data)
  se_data = SummarizedExperiment()
  metadata(se_data) = list(
    mwtab = processed_data,
    qc = list(
      id = processed_data$ID,
      failed = "create_se",
      why = why,
      succeeded = "none",
      do_next = "stop"
    )
  )
  se_data
}

grab_factor_types = function(
  factor_type,
  use_type = "subject_sample_factor"
) {
  keep_types = factor_type |>
    dplyr::filter(type %in% use_type) |>
    dplyr::pull(factor)
  keep_types
}

#' does two things:
#' 1 - checks for other factors that match the special
#'   "sample_id", and adds numbers to them to make them
#'   different.
#' 2 - checks for factors that become duplicates of each
#'   other after janitor::make_clean_names, for example
#'   "Raw file name" and "RAW_FILE_NAME", and makes them
#'   not duplicated.
clean_ssf_factors = function(ssf_data) {
  sample_factor = tibble::tibble(factor = "sample_id", type = "weird")
  distinct_ssf_factor = dplyr::bind_rows(sample_factor, ssf_data) |>
    dplyr::select(factor, type) |>
    dplyr::distinct() |>
    dplyr::arrange(dplyr::desc(type)) |>
    dplyr::mutate(
      clean_factor = janitor::make_clean_names(factor, allow_dupes = TRUE),
      clean_numbered = janitor::make_clean_names(factor, allow_dupes = FALSE)
    )

  if (any(distinct_ssf_factor$clean_numbered %in% "subject")) {
    loc_subject = which(distinct_ssf_factor$clean_numbered %in% "subject")
    distinct_ssf_factor$clean_numbered[loc_subject] = "subject_1"
  }

  for (irow in seq_len(nrow(distinct_ssf_factor))) {
    check_clean = distinct_ssf_factor$clean_numbered[irow]
    check_rows = rep(FALSE, nrow(distinct_ssf_factor))
    check_rows[
      (distinct_ssf_factor$type %in% distinct_ssf_factor$type[irow])
    ] = TRUE
    check_rows[irow] = FALSE
    if (any(check_clean %in% distinct_ssf_factor$factor[check_rows])) {
      split_under = stringr::str_split(check_clean, "_")[[1]]
      replace_val = "asd"
      split_rep = split_under
      split_rep[length(split_rep)] = replace_val

      new_clean = paste0(split_rep, collapse = "_")

      if (!any(new_clean %in% distinct_ssf_factor$factor[check_rows])) {
        distinct_ssf_factor$clean_numbered[irow] = new_clean
      } else {
        split_rep[length(split_rep)] = "z100"
        new_clean = paste0(split_rep, collapse = "_")
        if (!any(new_clean %in% distinct_ssf_factors$factor[check_rows])) {
          distinct_ssf_factors$clean_numbered[irow] = new_clean
        } else {
          warning("Coulnd't create a clean_factor that doesn't already exist.")
        }
      }
    }
  }

  for (iid in seq_len(nrow(distinct_ssf_factor))) {
    match_locs = (ssf_data$factor %in% distinct_ssf_factor$factor[iid]) &
      (ssf_data$type %in% distinct_ssf_factor$type[iid])
    ssf_data$factor[match_locs] = distinct_ssf_factor$clean_numbered[iid]
  }

  ssf_data
}

intersect_samples = function(abundance_data, sample_ids) {
  keep_samples = base::intersect(colnames(abundance_data), sample_ids)
  keep_cols = c(keep_samples, "feature_id")
  abundance_data[, keep_cols]
}
