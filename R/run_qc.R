sets_readrds_convert_to_se = function(set_df, processed_group) {
  if (inherits(set_df, "character")) {
    stop(
      "Run to processed_group first, and then the rest!  tar_make(processed_group)"
    )
  }
  all_files = set_df$file
  all_md5 = set_df$md5
  all_id = set_df$id
  all_se = furrr::future_map2(all_files, all_md5, readrds_convert_to_se)
  names(all_se) = all_id
  all_se
}

sets_run_dyn_neg = function(se_set) {
  # setup_set = tar_read(setup_g01); tar_source("R")
  # setup_set = tar_read(setup_g02); tar_source("R")

  # then we do dynamic range calculations for those that can.
  dyn_range_set = purrr::map(se_set, calculate_dynamic_range_check_neg)
}

sets_run_qc_neg = function(se_set) {
  # setup_set = tar_read(setup_g01); tar_source("R")
  # setup_set = tar_read(setup_g02); tar_source("R")

  # then we do dynamic range calculations for those that can.
  dyn_range_set = purrr::map(se_set, calculate_dynamic_range)
  qc_info = get_qc_outputs(dyn_range_set)

  neg_qc = qc_info |>
    dplyr::filter(why %in% "negative values present") |>
    dplyr::pull(id) |>
    unique()

  neg_list = dyn_range_set[neg_qc]
  neg_list
}

get_qc = function(se_data) {
  metadata(se_data)$qc
}


get_qc_outputs = function(run_set) {
  # run_set = tar_read(run_g01); tar_source("R")
  all_qc = purrr::map(run_set, \(x) {
    as.data.frame(get_qc(x))
  }) |>
    dplyr::bind_rows()
  all_qc
}

check_qc_se_1_3 = function() {
  qc_run1 = tar_read(neg_g01) |> get_qc_outputs()
  qc_run2 = tar_read(neg_g02) |> get_qc_outputs()
  qc_run3 = tar_read(neg_g03) |> get_qc_outputs()
  all_qc = dplyr::bind_rows(qc_run1, qc_run2, qc_run3)
  n_id = length(unique(all_qc$id))
  all_qc |>
    dplyr::summarise(N = dplyr::n(), .by = why) |>
    dplyr::mutate(perc = N / n_id * 100) |>
    paint::unpaint()
}

check_convert_character_matrix = function(raw_matrix) {
  if (is.character(raw_matrix[1, 1])) {
    safe_matrix = apply(
      raw_matrix,
      2,
      \(x) {
        stringr::str_replace(x, "\\<|\\>", "-")
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

readrds_convert_to_se = function(rds_file, md5_sum) {
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
  processed_data = readRDS(rds_file)
  check_md5(processed_data, md5_sum)

  ssf_data = processed_data$SUBJECT_SAMPLE_FACTORS
  abundance_data = processed_data$MEASUREMENTS
  metabolite_data = processed_data$METABOLITES

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

  low_ratio = 0.95
  hi_ratio = 1.05
  mean_median_ratio = calculate_mean_median_ratio(safe_matrix)
  out_matrix = safe_matrix
  processed_data$transform = "none"
  if (
    !(is.na(mean_median_ratio)) ||
      is.nan(mean_median_ratio) ||
      is.infinite(mean_median_ratio)
  ) {
    if ((mean_median_ratio >= low_ratio) && (mean_median_ratio <= hi_ratio)) {
      out_matrix = 2^safe_matrix
      processed_data$transform = "exponentiated"
    }
  }
  processed_data$mean_median_ratio = mean_median_ratio

  feature_data = metabolite_data |>
    dplyr::filter(feature_id %in% rownames(out_matrix)) |>
    DataFrame()
  rownames(feature_data) = feature_data$feature_id
  feature_data = feature_data[rownames(out_matrix), ]

  ssf = clean_ssf_factors(ssf_data)

  if (length(unique(ssf$subject)) == 1) {
    ssf$subject = NULL
    id_cols = "sample_id"
  } else {
    id_cols = c("subject", "sample_id")
  }
  factor_type = ssf |>
    dplyr::select(factor, type) |>
    dplyr::distinct() |>
    dplyr::mutate(factor = janitor::make_clean_names(factor))

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
    dplyr::filter(sample_id %in% colnames(out_matrix))
  out_matrix = out_matrix[, ssf_wide$sample_id, drop = FALSE]

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

  pooled_or_qc = determine_qc_samples(ssf)

  ssf_use = ssf_use |>
    dplyr::mutate(
      experiment_or_qc = dplyr::case_when(
        sample_id %in% pooled_or_qc ~ "qc",
        TRUE ~ "experiment"
      )
    )

  se_out = SummarizedExperiment(
    assays = SimpleList(raw = out_matrix),
    rowData = feature_data,
    colData = ssf_use,
    metadata = list(
      mwtab = clean_processed_data(processed_data),
      factor_types = factor_type,
      transformation = list(
        mean_median_ratio = mean_median_ratio,
        transform = transform
      ),
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
clean_ssf_factors = function(ssf) {
  sample_factor = tibble::tibble(factor = "sample_id", type = "weird")
  distinct_ssf_factor = dplyr::bind_rows(sample_factor, ssf) |>
    dplyr::select(factor, type) |>
    dplyr::distinct() |>
    dplyr::arrange(dplyr::desc(type)) |>
    dplyr::mutate(
      clean_factor = janitor::make_clean_names(factor, allow_dupes = TRUE),
      clean_numbered = janitor::make_clean_names(factor, allow_dupes = FALSE)
    )

  replace_now = distinct_ssf_factor |>
    dplyr::filter(clean_factor != clean_numbered)
  if (nrow(replace_now) > 0) {
    for (iid in seq_len(nrow(replace_now))) {
      ssf$factor = gsub(
        replace_now$factor[iid],
        replace_now$clean_numbered[iid],
        ssf$factor,
        fixed = TRUE
      )
    }
  }
  ssf
}

calculate_mean_median_ratio = function(safe_matrix) {
  safe_median = median(safe_matrix, na.rm = TRUE)
  safe_mean = mean(safe_matrix, na.rm = TRUE)
  abs(safe_mean / safe_median)
}

run_qc = function(in_data) {
  # the idea here is going to be:
  #   - check if we can create the SE object;
  #   - create the SE object;
  #   - set all the possible qc steps to FALSE;
  #   - try the first one, if it succeeds or returns minimum for next;
  #   - do the next one;
  #   - and so on.
  #   - if at any point you get a NO on the next one, we exit and return the object
  #     with what was actually done.
  #   - now, the NO will have to be really bad to stop it.
  #   If we have 10 features and 3 samples, we are going to go for it.
  #   But the first check is absolutely going to be, for each sample, are there negative values, or missing values, and what percentage, so we can exclude it right there and then.
}
