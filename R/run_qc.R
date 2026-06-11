sets_convert_to_se = function(setup_set) {
  se_set = purrr::map(setup_set, convert_to_se)
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

convert_to_se = function(setup_data) {
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

  if (!is.null(setup_data$MEASUREMENTS)) {
    raw_matrix = get_just_values(setup_data$MEASUREMENTS)

    n_feature = nrow(raw_matrix)
    n_sample = ncol(raw_matrix)
  } else {
    se_data = create_empty_se(setup_data, why = "no feature data")
    return(se_data)
  }

  if (n_feature == 0) {
    filtered_features = setup_data$FILTERED_FEATURES
    if (filtered_features$kept[1] == 0) {
      why = "no features; didn't match metabolites"
    } else {
      why = "no features; unknown"
    }
    se_data = create_empty_se(setup_data, why = why)
    return(se_data)
  }

  if (n_sample == 0) {
    filtered_samples = setup_data$FILTERED_SAMPLES
    if (filtered_samples$kept[2] == 0) {
      why = "no samples; didn't match ssf"
    } else {
      why = "no samples; unknown"
    }
    se_data = create_empty_se(setup_data, why = why)
    return(se_data)
  }

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
    na_vals = is.na(safe_matrix)
    if (all(na_vals)) {
      se_data = create_empty_se(
        setup_data,
        why = "all character values when converting to matrix"
      )
      return(se_data)
    }
  } else {
    safe_matrix = raw_matrix
  }

  feature_data = setup_data$METABOLITES |> DataFrame()
  rownames(feature_data) = feature_data$feature_id
  feature_data = feature_data[rownames(safe_matrix), ]

  ssf = setup_data$SUBJECT_SAMPLE_FACTORS
  if (is.null(ssf)) {
    se_data = create_empty_se(setup_data, why = "no subject sample factors")
    return(se_data)
  }
  ssf = clean_ssf_factors(ssf)

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
    se_data = create_empty_se(setup_data, why = "duplicate samples somewhere")
    return(se_data)
  }
  is_list = purrr::map_lgl(ssf_wide, is.list)
  has_p_n = grepl("\\+$|\\-$", ssf$item)
  if (any(is_list) && any(has_p_n)) {
    tmp_item = ssf$item
    tmp_item = gsub("\\+$", "_p", tmp_item)
    tmp_item = gsub("\\-$", "_n", tmp_item)
    ssf2 = ssf
    ssf2$item = tmp_item
    ssf2_wide = widen_ssf_json(ssf2, id_cols)
    is_list2 = purrr::map_lgl(ssf2_wide, is.list)
    if (any(is_list2)) {
      se_data = create_empty_se(setup_data, why = "duplicate sample attributes")
      return(se_data)
    }
  }
  if (any(is_list) && !any(has_p_n)) {
    se_data = create_empty_se(setup_data, why = "duplicate sample attributes")
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

  pooled_or_qc = determine_qc_samples(ssf)

  ssf_use = ssf_use |>
    dplyr::mutate(
      experiment_or_qc = dplyr::case_when(
        sample_id %in% pooled_or_qc ~ "qc",
        TRUE ~ "experiment"
      )
    )

  if (!is.null(setup_data$MS_METABOLITE_DATA)) {
    setup_data$MS_METABOLITE_DATA = NULL
  }
  if (!is.null(setup_data$NMR_METABOLITE_DATA)) {
    setup_data$NMR_METABOLITE_DATA = NULL
  }
  if (!is.null(setup_data$NMR_BINNED_DATA)) {
    setup_data$NMR_BINNED_DATA = NULL
  }
  setup_data$MEASUREMENTS = NULL
  setup_data$METABOLITES = NULL

  se_out = SummarizedExperiment(
    assays = SimpleList(raw = safe_matrix),
    rowData = feature_data,
    colData = ssf_use,
    metadata = list(
      mwtab = setup_data,
      factor_types = factor_type,
      qc = list(
        id = setup_data$ID,
        failed = "",
        succeeded = "create_se",
        do_next = "go"
      )
    )
  )
  se_out
}

create_empty_se = function(setup_data, why = "") {
  se_data = SummarizedExperiment()
  metadata(se_data) = list(
    mwtab = setup_data,
    qc = list(
      id = setup_data$ID,
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
