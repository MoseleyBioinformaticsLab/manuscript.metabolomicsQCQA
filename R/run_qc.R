sets_run_qc = function(setup_set) {
  # setup_set = tar_read(setup_g01); tar_source("R")

  se_set = sets_convert_to_se(setup_set)
}

sets_convert_to_se = function(setup_set) {
  se_set = furrr::future_map(setup_set, convert_to_se)
  se_set
}


get_qc_info = function(se_data) {
  metadata(se_data)$qc
}
convert_to_se = function(setup_data) {
  # setup_set = tar_read(setup_g01); setup_data = setup_set[[1]]; tar_source("R")
  # setup_data = setup_set[[4]]
  # setup_data = setup_set[[55]]
  # setup_data = setup_set[["AN004188"]]

  if (!is.null(setup_data$MEASUREMENTS)) {
    raw_matrix = get_just_values(setup_data$MEASUREMENTS)

    n_feature = nrow(raw_matrix)
    n_sample = ncol(raw_matrix)
  } else {
    n_feature = 0
    n_sample = 0
  }

  if ((n_feature == 0) || (n_sample == 0)) {
    se_data = SummarizedExperiment()
    metadata(se_data) = list(
      mwtab = setup_data,
      qc = list(
        id = setup_data$ID,
        failed = "create_se",
        why = "no features OR no samples",
        succeeded = "none",
        do_next = "stop"
      )
    )
    return(se_data)
  }

  feature_data = setup_data$METABOLITES |> DataFrame()
  rownames(feature_data) = feature_data$feature_id
  feature_data = feature_data[rownames(raw_matrix), ]

  ssf = setup_data$SUBJECT_SAMPLE_FACTORS

  factor_type = ssf |>
    dplyr::select(factor, type) |>
    dplyr::distinct() |>
    dplyr::mutate(factor = janitor::make_clean_names(factor))

  ssf_wide = widen_ssf_json(ssf)
  is_list = purrr::map_lgl(ssf_wide, is.list)
  has_p_n = grepl("\\+$|\\-$", ssf$item)
  if (any(is_list) && any(has_p_n)) {
    tmp_item = ssf$item
    tmp_item = gsub("\\+$", "_p", tmp_item)
    tmp_item = gsub("\\-$", "_n", tmp_item)
    ssf2 = ssf
    ssf2$item = tmp_item
    ssf2_wide = widen_ssf_json(ssf2)
    is_list2 = purrr::map_lgl(ssf2_wide, is.list)
    if (any(is_list2)) {
      se_data = SummarizedExperiment()
      metadata(se_data) = list(
        mwtab = setup_data,
        qc = list(
          failed = "create_se",
          why = "duplicate sample attributes",
          succeeded = "none",
          do_next = "stop"
        )
      )
      return(se_data)
    }
  }
  ssf_wide = ssf_wide |>
    dplyr::filter(sample_id %in% colnames(raw_matrix))
  raw_matrix = raw_matrix[, ssf_wide$sample_id, drop = FALSE]

  use_for_factors = c("subject", "sample_id", grab_factor_types(factor_type))
  ssf_concat = create_ssf_groups_json(
    ssf_wide |> dplyr::select(tidyselect::all_of(use_for_factors))
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
    assays = SimpleList(raw = raw_matrix),
    rowData = feature_data,
    colData = ssf_use,
    metadata = list(
      mwtab = setup_data,
      factor_types = factor_type,
      qc = list(failed = "", succeeded = "create_se", do_next = "go")
    )
  )
  se_out
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
