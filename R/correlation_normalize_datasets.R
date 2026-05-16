sets_run_group_correlations = function(in_setup) {
  # in_setup = tar_read(setup_g01)

  cor_data = purrr::map(
    in_setup,
    do_correlation
  )
}

get_just_values = function(in_df) {
  df_matrix = in_df |>
    dplyr::select(-feature_id) |>
    as.matrix()
  rownames(df_matrix) = in_df$feature_id
  df_matrix
}

median_normalize_matrix = function(in_matrix) {
  med_values = apply(in_matrix, 2, median, na.rm = TRUE)
  med_matrix = matrix(
    med_values,
    nrow = nrow(in_matrix),
    ncol = ncol(in_matrix),
    byrow = TRUE
  )
  norm_matrix = in_matrix / med_matrix
  return(norm_matrix)
}

count_missingness_testleftcensored = function(raw_matrix, ssf_wide) {
  perc_missing = sum(is.na(raw_matrix) | (raw_matrix == 0)) /
    (nrow(raw_matrix) * ncol(raw_matrix)) *
    100

  raw_matrix = raw_matrix[, ssf_wide$sample_id]

  if (perc_missing > 1) {
    lc_results = try(ICIKendallTau::test_left_censorship(
      raw_matrix,
      sample_classes = ssf_wide$factors
    ))
    if (inherits(lc_results, "try-error")) {
      out_res = list(lc_test = NULL, perc_missing = perc_missing)
    } else {
      out_res = list(lc_test = lc_results, perc_missing = perc_missing)
    }
  } else {
    out_res = list(lc_test = NULL, perc_missing = perc_missing)
  }
  out_res
}

do_correlation = function(setup_data) {
  # setup_data = in_group[[1]]
  raw_matrix = get_just_values(setup_data$MEASUREMENTS)

  n_cores = future::nbrOfWorkers()

  n_item = ncol(raw_matrix) * nrow(raw_matrix)
  # try halving number of workers to get around memory issues
  if (n_item > 1e6) {
    plan(multicore, workers = ceiling(n_cores / 2))
  }

  norm_matrix = median_normalize_matrix(raw_matrix)
  norm_df = as.data.frame(norm_matrix)
  norm_df$feature_id = setup_data$MEASUREMENTS$feature_id
  setup_ici = ICIKendallTau::ici_kendalltau(raw_matrix, return_matrix = FALSE)
  ssf_wide = widen_ssf_json(
    setup_data$SUBJECT_SAMPLE_FACTORS |>
      dplyr::filter(type %in% "subject_sample_factor")
  ) |>
    create_factor_groups_json()
  # make sure to put them back
  plan(multicore, workers = n_cores)

  cor_data = list(
    id = setup_data$ID,
    raw = setup_data$MEASUREMENTS,
    normalized = norm_df,
    icikt_cor = setup_ici$cor,
    icikt_runtime = setup_ici$run_time,
    ssf = ssf_wide
  )
  cor_data
}


sets_add_mwtab_do_left_censored = function(cor_set, setup_set) {
  # cor_set = tar_read(cor_g01); setup_set = tar_read(setup_g01)

  lc_set = purrr::map2(cor_set, setup_set, add_mwtab_do_left_censored)
  lc_set
}

add_mwtab_do_left_censored = function(cor_data, setup_data) {
  # cor_data = cor_set[[11]]; setup_data = setup_set[[11]]
  stopifnot(cor_data$id == setup_data$ID)
  raw_matrix = get_just_values(cor_data$raw)
  ssf_wide = cor_data$ssf

  feature_data = setup_data$METABOLITES

  matching_features = base::intersect(
    feature_data$feature_id,
    rownames(raw_matrix)
  )

  feature_data = feature_data |>
    dplyr::filter(feature_id %in% matching_features)

  raw_matrix = raw_matrix[feature_data$feature_id, ssf_wide$sample_id]

  norm_matrix = get_just_values(cor_data$normalized)[
    feature_data$feature_id,
    ssf_wide$sample_id
  ]

  qc_samples = determine_qc_samples(setup_data$SUBJECT_SAMPLE_FACTORS)
  ssf_wide = ssf_wide |>
    dplyr::mutate(
      experiment_or_qc = dplyr::case_when(
        sample_id %in% qc_samples ~ "qc",
        TRUE ~ "experiment"
      )
    )
  se_data = SummarizedExperiment(
    assays = SimpleList(raw = raw_matrix, normalized = norm_matrix),
    rowData = feature_data,
    colData = ssf_wide
  )

  lc_data = count_missingness_testleftcensored(raw_matrix, ssf_wide)

  cor_data$raw = NULL
  cor_data$normalized = NULL
  cor_data$ssf = NULL
  cor_data$abundance = se_data
  cor_data$mwtab = setup_data
  cor_data$left_censorship = lc_data
  cor_data
}

determine_qc_samples = function(ssf_long) {
  # ssf_long = in_lc$mwtab$SUBJECT_SAMPLE_FACTORS
  is_pooled_qc = ssf_long |>
    dplyr::filter(
      grepl(
        ".*pool.*|.*qc.*|.*blank.*",
        sample_id,
        ignore.case = TRUE
      ) |
        grepl(".*pool.*|.*qc.*|.*blank.*", item, ignore.case = TRUE)
    ) |>
    dplyr::pull(sample_id)
}
