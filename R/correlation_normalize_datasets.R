#' filters measurements and ssf to matching samples
#'
#' Often, the analyses have sample factors for samples
#' that aren't necessarily in the measurements, or vice versa.
#' So this code takes the analysis, and filters the subject
#' factors and measurements to the common samples between them.
filter_matching_samples = function(processed_data) {
  ssf_data = processed_data$SUBJECT_SAMPLE_FACTORS
  abundance_data = processed_data$MEASUREMENTS

  all_samples = unique(ssf_data$sample_id)
  all_abundance_cols = colnames(abundance_data)

  match_samples = base::intersect(all_samples, all_abundance_cols)
  ssf_diff = base::setdiff(all_samples, match_samples) |>
    paste0(collapse = "; ")
  meas_diff = base::setdiff(all_abundance_cols, match_samples) |>
    paste0(collapse = "; ")
  ssf_data = ssf_data |>
    dplyr::filter(sample_id %in% match_samples)

  abundance_cols = c(match_samples, "feature_id")
  abundance_data = abundance_data[, abundance_cols]

  processed_data$MEASUREMENTS = abundance_data
  processed_data$SUBJECT_SAMPLE_FACTORS = ssf_data

  n_filtered = tibble::tribble(
    ~source        , ~original                  , ~kept                 , ~filtered ,
    "ssf"          , length(all_samples)        , length(match_samples) , ssf_diff  ,
    "measurements" , length(all_abundance_cols) , length(match_samples) , meas_diff
  )
  processed_data$FILTERED = n_filtered
  processed_data
}


run_group_correlations = function(in_setup) {
  # in_setup = tar_read(setup_g01)

  cor_data = purrr::map(
    in_setup,
    do_missingness_correlation,
    .progress = TRUE
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

do_missingness_correlation = function(setup_data) {
  # setup_data = in_group[[1]]
  raw_matrix = get_just_values(setup_data$MEASUREMENTS)
  norm_matrix = median_normalize_matrix(raw_matrix)
  norm_df = as.data.frame(norm_matrix)
  norm_df$feature_id = setup_data$MEASUREMENTS$feature_id
  setup_ici = ICIKendallTau::ici_kendalltau(raw_matrix, return_matrix = FALSE)
  ssf_wide = widen_ssf_json(
    setup_data$SUBJECT_SAMPLE_FACTORS |>
      dplyr::filter(type %in% "subject_sample_factor")
  ) |>
    create_factor_groups_json()

  mwtab = setup_data
  mwtab$MEASUREMENTS = NULL
  if ("MS_METABOLITE_DATA" %in% names(mwtab)) {
    mwtab$MS_METABOLITE_DATA = NULL
  }
  if ("NMR_METABOLITE_DATA" %in% names(mwtab)) {
    mwtab$NMR_METABOLITE_DATA = NULL
  }

  lc_data = count_missingness_testleftcensored(raw_matrix, ssf_wide)

  qcqa_data = list(
    id = mwtab$ID,
    raw = setup_data$MEASUREMENTS,
    normalized = norm_df,
    icikt_cor = setup_ici$cor,
    icikt_runtime = setup_ici$run_time,
    leftcensored = lc_data,
    ssf = ssf_wide,
    mwtab = mwtab
  )
  qcqa_data
}
