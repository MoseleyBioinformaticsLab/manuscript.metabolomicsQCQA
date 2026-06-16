sets_run_setup = function(in_df, tmp_sets) {
  # group_sets = tar_read(processed_group); in_df = group_sets$set[[1]]; tar_source("R")
  # group_sets = tar_read(processed_group); in_df = group_sets$set[[29]]; tar_source("R")
  if (inherits(in_df, "character")) {
    stop("Run to tmp_sets first, and then the rest!  tar_make(tmp_sets)")
  }
  all_files = in_df$file
  all_md5 = in_df$md5
  all_id = in_df$id
  all_up = furrr::future_map2(all_files, all_md5, read_rds_setup_data)
  names(all_up) = all_id
  all_up
}

read_rds_setup_data = function(rds_file, md5_sum) {
  # rds_file = "parsed_rds/AN000028.rds"
  # md5_sum = "0e22eb16e4debdc1fc55a4c2b58767df"

  # rds_file = "parsed_rds/AN000112.rds"
  # rds_file = "parsed_rds/AN003604.rds"

  # id = "AN000076"; rds_file = "parsed_rds/AN000076.rds"; md5_sum = "13e4670918e908f47145b80e0669421c"; tar_source("R")

  # id = "AN000072"; rds_file = "parsed_rds/AN000072.rds"; md5_sum = "781a19811ab6cf3ca6524224ab2e3507"; tar_source("R")
  # rds_file = "parsed_rds/AN000163.rds"; id = "AN000163"; md5_sum = "92913d26c46f3e0f07a3dab936267e6c"; tar_source("R")
  # rds_file = "parsed_rds/AN001953.rds"; tar_source("R")
  # rds_file = "parsed_rds/AN003244.rds"; tar_source("R")
  # this next one had a really weird thing due to not handling errors correctly.
  # rds_file = "parsed_rds/AN004699.rds"; tar_source("R")
  processed_data = readRDS(rds_file)
  check_md5(processed_data, md5_sum)

  processed_data = filter_matching_samples_features(processed_data)
  processed_data
}

check_md5 = function(in_data, md5_sum) {
  new_md5 = digest::digest(in_data, algo = "md5")
  stopifnot(new_md5 == md5_sum)
  return(invisible(NULL))
}

transform_abundances = function(raw_matrix) {
  low_ratio = 0.95
  hi_ratio = 1.05
  mean_median_ratio = calculate_mean_median_ratio(raw_matrix)
  out_matrix = safe_matrix
  transformation = "none"
  if (
    !(is.na(mean_median_ratio) ||
      is.nan(mean_median_ratio) ||
      is.infinite(mean_median_ratio))
  ) {
    if ((mean_median_ratio >= low_ratio) && (mean_median_ratio <= hi_ratio)) {
      out_matrix = 2^raw_matrix
      transformation = "exponentiated"
    }
  }
  list(
    matrix = out_matrix,
    transformation = transformation,
    mean_median_ratio = mean_median_ratio
  )
}

intersect_samples = function(abundance_data, sample_ids) {
  keep_samples = base::intersect(colnames(abundance_data), sample_ids)
  keep_cols = c(keep_samples, "feature_id")
  abundance_data[, keep_cols]
}

#' filters measurements and ssf to matching samples
#'
#' Often, the analyses have sample factors for samples
#' that aren't necessarily in the measurements, or vice versa.
#' So this code takes the analysis, and filters the subject
#' factors and measurements to the common samples between them.
filter_matching_samples_features = function(processed_data) {
  ssf_data = processed_data$SUBJECT_SAMPLE_FACTORS
  abundance_data = processed_data$MEASUREMENTS
  metabolite_data = processed_data$METABOLITES

  if ((is.null(ssf_data) || is.null(abundance_data))) {
    return(processed_data)
  }

  n_feature_abund = nrow(abundance_data)
  n_feature_metabolite = nrow(metabolite_data)
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

  match_features = base::intersect(
    abundance_data$feature_id,
    metabolite_data$feature_id
  )
  abundance_data = abundance_data |>
    dplyr::filter(feature_id %in% match_features)
  metabolite_data = metabolite_data |>
    dplyr::filter(feature_id %in% match_features)

  if ((nrow(abundance_data) > 0) && (ncol(abundance_data) > 1)) {
    abundance_transformed = transform_abundances(abundance_data)
    abundance_data = abundance_transformed$measurements

    raw_matrix = get_just_values(abundance_data)
    keep_things = try(visualizationQualityControl::keep_non_missing_percentage(
      raw_matrix,
      keep_num = 1,
      missing_value = c(0, NA)
    ))
    if (!inherits(keep_things, "try-error")) {
      keep_features = rownames(raw_matrix)[keep_things]
    } else {
      keep_features = abundance_data$feature_id
    }
  } else {
    keep_features = abundance_data$feature_id
    abundance_transformed = list(transform = "none")
  }

  abundance_data = abundance_data |>
    dplyr::filter(feature_id %in% keep_features)
  metabolite_data = metabolite_data |>
    dplyr::filter(feature_id %in% keep_features)

  processed_data$MEASUREMENTS = abundance_data
  processed_data$METABOLITES = metabolite_data
  processed_data$SUBJECT_SAMPLE_FACTORS = ssf_data
  processed_data$TRANSFORM = abundance_transformed$transform

  n_filtered_samples = tibble::tribble(
    ~source        , ~original                  , ~kept                 , ~filtered ,
    "ssf"          , length(all_samples)        , length(match_samples) , ssf_diff  ,
    "measurements" , length(all_abundance_cols) , length(match_samples) , meas_diff
  )
  n_filtered_features = tibble::tribble(
    ~source       , ~original            , ~kept                 ,
    "abundance"   , n_feature_abund      , nrow(abundance_data)  ,
    "metabolites" , n_feature_metabolite , nrow(metabolite_data)
  )
  processed_data$FILTERED_SAMPLES = n_filtered_samples
  processed_data$FILTERED_FEATURES = n_filtered_features
  processed_data
}

sets_check_features = function(in_setup) {
  purrr::imap(in_setup, \(setup, id) {
    tibble::tibble(id = id, n_features = nrow(setup$MEASUREMENTS))
  }) |>
    purrr::list_rbind()
}
