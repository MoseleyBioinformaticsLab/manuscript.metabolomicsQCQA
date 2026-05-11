run_setup = function(in_df, tmp_sets) {
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

  processed_data = readRDS(rds_file)
  check_md5(processed_data, md5_sum)
  processed_data = filter_matching_samples(processed_data)
  transformed_data = transform_abundances(processed_data$MEASUREMENTS)
  processed_data$MEASUREMENTS = transformed_data$measurements
  processed_data$TRANSFORM = transformed_data$transform
  processed_data
}

check_md5 = function(in_data, md5_sum) {
  new_md5 = digest::digest(in_data, algo = "md5")
  stopifnot(new_md5 == md5_sum)
  return(invisible(NULL))
}

transform_abundances = function(in_abundances, max_min = 20) {
  transform = "none"
  abundance_matrix = as.matrix(in_abundances |> dplyr::select(-feature_id))
  value_range = range(abundance_matrix, na.rm = TRUE)
  if (max(value_range) < max_min) {
    sample_cols = colnames(abundance_matrix)
    for (icol in sample_cols) {
      in_abundances[[icol]] = 2^in_abundances[[icol]]
    }
    transform_matrix = as.matrix(in_abundances |> dplyr::select(-feature_id))
    na_org = is.na(abundance_matrix)
    na_new = is.na(transform_matrix)
    if (!all.equal(na_org, na_new)) {
      stop("Transforming made different missingness patterns!")
    }
    transform = "2-exponentiated"
  }
  list(measurements = in_abundances, transform = transform)
}

#' filters measurements and ssf to matching samples
#'
#' Often, the analyses have sample factors for samples
#' that aren't necessarily in the measurements, or vice versa.
#' So this code takes the analysis, and filters the subject
#' factors and measurements to the common samples between them.
filter_matching_samples = function(processed_data) {
  ssf_data = processed_data$SUBJECT_SAMPLE_FACTORS
  abundance_data = processed_data$MEASUREMENTS

  all_samples = ssf_data$sample_id
  all_abundance_cols = colnames(abundance_data)

  unique_samples = unique(ssf_data$sample_id)
  match_samples = purrr::map(unique_samples, \(x) {
    tmp = grep(x, colnames(abundance_data), value = TRUE)
    if (length(tmp) == 1) {
      return(tibble::tibble(ssf = x, measurement_id = tmp))
    } else {
      return(NULL)
    }
  }) |>
    purrr::list_rbind()

  ssf_data = dplyr::left_join(
    match_samples,
    ssf_data,
    by = c("ssf" = "sample_id")
  ) |>
    dplyr::mutate(sample_id = measurement_id)

  abundance_cols = c(match_samples$measurement_id, "feature_id")
  abundance_data = abundance_data[, abundance_cols]

  processed_data$MEASUREMENTS = abundance_data
  processed_data$SUBJECT_SAMPLE_FACTORS = ssf_data
  processed_data
}
