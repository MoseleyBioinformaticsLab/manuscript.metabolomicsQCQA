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

calculate_mean_median_ratio = function(safe_matrix) {
  safe_median = median(safe_matrix, na.rm = TRUE)
  safe_mean = mean(safe_matrix, na.rm = TRUE)
  abs(safe_mean / safe_median)
}

sets_check_mmr_negatives = function(in_set) {
  # in_set = tar_read(se_g01); tar_source("R");
  sets_mmr = furrr::future_map(in_set, check_mmr_negatives) |>
    purrr::list_rbind()
  sets_mmr
}

check_mmr_negatives = function(in_se) {
  # in_se = in_set[[1]]
  qc_info = get_qc(in_se) |> tibble::as_tibble()
  if (qc_info$do_next %in% "stop") {
    return(NULL)
  }
  raw_data = assays(in_se)$raw

  mmr_value = calculate_mean_median_ratio(raw_data)
  n_vals = nrow(raw_data) * ncol(raw_data)
  n_miss = sum(is.na(raw_data))
  n_neg = sum(raw_data < 0, na.rm = TRUE)
  max_value = max(raw_data, na.rm = TRUE)
  min_value = min(raw_data, na.rm = TRUE)

  dplyr::bind_cols(
    qc_info,
    tibble::tibble(
      mmr = mmr_value,
      N = n_vals,
      N_present = n_vals - n_miss,
      N_neg = n_neg,
      max = max_value,
      min = min_value
    )
  )
}


working_mmr_data = function(se_mmr_combine) {
  use_mmr = se_mmr_combine |>
    dplyr::filter(mmr > 0.01, mmr <= 200) |>
    dplyr::mutate(raw_max = 2^max)

  # if we think things might possibly be log transformed, then a max cutoff of
  # 100 for the max intensity is probably reasonable as well.
  use_mmr = use_mmr |>
    dplyr::filter(max <= 100)
}
