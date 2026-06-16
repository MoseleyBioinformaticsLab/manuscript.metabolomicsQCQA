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
