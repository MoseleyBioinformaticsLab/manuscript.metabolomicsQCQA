calculate_dynamic_range_check_neg = function(se_data) {
  # first step is always to get the qc metadata, and
  # check if we can even do this step.
  # if not, we just return the data as is.
  qc = metadata(se_data)$qc
  if (qc$do_next %in% "stop") {
    return(se_data)
  }

  # we should actually check for too much negative or missing values
  # right here, because if we need to exclude a sample, we want to
  # exclude it right now, before we try and do left-censorship, correlations and PCA,
  # and set up the missingness correctly.
  data_matrix = assays(se_data)$raw
  data_matrix[data_matrix == 0] = NA
  assays(se_data)$raw = data_matrix

  sample_mins = MatrixGenerics::colMins(data_matrix, na.rm = TRUE)
  sample_maxs = MatrixGenerics::colMaxs(data_matrix, na.rm = TRUE)
  se_data$value_min = sample_mins
  se_data$value_max = sample_maxs
  se_data$dynamic_range = se_data$value_max / se_data$value_min

  data_range = range(data_matrix, na.rm = TRUE)
  data_dyn_range = data_range[2] / data_range[1]
  col_na = MatrixGenerics::colCounts(data_matrix, value = NA)
  data_sign = sign(data_matrix)
  col_neg = MatrixGenerics::colCounts(data_matrix, value = -1, na.rm = TRUE)

  n_feature = nrow(data_matrix)
  perc_na = col_na / n_feature
  perc_neg = col_neg / n_feature

  se_data$removed = NA_character_
  se_data$removed[perc_na >= 0.9] = "> 90% missing"
  se_data$removed[perc_neg >= 0.1] = "> 10% negative values"

  metadata(se_data)$factor_types = rbind(
    metadata(se_data)$factor_types,
    data.frame(
      factor = c("value_min", "value_max", "dynamic_range", "removed"),
      type = "computed"
    )
  )

  # data_matrix[data_sign == -1] = NA
  data_matrix = data_matrix[, is.na(se_data$removed), drop = FALSE]
  qc$succeeded = c(qc$succeeded, "calculate_dynamic_range")
  if ((nrow(data_matrix) >= 10) && (ncol(data_matrix) >= 3)) {
    if (sum(data_sign == -1, na.rm = TRUE) >= 2) {
      qc$do_next = "stop"
      qc$why = "negative values present"
    } else {
      qc$do_next = "go"
    }
  } else {
    qc$do_next = "stop"
    qc$why = "too few features or samples"
  }
  metadata(se_data)$qc = qc
  return(se_data)
}


calculate_dynamic_range = function(se_data) {
  # first step is always to get the qc metadata, and
  # check if we can even do this step.
  # if not, we just return the data as is.
  qc = metadata(se_data)$qc
  if (qc$do_next %in% "stop") {
    return(se_data)
  }

  # we should actually check for too much negative or missing values
  # right here, because if we need to exclude a sample, we want to
  # exclude it right now, before we try and do left-censorship, correlations and PCA,
  # and set up the missingness correctly.
  data_matrix = assays(se_data)$raw
  data_matrix[data_matrix == 0] = NA
  assays(se_data)$raw = data_matrix

  sample_mins = MatrixGenerics::colMins(data_matrix, na.rm = TRUE)
  sample_maxs = MatrixGenerics::colMaxs(data_matrix, na.rm = TRUE)
  se_data$value_min = sample_mins
  se_data$value_max = sample_maxs
  se_data$dynamic_range = se_data$value_max / se_data$value_min

  data_range = range(data_matrix, na.rm = TRUE)
  data_dyn_range = data_range[2] / data_range[1]
  col_na = MatrixGenerics::colCounts(data_matrix, value = NA)
  data_sign = sign(data_matrix)
  col_neg = MatrixGenerics::colCounts(data_matrix, value = -1, na.rm = TRUE)

  n_feature = nrow(data_matrix)
  perc_na = col_na / n_feature
  perc_neg = col_neg / n_feature

  se_data$removed = NA_character_
  se_data$removed[perc_na >= 0.9] = "> 90% missing"
  se_data$removed[perc_neg >= 0.1] = "> 10% negative values"

  metadata(se_data)$factor_types = rbind(
    metadata(se_data)$factor_types,
    data.frame(
      factor = c("value_min", "value_max", "dynamic_range", "removed"),
      type = "computed"
    )
  )

  # data_matrix[data_sign == -1] = NA
  data_matrix = data_matrix[, is.na(se_data$removed), drop = FALSE]
  qc$succeeded = c(qc$succeeded, "calculate_dynamic_range")
  if ((nrow(data_matrix) >= 10) && (ncol(data_matrix) >= 3)) {
    if (sum(data_sign == -1, na.rm = TRUE) >= 2) {
      qc$do_next = "stop"
      qc$why = "negative values present"
    } else {
      qc$do_next = "go"
    }
  } else {
    qc$do_next = "stop"
    qc$why = "too few features or samples"
  }
  metadata(se_data)$qc = qc
  return(se_data)
}

calculate_rsd_matrix = function(
  abundance_matrix,
  sample_metadata,
  all_factors
) {
  sample_metadata = sample_metadata |>
    dplyr::arrange(factors, sample_id)
  abundance_matrix = abundance_matrix[, sample_metadata$sample_id]

  sample_metadata$factors = factor(
    sample_metadata$factors,
    levels = all_factors
  )

  abundance_df = abundance_matrix |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "feature_id") |>
    tidyr::pivot_longer(cols = -feature_id, names_to = "sample_id")
  abundance_df = dplyr::left_join(
    abundance_df,
    sample_metadata[, c("sample_id", "factors")],
    by = "sample_id"
  )

  abundance_summary = abundance_df |>
    dplyr::filter(!is.na(value)) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(value),
      sd = sd(value),
      rsd = sd / mean,
      .by = c(feature_id, factors)
    )
  abundance_summary
}

create_rsd_plot = function(all_rsd) {
  # all_rsd = use_pca$pca_groups$all$rsd$all_rsd

  rsd_plot = all_rsd |>
    ggplot(aes(x = rsd, y = factors_data, color = data, group = factors_data)) +
    geom_sina() +
    geom_boxplot(alpha = 0, color = "black") +
    theme(
      legend.position = "top",
      legend.location = "plot",
      legend.justification.top = "left",
      plot.title.position = "plot"
    )
  rsd_plot
}

create_rsd_comparison_table = function(all_rsd_median) {
  # all_rsd_median = use_pca$pca_groups$all$rsd$all_rsd_median
  rsd_compare_table = all_rsd_median |>
    tidyr::pivot_wider(
      id_cols = factors,
      names_from = data,
      values_from = median
    ) |>
    dplyr::transmute(
      factors = factors,
      raw = r,
      normalized = n,
      `normalized - raw` = normalized - raw
    ) |>
    gt::gt() |>
    gt::fmt_number(decimals = 2)
  rsd_compare_table
}
