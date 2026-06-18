sets_calculate_dynamic_range = function(sets_se) {
  # sets_se = tar_read(se_g01); tar_source("R")
  sets_dynamic_range = furrr::future_map(sets_se, calculate_dynamic_range)
  sets_dynamic_range
}


calculate_dynamic_range = function(se_data) {
  # first step is always to get the qc metadata, and
  # check if we can even do this step.
  # if not, we just return the data as is.
  qc = metadata(se_data)$qc
  if (qc$do_next %in% "stop") {
    return(se_data)
  }

  data_matrix = assays(se_data)$raw

  max_data = max(data_matrix, na.rm = TRUE)
  neg_data = sum(data_matrix < 0, na.rm = TRUE)
  n_original_features = nrow(data_matrix)
  n_original_samples = ncol(data_matrix)

  if (max_data <= 1000) {
    transform_info = transform_abundances(data_matrix)
    if (transform_info$transformation %in% "exponentiated") {
      data_matrix = transform_info$matrix
    }
    transform_info$matrix = NULL
  } else {
    transform_info = list(
      transformation = "none",
      mean_median_ratio = calculate_mean_median_ratio(data_matrix)
    )
  }

  transform_info$original_max = max_data
  transform_info$n_neg = neg_data

  data_matrix[data_matrix <= 0] = NA

  # now that we've set the negatives and zeros to missing proper, we can try
  # doing the filtering to being present in at least one sample again.
  data_keep = try(visualizationQualityControl::keep_non_missing_percentage(
    data_matrix,
    keep_num = 1
  ))

  if (!inherits(data_keep, "try-error")) {
    data_matrix = data_matrix[data_keep, , drop = FALSE]
    se_data = se_data[data_keep, ]
  }
  assays(se_data)$raw = data_matrix

  data_normalized = median_normalize_matrix(data_matrix)
  assays(se_data)$normalized = data_normalized
  sample_mins = MatrixGenerics::colMins(data_normalized, na.rm = TRUE)
  sample_maxs = MatrixGenerics::colMaxs(data_normalized, na.rm = TRUE)
  se_data$value_min = sample_mins
  se_data$value_max = sample_maxs
  se_data$dynamic_range = se_data$value_max / se_data$value_min

  data_range = range(data_normalized, na.rm = TRUE)
  data_dyn_range = data_range[2] / data_range[1]
  n_feature = nrow(data_matrix)
  col_not_na = n_feature - MatrixGenerics::colCounts(data_matrix, value = NA)

  se_data$removed = NA_character_
  se_data$removed[col_not_na < 10] = "< 10 present values"

  metadata(se_data)$factor_types = rbind(
    metadata(se_data)$factor_types,
    data.frame(
      factor = c("value_min", "value_max", "dynamic_range", "removed"),
      type = "computed"
    )
  )

  n_new_features = nrow(se_data)
  n_new_samples = sum(is.na(se_data$removed))

  filtering_df = tibble::tibble(
    n_features = c(n_original_features, n_new_features),
    n_samples = c(n_original_samples, n_new_samples),
    which = c("original", "filtered")
  )
  # data_matrix[data_sign == -1] = NA
  data_matrix = data_matrix[, is.na(se_data$removed), drop = FALSE]
  qc$succeeded = c(qc$succeeded, "calculate_dynamic_range / normalized")
  if ((nrow(data_matrix) < 10)) {
    qc$do_next = "stop"
    qc$why = "too few features"
  }
  if (ncol(data_matrix) < 3) {
    qc$do_next = "stop"
    qc$why = "too few samples"
  }
  metadata(se_data)$transformation = transform_info
  metadata(se_data)$dynamic_range_filtering = filtering_df
  metadata(se_data)$qc = qc
  return(se_data)
}

sets_get_negative_information = function(sets_se) {
  # sets_se = tar_read(se_dynamic_range_g14); tar_source("R")
  sets_info = purrr::map(sets_se, get_negative_information) |>
    purrr::list_rbind()
}

get_negative_information = function(se_data) {
  # se_data = sets_se[[1]]
  transformation_info = metadata(se_data)$transformation |>
    tibble::as_tibble()
  progress_info = metadata(se_data)$qc[c("id", "do_next")] |>
    tibble::as_tibble()

  dplyr::bind_cols(transformation_info, progress_info)
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

  max_data = max(data_matrix, na.rm = TRUE)
  neg_data = sum(data_matrix < 0, na.rm = TRUE)

  if (max_data <= 1000) {
    transform_info = transform_abundances(data_matrix)
    if (transform_info$transformation %in% "exponentiated") {
      data_matrix = transform_info$matrix
    }
  } else {
    transform_info = list(
      transformation = "none",
      mean_median_ratio = calculate_mean_median_ratio(data_matrix)
    )
  }

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

  metadata(se_data)$factor_types = rbind(
    metadata(se_data)$factor_types,
    data.frame(
      factor = c("value_min", "value_max", "dynamic_range", "removed"),
      type = "computed"
    )
  )

  # data_matrix[data_sign == -1] = NA
  # data_matrix = data_matrix[, is.na(se_data$removed), drop = FALSE]
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
  assays(se_data)$raw = data_matrix
  metadata(se_data)$qc = qc
  return(se_data)
}
