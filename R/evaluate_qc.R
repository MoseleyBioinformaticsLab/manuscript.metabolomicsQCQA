sets_evaluate_qc = function(sets_se, data_limits) {
  # sets_se = tar_read(se_pca_g01); tar_source("R"); tar_load(data_limits)
  # sets_se = tar_read(se_pca_g02); tar_source("R")
  # sets_se = tar_read(se_pca_g20); se_data = sets_se[[194]]; tar_source("R"); tar_load(data_limits)
  sets_qc = purrr::map(sets_se, \(x) {
    evaluate_qc(x, data_limits) |> tibble::as_tibble()
  }) |>
    purrr::list_rbind()
  sets_qc
}

evaluate_qc = function(se_data, data_limits) {
  # get all the metabolomics workbench info out first
  # the goal here is to have analysis attributes
  # that are searchable and filterable by consumers
  # of the qc-qa analysis

  # se_data = sets_se[[1]]
  # se_data = sets_se[[92]]
  # se_data = sets_se[[123]]
  # se_data = sets_se[[159]]

  use_defaults = FALSE
  mwtab_info = metadata(se_data)$mwtab

  measurement_list = list(
    measurement_type = mwtab_info$MEASUREMENT_TYPE,
    data_source = mwtab_info$DATASOURCE,
    raw_data_available = if (!is.null(mwtab_info$RAW_DATA)) {
      "yes"
    } else {
      "no"
    }
  )
  if (is.null(measurement_list$data_source)) {
    measurement_list$data_source = "none"
  }
  project_list = get_mwtab_project_info(mwtab_info)

  instrument_list = get_mwtab_instrument_chromatography_info(mwtab_info)
  if (is.null(instrument_list)) {
    return(NULL)
  }
  organism_list = get_mwtab_organism_sample(mwtab_info)

  report_list = check_qc_outputs(se_data)
  last_step = get_qc(se_data)
  internal_id = last_step$id

  analysis_list = list(
    completed_steps = paste0(report_list$steps, collapse = "; "),
    failed_steps = report_list$failed[["last"]],
    failed_why = report_list$failed[["why"]]
  )
  mwtab_list = c(
    list(internal_id = internal_id),
    measurement_list,
    project_list,
    instrument_list,
    organism_list
  )

  # next is to start looking for red flags in the data
  # using the various criteria, and fitting them into the
  # same list / data.frame format that allows us to easily
  # filter / group / exclude / include different analyses

  mwtab_data_list = get_mwtab_data_counts(se_data)

  mwtab_list = c(mwtab_list, mwtab_data_list)
  mwtab_features_zero_or_missing = is.na(
    mwtab_data_list$mwtab_matched_features
  ) ||
    (mwtab_data_list$mwtab_matched_features == 0)
  mwtab_samples_zero_or_missing = is.na(
    mwtab_data_list$mwtab_matched_samples
  ) ||
    (mwtab_data_list$mwtab_matched_samples == 0)
  if (
    mwtab_features_zero_or_missing ||
      mwtab_samples_zero_or_missing
  ) {
    use_defaults = TRUE
  }

  if (!("se" %in% report_list$steps)) {
    use_defaults = TRUE
  }

  se_data_list = get_se_data_counts_status(se_data, use_defaults)

  se_samples_removed_list = get_samples_removed(se_data, use_defaults)
  n_samples_list = get_n_samples_ssf(se_data, use_defaults)

  qc_status_list = check_qc_samples(se_data, use_defaults)

  if (!is.null(metadata(se_data)$left_censorship)) {
    percent_missing_list = check_percent_missing(
      se_data,
      data_limits$missing,
      FALSE
    )
  } else {
    percent_missing_list = check_percent_missing(
      se_data,
      data_limits$missing,
      TRUE
    )
  }

  if (!("correlation" %in% report_list$steps)) {
    use_defaults = TRUE
  }
  median_correlations_list = check_median_correlations(
    se_data,
    data_limits$median_correlation,
    use_defaults
  )

  if (!("pca" %in% report_list$steps)) {
    use_defaults = TRUE
  }

  rsd_list = check_rsd_values(
    se_data,
    data_limits$rsd,
    use_defaults
  )

  pca_list = check_pca(se_data, data_limits$pca, use_defaults)
  score_list = get_scores(se_data, use_defaults)

  data_list = c(
    se_data_list,
    se_samples_removed_list,
    n_samples_list,
    qc_status_list,
    percent_missing_list,
    median_correlations_list,
    rsd_list,
    pca_list,
    score_list
  )
  return(c(mwtab_list, data_list, analysis_list))
}

get_samples_removed = function(se_data, default = FALSE) {
  if (default) {
    return(list(removed_samples = "None"))
  }
  removed_samples = se_data$removed
  if (sum(!is.na(removed_samples)) == 0) {
    return(list(removed_samples = "None"))
  }
  removed_reasons = removed_samples[!is.na(removed_samples)]
  split_reasons = split(removed_reasons, removed_reasons)
  n_reasons = purrr::map(split_reasons, length)
  reasons_string = purrr::imap(n_reasons, \(n_x, id_x) {
    paste0(id_x, ": ", n_x)
  }) |>
    paste0(collapse = "; ")
  return(list(removed_samples = reasons_string))
}

get_scores = function(se_data, default = FALSE) {
  if (default) {
    return(list(score_diff_summary = "NA", score_diff_status = "red"))
  }
  pca_data = metadata(se_data)$pca_groups
  last_run = check_qc_outputs(se_data)

  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  ssf_scores = pca_data_use$scores

  ssf_noall = ssf_scores |>
    dplyr::filter(!ssf %in% "all")

  ssf_scores_wide = ssf_noall |>
    dplyr::select(-differences_fraction) |>
    tidyr::pivot_wider(names_from = "method", values_from = "n_differences")

  ssf_scores_wide = ssf_scores_wide |>
    dplyr::mutate(
      rle_frac = rle / group_size,
      rle_status = dplyr::case_when(
        (group_size <= 5) & (rle > 2) ~ "red",
        (group_size <= 5) & (rle > 1) ~ "yellow",
        (group_size <= 5) & (rle == 0) ~ "green",
        (group_size <= 10) & (rle > 4) ~ "red",
        (group_size <= 10) & (rle > 1) ~ "yellow",
        (group_size <= 10) & (rle <= 1) ~ "green",
        (group_size > 10) & (rle_frac > 0.5) ~ "red",
        (group_size > 10) & (rle_frac > 0.25) ~ "yellow",
        (group_size > 10) & (rle_frac <= 0.25) ~ "green"
      ),
      ici_frac = `median-ici` / group_size,
      ici_status = dplyr::case_when(
        (group_size <= 5) & (`median-ici` > 1) ~ "red",
        (group_size <= 5) & (`median-ici` > 0) ~ "yellow",
        (group_size <= 5) & (`median-ici` == 0) ~ "green",
        (group_size <= 10) & (`median-ici` > 3) ~ "red",
        (group_size <= 10) & (`median-ici` >= 2) ~ "yellow",
        (group_size <= 10) & (`median-ici` <= 1) ~ "green",
        (group_size > 10) & (ici_frac > 0.3) ~ "red",
        (group_size > 10) & (ici_frac > 0.1) ~ "yellow",
        (group_size > 10) & (ici_frac <= 0.1) ~ "green"
      ),
    )

  ssf_summaries = ssf_scores_wide |>
    dplyr::summarise(
      ici_red = sum(ici_status %in% "red") / nrow(ssf_scores_wide),
      ici_yellow = sum(ici_status %in% "yellow") / nrow(ssf_scores_wide),
      ici_green = sum(ici_status %in% "green") / nrow(ssf_scores_wide),
      rle_red = sum(rle_status %in% "red") / nrow(ssf_scores_wide),
      rle_yellow = sum(rle_status %in% "yellow") / nrow(ssf_scores_wide),
      rle_green = sum(rle_status %in% "green") / nrow(ssf_scores_wide)
    )

  ssf_status = dplyr::case_when(
    ssf_summaries$ici_red > 0.5 ~ "red",
    ssf_summaries$ici_yellow > 0.5 ~ "yellow",
    ssf_summaries$ici_green > 0.75 & ssf_summaries$rle_yellow < 0.5 ~ "green",
    ssf_summaries$ici_green > 0.75 & ssf_summaries$rle_green > 0.5 ~ "green",
    ssf_summaries$ici_green > 0.75 & ssf_summaries$rle_red > 0.5 ~ "yellow",
  )
  list(
    score_diff_summary = paste0(
      "icikt: ",
      paste0(
        ssf_scores_wide$`median-ici`,
        "/",
        ssf_scores_wide$group_size,
        collapse = ";"
      ),
      "; rle: ",
      paste0(
        ssf_scores_wide$rle,
        "/",
        ssf_scores_wide$group_size,
        collapse = ";"
      )
    ),
    score_diff_status = ssf_status
  )
}

get_n_samples_ssf = function(se_data, default = FALSE) {
  if (default) {
    out_n = tibble::tibble(ssf = "NA", `N-Samples` = NA)
    return(list(
      n_samples_ssf = "NA",
      n_ssf_n = NA,
      n_sample_status = "red",
      n_ssf_status = "red"
    ))
  }
  factor_types = metadata(se_data)$factor_types
  se_data = se_data[, is.na(se_data$removed), drop = FALSE]

  if (!is.null(se_data$outlier)) {
    se_data = se_data[, !se_data$outlier, drop = FALSE]
  }
  if (!is.null(se_data$experiment_or_qc)) {
    se_data = se_data[,
      se_data$experiment_or_qc %in% "experiment",
      drop = FALSE
    ]
  }
  ssf_info = get_vardata(se_data)
  ssf_n = determine_experimental_design(ssf_info, factor_types)
  out_n = ssf_n$n_by_factor_factors

  n_has_5 = out_n$`N-Samples` >= 5
  n_has_3 = out_n$`N-Samples` >= 3

  n_sample_status = dplyr::case_when(
    !all(n_has_3) ~ "red",
    !all(n_has_5) && all(n_has_3) ~ "yellow",
    all(n_has_5) ~ "green"
  )

  ssf_metadata_status = dplyr::case_when(
    nrow(out_n) >= 20 ~ "red",
    dplyr::between(nrow(out_n), 10, 20) ~ "yellow",
    nrow(out_n) == 1 ~ "yellow",
    nrow(out_n) < 10 ~ "green"
  )

  out_list = list(
    n_samples_ssf = paste0(out_n$`N-Samples`, collapse = ":"),
    n_ssf_n = nrow(out_n),
    n_sample_status = n_sample_status,
    ssf_metadata_status = ssf_metadata_status
  )
  out_list
}

check_pca = function(se_data, pca_limits, default = FALSE) {
  # pca_limits = data_limits$pca
  if (default) {
    return(list(
      pca_pc_index_90 = NA,
      pca_pc_index_90_status = "red",
      pca_anova_pvalue = NA,
      pca_anova_pvalue_PC = "NA",
      pca_anova_pvalue_status = "red"
    ))
  }
  pca_data = metadata(se_data)$pca_groups
  last_run = check_qc_outputs(se_data)

  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  pca_variance = pca_data_use$pca$variance

  pca_pc_index_90 = which(pca_variance$cumulative > 0.9) |> min()
  pca_pc_index_90_status = dplyr::case_when(
    pca_pc_index_90 <= pca_limits[1] ~ "green",
    pca_pc_index_90 > pca_limits[1] &&
      pca_pc_index_90 <= pca_limits[2] ~ "yellow",
    pca_pc_index_90 > pca_limits[2] ~ "red"
  )
  pca_anova = pca_data_use$pca$anova |>
    dplyr::mutate(p.adjust = p.adjust(p.value, method = "bonferroni"))

  pca_anova_pvalue = pca_anova |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = 1)

  pca_anova_perc = pca_variance |>
    dplyr::filter(pc %in% pca_anova_pvalue$PC) |>
    dplyr::pull(percent)
  pca_anova_pvalue_status = dplyr::case_when(
    pca_anova_pvalue$p.adjust > 0.05 ~ "red",
    pca_anova_pvalue$p.adjust <= 0.05 &&
      (pca_anova_perc >= 0.05) ~ "green",
    pca_anova_pvalue$p.adjust <= 0.05 &&
      (pca_anova_perc >= 0.01 && pca_anova_perc < 0.05) ~ "yellow",
    pca_anova_pvalue$p.adjust <= 0.05 &&
      (pca_anova_perc < 0.01) ~ "red",
  )
  return(list(
    pca_pc_index_90 = pca_pc_index_90,
    pca_pc_index_90_status = pca_pc_index_90_status,
    pca_anova_pvalue = pca_anova_pvalue$p.adjust,
    pca_anova_pvalue_PC = pca_anova_pvalue$PC,
    pca_anova_pvalue_status = pca_anova_pvalue_status
  ))
}

check_rsd_values = function(se_data, rsd_limits, default = FALSE) {
  # rsd_limits = data_limits$rsd; default = FALSE
  if (default) {
    return(list(
      rsd = NA,
      rsd_status = "red"
    ))
  }
  pca_data = metadata(se_data)$pca_groups
  last_run = check_qc_outputs(se_data)

  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  if (is.null(pca_data_use$rsd$all_rsd)) {
    return(list(rsd = NA, rsd_status = "red"))
  }
  rsd_use = pca_data_use$rsd$all_rsd_median |>
    dplyr::filter(data %in% "r") |>
    dplyr::pull(median) |>
    max(na.rm = TRUE)

  list(
    rsd = rsd_use,
    rsd_status = dplyr::case_when(
      rsd_use <= rsd_limits[1] ~ "green",
      rsd_use < rsd_limits[1] && rsd_use <= rsd_limits[2] ~ "yellow",
      rsd_use > rsd_limits[2] ~ "red"
    )
  )
}

check_median_correlations = function(
  se_data,
  correlation_limits,
  default = FALSE
) {
  if (default) {
    return(list(
      max_median_correlation = NA,
      max_median_correlation_status = "red"
    ))
  }
  # correlation_limits = data_limits$median_correlation; default = FALSE
  cor_info = se_data$median_ici_correlation
  max_median_correlation = max(cor_info, na.rm = TRUE)
  max_median_correlation_status = dplyr::case_when(
    max_median_correlation >= correlation_limits[1] ~ "green",
    max_median_correlation < correlation_limits[1] &&
      max_median_correlation >= correlation_limits[2] ~ "yellow",
    max_median_correlation < correlation_limits[2] ~ "red"
  )
  list(
    max_median_correlation = max_median_correlation,
    max_median_correlation_status = max_median_correlation_status
  )
}

check_percent_missing = function(se_data, percent_limits, default = FALSE) {
  # percent_limits = data_limits$missing; default = FALSE
  if (!default) {
    lc_info = metadata(se_data)$left_censorship
  } else {
    percent_list = list(percent_missing = NA, percent_missing_status = "red")
    lc_info = NULL
  }

  if (is.null(lc_info)) {
    default = TRUE
  }
  if (default) {
    percent_list = list(percent_missing = NA, percent_missing_status = "red")
    lc_test_list = list(
      left_censorship_test = NA,
      left_censorship_test_status = "red"
    )
    return(c(percent_list, lc_test_list))
  }
  percent_missing = lc_info$perc_missing
  percent_missing_status = dplyr::case_when(
    percent_missing <= percent_limits[1] ~ "green",
    percent_limits[1] < percent_missing &&
      percent_missing <= percent_limits[2] ~ "yellow",
    percent_missing > percent_limits[2] ~ "red"
  )

  percent_list = list(
    percent_missing = percent_missing,
    percent_missing_status = percent_missing_status
  )
  lc_test = lc_info$lc_test$binomial_test
  if (is.null(lc_test)) {
    default = TRUE
  }
  if (default) {
    lc_test_list = list(
      left_censorship_test = NA,
      left_censorship_test_status = "red"
    )
    return(c(percent_list, lc_test_list))
  }
  lc_test_list = list(
    left_censorship_test = lc_test$p.value,
    left_censorship_test_status = dplyr::case_when(
      percent_missing > 10 && lc_test$p.value <= 0.05 ~ "green",
      percent_missing > 10 && lc_test$p.value > 0.05 ~ "yellow",
      TRUE ~ "green"
    )
  )
  c(percent_list, lc_test_list)
}


get_ssf_group_sizes = function(ssf_data) {
  # ssf_data = mwtab_info$SUBJECT_SAMPLE_FACTORS
  ssf = clean_ssf_factors(ssf_data)

  if (length(unique(ssf$subject)) == 1) {
    ssf$subject = NULL
    id_cols = "sample_id"
  } else {
    id_cols = c("subject", "sample_id")
  }
  factor_type = ssf |>
    dplyr::select(factor, type) |>
    dplyr::distinct()

  ssf_wide = widen_ssf_json(ssf, id_cols)
  use_for_factors = c(id_cols, grab_factor_types(factor_type))
  ssf_concat = create_ssf_groups_json(
    ssf_wide |> dplyr::select(tidyselect::all_of(use_for_factors)),
    id_cols
  )

  ssf_design = determine_experimental_design(ssf_concat, factor_type)
  ssf_design
}

get_mwtab_data_counts = function(se_data) {
  mwtab_info = metadata(se_data)$mwtab

  sample_filtering = mwtab_info$FILTERED_SAMPLES
  feature_filtering = mwtab_info$FILTERED_FEATURES

  if (!is.null(mwtab_info$SUBJECT_SAMPLE_FACTORS)) {
    ssf_n_size = try(get_ssf_group_sizes(mwtab_info$SUBJECT_SAMPLE_FACTORS))
    if (inherits(ssf_n_size, "try-error")) {
      mwtab_n_ssf = NA
      mwtab_n_ssf_samples = NA
    } else {
      mwtab_n_ssf = ssf_n_size$summary$n_groups
      mwtab_n_ssf_samples = max(ssf_n_size$n_by_factor[["N-Samples"]])
    }
  } else {
    mwtab_n_ssf = NA
    mwtab_n_ssf_samples = NA
  }

  out_list =
    list(
      mwtab_measurement_features = feature_filtering$original[2],
      mwtab_matched_features = feature_filtering$kept[2],
      mwtab_measurement_samples = sample_filtering$original[2],
      mwtab_matched_samples = sample_filtering$kept[2],
      mwtab_n_ssf = mwtab_n_ssf,
      mwtab_n_ssf_samples = mwtab_n_ssf_samples
    )
  out_list2 = purrr::map(out_list, \(x) {
    if (is.null(x)) {
      return(NA)
    } else {
      (return(x))
    }
  })
  out_list2
}

get_se_data_counts_status = function(se_data, default = FALSE) {
  if (default) {
    return(list(
      summarized_features = 0,
      summarized_features_status = "none",
      summarized_samples = 0,
      summarized_samples_status = "none"
    ))
  }

  dyn_range_filtering = metadata(se_data)$dynamic_range_filtering

  if (is.null(dyn_range_filtering)) {
    summarized_features = 0
    summarized_samples = 0
  } else {
    summarized_features = dyn_range_filtering$n_features[2]
    summarized_samples = dyn_range_filtering$n_samples[2]
  }

  n_samples_status = dplyr::case_when(
    summarized_samples <= 3 ~ "low",
    dplyr::between(summarized_samples, 4, 9) ~ "medium",
    summarized_samples >= 10 ~ "high"
  )

  n_features_status = dplyr::case_when(
    summarized_features <= 9 ~ "low",
    dplyr::between(summarized_features, 10, 29) ~ "medium",
    summarized_features >= 30 ~ "high"
  )
  list(
    summarized_features = summarized_features,
    summarized_features_status = n_features_status,
    summarized_samples = summarized_samples,
    summarized_samples_status = n_samples_status
  )
}


check_qc_samples = function(se_data, default = FALSE) {
  if (default) {
    out_list = list(has_qc_samples = "no", has_qc_samples_status = "red")
    return(out_list)
  }
  has_qc = if (sum(se_data$experiment_or_qc %in% "qc") > 0) {
    "yes"
  } else {
    "no"
  }
  list(
    has_qc_samples = has_qc,
    has_qc_samples_status = if (has_qc %in% "yes") {
      "green"
    } else {
      "yellow"
    }
  )
}

get_mwtab_project_info = function(mwtab_info) {
  mw_fields = tibble::tibble(
    field = c("PROJECT_ID", "STUDY_ID", "ANALYSIS_ID", "CREATED_ON")
  )
  project_fields = tibble::tibble(
    field = c("PROJECT_TITLE", "INSTITUTE", "LABORATORY", "DOI")
  )
  mw_list = dplyr::left_join(
    mw_fields,
    mwtab_info$`METABOLOMICS WORKBENCH`,
    by = "field"
  ) |>
    tibble_to_list()
  project_list = dplyr::left_join(
    project_fields,
    mwtab_info$PROJECT,
    by = "field"
  ) |>
    tibble_to_list()
  c(mw_list, project_list)
}

get_mwtab_instrument_chromatography_info = function(mwtab_info) {
  nmr_ms_info_fields = tibble::tibble(
    field = c(
      "INSTRUMENT_NAME",
      "INSTRUMENT_TYPE",
      "MS_TYPE",
      "ION_MODE",
      "NMR_EXPERIMENT_TYPE",
      "SPECTROMETER_FREQUENCY"
    )
  )
  if (!is.null(mwtab_info$MS)) {
    use_info = mwtab_info$MS
  }
  if (!is.null(mwtab_info$NM)) {
    use_info = mwtab_info$NM
  }

  if (!exists("use_info")) {
    return(NULL)
  }
  instr_info = dplyr::left_join(nmr_ms_info_fields, use_info, by = "field")

  instr_list = tibble_to_list(instr_info)

  chrom_fields = tibble::tibble(
    field = c("CHROMATOGRAPHY_TYPE", "INSTRUMENT_NAME", "COLUMN_NAME"),
    use_field = c(
      "chromatography_type",
      "chromatography_instrument",
      "column_name"
    )
  )
  chrom_data = mwtab_info$CHROMATOGRAPHY
  if (is.null(chrom_data) || nrow(chrom_data) == 0) {
    chrom_info = chrom_fields
    chrom_info$value = "NA"
  } else {
    chrom_info = dplyr::left_join(
      chrom_fields,
      mwtab_info$CHROMATOGRAPHY,
      by = "field"
    )
  }

  chrom_list = tibble_to_list(chrom_info, "use_field")

  out_list = c(instr_list, chrom_list)
  out_list
}

get_mwtab_organism_sample = function(mwtab_info) {
  subject_fields = tibble::tibble(
    field = c(
      "SUBJECT_TYPE",
      "SUBJECT_SPECIES",
      "TAXONOMY_ID",
      "SPECIES_GROUP"
    )
  )
  subject_info = dplyr::left_join(
    subject_fields,
    mwtab_info$SUBJECT,
    by = "field"
  )

  subject_list = tibble_to_list(subject_info)

  sampleprep_fields = tibble::tibble(
    field = c("ORGAN", "ORGAN_SPECIFICATION")
  )
  sampleprep_info = dplyr::left_join(
    sampleprep_fields,
    mwtab_info$SAMPLEPREP,
    by = "field"
  )
  sampleprep_list = tibble_to_list(sampleprep_info)

  collection_fields = tibble::tibble(
    field = c("SAMPLE_TYPE", "COLLECTION_METHOD", "COLLECTION_LOCATION")
  )
  collection_info = dplyr::left_join(
    collection_fields,
    mwtab_info$COLLECTION,
    by = "field"
  )
  collection_list = tibble_to_list(collection_info)

  c(subject_list, sampleprep_list, collection_list)
}

tibble_to_list = function(in_tbl, name_field = "field") {
  if (is.null(in_tbl$value)) {
    out_list = rep("NA", nrow(in_tbl)) |> as.list()
  } else {
    out_list = as.list(in_tbl$value)
  }

  names(out_list) = tolower(in_tbl[[name_field]])
  out_list
}


sets_percent_missing = function(se_pca) {
  # se_pca = tar_read(se_pca_g01); tar_source("R");
  se_missing = purrr::map(se_pca, get_percent_missing) |>
    purrr::list_rbind()
  se_missing
}

get_percent_missing = function(se_data) {
  id = metadata(se_data)$qc$id
  lc_info = metadata(se_data)$left_censorship
  if (!is.null(lc_info$lc_test)) {
    p_value = lc_info$lc_test$binomial_test$p.value
  } else {
    p_value = NA
  }
  if (is.null(lc_info)) {
    return(tibble::tibble(id = id, missing = NA, p_value = p_value))
  } else {
    return(tibble::tibble(
      id = id,
      missing = lc_info$perc_missing,
      p_value = p_value
    ))
  }
}

sets_median_correlations = function(se_pca) {
  se_median_cor = purrr::map(se_pca, get_median_correlations) |>
    purrr::list_rbind()
  se_median_cor
}

get_median_correlations = function(se_data) {
  id = metadata(se_data)$qc$id
  if (!is.null(se_data[["median_ici_correlation"]])) {
    get_vardata(se_data) |>
      dplyr::select(sample_id, median_ici_correlation) |>
      dplyr::mutate(id = id)
  } else {
    return(NULL)
  }
}

sets_rsd = function(se_pca) {
  se_rsd = purrr::map(se_pca, get_rsd) |>
    purrr::list_rbind()
  se_rsd
}

get_rsd = function(se_data) {
  id = metadata(se_data)$qc$id

  pca_groups = metadata(se_data)$pca_groups
  if (is.null(pca_groups)) {
    return(NULL)
  }
  use_rsd = pca_groups$all$rsd$all_rsd_median
  if (is.null(use_rsd)) {
    return(NULL)
  }
  out_rsd = use_rsd |>
    dplyr::filter(data %in% "r") |>
    dplyr::mutate(id = id)
  out_rsd
}

sets_get_pc = function(se_pca) {
  sets_pc = purrr::map(se_pca, get_pc_90) |>
    purrr::list_rbind()
  sets_pc
}

get_pc_90 = function(se_data) {
  id = metadata(se_data)$qc$id

  pca_groups = metadata(se_data)$pca_groups

  if (is.null(pca_groups)) {
    return(NULL)
  }
  last_id = check_qc_outputs(se_data)
  pc_index = last_id$pca[length(last_id$pca)]
  use_pca = pca_groups[[pc_index]]$pca$variance
  pc_90 = min(which(use_pca$cumulative > 0.9))
  tibble::tibble(
    analysis_id = id,
    pc_90 = pc_90,
    n_pc = nrow(use_pca),
    fraction_90 = pc_90 / nrow(use_pca)
  )
}


derive_data_limits = function(
  sets_percent_missing_result,
  sets_median_correlation_result,
  sets_rsd_result,
  sets_pc_number_result
) {
  # tar_load(c(sets_percent_missing_result, sets_median_correlation_result, sets_rsd_result))

  use_percent_missing = sets_percent_missing_result |>
    dplyr::filter(!is.na(missing), missing > 0)
  p_missing_sd = sd(use_percent_missing$missing)
  p_missing_limits = c(2, 4) * p_missing_sd

  med_cor_summary = sets_median_correlation_result |>
    dplyr::summarise(
      max = max(median_ici_correlation, na.rm = TRUE),
      .by = id
    ) |>
    dplyr::filter(!is.infinite(max)) |>
    dplyr::summarise(median = median(max), sd = sd(max))
  med_cor_limits = med_cor_summary$median - (c(2, 4) * med_cor_summary$sd)

  rsd_summary = sets_rsd_result |>
    dplyr::summarise(
      median_value = median(median, na.rm = TRUE),
      sd = sd(median, na.rm = TRUE)
    )
  rsd_limits = rsd_summary$median_value + (c(2, 4) * rsd_summary$sd)

  use_pc_number = sets_pc_number_result |>
    dplyr::mutate(pc_log = log10(pc_90)) |>
    dplyr::filter(!is.infinite(pc_90), pc_log >= 1)

  pc_summary = use_pc_number |>
    dplyr::summarise(median_value = 1, sd = sd(pc_log))

  pc_limits = 10^(pc_summary$median_value + (c(2, 4) * pc_summary$sd)) |>
    floor()
  list(
    missing = p_missing_limits,
    median_correlation = med_cor_limits,
    rsd = rsd_limits,
    pca = pc_limits
  )
}

create_review_list = function(sets_qc_evaluation, mwtab_datasets) {
  # tar_load(sets_qc_evaluation)
  interesting_failures = list(failed_steps = c("create_se", "correlation"))

  review_data = purrr::imap(interesting_failures, \(x, icol) {
    message(x)
    message(icol)
    tmp_column = sets_qc_evaluation[[icol]]
    out_match = purrr::map(x, \(in_step) {
      sets_qc_evaluation[["internal_id"]][tmp_column %in% in_step]
    }) |>
      purrr::list_c() |>
      unique()
  }) |>
    purrr::list_c() |>
    unique()
  out_qc = sets_qc_evaluation |>
    dplyr::filter(internal_id %in% review_data) |>
    dplyr::select(
      internal_id,
      analysis_id,
      study_id,
      project_id,
      failed_steps,
      failed_why
    ) |>
    dplyr::arrange(failed_steps, failed_why)
  mwtab_datasets$json = fs::path_file(mwtab_datasets$file)
  out_qc = dplyr::left_join(
    out_qc,
    mwtab_datasets[, c("id", "json")],
    by = c("internal_id" = "id")
  )
  write.table(
    out_qc,
    file = "outputs/weird_issues_2026-07-09.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}

summarize_qc_results_table = function(sets_qc_evaluation, group_none = FALSE) {
  n_qc = nrow(sets_qc_evaluation)
  qc_summary = sets_qc_evaluation |>
    dplyr::select(completed_steps, failed_steps, failed_why) |>
    dplyr::summarize(
      n = dplyr::n(),
      perc = n / n_qc,
      perc_display = paste0(format(perc * 100, digits = 2), "%"),
      .by = c(completed_steps, failed_steps, failed_why)
    ) |>
    dplyr::arrange(completed_steps, failed_steps, failed_why)

  qc_none = sets_qc_evaluation |>
    dplyr::filter(failed_steps %in% "create_se") |>
    dplyr::select(completed_steps, failed_steps) |>
    dplyr::summarise(
      n = dplyr::n(),
      perc = n / n_qc,
      perc_display = paste0(format(perc * 100, digits = 2), "%"),
      .by = c(completed_steps)
    )

  dplyr::bind_rows(qc_summary, qc_none)
}

find_green_analyses = function(sets_qc_evaluation) {
  status_cols = grep("status$", names(sets_qc_evaluation), value = TRUE)
  status_cols = grep(
    "^summarized|^has_qc",
    status_cols,
    invert = TRUE,
    value = TRUE
  )

  green_status = yellow_status = rep(TRUE, nrow(sets_qc_evaluation))

  for (istatus in status_cols) {
    green_status = green_status &
      (sets_qc_evaluation[[istatus]] %in% c("green"))
    yellow_status = yellow_status &
      (sets_qc_evaluation[[istatus]] %in% c("yellow", "green"))
  }
  overall_status = rep("red", nrow(sets_qc_evaluation))
  overall_status[yellow_status] = "yellow"
  overall_status[green_status] = "green"
  sets_qc_evaluation$overall_status = overall_status
  sets_qc_evaluation
}

summarize_status = function(overall_status) {
  n_qc = nrow(overall_status)
  color_summary = overall_status |>
    dplyr::select(internal_id, overall_status) |>
    dplyr::summarize(
      n = dplyr::n(),
      perc = n / n_qc,
      perc_display = paste0(format(perc * 100, digits = 2), "%"),
      .by = c(overall_status)
    ) |>
    dplyr::arrange(overall_status)
}

format_in_sentence = function(n_tibble) {
  n_tibble |>
    glue::glue_data("{n} ({perc_display})")
}

calculate_mode = function(values) {
  density_estimate = stats::density(values)
  mode_value = density_estimate$x[which.max(density_estimate$y)]
  mode_value
}
