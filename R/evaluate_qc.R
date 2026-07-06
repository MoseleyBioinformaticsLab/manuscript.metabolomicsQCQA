sets_evaluate_qc = function(sets_se) {
  # sets_se = tar_read(se_pca_g01); tar_source("R")
  # sets_se = tar_read(se_pca_g02); tar_source("R")
  # sets_se = tar_read(se_pca_g20); se_data = sets_se[[194]]; tar_source("R")
  sets_qc = purrr::map(sets_se, \(x) {
    evaluate_qc(x) |> tibble::as_tibble()
  }) |>
    purrr::list_rbind()
  sets_qc
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

get_se_data_counts = function(se_data) {
  dyn_range_filtering = metadata(se_data)$dynamic_range_filtering

  if (is.null(dyn_range_filtering)) {
    out_list = list(
      summarized_features = 0,
      summarized_samples = 0
    )
  } else {
    out_list = list(
      summarized_features = dyn_range_filtering$n_features[2],
      summarized_samples = dyn_range_filtering$n_samples[2]
    )
  }
  out_list
}

evaluate_qc = function(se_data) {
  # get all the metabolomics workbench info out first
  # the goal here is to have analysis attributes
  # that are searchable and filterable by consumers
  # of the qc-qa analysis

  # se_data = sets_se[[92]]
  # se_data = sets_se[[123]]
  # se_data = sets_se[[159]]

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

  c(mwtab_list, mwtab_data_list)
  se_data_list = get_se_data_counts(se_data)

  # has_qc = if (sum(se_data$experiment_or_qc %in% "qc") > 0) {
  #   "yes"
  # } else {
  #   "no"
  # }
  # # how many experimental samples in each SSF after removing
  # # qc and outliers?

  # f_types = get_specific_metadata(se_data, "factor_types")

  # # ANOVA p-value after removing qc and outliers?
  # if (last_step$succeeded %in% "pca") {
  #   last_pca = report_list$pca[length(report_list$pca)]
  #   pca_res = metadata(se_data)$pca_groups[[last_pca]]
  #   best_p = pca_res$pca$anova |>
  #     dplyr::mutate(p.adjust = p.adjust(p.value, method = "bonferroni")) |>
  #     dplyr::arrange(p.adjust) |>
  #     dplyr::slice_head(n = 1) |>
  #     dplyr::pull(p.adjust)
  #   variances = pca_res$pca$variance
  #   n_pc = nrow(variances)
  #   if (n_pc < 5) {
  #     slice_pc = n_pc
  #   } else {
  #     slice_pc = 5
  #   }
  #   var_first_5 = variances |>
  #     dplyr::slice_head(n = slice_pc) |>
  #     dplyr::slice_tail(n = 1) |>
  #     dplyr::pull(cumulative)
  # }
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
