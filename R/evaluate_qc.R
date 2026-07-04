evaluate_qc = function(se_data) {
  # get all the metabolomics workbench info out first
  # the goal here is to have analysis attributes
  # that are searchable and filterable by consumers
  # of the qc-qa analysis
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
  project_list = get_mwtab_project_info(mwtab_info)
  instrument_list = get_mwtab_instrument_chromatography_info(mwtab_info)
  organism_list = get_mwtab_organism_sample(mwtab_info)

  mwtab_list = c(measurement_list, project_list, instrument_list, organism_list)

  # next is to start looking for red flags in the data
  # using the various criteria, and fitting them into the
  # same list / data.frame format that allows us to easily
  # filter / group / exclude / include different analyses
  report_list = check_qc_outputs(se_data)
  last_step = get_qc(se_data)
  if (nrow(se_data) == 0 || ncol(se_data) == 0) {
    data_list = list(
      n_initial_samples = 0,
      n_initial_metabolites = 0,
      n_initial_factor_groups = 0,
      qc_samples = "no",
      why = last_step$why
    )
  } else {}

  has_qc = if (sum(se_data$experiment_or_qc %in% "qc") > 0) {
    "yes"
  } else {
    "no"
  }
  # how many experimental samples in each SSF after removing
  # qc and outliers?

  f_types = get_specific_metadata(se_data, "factor_types")

  # ANOVA p-value after removing qc and outliers?
  if (last_step$succeeded %in% "pca") {
    last_pca = report_list$pca[length(report_list$pca)]
    pca_res = metadata(se_data)$pca_groups[[last_pca]]
    best_p = pca_res$pca$anova |>
      dplyr::mutate(p.adjust = p.adjust(p.value, method = "bonferroni")) |>
      dplyr::arrange(p.adjust) |>
      dplyr::slice_head(n = 1) |>
      dplyr::pull(p.adjust)
    variances = pca_res$pca$variance
    n_pc = nrow(variances)
    if (n_pc < 5) {
      slice_pc = n_pc
    } else {
      slice_pc = 5
    }
    var_first_5 = variances |>
      dplyr::slice_head(n = slice_pc) |>
      dplyr::slice_tail(n = 1) |>
      dplyr::pull(cumulative)
  }
}


get_mwtab_project_info = function(mwtab_info) {
  mw_fields = tibble::tibble(field = c("STUDY_ID", "ANALYSIS_ID", "CREATED_ON"))
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
  if (nrow(chrom_data) == 0) {
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
  out_list = as.list(in_tbl$value)
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
