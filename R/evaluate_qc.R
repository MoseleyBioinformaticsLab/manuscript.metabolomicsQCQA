evaluate_qc = function(se_data) {
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
  mwtab_list = c(mwtab_list, instrument_list)
}

get_mwtab_project_info = function(mwtab_info) {
  mw_fields = tibble::tibble(field = c("STUDY_ID", "ANALYSIS_ID", "CREATED"))
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
