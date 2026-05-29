sets_create_qcqa_reports = function(
  pca_sets,
  save_loc = "qcqa_reports",
  report_file = "resources/example_report.qmd"
) {
  # qcqa_plots_sets = tar_read(qcqa_plots_g01); pca_sets = tar_read(pca_g01)
  out_reports = purrr::map(
    pca_sets,
    create_report_from_example,
    save_loc,
    report_file
  )
  out_reports
}
create_report_from_example = function(in_pca, save_loc, report_file) {
  analysis_id = in_pca$id
  if (is.null(in_pca[["pca_groups"]])) {
    return(list(id = analysis_id, file = NULL, md5 = NULL))
  }
  rds_line = 59
  study_link_line = 17

  report_text = suppressWarnings(readLines(report_file))
  title_text = glue::glue("title: QC/QA Report for MW {analysis_id}")
  report_text[2] = title_text

  study_template_text = report_text[study_link_line]
  study_info = get_study_id_link(in_pca$mwtab[["METABOLOMICS WORKBENCH"]])
  study_id = study_info$study
  study_link = study_info$study_link
  analysis_link = study_info$analysis_link

  study_template_out = glue::glue(study_template_text)
  report_text[study_link_line] = study_template_out

  tmp_files = paste0(tempfile(pattern = in_pca$id), c(".rds", ".qmd"))
  saveRDS(in_pca, file = tmp_files[1])
  read_files = glue::glue('use_pca = readRDS("{tmp_files[1]}")')
  report_text[rds_line] = read_files
  cat(report_text, file = tmp_files[2], sep = "\n")
  final_html = fs::path(save_loc, glue::glue("{in_pca$id}.html"))
  q_res = quarto::quarto_render(
    tmp_files[2],
    execute_dir = here::here(),
    quiet = FALSE
  )
  tmp_qmd = tmp_files[2]
  tmp_html = tmp_qmd
  fs::path_ext(tmp_html) = "html"

  fs::file_copy(tmp_html, final_html, overwrite = TRUE)
  return(list(
    id = analysis_id,
    file = final_html,
    md5 = digest::digest(final_html, algo = "md5")
  ))
}

get_study_id_link = function(mwtab_data) {
  #mwtab_data = in_pca$mwtab[["METABOLOMICS WORKBENCH"]]
  study_id = mwtab_data |>
    dplyr::filter(field %in% "STUDY_ID") |>
    dplyr::pull(value)
  analysis_id = mwtab_data |>
    dplyr::filter(field %in% "ANALYSIS_ID") |>
    dplyr::pull(value)
  study_link = glue::glue(
    "https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={study_id}"
  )
  analysis_link = glue::glue(
    "https://metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID={study_id}&ANALYSIS_ID={analysis_id}"
  )
  return(list(
    study = study_id,
    analysis_link = analysis_link,
    study_link = study_link
  ))
}

create_pca_anova_table = function(anova_results) {
  tmp_anova = anova_results |>
    dplyr::filter(!is.na(df))
  use_pcs = c("PC1", "PC2", "PC3", "PC4", "PC5")
  tmp_anova |>
    dplyr::select(PC, variable, p.value, df, sumsq, meansq, statistic) |>
    dplyr::filter(PC %in% use_pcs, p.value <= 0.05) |>
    dplyr::arrange(PC, p.value) |>
    gt::gt() |>
    gt::fmt_integer(columns = df) |>
    gt::fmt_number(columns = c(sumsq, meansq, statistic), decimals = 1) |>
    gt::fmt_engineering(columns = p.value, decimals = 2)
}

calculate_variance_stuff = function(in_abundance) {
  sample_metadata = get_vardata(in_abundance)
  all_factors = unique(sample_metadata$factors)

  raw_abundance = assays(in_abundance)$raw
  norm_abundance = assays(in_abundance)$normalized

  raw_rsd = calculate_rsd_matrix(raw_abundance, sample_metadata, all_factors) |>
    dplyr::filter(n >= 3) |>
    dplyr::mutate(data = "r")
  norm_rsd = calculate_rsd_matrix(
    norm_abundance,
    sample_metadata,
    all_factors
  ) |>
    dplyr::filter(n >= 3) |>
    dplyr::mutate(data = "n")

  all_rsd = dplyr::bind_rows(raw_rsd, norm_rsd) |>
    dplyr::mutate(factors_data = paste0(factors, "-", data))
  all_rsd_median = all_rsd |>
    dplyr::summarise(median = median(rsd, na.rm = TRUE), .by = c(factors, data))
  raw_medians = matrixStats::colMedians(raw_abundance, na.rm = TRUE)
  raw_medians_rsd = sd(raw_medians) / mean(raw_medians)
  norm_medians = matrixStats::colMedians(norm_abundance, na.rm = TRUE)
  norm_medians_rsd = sd(norm_medians) / mean(norm_medians)

  list(
    all_rsd = all_rsd,
    all_rsd_median = all_rsd_median,
    raw_medians_rsd = raw_median_rsd,
    norm_medians_rsd = norm_medians_rsd
  )
}

setup_missing_matrix = function(data_matrix, global_na) {
  exclude_loc = matrix(
    FALSE,
    nrow = nrow(data_matrix),
    ncol = ncol(data_matrix)
  )
  rownames(exclude_loc) = rownames(data_matrix)
  colnames(exclude_loc) = colnames(data_matrix)
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
  }
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }

  return(exclude_loc)
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
