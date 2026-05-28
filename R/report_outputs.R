sets_create_qcqa_reports = function(
  pca_sets,
  save_loc = "qcqa_reports",
  report_file = "R/example_report.qmd"
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
  rds_line = 53
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
    quiet = TRUE
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
