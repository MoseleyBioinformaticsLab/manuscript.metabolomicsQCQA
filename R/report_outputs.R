sets_create_qcqa_reports = function(
  pca_sets,
  save_loc = "qcqa_reports",
  report_file = "resources/example_report.qmd"
) {
  # pca_sets = tar_read(pca_g01); tar_source("R"); save_loc = "qcqa_reports"; report_file = "resources/example_report.qmd"
  # pca_sets = tar_read(pca_g13); tar_source("R"); save_loc = "qcqa_reports";   report_file = "resources/example_report.qmd"

  sets_tmp_files = purrr::map(names(pca_sets), \(in_id) {
    paste0(tempfile(paste0(in_id, "_")), c(".rds", ".qmd"))
  })
  names(sets_tmp_files) = names(pca_sets)
  out_reports = furrr::future_map2(
    pca_sets,
    sets_tmp_files,
    create_report_from_example,
    save_loc,
    report_file
  )
  out_reports
}
create_report_from_example = function(
  in_pca,
  in_tmp_files,
  save_loc = "qcqa_reports",
  report_file = "resources/example_report.qmd"
) {
  analysis_id = in_pca$id
  if (is.null(in_pca[["pca_groups"]])) {
    return(list(id = analysis_id, file = NULL, md5 = NULL))
  }
  rds_line = 63
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
  if (length(study_template_out) == 0) {
    study_template_out = "No study information found."
  }
  report_text[study_link_line] = study_template_out

  saveRDS(in_pca, file = in_tmp_files[1])
  read_files = glue::glue('use_pca = readRDS("{in_tmp_files[1]}")')
  report_text[rds_line] = read_files
  cat(report_text, file = in_tmp_files[2], sep = "\n")
  final_html = fs::path(save_loc, glue::glue("{in_pca$id}.html"))
  q_res = quarto::quarto_render(
    in_tmp_files[2],
    execute_dir = here::here(),
    quiet = TRUE
  )
  tmp_qmd = in_tmp_files[2]
  tmp_html = fs::path_ext_set(tmp_qmd, "html")

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
    raw_medians_rsd = raw_medians_rsd,
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

create_grouping_scores_table = function(
  score_table,
  score_method = "median-ici"
) {
  score_table |>
    dplyr::filter(method %in% score_method) |>
    dplyr::select(-method) |>
    dplyr::select(factors, differences_fraction, n_differences, group_size) |>
    gt::gt() |>
    gt::fmt_number(columns = differences_fraction, decimals = 2)
}

create_variance_values = function(variance_table, ind = "PC2", group = "PC5") {
  # variance_table = use_pca$pca_groups$all$pca$variance; ind = "PC2"; group = "PC5"
  ind_indices = seq_len(which(variance_table$pc %in% ind))
  ind_values = variance_table$percent[ind_indices]
  names(ind_values) = variance_table$pc[ind_indices]

  ind_values = c(
    ind_values,
    "cumulative" = variance_table$cumulative[ind_indices[length(ind_indices)]]
  )

  ind_print = paste0(format(ind_values * 100, digits = 2), "%")

  plot_result = glue::glue(
    "PCA sample scores plot, with samples colored by combination of factors. Variance explained by PC1: {ind_print[1]}; PC2: {ind_print[2]}; combined: {ind_print[3]}."
  )

  if (!(group %in% variance_table$pc)) {
    group = variance_table$pc[nrow(variance_table)]
  }
  group_value = c(
    group = variance_table |>
      dplyr::filter(pc %in% group) |>
      dplyr::pull(cumulative)
  )
  n_group = which(variance_table$pc %in% group)
  group_print = paste0(format(group_value * 100, digits = 2), "%")
  group_result = glue::glue(
    "PC-variable ANOVA results for PCs 1 - {n_group}, and p-value <= 0.05, arranged by PC and p-value. Cumulative variance accounted for is {group_print}."
  )
  names(group_value) = group

  list(plot = plot_result, table = group_result)
}

add_total_samples_by_factor = function(n_by_factor_factors) {
  # n_by_factor = use_pca$pca_groups$all$design$n_by_factor
  tot = n_by_factor_factors |>
    dplyr::summarise(
      factors = "All",
      `N-Samples` = sum(`N-Samples`)
    )
  dplyr::bind_rows(n_by_factor_factors, tot)
}
