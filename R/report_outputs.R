sets_create_qcqa_reports = function(
  pca_sets,
  save_loc = "qcqa_reports"
) {
  # pca_sets = tar_read(se_pca_g01); tar_source("R"); save_loc = "qcqa_reports"; report_file = "resources/example_report.qmd"; quiet = FALSE

  sets_tmp_files = purrr::map(names(pca_sets), \(in_id) {
    paste0(tempfile(paste0(in_id, "_")), c(".rds", ".qmd"))
  })
  names(sets_tmp_files) = names(pca_sets)
  out_reports = furrr::future_map2(
    pca_sets,
    sets_tmp_files,
    create_report_from_templates,
    save_loc,
    quiet = TRUE
  )
  out_reports
}
create_report_from_templates = function(
  se_data,
  in_tmp_files,
  save_loc = "qcqa_reports",
  quiet = TRUE
) {
  # se_data = pca_sets[[1]]; in_tmp_files = sets_tmp_files[[1]]; quiet = FALSE; tar_source("R")
  analysis_id = get_specific_metadata(se_data, "qc")$id
  # message(analysis_id)
  mwtab_info = get_specific_metadata(se_data, "mwtab")
  report_text = suppressWarnings(readLines("resources/has_mwtab.qmd"))
  title_text = glue::glue("title: QC/QA Report for MW {analysis_id}")
  report_text[2] = title_text

  study_link_line = grep("This is a quality-control", report_text)[1]
  rds_line = grep("# se_data = readRDS", report_text)[1]

  study_template_text = report_text[study_link_line]
  study_info = get_study_id_link(mwtab_info[["METABOLOMICS WORKBENCH"]])
  study_id = study_info$study
  study_link = study_info$study_link
  analysis_link = study_info$analysis_link

  study_template_out = glue::glue(study_template_text)
  if (length(study_template_out) == 0) {
    study_template_out = "No study information found."
  }
  report_text[study_link_line] = study_template_out

  saveRDS(se_data, file = in_tmp_files[1])
  read_files = glue::glue('se_data = readRDS("{in_tmp_files[1]}")')
  report_text[rds_line] = read_files

  # the way to deal with the conditional information, is not to try and have
  # conditional includes.
  # its to do a wee little bit of computation out here, and then based on
  # that, read the correct include file here, and append it to
  # the basic file.
  # I don't like it, but I think it's the way to go for this use case.
  report_list = check_qc_outputs(se_data)

  if (report_list$steps[1] %in% "none") {
    dynamic_range_report = "resources/has_no_se.qmd"
  } else {
    dynamic_range_report = "resources/has_dynamic_range.qmd"
  }

  dyn_text = suppressWarnings(readLines(dynamic_range_report))

  all_text = c(report_text, dyn_text)

  n_kept = sum(is.na(se_data$removed))

  if (n_kept > 0) {
    boxplot_text = suppressWarnings(readLines(
      "resources/has_samples_boxplots.qmd"
    ))
    all_text = c(all_text, boxplot_text)
  }

  if ("pca" %in% report_list$steps) {
    pca_files = c(
      "all" = "resources/has_pca_all.qmd",
      "not_outlier" = "resources/has_pca_nooutliers.qmd",
      "not_qc" = "resources/has_pca_noqc.qmd",
      "not_outlier_qc" = "resources/has_pca_nooutliers_noqc.qmd"
    )
    keep_pca = pca_files[report_list$pca]

    pca_text = suppressWarnings(readLines(pca_files[1]))

    if (length(keep_pca) > 1) {
      for (ipca in seq(2, length(keep_pca))) {
        pca_text = c(
          pca_text,
          "\n\n",
          suppressWarnings(readLines(pca_files[ipca]))
        )
      }
    }
    all_text = c(all_text, pca_text)
  }

  cat(all_text, file = in_tmp_files[2], sep = "\n")
  final_html = fs::path(save_loc, glue::glue("{analysis_id}.html"))
  q_res = try(quarto::quarto_render(
    in_tmp_files[2],
    execute_dir = here::here(),
    quiet = quiet
  ))

  if (inherits(q_res, "try-error")) {
    return(tibble::tibble(
      id = analysis_id,
      file = "none",
      md5 = "quarto failure"
    ))
  }

  tmp_qmd = in_tmp_files[2]
  tmp_html = fs::path_ext_set(tmp_qmd, "html")

  fs::file_copy(tmp_html, final_html, overwrite = TRUE)
  unlink(in_tmp_files[1])
  unlink(in_tmp_files[2])
  unlink(tmp_html)
  return(tibble::tibble(
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
    dplyr::filter(!is.na(df)) |>
    dplyr::mutate(p.adjust = p.adjust(p.value, method = "bonferroni"))
  use_pcs = c("PC1", "PC2", "PC3", "PC4", "PC5")
  tmp_anova |>
    dplyr::select(
      PC,
      variable,
      p.adjust,
      p.value,
      df,
      sumsq,
      meansq,
      statistic
    ) |>
    dplyr::filter(PC %in% use_pcs) |>
    dplyr::arrange(PC, p.value) |>
    gt::gt() |>
    gt::fmt_integer(columns = df) |>
    gt::fmt_number(columns = c(sumsq, meansq, statistic), decimals = 1) |>
    gt::fmt_engineering(columns = c(p.adjust, p.value), decimals = 2) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        rows = p.adjust <= 0.05
      )
    )
}

calculate_variance_stuff = function(in_abundance) {
  sample_metadata = get_vardata(in_abundance)
  all_factors = unique(sample_metadata$ssf)

  raw_abundance = assays(in_abundance)$raw
  norm_abundance = assays(in_abundance)$normalized

  all_rsd = NULL
  all_rsd_median = NULL
  all_rsd_median_diff = NULL
  raw_medians_rsd = NULL
  norm_medians_rsd = NULL

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

  if ((nrow(raw_rsd) >= 5) && (nrow(raw_rsd) == nrow(norm_rsd))) {
    all_rsd = dplyr::bind_rows(raw_rsd, norm_rsd) |>
      dplyr::mutate(ssf_data = paste0(ssf, "-", data))
    all_rsd_wide_rn = all_rsd |>
      tidyr::pivot_wider(
        id_cols = c(feature_id, ssf),
        names_from = data,
        values_from = rsd
      ) |>
      dplyr::mutate(r_n = r - n)
    all_rsd_median = all_rsd |>
      dplyr::summarise(median = median(rsd, na.rm = TRUE), .by = c(ssf, data))
    all_rsd_median_diff = all_rsd_wide_rn |>
      dplyr::summarise(median_rn = median(r_n, na.rm = TRUE), .by = c(ssf))
    raw_medians = matrixStats::colMedians(raw_abundance, na.rm = TRUE)
    raw_medians_rsd = sd(raw_medians) / mean(raw_medians)
    norm_medians = matrixStats::colMedians(norm_abundance, na.rm = TRUE)
    norm_medians_rsd = sd(norm_medians) / mean(norm_medians)
  }

  list(
    all_rsd = all_rsd,
    all_rsd_median = all_rsd_median,
    all_rsd_median_diff = all_rsd_median_diff,
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


create_grouping_scores_table = function(
  score_table,
  score_method = "median-ici"
) {
  score_table |>
    dplyr::filter(method %in% score_method) |>
    dplyr::select(-method) |>
    dplyr::select(ssf, differences_fraction, n_differences, group_size) |>
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
    "PC-variable ANOVA results for PCs 1 - {n_group}, arranged by PC and p-value. Adjusted p-values are bolded. Cumulative variance accounted for is {group_print}."
  )
  names(group_value) = group

  list(plot = plot_result, table = group_result)
}

add_total_samples_by_factor = function(n_by_factor_factors) {
  # n_by_factor = use_pca$pca_groups$all$design$n_by_factor
  tot = n_by_factor_factors |>
    dplyr::summarise(
      ssf = "All",
      `N-Samples` = sum(`N-Samples`)
    )
  dplyr::bind_rows(n_by_factor_factors, tot)
}

check_for_outliers_qc = function(in_abundance) {
  # use_pca = pca_sets[[1]]; in_abundance = use_pca$abundance; tar_source("R")
  sample_info = in_abundance |> get_vardata()

  out_todo = list(
    outliers = any(sample_info$outlier),
    qc = any(sample_info$experiment_or_qc %in% "qc")
  )
  out_todo$both = out_todo$outliers && out_todo$qc
  out_todo
}

create_lc_value_table = function(lc_info) {
  if (is.null(lc_info$lc_test)) {
    tibble::tibble(ssf = NA, trials = NA, success = NA) |>
      gt::gt()
  } else {
    lc_info$lc_test$values |>
      dplyr::select(class, trials, success) |>
      dplyr::transmute(ssf = class, trials = trials, success = success) |>
      gt::gt()
  }
}

create_lc_test_table = function(lc_info) {
  if (is.null(lc_info$lc_test)) {
    tibble::tibble(`P-Value` = NA, `Less-Median` = NA, `Missing` = NA) |>
      gt::gt()
  } else {
    lc_info$lc_test$binomial_test |>
      dplyr::transmute(
        `P-Value` = p.value,
        `Less-Median` = statistic,
        `Missing` = parameter
      ) |>
      gt::gt() |>
      gt::fmt_scientific(columns = `P-Value`, n_sigfig = 3)
  }
}


bind_analysis_related_tables = function(mwtab_data) {
  # mwtab_data = metadata(se_data)$mwtab
  analysis_data = mwtab_data$ANALYSIS

  if (!is.null(mwtab_data$NM)) {
    analysis_data = dplyr::bind_rows(analysis_data, mwtab_data$NM)
  }
  if (!is.null(mwtab_data$MS)) {
    analysis_data = dplyr::bind_rows(analysis_data, mwtab_data$MS)
  }
  if (!is.null(mwtab_data$CHROMATOGRAPHY)) {
    analysis_data = dplyr::bind_rows(analysis_data, mwtab_data$CHROMATOGRAPHY)
  }
  analysis_data
}
