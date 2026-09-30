run_pc_ssfgroup_diptest = function(se_data) {
  # sets_pca = tar_read(se_pca_g34); se_data = sets_pca[["AN000923"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g05); se_data = sets_pca[["AN000262"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g25); se_data = sets_pca[["AN005412"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g25); se_data = sets_pca[["AN000118"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g02); se_data = sets_pca[["AN002914"]]; tar_source("R")

  report_list = check_qc_outputs(se_data)

  pca_groups = get_specific_metadata(se_data, "pca_groups")
  if (!is.null(report_list$pca)) {
    use_group = pca_groups[report_list$pca[length(report_list$pca)]]
  } else {
    return(NULL)
  }

  pca_scores = use_group[[1]]$pca$pca
  pca_variance = use_group[[1]]$pca$variance
  pca_anova = use_group[[1]]$pca$anova

  if (nrow(pca_scores) < 4) {
    return(NULL)
  }

  test_pcs = unique(pca_anova$PC)

  ssf_data = pca_scores[, c("sample_id", "ssf")]

  samples_ssf = split(ssf_data$sample_id, ssf_data$ssf)
  dip_pcs = purrr::map(test_pcs, \(in_pc) {
    # in_pc = "PC1"
    dip_ssf = purrr::imap(samples_ssf, \(in_samples, ssf_id) {
      # in_samples = samples_ssf[[1]]
      if (length(in_samples) < 2) {
        return(NULL)
      }
      use_scores = pca_scores |>
        dplyr::filter(sample_id %in% in_samples) |>
        dplyr::pull(tidyselect::all_of(in_pc))
      out_dip = broom::tidy(diptest::dip.test(use_scores))
      out_dip$ssf = ssf_id
      out_dip
    }) |>
      purrr::list_rbind()
    if (nrow(dip_ssf) > 0) {
      dip_ssf$PC = in_pc
    }

    dip_ssf
  }) |>
    purrr::list_rbind()

  if (nrow(dip_pcs) > 1) {
    dip_pcs = dip_pcs |>
      dplyr::mutate(
        padjust = p.adjust(p.value, method = "bonferroni"),
        .by = ssf
      )
  } else if (nrow(dip_pcs) == 1) {
    dip_pcs = dip_pcs |>
      dplyr::mutate(padjust = p.value)
  } else {
    return(NULL)
  }

  list(
    analysis_id = metadata(se_data)$mwtab$ID,
    variance = pca_variance,
    anova = pca_anova,
    dip = dip_pcs,
    scores = pca_scores,
    factor_types = get_specific_metadata(se_data, "factor_types")
  )
}

sets_run_pc_ssfgroup_diptest = function(sets_pca) {
  out_info = purrr::map(sets_pca, run_pc_ssfgroup_diptest)
  out_info
}


add_confounding_factor_status = function(
  sets_qc_color,
  sets_confounding_factors
) {
  # tar_load(c(sets_qc_color, sets_confounding_factors))
  sets_yellow_red = sets_qc_color |>
    dplyr::filter(pca_anova_pvalue_status %in% c("yellow", "red"))
  possible_confound = sets_confounding_factors |>
    dplyr::filter(
      analysis_id %in% sets_yellow_red$analysis_id,
      p.value <= 1e-10,
      statistic >= 300,
      percent >= 0.3
    ) |>
    dplyr::pull(analysis_id)
  sets_qc_confounds = sets_qc_color |>
    dplyr::mutate(
      confounding_factor_detected = dplyr::case_when(
        analysis_id %in% possible_confound ~ "yes",
        TRUE ~ "no"
      )
    )
  sets_qc_confounds
}

confound_info_analysis = function(not_gray_pca) {
  has_confound = not_gray_pca |>
    dplyr::filter(confounding_factor_detected %in% "yes")

  confound_yellow = has_confound |>
    dplyr::filter(pca_anova_pvalue_status %in% "yellow")
  confound_red = has_confound |>
    dplyr::filter(pca_anova_pvalue_status %in% "red")

  list(
    n_overall = nrow(has_confound),
    n_yellow = nrow(confound_yellow),
    n_red = nrow(confound_red)
  )
}

examine_confounding_pvalues = function(sets_confounding_factors) {
  # tar_load(sets_confounding_factors)
  null_factors = purrr::map_lgl(sets_confounding_factors, is.null)
  keep_confounding = sets_confounding_factors[!null_factors]

  confounding_compare = purrr::map(keep_confounding, \(in_confounding) {
    # in_confounding = keep_confounding[[1]]
    # in_confounding = keep_confounding[["se_confounds_g05.AN000262"]]
    cf_vals = in_confounding |>
      dplyr::filter(variable %in% "confounding") |>
      dplyr::mutate(
        sig = dplyr::case_when(
          p.value <= 0.05 ~ "yes",
          TRUE ~ "no"
        )
      )

    other_vals = in_confounding |>
      dplyr::filter_out(variable %in% "confounding") |>
      dplyr::mutate(
        sig = dplyr::case_when(
          p.value <= 0.05 ~ "yes",
          TRUE ~ "no"
        )
      )

    together_vals = dplyr::full_join(
      other_vals[, c(
        "statistic",
        "p.value",
        "padjust",
        "percent",
        "cumulative",
        "sig",
        "variable",
        "PC"
      )],
      cf_vals[, c(
        "statistic",
        "p.value",
        "padjust",
        "sig",
        "variable",
        "PC"
      )],
      by = "PC",
      suffix = c(".variable", ".confound")
    )

    together_vals$analysis_id = cf_vals$analysis_id[1]
    together_vals
  }) |>
    purrr::list_rbind()
}
