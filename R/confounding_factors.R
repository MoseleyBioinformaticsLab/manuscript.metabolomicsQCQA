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

do_fishers_justp = function(p_vals) {
  out_res = poolr::fisher(p_vals)
  out_res$p
}


add_confounding_factor_status = function(
  sets_qc_color,
  sets_diptest
) {
  # tar_load(c(sets_qc_color, sets_diptest)); tar_source("R")

  null_dipres = purrr::map_lgl(sets_diptest, is.null)
  keep_dipres = sets_diptest[!null_dipres]
  keep_sets_diptest = purrr::map(keep_dipres, \(in_set) {
    # in_factor = keep_sets[[1]]
    # in_factor = keep_sets[[162]]
    # in_factor = keep_sets[[415]]
    # in_factor = keep_sets[[1725]]

    if (all(is.na(in_set$anova$p.value)) || nrow(in_set$dip) == 0) {
      return(NULL)
    }
    if (nrow(in_set$dip) == 1) {
      summarized_dip = tibble::tibble(
        PC = in_set$dip$PC,
        dippooled = in_set$dip$p.value
      )
    } else {
      summarized_dip = in_set$dip |>
        dplyr::summarise(dippooled = do_fishers_justp(padjust), .by = PC)
    }

    in_set$summarized_dip = summarized_dip
    in_set
  })

  null_sets2 = purrr::map_lgl(keep_sets_diptest, is.null)
  keep_sets_diptest2 = keep_sets_diptest[!null_sets2]

  dip_stats = purrr::map(keep_sets_diptest2, \(in_set) {
    in_set$summarized_dip |>
      dplyr::mutate(internal_id = in_set$analysis_id)
  }) |>
    purrr::list_rbind()

  dip_min = dip_stats |>
    dplyr::summarise(min_p = min(dippooled), .by = internal_id)
  dip_pc1 = dip_stats |>
    dplyr::filter(PC %in% "PC1") |>
    dplyr::mutate(
      confounding_factor_detected = dplyr::case_when(
        dippooled <= 0.05 ~ "yes",
        TRUE ~ "no"
      )
    )

  sets_qc_color_confounds = dplyr::left_join(
    sets_qc_color,
    dip_pc1,
    by = "internal_id"
  )

  return(sets_qc_color_confounds)
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

sets_grab_confounded = function(se_pca, not_gray_pca) {
  # se_pca = tar_read(se_pca_g01); tar_load(not_gray_pca); tar_source("R")
  has_confounded = not_gray_pca |>
    dplyr::filter(confounding_factor_detected %in% "yes")

  match_pca = which(names(se_pca) %in% has_confounded$internal_id)

  if (length(match_pca) > 0) {
    return(se_pca[match_pca])
  } else {
    return(NULL)
  }
}
