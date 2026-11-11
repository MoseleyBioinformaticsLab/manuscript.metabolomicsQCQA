run_pc_ssfgroup_diptest = function(se_data, just_ssf = "yes") {
  # sets_pca = tar_read(se_pca_g34); se_data = sets_pca[["AN000923"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g05); se_data = sets_pca[["AN000262"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g25); se_data = sets_pca[["AN005412"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g25); se_data = sets_pca[["AN000118"]]; tar_source("R")
  # sets_pca = tar_read(se_pca_g02); se_data = sets_pca[["AN002914"]]; tar_source("R")

  # sets_pca = tar_read(se_pca_g01); sets_pca = sets_pca[1:100]; tar_source("R")
  # sets_pca = readRDS("tests_problems/dip_issue.rds");  se_data = sets_pca[[97]]; just_ssf = "yes"; tar_source("R")

  if (metadata(se_data) |> check_isstop()) {
    return(se_data)
  }
  in_pca = get_last_pca(se_data)$pca

  pca_scores = in_pca$pca
  pca_variance = in_pca$variance
  pca_anova = in_pca$anova

  if (nrow(pca_scores) < 4) {
    dip_pc_factor = NULL
    metadata(se_data)$qc = add_step(
      metadata(se_data)$qc,
      "failed",
      "diptest"
    )
    return(se_data)
  }

  test_pcs = unique(pca_anova$PC)

  if (just_ssf %in% "yes") {
    test_factors = "ssf"
  } else {
    factor_types = get_specific_metadata(se_data, "factor_types")
    test_factors = factor_types |>
      dplyr::filter(type %in% c("additional_sample_data", "additional_ssf")) |>
      dplyr::pull(factor)
    test_factors = c("ssf", test_factors)
  }

  dip_pc_factor = purrr::map(test_pcs, \(in_pc) {
    # in_pc = "PC1"
    out_factor = purrr::map(test_factors, \(in_factor) {
      # in_factor = "cohort"
      if (in_factor %in% names(pca_scores)) {
        tmp_df = data.frame(
          PC = pca_scores[[in_pc]],
          group = pca_scores[[in_factor]]
        )
        dip_df = tmp_df |>
          dplyr::summarise(dip = calculate_dip(PC), .by = group) |>
          tidyr::unnest(cols = c(dip)) |>
          dplyr::mutate(
            padjust = p.adjust(p.value, method = "bonferroni"),
            PC = in_pc,
            variable = in_factor
          )
        return(dip_df)
      } else {
        return(NULL)
      }
    }) |>
      purrr::list_rbind()
    out_factor = out_factor |>
      dplyr::filter_out(is.na(p.value))

    if (nrow(out_factor) > 0) {
      pooled_p = out_factor |>
        dplyr::summarise(
          dippooled = do_fishers_justp(padjust),
          .by = variable
        ) |>
        dplyr::mutate(PC = out_factor$PC[1])
    } else {
      pooled_p = tibble::tibble(
        variable = NA_character_,
        dippooled = NA,
        PC = in_pc
      )
    }

    return(pooled_p)
  }) |>
    purrr::list_rbind()

  if (all(is.na(dip_pc_factor$dippooled))) {
    metadata(se_data)$qc = add_step(
      metadata(se_data)$qc,
      "failed",
      "diptest"
    )
    return(se_data)
  }
  metadata(se_data)$diptest = dip_pc_factor
  metadata(se_data)$qc = add_step(
    metadata(se_data)$qc,
    "succeeded",
    "diptest"
  )
  return(se_data)
}

sets_run_pc_ssfgroup_diptest = function(sets_pca) {
  out_info = purrr::map(sets_pca, run_pc_ssfgroup_diptest)
  out_info
}

do_fishers_justp = function(p_vals) {
  out_res = poolr::fisher(p_vals)
  out_res$p
}


create_final_dipresults = function(sets_diptest) {
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

  dip_stats
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

sets_determine_confound_is_missing_ssf = function(
  sets_confounded
) {
  # tar_load(sets_confounded); tar_source("R")
  p_do = progressr::progressor(length(sets_confounded))
  n_do = seq_along(sets_confounded)
  is_missing_ssf = furrr::future_map2(
    sets_confounded,
    n_do,
    \(in_set, in_n) {
      p_do(sprintf("Confounding N=%g", in_n))
      determine_confound_is_missing_ssf(in_set)
    }
  ) |>
    purrr::list_rbind()
  is_missing_ssf
}

determine_confound_is_missing_ssf = function(se_data) {
  # se_data = sets_confounded[[1]]
  # se_data = sets_confounded[["se_confounded_g20.AN006863"]]
  # se_data = sets_confounded[[5]]
  # se_data = sets_confounded[[12]]
  # se_data = sets_confounded[[15]]

  se_tmp = run_pc_ssfgroup_diptest(se_data, just_ssf = "no")

  dip_all = get_specific_metadata(se_tmp, "diptest")

  use_pca = get_last_pca(se_data)
  use_anova = use_pca$pca$anova

  use_pcs = c("PC1", "PC2")

  dip_ssf = dip_all |>
    dplyr::filter(variable %in% "ssf", PC %in% use_pcs) |>
    dplyr::arrange(PC)
  dip_other = dip_all |>
    dplyr::filter(PC %in% use_pcs) |>
    dplyr::filter_out(variable %in% "ssf") |>
    dplyr::arrange(PC)
  dip_compare = dplyr::left_join(
    dip_other,
    dip_ssf,
    by = "PC",
    suffix = c(".other", ".ssf")
  )
  dip_compare_nozero = replace_zero_pvalue(dip_compare)

  anova_ssf = use_anova |>
    dplyr::filter(
      variable %in% dip_all$variable,
      variable %in% "ssf",
      PC %in% use_pcs
    ) |>
    dplyr::arrange(PC) |>
    dplyr::select(PC, variable, padjust)
  anova_other = use_anova |>
    dplyr::filter(variable %in% dip_all$variable, PC %in% use_pcs) |>
    dplyr::filter_out(variable %in% "ssf") |>
    dplyr::arrange(PC) |>
    dplyr::select(PC, variable, padjust)
  anova_compare = dplyr::left_join(
    anova_other,
    anova_ssf,
    by = "PC",
    suffix = c(".other", ".ssf")
  )
  anova_compare_nozero = anova_compare |>
    dplyr::filter_out(is.na(padjust.other) | is.na(padjust.ssf)) |>
    replace_zero_pvalue()

  dip_compare_nozero = dip_compare_nozero |>
    dplyr::mutate(
      diff = -1 * log10(dippooled.other / dippooled.ssf)
    )

  analysis_id = metadata(se_data)$mwtab$ID

  if (any(dip_compare_nozero$dippooled.ssf <= 0.05)) {
    verify_confound = "yes"
  } else {
    verify_confound = "not confounded"
  }

  less_dip = dip_compare_nozero |>
    dplyr::filter(diff < 0)
  if (nrow(less_dip) > 0) {
    notsig_dip = less_dip |>
      dplyr::filter(dippooled.other > 0.2)
    if (nrow(notsig_dip) > 0) {
      dip_possible_missing_ssf = gsub("ssf_", "", notsig_dip$variable.other) |>
        unique() |>
        paste0(collapse = "; ")
      dip_pdiff = min(notsig_dip$diff)
      dip_possible_missing_evidence = "high"
    } else {
      dip_pdiff = min(less_dip$diff)
      dip_possible_missing_ssf = "unknown"
      dip_possible_missing_evidence = "medium"
    }
  } else {
    dip_pdiff = NA
    dip_possible_missing_ssf = "unknown"
    dip_possible_missing_evidence = "low"
  }

  if (nrow(anova_compare_nozero) == 0) {
    anova_pdiff = NA
    anova_possible_missing_ssf = "unknown"
    anova_possible_missing_evidence = "none"
  } else {
    anova_compare_nozero = anova_compare_nozero |>
      dplyr::mutate(diff = -1 * log10(padjust.other / padjust.ssf))

    more_anova = anova_compare_nozero |>
      dplyr::filter(diff > 0)
    if (nrow(more_anova) > 0) {
      moresig_anova = more_anova |>
        dplyr::filter(padjust.other < 0.1)
      if (nrow(moresig_anova) > 0) {
        anova_possible_missing_ssf = gsub(
          "ssf_",
          "",
          moresig_anova$variable.other
        ) |>
          unique() |>
          paste0(collapse = "; ")
        anova_pdiff = max(more_anova$diff)
        anova_possible_missing_evidence = "high"
      } else {
        anova_pdiff = max(more_anova$diff)
        anova_possible_missing_ssf = "unknown"
        anova_possible_missing_evidence = "medium"
      }
    } else {
      anova_pdiff = NA
      anova_possible_missing_ssf = "unknown"
      anova_possible_missing_evidence = "low"
    }
  }

  out_res = tibble::tibble(
    internal_id = analysis_id,
    verified_confounded = verify_confound,
    dip_pdiff = dip_pdiff,
    dip_missing_ssf = dip_possible_missing_ssf,
    dip_missing_evidence = dip_possible_missing_evidence,
    anova_pdiff = anova_pdiff,
    anova_missing_ssf = anova_possible_missing_ssf,
    anova_missing_evidence = anova_possible_missing_evidence
  )
  out_res
}

replace_zero_pvalue = function(in_df) {
  # in_df = anova_compare
  if (nrow(in_df) == 0) {
    return(in_df)
  }
  pval_cols = c(
    "padjust.ssf",
    "padjust.other",
    "dippooled.ssf",
    "dippooled.other"
  )

  all_pval = in_df |>
    dplyr::select(tidyselect::any_of(pval_cols)) |>
    unlist() |>
    unique()

  min_nonzero = min(all_pval[all_pval > 0])
  if (any(all_pval == 0)) {
    if (min_nonzero > 1e-16) {
      replace_pval = 1e-16
    } else if (min_nonzero > 1e-32) {
      replace_pval = 1e-32
    } else if (min_nonzero > 1e-99) {
      replace_pval = 1e-100
    } else if (min_nonzero > 1e-199) {
      replace_pval = 1e-200
    }

    plocs = names(in_df)[names(in_df) %in% pval_cols]

    for (iloc in plocs) {
      zero_loc = in_df[[iloc]] == 0
      in_df[[iloc]][zero_loc] = replace_pval
    }
  }

  return(in_df)
}

calculate_dip = function(values) {
  if (length(values) < 2) {
    out_dip = structure(
      list(
        statistic = c(D = NA),
        p.value = NA,
        nobs = NA_integer_,
        alternative = NA_character_,
        method = "Hartigans' dip test for unimodality / multimodality",
        data.name = NA_character_
      ),
      class = "htest"
    )
  } else {
    out_dip = diptest::dip.test(values)
  }
  tidy_dip = broom::tidy(out_dip)
  tidy_dip
}
