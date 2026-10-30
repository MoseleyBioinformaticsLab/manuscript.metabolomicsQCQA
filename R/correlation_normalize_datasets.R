sets_calculate_correlations_pca = function(sets_se) {
  # sets_se = tar_read(se_dynamic_range_g01); tar_source("R"); sets_se = sets_se[1:60];
  # sets_se = tar_read(se_dynamic_range_g02); tar_source("R"); sets_se = sets_se[1:60];
  # sets_se = tar_read(se_dynamic_range_g09); tar_source("R");
  # message("Running correlations ...")
  # sets_se = tar_read(se_dynamic_range_g20); tar_source("R"); se_data = sets_se[["AN002453"]]
  p_cor = progressr::progressor(length(sets_se))
  n_set = seq_along(sets_se)
  cor_sets = purrr::map2(
    sets_se,
    n_set,
    \(in_set, in_n) {
      p_cor(sprintf("Correlating N=%g", in_n))
      do_correlation(in_set)
    }
  )

  p_pca = progressr::progressor(length(cor_sets))
  n_set = seq_along(cor_sets)
  pca_sets = furrr::future_map2(
    cor_sets,
    n_set,
    \(in_set, in_n) {
      p_pca(sprintf("PCA'ing N=%g", in_n))
      run_pca_experimental_design(in_set)
    }
  )

  p_dip = progressr::progressor(length(pca_sets))
  n_set = seq_along(pca_sets)
  dip_sets = furrr::future_map2(
    pca_sets,
    n_set,
    \(in_set, in_n) {
      p_dip(sprintf("Dip'ing N=%g", in_n))
      run_pc_ssfgroup_diptest(in_set)
    }
  )
  dip_sets
}

get_just_values = function(in_df) {
  df_matrix = in_df |>
    dplyr::select(-feature_id) |>
    as.matrix()
  rownames(df_matrix) = in_df$feature_id
  df_matrix
}

median_normalize_matrix = function(in_matrix) {
  med_values = apply(in_matrix, 2, median, na.rm = TRUE)
  med_matrix = matrix(
    med_values,
    nrow = nrow(in_matrix),
    ncol = ncol(in_matrix),
    byrow = TRUE
  )
  norm_matrix = in_matrix / med_matrix
  return(norm_matrix)
}

check_isstop = function(in_metadata) {
  in_metadata$qc$do_next %in% "stop"
}

count_missingness_testleftcensored = function(raw_matrix, ssf_wide) {
  perc_missing = sum(is.na(raw_matrix)) /
    (nrow(raw_matrix) * ncol(raw_matrix)) *
    100

  raw_matrix = raw_matrix[, ssf_wide$sample_id]

  if (!is.na(perc_missing) && (perc_missing > 0)) {
    lc_results = try(ICIKendallTau::test_left_censorship(
      raw_matrix,
      sample_classes = ssf_wide$ssf
    ))
    if (inherits(lc_results, "try-error")) {
      out_res = list(lc_test = NULL, perc_missing = perc_missing)
    } else {
      tibble_lc2 = broom::tidy(lc_results$binomial_test)
      lc_results$binomial_test = tibble_lc2
      out_res = list(lc_test = lc_results, perc_missing = perc_missing)
    }
  } else {
    out_res = list(lc_test = NULL, perc_missing = perc_missing)
  }
  out_res
}

check_na_cor_single_value = function(se_data) {
  se_cor = metadata(se_data)$icikt
  se_cor_df = se_cor$cor
  if (!any(is.na(se_cor_df$raw))) {
    return(se_data)
  }

  cor_matrix = ICIKendallTau::long_df_2_cor_matrix(se_cor_df)
  na_cor = se_cor_df |>
    dplyr::filter(is.na(cor))
  all_na_cor = MatrixGenerics::colSums(is.na(cor_matrix))
  n_samples = ncol(cor_matrix)
  problem_samples = which(all_na_cor >= (n_samples - 1)) |> names()

  verify_na = se_cor_df |>
    dplyr::filter((s1 %in% problem_samples) | (s2 %in% problem_samples)) |>
    dplyr::filter(!(s1 == s2))

  se_cor_df = se_cor_df |>
    dplyr::filter(!((s1 %in% problem_samples) | (s2 %in% problem_samples)))
  metadata(se_data)$icikt$cor = se_cor_df

  expect_na = length(problem_samples) * (n_samples - 1)

  if (!(nrow(verify_na) == expect_na)) {
    metadata(
      se_data
    )$qc$warnings = "Number of correlation NA doesn't seem quite right."
  }

  raw_data = assays(se_data)$raw[, problem_samples, drop = FALSE]
  n_unique = apply(raw_data, 2, \(x) {
    tmp_x = x[!is.na(x)]
    length(unique(tmp_x))
  })
  removal_message = purrr::map(n_unique, \(n_values) {
    if (n_values == 1) {
      "NA ICI-Kt correlations, 1 unique value in sample."
    } else {
      "NA ICI-Kt correlations, multiple values in sample."
    }
  }) |>
    purrr::list_c()
  names(removal_message) = names(n_unique)

  removed_vector = se_data$removed
  names(removed_vector) = se_data$sample_id
  removed_vector[names(removal_message)] = removal_message
  se_data$removed = removed_vector

  n_keep = sum(is.na(se_data$removed))

  if (n_keep < 3) {
    metadata(se_data)$qc = add_step(
      metadata(se_data)$qc,
      "failed",
      "correlation"
    )
    metadata(se_data)$qc$do_next = "stop"
    metadata(se_data)$qc$why = "too few samples due to NA correlationsq"
  } else {
    metadata(se_data)$qc = add_step(
      metadata(se_data)$qc,
      "succeeded",
      "correlation"
    )
  }
  return(se_data)
}

do_correlation = function(se_data) {
  if (metadata(se_data) |> check_isstop()) {
    return(se_data)
  }
  # se_data = sets_se[[1]]
  use_cols = is.na(se_data$removed)

  raw_matrix = assays(se_data)$raw[, use_cols, drop = FALSE]

  # this actually should never happen, but it's a failsafe
  # and can point out possibly problematic datasets.
  min_val = min(raw_matrix, na.rm = TRUE)
  if (min_val < 0) {
    warning("There are negative values!")
    metadata(se_data)$qc = add_step(
      metadata(se_data)$qc,
      "failed",
      "correlation"
    )
    return(se_data)
  }
  n_cores = future::nbrOfWorkers()

  n_item = ncol(raw_matrix) * nrow(raw_matrix)
  n_comparison = choose(ncol(raw_matrix), 2)
  # try halving number of workers to get around memory issues
  if (n_item > 1e6 || n_comparison > 1e7) {
    plan(multicore, workers = ceiling(n_cores / 2))
  }
  se_ici = ICIKendallTau::ici_kendalltau(raw_matrix, return_matrix = FALSE)
  # make sure to put them back
  plan(multicore, workers = n_cores)

  metadata(se_data)$icikt = se_ici

  # here we check for NA correlations, why they may have happened, and
  # decide if we are continuing or stopping due to having < 3 samples
  # that can still be used, instead of just stopping due to blanket
  # NA correlation.
  se_data = check_na_cor_single_value(se_data)
  se_data
}


add_step = function(qc, id, message) {
  tmp_id = qc[[id]]
  keep_id = tmp_id[nchar(tmp_id) > 0]
  out_id = c(keep_id, message)
  qc[[id]] = out_id
  qc
}

count_features_samples_factors = function(lc_sets) {
  # lc_sets = tar_read(lc_g04)
  out_n = purrr::map(lc_sets, \(in_set) {
    # in_set = lc_sets[[1]]
    tibble::tibble(
      id = in_set$id,
      n_features = nrow(in_set$abundance),
      n_samples = ncol(in_set$abundance),
      n_factors = ncol(colData(in_set$abundance)) - 3
    )
  }) |>
    purrr::list_rbind()
}

sets_add_mwtab_do_left_censored = function(cor_set, setup_set) {
  # cor_set = tar_read(cor_g01); setup_set = tar_read(setup_g01)

  lc_set = purrr::map2(cor_set, setup_set, add_mwtab_do_left_censored)
  lc_set
}

add_mwtab_do_left_censored = function(cor_data, setup_data) {
  # cor_data = cor_set[[11]]; setup_data = setup_set[[11]]
  # cor_data = cor_set$AN000123; setup_data = setup_set$AN000123
  stopifnot(cor_data$id == setup_data$ID)
  raw_matrix = get_just_values(cor_data$raw)
  ssf_wide = cor_data$ssf

  feature_data = setup_data$METABOLITES

  matching_features = base::intersect(
    feature_data$feature_id,
    rownames(raw_matrix)
  )

  feature_data = feature_data |>
    dplyr::filter(feature_id %in% matching_features)

  raw_matrix = raw_matrix[feature_data$feature_id, ssf_wide$sample_id]

  norm_matrix = get_just_values(cor_data$normalized)[
    feature_data$feature_id,
    ssf_wide$sample_id
  ]

  qc_samples = determine_qc_samples(setup_data$SUBJECT_SAMPLE_FACTORS)
  ssf_wide = ssf_wide |>
    dplyr::mutate(
      experiment_or_qc = dplyr::case_when(
        sample_id %in% qc_samples ~ "qc",
        TRUE ~ "experiment"
      )
    )
  se_data = SummarizedExperiment(
    assays = SimpleList(raw = raw_matrix, normalized = norm_matrix),
    rowData = feature_data,
    colData = ssf_wide
  )

  lc_data = count_missingness_testleftcensored(raw_matrix, ssf_wide)

  # removing pieces we don't need because we have
  # the SE data object now.
  if (!is.null(setup_data$MS_METABOLITE_DATA)) {
    setup_data$MS_METABOLITE_DATA = NULL
  }
  if (!is.null(setup_data$NMR_METABOLITE_DATA)) {
    setup_data$NMR_METABOLITE_DATA = NULL
  }
  if (!is.null(setup_data$NMR_BINNED_DATA)) {
    setup_data$NMR_BINNED_DATA = NULL
  }
  setup_data$MEASUREMENTS = NULL
  setup_data$METABOLITES = NULL
  cor_data$raw = NULL
  cor_data$normalized = NULL
  cor_data$ssf = NULL
  cor_data$abundance = se_data
  cor_data$mwtab = setup_data
  cor_data$left_censorship = lc_data
  cor_data
}

determine_qc_samples = function(ssf_long) {
  if (!is.null(ssf_long[["sample_id"]])) {
    sample_qc = grepl(
      ".*pool.*|.*qc.*|.*blank.*|.*quality.*",
      ssf_long$sample_id,
      ignore.case = TRUE
    )
  } else {
    sample_qc = rep(FALSE, nrow(ssf_long))
  }

  if (!is.null(ssf_long[["subject"]])) {
    subject_qc = grepl(
      ".*pool.*|.*qc.*|.*blank.*|.*quality.*",
      ssf_long$subject,
      ignore.case = TRUE
    )
  } else {
    subject_qc = rep(FALSE, nrow(ssf_long))
  }

  if (!is.null(ssf_long[["item"]])) {
    item_qc = grepl(
      ".*pool.*|.*qc.*|.*blank.*|.*quality.*",
      ssf_long$item,
      ignore.case = TRUE
    )
  } else {
    item_qc = rep(FALSE, nrow(ssf_long))
  }

  all_locs = sample_qc | subject_qc | item_qc
  unique(ssf_long$sample_id[all_locs])
}
