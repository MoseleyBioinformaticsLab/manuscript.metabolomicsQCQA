detect_technical_replicate_samples = function(se_data) {
  # possible replicate datasets to work with:
  #  AN001764 - no replicate in the metadata; THIS IS NOT ACTUALLY REPLICATES GROUPING TOGETHER!
  #  AN004485 - has biological replicate in metadata;
  #  AN003552 - no replicate in the metadata, maybe replicates? Hard to tell
  #  AN006191 - has triplicate technical replicates
  #
  # how to detect these?
  # Can we correlate leveshtein distances with correlation based distances?
  # se_pca = tar_read(se_pca_g01); se_data = se_pca[["AN006191"]]; tar_source("R")
  # se_pca = tar_read(se_pca_g01); se_data = se_pca[["AN006248"]]; tar_source("R")

  # need the ici data here too!
  pca_data = get_specific_metadata(se_data, "pca_groups")
  se_ici = get_specific_metadata(se_data, "icikt")$cor
  ici_square = ICIKendallTau::long_df_2_cor_matrix(se_ici)

  last_run = check_qc_outputs(se_data)

  if (is.null(last_run$pca)) {
    return(NULL)
  }
  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  ici_order = pca_data_use$order$ici_order
  ici_square = ici_square[ici_order, ici_order]
  string_dist = calculate_stringdist_distances(ici_order, "cosine")
  ici_dist = 1 - ici_square

  order_data = pca_data_use$order$dendrogram$indices
  ici_dist = ici_dist[order_data, order_data]
  string_dist = string_dist[order_data, order_data]

  ici_df = ICIKendallTau::cor_matrix_2_long_df(ici_dist) |>
    dplyr::mutate(cor_dist = cor)
  string_df = ICIKendallTau::cor_matrix_2_long_df(string_dist) |>
    dplyr::mutate(string_dist = cor)

  dist_df = dplyr::left_join(
    ici_df[, c("s1", "s2", "cor_dist")],
    string_df[, c("s1", "s2", "string_dist")],
    by = c("s1", "s2")
  )

  dist_df = dist_df |>
    dplyr::filter(!(s1 == s2))

  dist_summary = dist_df |>
    dplyr::summarise(
      sd_cor = sd(cor_dist),
      mean_cor = mean(cor_dist),
      sd_string = sd(string_dist),
      mean_string = mean(string_dist),
      limit_cor = mean_cor - (2.5 * sd_cor),
      limit_string = mean_string - (1 * sd_string)
    )

  test_plot = dist_df |>
    ggplot(aes(x = cor_dist, y = string_dist)) +
    geom_jitter(height = 0.01, width = 0.01, alpha = 0.2) +
    geom_vline(data = dist_summary, aes(xintercept = limit_cor)) +
    geom_hline(
      data = dist_summary,
      aes(yintercept = limit_string)
    )

  dist_df = dist_df |>
    dplyr::mutate(
      replicate = dplyr::case_when(
        cor_dist <= dist_summary$limit_cor &
          string_dist <= dist_summary$limit_string ~ 1,
        TRUE ~ 0
      )
    )
  dist_rle = rle_2_df(rle(dist_df$replicate))

  dist_rle$locs = vector(mode = "list", length = nrow(dist_rle))

  if (!any(dist_rle$value == 1)) {
    return(NULL)
  }

  for (irow in seq_len(nrow(dist_rle))) {
    if (irow == 1) {
      dist_rle$locs[[irow]] = seq(1, dist_rle$length[irow])
    } else {
      prev_loc = dist_rle[irow - 1, ]
      start_loc = prev_loc$locs[[1]][prev_loc$length] + 1
      dist_rle$locs[[irow]] = seq(
        start_loc,
        start_loc + dist_rle$length[irow] - 1
      )
    }
  }
  dist_rle = dist_rle |>
    dplyr::filter(value == 1)
  dist_rle$replicate_group = paste0("r", seq_len(nrow(dist_rle)))

  replicate_samples = purrr::map(
    dist_rle$locs,
    \(in_loc) {
      tmp_df = dist_df[in_loc, ]
      sort(unique(unlist(tmp_df[, c("s1", "s2")])))
    }
  ) |>
    unique()
  names(replicate_samples) = paste0(
    "technical_replicate_",
    seq_len(length(replicate_samples))
  )
  replicate_df = tibble::enframe(
    replicate_samples,
    name = "replicate_group",
    value = "sample_id"
  ) |>
    tidyr::unnest(cols = sample_id)

  if (sum(duplicated(replicate_df$sample_id)) > 0) {
    warning("Samples in multiple replicates, returning NULL.")
    return(NULL)
  }

  n_rep = nrow(replicate_df)
  n_samples = nrow(ici_data)
  if ((n_rep / n_samples) < 0.5) {
    warning(
      "Number of samples in replicates are not at least 50% of samples, returning NULL."
    )
    return(NULL)
  }

  return(replicate_df)
}


calculate_stringdist_distances = function(in_strings, method = "jaccard") {
  n_strings = length(in_strings)
  out_matrix = matrix(NA, nrow = n_strings, ncol = n_strings)
  rownames(out_matrix) = colnames(out_matrix) = in_strings

  for (i_string in seq(1, n_strings)) {
    rem_strings = seq(i_string, n_strings)
    rem_dists = stringdist::stringdist(
      in_strings[i_string],
      in_strings[rem_strings],
      method = method
    )
    out_matrix[i_string, rem_strings] = rem_dists
    out_matrix[rem_strings, i_string] = rem_dists
  }
  out_matrix
}

finding_n_per_ssf_1 = function(sets_qc_color) {
  n_per_ssf = sets_qc_color$n_samples_ssf

  split_n = strsplit(n_per_ssf, ":", fixed = TRUE)
  what_n = purrr::map_chr(split_n, \(in_n) {
    uniq_n = unique(as.numeric(in_n))
    if (length(uniq_n) == 1) {
      if (!is.na(uniq_n) && (uniq_n == 1)) {
        return("hit")
      } else {
        return("miss")
      }
    } else {
      return("miss")
    }
  })
}


finding_replicate_metadata = function() {
  pca_to_load = c(
    "se_pca_g01",
    "se_pca_g02",
    "se_pca_g03",
    "se_pca_g04",
    "se_pca_g05",
    "se_pca_g06",
    "se_pca_g07",
    "se_pca_g08"
  )

  possible_replicates = purrr::map(pca_to_load, \(in_pca) {
    se_pca = tar_read_raw(in_pca)

    rep_info = purrr::imap(se_pca, \(in_data, id) {
      if (nrow(in_data) == 0) {
        return(NULL)
      }
      f_types = get_specific_metadata(in_data, "factor_types")
      tibble::tibble(
        id = id,
        factors = f_types |>
          dplyr::filter(
            grepl("replicate", factor, ignore.case = TRUE),
            type %in% c("subject_sample_factor", "additional_sample_data")
          ) |>
          dplyr::pull(factor)
      )
    }) |>
      purrr::list_rbind()
    rep_info
  }) |>
    purrr::list_rbind()
}
