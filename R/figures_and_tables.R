calculate_losses_at_steps = function(sets_qc_color) {
  # tar_load(sets_qc_color); tar_source("R")
  n_each = summarize_qc_results_table(sets_qc_color, group_none = "yes")
  n_total = nrow(sets_qc_color)

  stopped_se = n_each |>
    dplyr::filter(completed_steps %in% "none", is.na(failed_steps)) |>
    dplyr::pull(n)
  stopped_n3 = n_each |>
    dplyr::filter(
      completed_steps %in% "se; dynamic_range",
      failed_steps %in% "none"
    ) |>
    dplyr::pull(n)
  stopped_correlation_n3 = n_each |>
    dplyr::filter(failed_steps %in% "correlation") |>
    dplyr::pull(n)
  gray_pca = sum(sets_qc_color$pca_anova_pvalue_status %in% "gray")

  results_df = tibble::tribble(
    ~step                            , ~n                          , ~type   ,
    "Input"                          , n_total                     , "total" ,
    "No Data"                        , -1 * stopped_se             , "loss"  ,
    "< 3 Samples"                    , -1 * stopped_n3             , "loss"  ,
    "< 3 Samples After\nCorrelation" , -1 * stopped_correlation_n3 , "loss"  ,
    "Gray PCA"                       , -1 * gray_pca               , "loss"
  )
  results_df
}

create_loss_waterfall_plot = function(waterfall_input) {
  library(waterfalls)
  loss_waterfall_plot = waterfall_input |>
    waterfall(
      calc_total = TRUE,
      total_axis_text = "Analyzed",
      total_rect_color = "lightgreen",
      total_rect_text_color = "black",
      put_rect_text_outside_when_value_below = 300,
      rect_text_size = 1.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Step", y = "N-Datasets")
  loss_waterfall_plot
}

create_nose_reasons_bar_plot = function(sets_qc_color) {
  # tar_load(sets_qc_color); tar_source("R");
  n_each = summarize_qc_results_table(sets_qc_color, group_none = "no")
  none_steps = n_each |>
    dplyr::filter(completed_steps %in% "none") |>
    dplyr::arrange(n) |>
    dplyr::mutate(reason = factor(failed_why, levels = failed_why))
  nose_histogram_plot = none_steps |>
    ggplot(aes(x = n, y = reason)) +
    geom_col() +
    labs(x = "N-Datasets", y = "Reason") +
    scale_x_continuous(expan = expansion(mult = c(0, 0.05)))
  nose_histogram_plot
}

create_dataset_sizes_plot = function(sets_qc_color) {
  # tar_load(sets_qc_color); tar_source("R")
  all_sizes = purrr::map(seq_len(nrow(sets_qc_color)), \(in_row) {
    tmp_df = sets_qc_color$data_sizes_df[[in_row]]
    tmp_df$analysis_id = sets_qc_color$analysis_id[in_row]
    tmp_df
  }) |>
    purrr::list_rbind()

  all_sizes = all_sizes |>
    dplyr::mutate(
      Source = factor(
        Source,
        levels = c("Metabolomics Workbench", "Initial-SE", "Kept-SE")
      )
    )
  feature_plot = all_sizes |>
    ggplot(aes(x = log10(`N-Features`))) +
    geom_histogram(bins = 100) +
    facet_wrap(~Source, ncol = 1) +
    labs(x = "Log10(N-Features")

  sample_plot = all_sizes |>
    ggplot(aes(x = log10(`N-Samples`))) +
    geom_histogram(bins = 100) +
    facet_wrap(~Source, ncol = 1) +
    labs(x = "Log10(N-Samples)")

  (feature_plot | sample_plot) + plot_annotation(tag_levels = "A")
}


create_technical_replicate_plot = function(se_data) {
  # se_pca = tar_read(se_pca_g01);  se_data = se_pca[["AN006191"]];  tar_source("R")

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

  dist_df = dist_df |>
    dplyr::mutate(
      replicate = dplyr::case_when(
        cor_dist <= dist_summary$limit_cor &
          string_dist <= dist_summary$limit_string ~ 1,
        TRUE ~ 0
      )
    )

  rep_df = dist_df |>
    dplyr::filter(replicate == 1)

  dist_plot = dist_df |>
    dplyr::filter(replicate == 0) |>
    ggplot(aes(x = cor_dist, y = string_dist)) +
    geom_point(alpha = 0.2) +
    geom_vline(xintercept = dist_summary$limit_cor) +
    geom_hline(yintercept = dist_summary$limit_string) +
    geom_point(data = rep_df, color = "red", alpha = 0.2)

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

  uniq_rep = tibble::tibble(
    replicate_group = unique(replicate_df$replicate_group)
  ) |>
    dplyr::mutate(
      group = paste0(
        "g_",
        stringr::str_pad(
          seq_len(length(replicate_group)),
          2,
          side = "left",
          pad = "0"
        )
      )
    )

  tmp_df = data.frame(sample_id = ici_order)
  replicate_df = dplyr::left_join(tmp_df, replicate_df, by = "sample_id")
  replicate_df = dplyr::left_join(
    replicate_df,
    uniq_rep,
    by = "replicate_group"
  )
  rep_groups = unique(replicate_df$group)
  n_reps = length(rep_groups)

  rep_colors = scales::viridis_pal()(n_reps)
  names(rep_colors) = rep_groups

  heatmap_colors = list(
    group = rep_colors
  )

  left_var = ComplexHeatmap::HeatmapAnnotation(
    df = replicate_df[, "group", drop = FALSE],
    col = heatmap_colors,
    which = "row",
    show_annotation_name = FALSE
  )
  bottom_var = ComplexHeatmap::HeatmapAnnotation(
    df = replicate_df[, "group", drop = FALSE],
    col = heatmap_colors,
    which = "column",
    show_legend = FALSE
  )

  n_value = 40
  cor_colors = circlize::colorRamp2(
    seq(-1, 1, length.out = n_value),
    scico::scico(n_value, palette = "berlin")
  )
  out_heatmap = ComplexHeatmap::Heatmap(
    matrix = ici_square,
    col = cor_colors,
    name = "ICI-Kt",
    show_column_names = FALSE,
    row_labels = replicate_df$group,
    show_row_names = TRUE,
    cluster_columns = pca_data_use$order$dendrogram$dendrogram,
    cluster_rows = pca_data_use$order$dendrogram$dendrogram
  )

  return(list(dotplot = dist_plot, heatmap = out_heatmap))
}

create_example_confound_plot = function(se_data) {
  # se_pca = tar_read(se_pca_g34); se_data = se_pca[["AN000454"]]; tar_source("R")
  pca_data = get_specific_metadata(se_data, "pca_groups")
  last_run = check_qc_outputs(se_data)
  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  pca_scores = pca_data_use$pca$pca

  out_plot = pca_scores |>
    ggplot(aes(x = PC1, y = PC2, color = ssf)) +
    geom_point(size = 2) +
    theme(legend.position = "top") +
    labs(subtitle = metadata(se_data)$mwtab$ID)
  return(out_plot)
}

create_status_heatmap = function(
  parsed_percentage_overall,
  order_n = "yes",
  trim_90 = "yes",
  use_log = "no"
) {
  color_list = c("red" = "#D81B60", "green" = "#004D40", "yellow" = "#FFC107")

  if (order_n %in% "yes") {
    parsed_percentage_overall = parsed_percentage_overall |>
      dplyr::arrange(dplyr::desc(N))
  }

  if (trim_90 %in% "yes") {
    perc_n = cumsum(parsed_percentage_overall$N) /
      sum(parsed_percentage_overall$N)
    which_90 = which(perc_n >= 0.9) |> min()
    parsed_percentage_overall = parsed_percentage_overall[1:which_90, ]
  }
  color_df = parsed_percentage_overall[, 1:7]
  factor_df = purrr::map(color_df, \(in_status) {
    as.numeric(factor(in_status, levels = c("green", "yellow", "red")))
  }) |>
    dplyr::bind_cols()
  color_factor = c(
    "1" = color_list[["green"]],
    "2" = color_list[["yellow"]],
    "3" = color_list[["red"]]
  )
  factor_matrix = as.matrix(factor_df)

  row_counts = parsed_percentage_overall$N
  if (use_log %in% "yes") {
    row_ha = ComplexHeatmap::rowAnnotation(
      `Log10(N)` = ComplexHeatmap::anno_barplot(log10(row_counts))
    )
  } else {
    row_ha = ComplexHeatmap::rowAnnotation(
      N = ComplexHeatmap::anno_barplot(row_counts)
    )
  }

  show_colnames = gsub("_status", "", colnames(factor_matrix))
  show_colnames = purrr::map_chr(show_colnames, convert_description_to_sentence)
  out_heatmap = ComplexHeatmap::Heatmap(
    factor_matrix,
    col = color_factor,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_ha,
    show_heatmap_legend = FALSE,
    column_labels = show_colnames,
    column_names_side = "top",
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_names_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )

  out_heatmap
  # draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
}

create_mmr_logtransform_plots = function(se_mmr_combine) {
  # tar_load(se_mmr_combine)
  se_mmr_combine = se_mmr_combine |>
    dplyr::mutate(neg_fraction = N_neg / N_present)
  se_mmr_lowmax = se_mmr_combine |>
    dplyr::filter(max <= 10^12, dplyr::between(mmr, 0.2, 200))

  mmr_histogram_plot = se_mmr_lowmax |>
    ggplot(aes(x = log10(mmr))) +
    geom_histogram(bins = 100) +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    labs(x = "Log10(MMR)")

  mmr_vs_neg_plot = se_mmr_lowmax |>
    dplyr::filter(neg_fraction > 0) |>
    ggplot(aes(x = log10(mmr), y = neg_fraction)) +
    geom_point() +
    labs(x = "Log10(MMR)", y = "Negative-Fraction")

  mmr_histogram_less1000 = se_mmr_lowmax |>
    dplyr::filter(max <= 1000) |>
    ggplot(aes(x = log10(mmr))) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = log10(c(0.8, 1.2))) +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    labs(x = "Log10(MMR)")
  return(list(
    mmr_distribution = mmr_histogram_plot,
    mmr_vs_neg = mmr_vs_neg_plot,
    less1000 = mmr_histogram_less1000
  ))
}

create_rsd_related_plots = function(sets_qc_confounds_rsd) {
  rsd_measurement = sets_qc_confounds_rsd |>
    dplyr::select(
      internal_id,
      measurement_type,
      biological_rsd_median,
      technical_replicate_rsd_median
    ) |>
    tidyr::pivot_longer(
      c(-internal_id, -measurement_type),
      values_to = "rsd",
      names_to = "type"
    )
  rsd_plot = rsd_measurement |>
    dplyr::filter(rsd <= 1) |>
    ggplot(aes(x = rsd, color = type)) +
    geom_density(size = 2) +
    facet_wrap(~measurement_type, ncol = 1) +
    theme(legend.position = "top")

  rsd_modes_measurement = rsd_measurement |>
    dplyr::group_by(measurement_type, type) |>
    dplyr::summarise(mode = calculate_mode(rsd), .groups = "drop")
  rsd_modes_all = rsd_measurement |>
    dplyr::group_by(type) |>
    dplyr::summarise(mode = calculate_mode(rsd), .groups = "drop")

  rsd_summary = rsd_measurement |>
    dplyr::filter(type %in% "technical_replicate_rsd_median") |>
    dplyr::filter(rsd <= 1) |>
    dplyr::summarise(mean = mean(rsd, na.rm = TRUE), sd = sd(rsd, na.rm = TRUE))
}

calculate_mode = function(values) {
  values = values[!is.na(values)]
  density_estimate = stats::density(values)
  mode_value = density_estimate$x[which.max(density_estimate$y)]
  mode_value
}
