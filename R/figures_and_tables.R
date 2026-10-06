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
    labs(x = "Log10(N-Features)")

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
    cluster_rows = pca_data_use$order$dendrogram$dendrogram,
    row_names_gp = grid::gpar(fontsize = 5)
  )

  return(list(dotplot = dist_plot, heatmap = out_heatmap))
}

create_example_confound_plot = function(se_data, confounding_var = NULL) {
  # se_pca = tar_read(se_pca_g34); se_data = se_pca[["AN000454"]]; tar_source("R")
  # se_pca = tar_read(se_pca_g20); se_data = se_pca[["AN006863"]]; tar_source("R")
  pca_data = get_specific_metadata(se_data, "pca_groups")
  last_run = check_qc_outputs(se_data)
  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  pca_scores = pca_data_use$pca$pca

  pc1_dist_plot = pca_scores |>
    ggplot(aes(x = PC1)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ssf, ncol = 1) +
    theme(strip.background = element_blank())
  pc2_dist_plot = pca_scores |>
    ggplot(aes(x = PC2)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ssf, ncol = 1) +
    theme(strip.background = element_blank())

  scores_plot = pca_scores |>
    ggplot(aes(x = PC1, y = PC2, color = ssf)) +
    geom_point(size = 2) +
    theme(legend.position = "top") +
    labs(subtitle = metadata(se_data)$mwtab$ID)

  if (!is.null(confounding_var)) {
    confound_plot = pca_scores |>
      ggplot(aes(x = PC1, y = .data[[confounding_var]])) +
      geom_sina() +
      geom_boxplot(alpha = 0)
    return(list(
      pc1_histogram = pc1_dist_plot,
      pc2_histogram = pc2_dist_plot,
      scores = scores_plot,
      confound_sinaplot = confound_plot
    ))
  } else {
    return(list(
      pc1_histogram = pc1_dist_plot,
      pc2_histogram = pc2_dist_plot,
      scores = scores_plot
    ))
  }
}

create_color_list = function() {
  color_list = c("red" = "#D81B60", "green" = "#004D40", "yellow" = "#FFC107")
  color_list
}


create_status_heatmap = function(
  parsed_percentage_overall,
  order_n = "yes",
  trim_90 = "yes",
  use_log = "no",
  ...
) {
  # tar_load(parsed_percentage_overall); order_n = "yes"; trim_90 = "yes"; use_log = "no"
  color_list = create_color_list()
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
      `Log10(N)` = ComplexHeatmap::anno_barplot(log10(row_counts)),
      annotation_width = unit(3, "cm")
    )
  } else {
    row_ha = ComplexHeatmap::rowAnnotation(
      N = ComplexHeatmap::anno_barplot(row_counts),
      annotation_width = unit(3, "cm")
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
    column_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    ...
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

create_total_vs_tr_rsd_plot = function(
  sets_qc_confounds_rsd,
  by_measurement = "no"
) {
  rsd_measurement = sets_qc_confounds_rsd |>
    dplyr::select(
      internal_id,
      measurement_type,
      total_rsd_median,
      technical_replicate_rsd_median
    ) |>
    dplyr::transmute(
      `Total RSD` = total_rsd_median,
      `Technical Replicate RSD` = technical_replicate_rsd_median,
      internal_id = internal_id,
      measurement_type = measurement_type
    ) |>
    tidyr::pivot_longer(
      c(-internal_id, -measurement_type),
      values_to = "RSD",
      names_to = "Type"
    )
  if (by_measurement %in% "no") {
    rsd_plot = rsd_measurement |>
      dplyr::filter(RSD <= 1) |>
      ggplot(aes(x = RSD, color = Type)) +
      geom_density(position = "identity", linewidth = 2) +
      theme(legend.position = "top")
  } else {
    rsd_plot = rsd_measurement |>
      dplyr::filter(RSD <= 1) |>
      ggplot(aes(x = RSD, color = Type)) +
      geom_density(linewidth = 2) +
      theme(legend.position = "top") +
      facet_wrap(~measurement_type, ncol = 1)
  }

  rsd_plot
}

calculate_mode = function(values) {
  values = values[!is.na(values)]
  density_estimate = stats::density(values)
  mode_value = density_estimate$x[which.max(density_estimate$y)]
  mode_value
}

calculate_rsd_modes = function(sets_qc_confounds_rsd, measurement = "no") {
  if (measurement %in% "no") {
    rsd_modes = sets_qc_confounds_rsd |>
      dplyr::summarise(
        tr_mode = calculate_mode(technical_replicate_rsd_median),
        total_mode = calculate_mode(total_rsd_median)
      )
  } else if (measurement %in% "yes") {
    rsd_modes = sets_qc_confounds_rsd |>
      dplyr::summarise(
        tr_mode = calculate_mode(technical_replicate_rsd_median),
        total_mode = calculate_mode(total_rsd_median),
        .by = measurement_type
      )
  }

  rsd_modes
}

create_non_technical_rsd_plot = function(non_technical_rsd) {
  # tar_load(non_technical_rsd)
  wide_rsd = non_technical_rsd |>
    dplyr::filter_out(is.na(non_technical_rsd) | total_rsd_median > 1) |>
    dplyr::transmute(
      internal_id = internal_id,
      `Non Technical RSD` = non_technical_rsd,
      `Total RSD` = total_rsd_median,
      `Technical Replicate RSD` = technical_replicate_rsd_median
    ) |>
    tidyr::pivot_longer(-internal_id, values_to = "RSD", names_to = "Type") |>
    dplyr::mutate(
      Type = factor(
        Type,
        levels = c("Non Technical RSD", "Technical Replicate RSD", "Total RSD")
      )
    )
  wide_rsd_plot = wide_rsd |>
    ggplot(aes(x = RSD, y = Type)) +
    geom_sina() +
    geom_boxplot(alpha = 0)

  wide_rsd_plot
}


create_percent_missingness_plot = function(not_gray_pca) {
  # tar_load(not_gray_pca)
  nmr_max = not_gray_pca |>
    dplyr::filter(measurement_type %in% "NMR") |>
    dplyr::pull(percent_missing) |>
    max()
  ms_max = not_gray_pca |>
    dplyr::filter(measurement_type %in% "MS") |>
    dplyr::pull(percent_missing) |>
    max()

  max_subtitle = glue::glue(
    "Max % Missing, NMR: {round(nmr_max, digits = 0)}%, MS: {round(ms_max, digits = 0)}%"
  )
  p_miss_plot = not_gray_pca |>
    ggplot(aes(x = percent_missing, y = measurement_type)) +
    geom_sina() +
    labs(x = "% Missingness", y = "Measurement Type", subtitle = max_subtitle)
  p_miss_plot
}

create_example_raw_log_plot = function() {
  # se_data = tar_read(se_pca_g01)[["AN002054"]]
  values = tibble::tibble(
    id = seq(1, 10000),
    raw = rlnorm(10000, meanlog = 2, sdlog = 1.5)
  ) |>
    dplyr::mutate(log = log10(raw)) |>
    tidyr::pivot_longer(
      cols = c("raw", "log"),
      names_to = "transform",
      values_to = "value"
    ) |>
    dplyr::mutate(transform = factor(transform, levels = c("raw", "log")))
  summary = values |>
    dplyr::summarise(
      mean = mean(value),
      median = median(value),
      .by = transform
    ) |>
    tidyr::pivot_longer(-transform, names_to = "statistic", values_to = "value")

  out_plot = values |>
    dplyr::filter(value <= 250) |>
    ggplot(aes(x = value)) +
    geom_histogram(bins = 100) +
    geom_vline(data = summary, aes(xintercept = value, color = statistic)) +
    facet_wrap(~transform, nrow = 1, scales = "free") +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    theme(
      strip.background = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.2, 0.6)
    )
  out_plot
}


create_databased_cutoffs_plot = function(
  sets_percent_missing_result,
  sets_median_correlation_result,
  sets_rsd_result,
  sets_pc_number_result,
  data_limits
) {
  # tar_load(c(
  #   sets_percent_missing_result,
  #   sets_median_correlation_result,
  #   sets_rsd_result,
  #   sets_pc_number_result,
  #   data_limits
  # ))
  color_list = c("red" = "#D81B60", "green" = "#004D40", "yellow" = "#FFC107")
  pm_plot = sets_percent_missing_result |>
    dplyr::filter_out(is.na(missing) | missing == 0) |>
    ggplot(aes(x = missing)) +
    geom_histogram(bins = 100) +
    geom_vline(
      xintercept = data_limits$missing,
      color = color_list[c("yellow", "red")],
      linewidth = 2
    ) +
    geom_vline(xintercept = 0, color = color_list[c("green")], linewidth = 2) +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    scale_x_continuous(expand = expansion(mult = c(0, NA))) +
    labs(x = "% Missingness")

  rsd_plot = sets_rsd_result |>
    dplyr::filter_out(data %in% "n", (is.na(median) | median == 0)) |>
    ggplot(aes(x = median)) +
    geom_histogram(bins = 100) +
    geom_vline(
      xintercept = data_limits$rsd,
      color = color_list[c("yellow", "red")],
      linewidth = 2
    ) +
    geom_vline(
      aes(xintercept = median(median)),
      color = color_list["green"],
      linewidth = 2
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    scale_x_continuous(expand = expansion(mult = c(0, NA)), limits = c(0, 3)) +
    labs(x = "Median(RSD)")
  med_cor_plot = sets_median_correlation_result |>
    dplyr::summarise(
      max = max(median_ici_correlation, na.rm = TRUE),
      .by = id
    ) |>
    dplyr::filter(!is.infinite(max)) |>
    ggplot(aes(x = max)) +
    geom_histogram(bins = 100) +
    geom_vline(
      xintercept = data_limits$median_correlation,
      color = color_list[c("yellow", "red")],
      linewidth = 2
    ) +
    geom_vline(
      aes(xintercept = median(max)),
      color = color_list["green"],
      linewidth = 2
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    scale_x_continuous(expand = expansion(mult = c(0, NA))) +
    labs(x = "Max(Median(Sample-Sample Correlation))")
  pc_90_plot = sets_pc_number_result |>
    dplyr::mutate(pc_log = log10(pc_90)) |>
    dplyr::filter(!is.infinite(pc_90), pc_log >= 1) |>
    ggplot(aes(x = pc_log)) +
    geom_histogram(bins = 30) +
    geom_vline(
      xintercept = log10(data_limits$pca),
      color = color_list[c("yellow", "red")],
      linewidth = 2
    ) +
    geom_vline(xintercept = 1, color = color_list["green"], linewidth = 2) +
    scale_y_continuous(expand = expansion(mult = c(0, NA))) +
    scale_x_continuous(expand = expansion(mult = c(0, NA))) +
    labs(x = "Log10(PC-90%)")
  return(list(
    pm = pm_plot,
    rsd = rsd_plot,
    med_cor = med_cor_plot,
    pc_90 = pc_90_plot
  ))
}


create_count_nmostcommon_barplot = function(sets_qc_color_confounds) {
  # tar_load(sets_qc_color_confounds); tar_source("R")
  n_sample_statuses = sets_qc_color_confounds |>
    dplyr::filter_out(overall_status %in% "gray") |>
    dplyr::filter_out(failed_steps %in% "create_se") |>
    dplyr::filter_out(is.na(most_common_n)) |>
    dplyr::select(n_sample_status, most_common_n)
  count_most_common = rle_2_df(sort(n_sample_statuses$most_common_n) |> rle())
  count_most_common = count_most_common |>
    dplyr::filter_out(value > 25)

  color_list = create_color_list()
  color_df = tibble::tibble(status = names(color_list), color = color_list)

  count_most_common = count_most_common |>
    dplyr::mutate(
      status = dplyr::case_when(
        value < 3 ~ "red",
        value >= 3 & value < 5 ~ "yellow",
        TRUE ~ "green"
      )
    )

  out_plot = count_most_common |>
    ggplot(aes(x = value, y = length, fill = status)) +
    geom_col() +
    scale_fill_manual(values = color_list) +
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(c(0, NA))) +
    labs(x = "Most Frequent N-Sample / SSF Group", y = "Count") +
    scale_x_continuous(expand = expansion(c(0, NA)))
  out_plot
}


find_pca_examples = function(sets_qc_color_confounds) {
  pca_green = sets_qc_color_confounds |>
    dplyr::filter(overall_status %in% "green") |>
    dplyr::filter(most_common_n > 10, pca_anova_pvalue <= 1e-8, n_ssf_n < 5)
  pca_red = sets_qc_color_confounds |>
    dplyr::filter(overall_status %in% "red") |>
    dplyr::filter(most_common_n > 10, pca_anova_pvalue >= 0.5, n_ssf_n < 5)
}


create_pca_example_plots = function(se_data, out_file = NULL, ...) {
  # se_data = tar_read(se_pca_g01)[["AN002342"]]; tar_source("R")
  pca_data_use = get_last_pca(se_data)
  id = metadata(se_data)$mwtab$ID

  pca_scores = pca_data_use$pca$pca
  pca_pvalue = pca_data_use$pca$anova |>
    dplyr::filter(variable %in% "ssf", PC %in% "PC1") |>
    dplyr::pull(padjust) |>
    format(digits = 3)

  out_plot = pca_scores |>
    ggplot(aes(x = PC1, y = PC2, color = ssf)) +
    geom_point(size = 2) +
    labs(subtitle = paste0(id, ", SSF Adjusted P-Value: ", pca_pvalue)) +
    theme(legend.position = "top")

  if (is.null(out_file)) {
    return(out_plot)
  } else {
    ragg::agg_png(filename = out_file, ...)
    print(out_plot)
    dev.off()
    return(out_file)
  }
}


get_last_pca = function(se_data) {
  pca_data = get_specific_metadata(se_data, "pca_groups")
  last_run = check_qc_outputs(se_data)
  pca_index = last_run$pca[length(last_run$pca)]
  pca_data_use = pca_data[[pca_index]]

  pca_data_use
}
