calculate_losses_at_steps = function(sets_qc_color) {
  # tar_load(sets_qc_evaluation); tar_source("R")
  n_each = summarize_qc_results_table(sets_qc_color, group_none = "yes")
  n_total = nrow(sets_qc_color)

  stopped_se = n_each |>
    dplyr::filter(completed_steps %in% "none", is.na(failed_steps)) |>
    dplyr::pull(n)
  stopped_dynamic_range = n_each |>
    dplyr::filter(
      completed_steps %in% "se; dynamic_range",
      failed_steps %in% "none"
    ) |>
    dplyr::pull(n)
  stopped_correlation = n_each |>
    dplyr::filter(failed_steps %in% "correlation") |>
    dplyr::pull(n)
  gray_pca = sum(sets_qc_color$pca_anova_pvalue_status %in% "gray")

  results_df = tibble::tribble(
    ~step           , ~n                         , ~type   ,
    "Start"         , n_total                    , "total" ,
    "No SE"         , -1 * stopped_se            , "loss"  ,
    "Dynamic Range" , -1 * stopped_dynamic_range , "loss"  ,
    "Correlation"   , -1 * stopped_correlation   , "loss"  ,
    "PCA"           , -1 * gray_pca              , "loss"
  )
  results_df
}

create_loss_waterfall_plot = function(waterfall_input) {
  library(waterfall)
  loss_waterfall_plot = waterfall_input |>
    waterfall(
      calc_total = TRUE,
      total_rect_color = "lightgreen",
      total_rect_text_color = "black"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  loss_waterfall_plot
}

create_nose_reasons_bar_plot = function(sets_qc_color) {
  # tar_load(sets_q_color); tar_source("R");
  n_each = summarize_qc_results_table(sets_qc_color, group_none = "no")
  none_steps = n_each |>
    dplyr::filter(completed_steps %in% "none") |>
    dplyr::arrange(n) |>
    dplyr::mutate(reason = factor(failed_why, levels = failed_why))
  nose_histogram_plot = none_steps |>
    ggplot(aes(x = n, y = reason)) +
    geom_col() +
    labs(x = "N-Datasets", y = "Reason")
  nose_histogram_plot
}
