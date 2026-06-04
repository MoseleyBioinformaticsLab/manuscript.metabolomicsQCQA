create_group_scores_comparison_plot = function(
  group_scores,
  use_sample_set = "all"
) {
  # group_scores = tar_read(group_scores_all); use_sample_set = "all"
  use_scores = group_scores |>
    dplyr::filter(factors %in% "all", sample_set %in% use_sample_set)

  compare_scores = use_scores |>
    tidyr::pivot_wider(
      id_cols = "analysis_id",
      values_from = "differences_fraction",
      names_from = "method"
    )
  compare_scores |>
    ggplot(aes(x = `median-ici`, y = rle)) +
    geom_abline(slope = 1, color = "blue") +
    geom_point(size = 2, alpha = 0.5)
}

create_group_scores_histogram = function(
  group_scores,
  use_sample_set = "all",
  which = "median-ici"
) {
  # group_scores = tar_read(group_scores_all); use_sample_set = "all"; which = "median-ici"
  use_scores = group_scores |>
    dplyr::filter(factors %in% "all", sample_set %in% use_sample_set)

  compare_scores = use_scores |>
    tidyr::pivot_wider(
      id_cols = "analysis_id",
      values_from = "differences_fraction",
      names_from = "method"
    )
  compare_scores |>
    ggplot(aes(x = {{ which }})) +
    geom_histogram(bins = 100)
}
