## Load your packages and source files e.g. library(targets).
targets::tar_source(c("./packages.R", "R"))

group_sets_file = "resources/group_sets.rds"
# todo:
#  - in check_features, it should check for 95% or more missing values instead of just no data.

if (!fs::file_exists(group_sets_file)) {
  group_sets = tibble::tibble(
    group = paste0(
      "g",
      stringr::str_pad(seq(1, 20), width = 2, side = "left", pad = "0")
    ),
    set = ""
  )
} else {
  group_sets = readRDS(group_sets_file)
  if (nrow(group_sets) == 0) {
    stop("something is wrong with group_sets.rds, delete it!")
  }
}

repaired_files = fs::dir_ls(
  "mwtab/2026-07-21/mwtab_repaired_json",
  regexp = "json"
)

repaired_datasets_t = tar_target(
  repaired_datasets,
  extract_mwtab_ids(repaired_files)
)

ancillary_files = fs::dir_ls(
  "mwtab/2026-07-21/mwtab_results_files/",
  regexp = "txt"
)
ancillary_datasets_t = tar_target(
  ancillary_datasets,
  extract_mwtab_ids(ancillary_files)
)
save_path = "parsed_rds_2026-07-21"

library(furrr)
plan(multicore, workers = 10)

# first time running, go to **processed_group** so the group file
# gets created properly.
processed_t = tar_target(
  processed,
  process_mwtab(repaired_datasets, ancillary_datasets, save_path)
)
processed_group_t = tar_target(
  processed_group,
  bind_group_processed(processed)
)

report_file_t = tar_file(report_files, fs::dir_ls("resources", regexp = "qmd"))

mwtab_ftp_file_t = tar_file(
  mwtab_ftp_file,
  "resources/mwbench_ftp_studies_files_2026-06-25.txt"
)

mwtab_raw_files_t = tar_target(
  mwtab_raw_files,
  get_mwtab_ftp_files(mwtab_ftp_file)
)

sets_map = tar_map(
  group_sets,
  names = group,
  # first we convert everything to summarized-experiment objects
  tar_target(
    se,
    sets_readrds_convert_to_se(set, processed_group, mwtab_raw_files)
  ), #1
  # we do the analysis to figure out how to detect things that should be
  # log-transformed
  tar_target(se_mmr, sets_check_mmr_negatives(se)), #2
  # this is for getting some qc data out
  tar_target(se_qc, get_qc_outputs(se)), #3
  # checking if we treat negative as missing,
  # how does that effect the left-censorship test
  tar_target(neg, sets_run_qc_neg(se)), #4
  # transforming if necessary, normalizing, calculating dynamic range,
  # and checking if any further steps can be done.
  tar_target(se_dynamic_range, sets_calculate_dynamic_range(se)), #5
  tar_target(se_pca, sets_calculate_correlations_pca(se_dynamic_range)), #6
  tar_target(
    se_negative_info,
    sets_get_negative_information(se_dynamic_range) #7
  ),
  tar_target(se_percent_missing, sets_percent_missing(se_pca)), #8
  tar_target(se_median_correlation, sets_median_correlations(se_pca)), #9,
  tar_target(se_qc_eval, sets_evaluate_qc(se_pca, data_limits)), #10
  tar_target(se_rsd, sets_rsd(se_pca)), #11
  tar_target(se_pc_number, sets_get_pc(se_pca)), #12
  tar_target(se_confounds, sets_detect_confounding_factors(se_pca)), #13
  tar_target(
    qcqa_reports,
    sets_create_qcqa_reports(se_pca, sets_qc_confounds_rsd, report_files) #14
  ),
  tar_target(
    tr_pca,
    sets_calculate_tr_experimental_rsd(se_pca, sets_qc_color_confounds) #15
  )
)

se_negative_information_t = tar_combine(
  se_negative_information,
  sets_map[[7]],
  command = dplyr::bind_rows(!!!.x)
)

se_mmr_combine_t = tar_combine(
  se_mmr_combine,
  sets_map[[2]],
  command = dplyr::bind_rows(!!!.x)
)

se_qc_combine_t = tar_combine(
  se_qc_combine,
  sets_map[[3]],
  command = dplyr::bind_rows(!!!.x)
)

sets_neg_t = tar_combine(
  sets_neg,
  sets_map[[4]],
  command = c(!!!.x)
)

sets_neg_test_t = tar_target(
  sets_neg_test,
  test_negative_value_behavior(sets_neg)
)

sets_neg_test_result_t = tar_target(
  sets_neg_test_result,
  plot_neg_tests(sets_neg_test, 100)
)

sets_percent_missing_t = tar_combine(
  sets_percent_missing_result,
  sets_map[[8]],
  command = dplyr::bind_rows(!!!.x)
)

sets_median_correlation_t = tar_combine(
  sets_median_correlation_result,
  sets_map[[9]],
  command = dplyr::bind_rows(!!!.x)
)

sets_rsd_t = tar_combine(
  sets_rsd_result,
  sets_map[[11]],
  command = dplyr::bind_rows(!!!.x)
)

sets_pc_number_t = tar_combine(
  sets_pc_number_result,
  sets_map[[12]],
  command = dplyr::bind_rows(!!!.x)
)

sets_qc_evaluation_t = tar_combine(
  sets_qc_evaluation,
  sets_map[[10]],
  command = dplyr::bind_rows(!!!.x)
)

data_limits_t = tar_target(
  data_limits,
  derive_data_limits(
    sets_percent_missing_result,
    sets_median_correlation_result,
    sets_rsd_result,
    sets_pc_number_result
  )
)

sets_qc_reports_t = tar_combine(
  sets_qc_reports,
  sets_map[[14]],
  command = dplyr::bind_rows(!!!.x)
)

sets_tr_experimental_rsd_t = tar_combine(
  sets_tr_experimental_rsd,
  sets_map[[15]],
  command = dplyr::bind_rows(!!!.x)
)

confounding_factors_t = tar_combine(
  sets_confounding_factors,
  sets_map[[13]],
  command = dplyr::bind_rows(!!!.x)
)

actual_analysis_t = tar_assign({
  sets_qc_color = set_overall_color_status(sets_qc_evaluation) |> tar_target()
  sets_qc_color_confounds = add_confounding_factor_status(
    sets_qc_color,
    sets_confounding_factors
  ) |>
    tar_target()
  sets_qc_confounds_rsd = dplyr::left_join(
    sets_qc_color_confounds,
    sets_tr_experimental_rsd,
    by = "internal_id"
  ) |>
    tar_target()
  not_gray_pca = find_pca_remove_gray(sets_qc_color_confounds) |> tar_target()
  parsed_grouped_statuses = parse_statuses(not_gray_pca) |> tar_target()
  parsed_percentage_overall = calculate_percentage_overall(
    parsed_grouped_statuses
  ) |>
    tar_target()
  technical_replicate_info = technical_replicate_info_analysis(not_gray_pca) |>
    tar_target()
  potential_confound_info = confound_info_analysis(not_gray_pca) |>
    tar_target()
})

figures_tables_t = tar_assign({
  waterfall_input = calculate_losses_at_steps(sets_qc_color) |> tar_target()
  loss_waterfall_plot = create_loss_waterfall_plot(waterfall_input) |>
    tar_target()
  nose_reasons_bar_plot = create_nose_reasons_bar_plot(sets_qc_color) |>
    tar_target()
  datasets_sizes_plot = create_dataset_sizes_plot(sets_qc_color) |>
    tar_target()
  mmr_logtransform_plot = create_mmr_logtransform_plots(se_mmr_combine) |>
    tar_target()
  technical_replicates_plot = create_technical_replicate_plot(se_pca_g01[[
    "AN006191"
  ]]) |>
    tar_target()
  example_confound_plot = create_example_confound_plot(se_pca_g34[[
    "AN000454"
  ]]) |>
    tar_target()
  status_heatmap_top90_plot = create_status_heatmap(
    parsed_percentage_overall,
    order_n = "yes",
    trim_90 = "yes",
    use_log = "no"
  ) |>
    tar_target()
  status_heatmap_all_plot = create_status_heatmap(
    parsed_percentage_overall,
    order_n = "yes",
    trim_90 = "no",
    use_log = "yes"
  ) |>
    tar_target()
  negatives_as_missing_plot = plot_neg_tests(
    sets_neg_test,
    nbin = 100,
    use_subset = "all"
  ) |>
    tar_target()

  parsed_percentage_csv = create_csv_file(
    parsed_percentage_overall,
    out_file = "outputs/parsed_percentage_overall.csv"
  ) |>
    tar_target()
  sets_qc_color_confounds_csv = create_csv_file(
    sets_qc_color_confounds,
    out_file = "outputs/sets_qc_color_status.csv"
  ) |>
    tar_target()
})

docs_t = tar_assign({
  manuscript = tar_quarto("docs/manuscript.qmd")
  reproducibility_report = tar_quarto("docs/reproducibility.qmd")
  contains_negatives_report = tar_quarto("docs/contains_negatives.qmd")
  negatives_testing_report = tar_quarto("docs/negatives_testing_as_missing.qmd")
  mean_median_ratios_report = tar_quarto(
    "docs/mean_median_ratios_for_logdetection.qmd"
  )
  setting_databased_cutoffs_report = tar_quarto(
    "docs/setting_databased_cutoffs.qmd"
  )
  iecm2026_abstract = tar_quarto("docs/iecm2026_abstract.qmd")
})

# manuscript_t = tar_quarto(
#   manuscript,
#   "docs/manuscript.qmd"
# )

# reproducibility_t = tar_quarto(
#   reproducibility_report,
#   "docs/reproducibility.qmd"
# )

# contains_negatives_report_t = tar_quarto(
#   contains_negatives_report,
#   "docs/contains_negatives.qmd"
# )

# negatives_testing_report_t = tar_quarto(
#   negatives_testing_report,
#   "docs/negatives_testing_as_missing.qmd"
# )

# mean_median_ratios_report_t = tar_quarto(
#   mean_median_ratios_report,
#   "docs/mean_median_ratios_for_logdetection.qmd"
# )

# setting_databased_cutoffs_report_t =
#   tar_quarto(setting_databased_cutoffs, "docs/setting_databased_cutoffs.qmd")

# iecm2026_abstract_t = tar_quarto(
#   iecm2026_abstract,
#   "docs/iecm2026_abstract.qmd"
# )
# setup_map = tar_map(
#   json_all_map[[4]],
#   tar_target(setup, run_setup)
# )

list(
  repaired_datasets_t,
  ancillary_datasets_t,
  mwtab_ftp_file_t,
  mwtab_raw_files_t,
  report_file_t,
  processed_t,
  processed_group_t,
  sets_map,
  se_mmr_combine_t,
  se_negative_information_t,
  se_qc_combine_t,
  sets_neg_t,
  sets_neg_test_t,
  sets_neg_test_result_t,
  sets_percent_missing_t,
  sets_median_correlation_t,
  sets_qc_evaluation_t,
  sets_rsd_t,
  sets_pc_number_t,
  data_limits_t,
  sets_qc_reports_t,
  sets_tr_experimental_rsd_t,
  confounding_factors_t,
  actual_analysis_t,
  figures_tables_t,
  docs_t
)
