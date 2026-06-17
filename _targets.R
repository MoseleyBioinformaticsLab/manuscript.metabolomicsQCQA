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

json_files = fs::dir_ls(
  "mwtab/repaired/json_2025-11-12"
)

mwtab_datasets = extract_mwtab_ids(json_files)

ancillary_path = "mwtab/ancillary/"
save_path = "parsed_rds"

# library(mirai)
# daemons(3)
library(furrr)
plan(multicore, workers = 10)


# first time running, go to **processed_group** so the group file
# gets created properly.
processed_t = tar_target(
  processed,
  process_mwtab(mwtab_datasets, ancillary_path, save_path)
)
processed_group_t = tar_target(
  processed_group,
  bind_group_processed(processed)
)

report_file_t = tar_file(report_file, "resources/example_report.qmd")

sets_map = tar_map(
  group_sets,
  names = group,
  tar_target(se, sets_readrds_convert_to_se(set, processed_group)),
  tar_target(se_mmr, sets_check_mmr_negatives(se)),
  tar_target(se_qc, get_qc_outputs(se)),
  tar_target(neg, sets_run_qc_neg(se)),
  tar_target(run, sets_run_qc(se))
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
  plot_neg_tests(sets_neg_test)
)

manuscript_t = tar_quarto(
  manuscript,
  "docs/manuscript.qmd",
  output_file = "docs/qcqa_manuscript.docx"
)

reproducibility_t = tar_quarto(
  reproducibility_report,
  "docs/reproducibility.qmd"
)

contains_negatives_report_t = tar_quarto(
  contains_negatives_report,
  "docs/contains_negatives.qmd"
)

negatives_testing_report_t = tar_quarto(
  negatives_testing_report,
  "docs/effect_of_negatives_as_missing.qmd"
)

mean_median_ratios_report_t = tar_quarto(
  mean_median_ratios_report,
  "docs/mean_median_ratios_for_logdetection.qmd"
)
# setup_map = tar_map(
#   json_all_map[[4]],
#   tar_target(setup, run_setup)
# )

list(
  processed_t,
  processed_group_t,
  sets_map,
  se_mmr_combine_t,
  se_qc_combine_t,
  sets_neg_t,
  sets_neg_test_t,
  sets_neg_test_result_t,
  contains_negatives_report_t,
  negatives_testing_report_t,
  mean_median_ratios_report_t
)
