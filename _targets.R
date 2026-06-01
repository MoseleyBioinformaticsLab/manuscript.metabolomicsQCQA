## Load your packages and source files e.g. library(targets).
targets::tar_source(c("./packages.R", "R"))

# todo:
#  - in check_dynamic_range, it should check for 99% or more missing values

if (!fs::file_exists("parsed_rds/group_sets.rds")) {
  group_sets = tibble::tibble(
    group = paste0(
      "g",
      stringr::str_pad(seq(1, 20), width = 2, side = "left", pad = "0")
    ),
    set = ""
  )
} else {
  group_sets = readRDS("parsed_rds/group_sets.rds")
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


## tar_plan supports drake-style targets and also tar_target()
processed_t = tar_target(processed, process_mwtab(mwtab_datasets))
checked_t = tar_target(checked, check_mwtab(processed))

check_result_t = tar_target(check_result, bind_check_results(checked))

# first time running, only go to here so the right file is saved.
# then run the other targets below.
group_sets_t = tar_target(
  tmp_sets,
  filter_good_experiments_splitup(check_result, "parsed_rds/group_sets.rds")
)

report_file_t = tar_file(report_file, "resources/example_report.qmd")

sets_map = tar_map(
  group_sets,
  names = group,
  tar_target(setup, sets_run_setup(set, tmp_sets)),
  tar_target(check_features, sets_check_features(setup)),
  tar_target(cor, sets_run_group_correlations(setup)),
  tar_target(lc, sets_add_mwtab_do_left_censored(cor, setup)),
  tar_target(pca, sets_run_pca_experimental_design(lc)),
  tar_target(group_scores, sets_extract_grouping_scores(pca)),
  tar_target(
    reports,
    sets_create_qcqa_reports(pca, "qcqa_reports", report_file)
  )
)

check_features_all_t = tar_combine(
  check_features_all,
  sets_map[[2]],
  command = dplyr::bind_rows(!!!.x)
)

group_scores_all_t = tar_combine(
  group_scores_all,
  sets_map[[6]],
  command = dplyr::bind_rows(!!!.x)
)

group_scores_plot_t = tar_target(
  group_scores_plots,
  plot_group_scores(group_scores_all)
)

manuscript_t = tar_quarto(
  manuscript,
  "docs/manuscript.qmd",
  output_file = "docs/qcqa_manuscript.docx"
)
# setup_map = tar_map(
#   json_all_map[[4]],
#   tar_target(setup, run_setup)
# )

list(
  processed_t,
  checked_t,
  check_result_t,
  group_sets_t,
  sets_map,
  check_features_all_t,
  group_scores_all_t,
  group_scores_plot_t,
  report_file_t
)
