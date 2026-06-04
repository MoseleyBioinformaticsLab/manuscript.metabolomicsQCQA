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
processed_t = tar_target(processed, process_mwtab(mwtab_datasets))
processed_group_t = tar_target(
  processed_group,
  bind_group_processed(processed)
)

report_file_t = tar_file(report_file, "resources/example_report.qmd")

sets_map = tar_map(
  group_sets,
  names = group,
  tar_target(setup, sets_run_setup(set, processed_group))
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

reproducibility_t = tar_quarto(
  reproducibility_report,
  "docs/reproducibility.qmd"
)
# setup_map = tar_map(
#   json_all_map[[4]],
#   tar_target(setup, run_setup)
# )

list(
  processed_t,
  processed_group_t,
  sets_map
)
