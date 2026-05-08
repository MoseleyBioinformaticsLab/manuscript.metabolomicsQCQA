## Load your packages and source files e.g. library(targets).
targets::tar_source(c("./packages.R", "R"))

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
processed_t = tar_target(processed, command = {
  furrr::future_map2(
    mwtab_datasets$file,
    mwtab_datasets$id,
    \(x, y) {
      parse_json_save(x, y, ancillary_path, save_path)
    },
    .progress = TRUE
  )
})
checked_t = tar_target(checked, command = {
  furrr::future_map(
    processed,
    run_mwtab_checks_json,
    min_n = 3,
    min_metabolites = 50,
    max_min_value = 20,
    min_ssf = 1,
    use_ssf_only = "yes",
    remove_pooled = "no",
    .progress = TRUE
  )
})

check_result_t = tar_target(check_result, bind_check_results(checked))

group_sets_t = tar_target(
  tmp_sets,
  filter_good_experiments_splitup(check_result, "parsed_rds/group_sets.rds")
)

setup_map = tar_map(
  group_sets,
  names = group,
  tar_target(setup, run_setup(set, check_result))
)

# setup_map = tar_map(
#   json_all_map[[4]],
#   tar_target(setup, run_setup)
# )

list(processed_t, checked_t, check_result_t, group_sets_t, setup_map)
