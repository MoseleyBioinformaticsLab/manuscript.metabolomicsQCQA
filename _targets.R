## Load your packages and source files e.g. library(targets).
targets::tar_source(c("./packages.R", "R"))

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
json_all_map = tar_assign({
  processed = furrr::future_map2(
    mwtab_datasets$file,
    mwtab_datasets$id,
    \(x, y) {
      parse_json_save(x, y, ancillary_path, save_path)
    },
    .progress = TRUE
  ) |>
    tar_target()

  checked = furrr::future_map(
    processed,
    run_mwtab_checks_json,
    min_n = 3,
    min_metabolites = 50,
    max_min_value = 20,
    min_ssf = 1,
    use_ssf_only = "yes",
    remove_pooled = "no",
    .progress = TRUE
  ) |>
    tar_target()

  check_result = purrr::map(checked, \(x) {
    tibble::as_tibble(x)
  }) |>
    purrr::list_rbind() |>
    tar_target()
  # we filter to what's useable, and then
  # split them up into 20 different items to iterate over
  # this is to try and control the number of objects
  # targets tries to iterate over when we query
  # anything about the project.
  split_sets = filter_good_experiments_splitup(check_result) |>
    tar_target()

  setup_sets = tar_target(
    name = setup,
    command = run_setup(split_sets),
    pattern = map(split_sets)
  )
})

# setup_map = tar_map(
#   json_all_map[[4]],
#   tar_target(setup, run_setup)
# )

list(json_all_map)
