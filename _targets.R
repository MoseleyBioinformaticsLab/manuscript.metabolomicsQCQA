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
    .progress = TRUE
  ) |>
    tar_target()

  check_result = purrr::map(checked, \(x) {
    tibble::as_tibble(x)
  }) |>
    purrr::list_rbind() |>
    tar_target()
})
# json_all_map = tar_map(
#   mwtab_datasets[1:200, ],
#   names = id,
#   tar_target(dataset, file, format = "file"), #1
#   tar_target(processed, parse_json(dataset, id, ancillary_path)), #2
#   tar_target(
#     checked,
#     run_mwtab_checks_json(
#       processed,
#       min_n = 3,
#       min_metabolites = 50,
#       max_min_value = 20,
#       min_ssf = 1,
#       use_ssf_only = "yes"
#     )
#   ), #3
#   tar_target(check_result, get_check(checked)) #4
# )

# check_comb = tar_combine(
#   check_data,
#   json_all_map[[4]],
#   command = bind_rows(!!!.x)
# )

list(json_all_map)
