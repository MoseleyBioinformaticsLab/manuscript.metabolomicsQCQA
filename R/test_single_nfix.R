test_single_n_fix = function() {
  sets_n1_hit = readRDS("resources/sets_qc_color_n1_ssf_2026-08-24.rds")
  group_sets = readRDS("resources/group_sets.rds")
  tar_load(mwtab_raw_files)
  tar_load(processed_group)
  tar_source("R")
  library(furrr)
  plan(multicore, workers = 10)

  all_sets = purrr::map(group_sets$set, \(x) {
    x
  }) |>
    purrr::list_rbind()
  test_datasets = all_sets |>
    dplyr::filter(id %in% sets_n1_hit$internal_id) |>
    dplyr::mutate(type = "n1")
  other_datasets = all_sets |>
    dplyr::filter(!id %in% sets_n1_hit$internal_id) |>
    dplyr::slice_sample(n = 20) |>
    dplyr::mutate(type = "other")

  run_datasets = dplyr::bind_rows(test_datasets, other_datasets)

  all_files = run_datasets$file
  all_md5 = run_datasets$md5
  run_out = purrr::map2(
    all_files,
    all_md5,
    readrds_convert_to_se,
    mwtab_raw_files
  )
}
