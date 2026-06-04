find_analysis_in_group = function(analysis_id) {
  group_sets = readRDS("resources/group_sets.rds")
  all_ids = purrr::map(group_sets$set, \(in_set) {
    in_set
  }) |>
    purrr::list_rbind()

  all_ids |>
    dplyr::filter(id %in% analysis_id)
}
