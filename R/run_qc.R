sets_run_dyn_neg = function(se_set) {
  # setup_set = tar_read(setup_g01); tar_source("R")
  # setup_set = tar_read(setup_g02); tar_source("R")

  # then we do dynamic range calculations for those that can.
  dyn_range_set = purrr::map(se_set, calculate_dynamic_range_check_neg)
}

sets_run_qc_neg = function(se_set) {
  # setup_set = tar_read(setup_g01); tar_source("R")
  # setup_set = tar_read(setup_g02); tar_source("R")

  # then we do dynamic range calculations for those that can.
  dyn_range_set = purrr::map(se_set, calculate_dynamic_range)
  qc_info = get_qc_outputs(dyn_range_set)

  neg_qc = qc_info |>
    dplyr::filter(why %in% "negative values present") |>
    dplyr::pull(id) |>
    unique()

  neg_list = dyn_range_set[neg_qc]
  neg_list
}

get_qc = function(se_data) {
  qc_metadata = metadata(se_data)$qc
  transform_metadata = metadata(se_data)$transform
  if (inherits(transform_metadata$transform, "function")) {
    transform_metadata$transform = NULL
  }
  c(qc_metadata, transform_metadata)
}


get_qc_outputs = function(run_set) {
  # run_set = tar_read(run_g01); tar_source("R")
  all_qc = purrr::map(run_set, \(x) {
    as.data.frame(get_qc(x))
  }) |>
    dplyr::bind_rows()
  all_qc
}

check_qc_se_1_3 = function() {
  qc_run1 = tar_read(neg_g01) |> get_qc_outputs()
  qc_run2 = tar_read(neg_g02) |> get_qc_outputs()
  qc_run3 = tar_read(neg_g03) |> get_qc_outputs()
  all_qc = dplyr::bind_rows(qc_run1, qc_run2, qc_run3)
  n_id = length(unique(all_qc$id))
  all_qc |>
    dplyr::summarise(N = dplyr::n(), .by = why) |>
    dplyr::mutate(perc = N / n_id * 100) |>
    paint::unpaint()
}


run_qc = function(in_data) {
  # the idea here is going to be:
  #   - check if we can create the SE object;
  #   - create the SE object;
  #   - set all the possible qc steps to FALSE;
  #   - try the first one, if it succeeds or returns minimum for next;
  #   - do the next one;
  #   - and so on.
  #   - if at any point you get a NO on the next one, we exit and return the object
  #     with what was actually done.
  #   - now, the NO will have to be really bad to stop it.
  #   If we have 10 features and 3 samples, we are going to go for it.
  #   But the first check is absolutely going to be, for each sample, are there negative values, or missing values, and what percentage, so we can exclude it right there and then.
}
