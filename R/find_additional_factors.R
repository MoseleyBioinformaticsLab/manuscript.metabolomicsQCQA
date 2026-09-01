find_additional_factors = function(parsed_grouped_statuses) {
  # tar_load(parsed_grouped_statuses); tar_source("R")
  all_status_cols = setup_statuses(names(parsed_grouped_statuses))

  interesting_status_df = tibble::tribble(
    ~n_sample_status , ~percent_missing_status , ~left_censorship_test_status , ~max_median_correlation_status , ~rsd_status , ~pca_pc_index_90_status , ~pca_anova_pvalue_status ,
    "green"          , "green"                 , "green"                      , "green"                        , "green"     , "green"                 , "green"                  ,
    "yellow"         , "green"                 , "green"                      , "green"                        , "green"     , "green"                 , "green"                  ,
    "green"          , "green"                 , "green"                      , "green"                        , "green"     , "green"                 , "red"                    ,
    "yellow"         , "green"                 , "green"                      , "green"                        , "green"     , "green"                 , "red"                    ,
    "green"          , "green"                 , "green"                      , "green"                        , "green"     , "yellow"                , "yellow"
  )

  interesting_status_indices = vector("list", nrow(interesting_status_df))
  for (irow in seq_len(nrow(interesting_status_df))) {
    interesting_status_indices[[irow]] = find_matching_status_index(
      parsed_grouped_statuses,
      interesting_status_df[irow, ]
    )
  }
  all_indices = unique(unlist(interesting_status_indices))

  all_matches = parsed_grouped_statuses[all_indices, "match"] |>
    tidyr::unnest(cols = "match") |>
    dplyr::filter(most_common_n >= 10) |>
    dplyr::arrange(dplyr::desc(most_common_n))
  all_matches
}

find_matching_status_index = function(keep_status, status_definition) {
  # keep_status = tar_read(not_gray_pca); status_definition = count_statuses[1, ]; tar_source("R")
  match_status = rep(TRUE, nrow(keep_status))

  use_cols = grep(
    "N",
    names(status_definition),
    ignore.case = FALSE,
    invert = TRUE,
    value = TRUE
  )
  for (icol in use_cols) {
    match_status = match_status &
      (keep_status[[icol]] %in% status_definition[[icol]])
  }

  which(match_status)
}


setup_statuses = function(status_names) {
  status_cols = grep("status$", status_names, value = TRUE)
  status_cols = grep(
    "^summarized|^has_qc|^score_icikt|^score_rle|^ssf_metadata|^overall|^technical_replicate",
    status_cols,
    invert = TRUE,
    value = TRUE
  )
  status_cols
}
