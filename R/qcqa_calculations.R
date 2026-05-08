read_rds_setup_data = function(rds_file, md5_sum) {
  processed_data = readRDS(rds_file)
  check_md5(processed_data, md5_sum)
}

check_md5 = function(in_data, md5_sum) {
  new_md5 = digest::digest(in_data, algo = "md5")
  stopifnot(new_md5 == md5_sum)
  return(invisible(NULL))
}
