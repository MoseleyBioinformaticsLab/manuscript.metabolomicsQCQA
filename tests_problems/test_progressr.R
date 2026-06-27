library(progressr)
handlers("cli")
handlers(global = TRUE)
run_something = function(x) {
  Sys.sleep(0.1)
  x
}


run_a_lot = function(a_lot = seq(1, 20)) {
  p = progressr::progressor(along = a_lot)
  out_res = furrr::future_map(a_lot, \(x) {
    run_something(x)
    p(message = c("Correlating ", sprintf("N=%g", x)))
  })
}

tmp = run_a_lot()
