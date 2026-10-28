
d <- read.csv(
  "data/Syrian_Inflation_1961_2025.csv",
  check.names = FALSE
)

stopifnot(
  nrow(d) == 65L,
  identical(names(d),
    c("YEAR", "EXR", "TRDGDP", "INF", "M2GDP")),
  identical(as.integer(d$YEAR), 1961:2025),
  all(is.finite(as.matrix(d))),
  all(d$EXR > 0)
)

source("R/agbqr_core.R")

stopifnot(
  is.function(find_qr),
  is.function(hac_components),
  is.function(calibrate_eta),
  is.function(sample_gibbs)
)

message("Smoke test passed")
