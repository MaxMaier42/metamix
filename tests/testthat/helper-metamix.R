# The tests are organised in three tiers:
#   test-fast_*.R   always run, also on CRAN (prior-only draws and tiny fits, < 1 min)
#   test-quick_*.R  skipped on CRAN; quick recovery and bridge-sampling checks (a few min)
#   test-slow_*.R   only run when the environment variable METAMIX_SLOW_TESTS is "true"

skip_if_not_slow <- function() {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("METAMIX_SLOW_TESTS"), "true"),
              "set METAMIX_SLOW_TESTS = \"true\" to run the slow tests")
}

# posterior model probabilities (equal prior probabilities) from bridge sampling
model_probs <- function(fits) {
  logml <- sapply(fits, function(f) bridgesampling::bridge_sampler(f, silent = TRUE)$logml)
  p <- exp(logml - max(logml))
  p / sum(p)
}
