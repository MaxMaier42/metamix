# All checks below fail before any Stan model is run, so these tests are instant.

y  <- c(0.1, 0.3, -0.2, 0.5, 0.0)
sd <- c(0.2, 0.2, 0.3, 0.1, 0.2)

test_that("re_mix rejects invalid data and settings", {
  expect_error(re_mix(y, sd[1:4], M = 1), "same length")
  expect_error(re_mix(c(y[-1], NA), sd, M = 1), "missing or infinite")
  expect_error(re_mix(y, c(sd[-1], 0), M = 1), "must be > 0")
  expect_error(re_mix(as.character(y), sd, M = 1), "numeric")
  expect_error(re_mix(y, sd, M = 0), "positive integer")
  expect_error(re_mix(y, sd, M = 1.5), "positive integer")
  expect_error(re_mix(y, sd, M = 1, mu_sd = 0), "mu_sd must be > 0")
  expect_error(re_mix(y, sd, M = 1, tau_sd = -1), "tau_sd must be > 0")
  expect_error(re_mix(y, sd, M = 1, tau_prior = "inv_gamma", tau_alpha = 0), "tau_alpha and tau_beta")
  expect_error(re_mix(y, sd, M = 2, mean_prior = "gap", gap_sdlog = 0), "gap_sdlog must be > 0")
  expect_error(re_mix(y, sd, M = 2, gap_min = -0.1), "gap_min must be >= 0")
  expect_error(re_mix(y, sd, M = 1, mean_prior = "other"), "arg")
})

test_that("sel_mix rejects invalid steps and data", {
  expect_error(sel_mix(y, sd, M = 1, steps = c(0.9, 0.95, 0.975)), "one or two steps")
  expect_error(sel_mix(y, sd, M = 1, steps = c(0.95, 0.9)), "strictly increasing")
  expect_error(sel_mix(y, sd, M = 1, steps = c(0.95, 1)), "between 0 and 1")
  expect_error(sel_mix(y, sd, M = 1, steps = c(0.4, 0.95), one_sided = FALSE), "larger .5")
  expect_error(sel_mix(y, sd[1:4], M = 1), "same length")
  expect_error(sel_mix(y, c(sd[-1], 0), M = 1), "must be > 0")
  expect_error(sel_mix(y, sd, M = 0), "positive integer")
  expect_error(sel_mix(y, sd, M = 1, mu_sd = 0), "mu_sd must be > 0")
})

test_that("sel_flexpb rejects invalid groups and data", {
  grp <- c(1, 1, 2, 2, 2)
  expect_error(sel_flexpb(y, sd, grp[1:4], M = 1), "same length as y")
  expect_error(sel_flexpb(y, sd, c(1, 1, 2, 2, NA), M = 1), "positive integers")
  expect_error(sel_flexpb(y, sd, c(0, 1, 2, 2, 2), M = 1), "positive integers")
  expect_error(sel_flexpb(y, sd, grp, M = 1, unselected = 3), "between 1 and G")
  expect_error(sel_flexpb(y, sd, grp, M = 1, unselected = 1:2), "All groups are marked unselected")
  expect_error(sel_flexpb(y, sd, grp, M = 1, steps = c(0.95, 0.9)), "strictly increasing")
  expect_error(sel_flexpb(y, sd[1:4], grp, M = 1), "same length")
  expect_error(sel_flexpb(y, c(sd[-1], 0), grp, M = 1), "must be > 0")
})
