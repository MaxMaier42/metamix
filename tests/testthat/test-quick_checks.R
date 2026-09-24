# Quick recovery and model-selection checks with small K and short chains.
# They catch gross errors (e.g. a wrong selection likelihood) within a few
# minutes; the precise recovery tests are in test-slow_*.R. Skipped on CRAN.

test_that("sel_mix corrects strong one-sided publication bias", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")

  # true mean 0, but almost only significant studies are published
  set.seed(1)
  dat <- sim_mix(K = 60, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(0.05, 0.2, 1), one_sided = TRUE)
  fit_re  <- suppressWarnings(re_mix(dat$y, dat$sds, M = 1,
                                     chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_sel <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975),
                                      chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  mu_re  <- rstan::summary(fit_re,  pars = "mu", probs = c(0.025, 0.975))$summary
  mu_sel <- rstan::summary(fit_sel, pars = "mu", probs = c(0.025, 0.975))$summary

  # the uncorrected estimate is clearly biased upwards
  expect_gt(mu_re[, "2.5%"], 0.1)
  # the corrected estimate covers the truth and is within 3 posterior SDs of it
  expect_lt(mu_sel[, "2.5%"], 0)
  expect_gt(mu_sel[, "97.5%"], 0)
  expect_lt(abs(mu_sel[, "mean"]), 3 * mu_sel[, "sd"])
  # non-significant studies are estimated to be much less likely published
  omega <- rstan::summary(fit_sel, pars = "omega")$summary[, "mean"]
  expect_lt(omega[1], 0.5)

  # bridge sampling prefers the selection model
  p <- model_probs(list(re_mix = fit_re, sel_mix = fit_sel))
  expect_gt(p[["sel_mix"]], 0.9)
})

test_that("sel_mix corrects two-sided publication bias", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")

  # true mean 0 and tau 0.1; two-sided selection inflates the heterogeneity
  set.seed(2)
  dat <- sim_mix(K = 60, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(0.05, 0.2, 1), one_sided = FALSE)
  fit_re  <- suppressWarnings(re_mix(dat$y, dat$sds, M = 1,
                                     chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_sel <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975), one_sided = FALSE,
                                      chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  tau_re  <- rstan::summary(fit_re,  pars = "tau", probs = c(0.025, 0.975))$summary
  tau_sel <- rstan::summary(fit_sel, pars = "tau", probs = c(0.025, 0.975))$summary
  mu_sel  <- rstan::summary(fit_sel, pars = "mu",  probs = c(0.025, 0.975))$summary

  expect_gt(tau_re[, "2.5%"], 0.1)
  expect_lt(tau_sel[, "2.5%"], 0.1)
  expect_gt(tau_sel[, "97.5%"], 0.1)
  expect_lt(abs(mu_sel[, "mean"]), 3 * mu_sel[, "sd"])

  p <- model_probs(list(re_mix = fit_re, sel_mix = fit_sel))
  expect_gt(p[["sel_mix"]], 0.9)
})

test_that("a two-component mixture is recovered under publication bias", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")

  set.seed(3)
  dat <- sim_mix(K = 60, M = 2, mu = c(0, 1), tau = c(0.1, 0.1), steps = c(0.95, 0.975),
                 weights = c(0.05, 0.2, 1), one_sided = TRUE)
  fit_re2  <- suppressWarnings(re_mix(dat$y, dat$sds, M = 2,
                                      chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_sel1 <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975),
                                       chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_sel2 <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 2, steps = c(0.95, 0.975),
                                       chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  mu <- rstan::summary(fit_sel2, pars = "mu", probs = c(0.025, 0.975))$summary

  # both component means are within 3 posterior SDs of the truth
  expect_lt(abs(mu[1, "mean"] - 0), 3 * mu[1, "sd"])
  expect_lt(abs(mu[2, "mean"] - 1), 3 * mu[2, "sd"])
  # the null component is under-represented among published studies, so the
  # selection-corrected mixing weight of that component is larger
  theta <- rstan::summary(fit_sel2, pars = c("theta", "theta_preselection"))$summary[, "mean"]
  expect_gt(theta[["theta_preselection[1]"]], theta[["theta[1]"]])

  p <- model_probs(list(re_mix2 = fit_re2, sel_mix1 = fit_sel1, sel_mix2 = fit_sel2))
  expect_equal(names(which.max(p)), "sel_mix2")
})

test_that("without publication bias the simplest model is preferred", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")

  set.seed(4)
  dat <- sim_mix(K = 60, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(1, 1, 1), one_sided = TRUE)
  fit_re1  <- suppressWarnings(re_mix(dat$y, dat$sds, M = 1,
                                      chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_re2  <- suppressWarnings(re_mix(dat$y, dat$sds, M = 2,
                                      chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_sel1 <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975),
                                       chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  mu <- rstan::summary(fit_re1, pars = "mu", probs = c(0.025, 0.975))$summary
  expect_lt(mu[, "2.5%"], 0.3)
  expect_gt(mu[, "97.5%"], 0.3)

  p <- model_probs(list(re_mix1 = fit_re1, re_mix2 = fit_re2, sel_mix1 = fit_sel1))
  expect_equal(names(which.max(p)), "re_mix1")
})

test_that("sel_flexpb with an unselected group corrects bias in the other group", {
  skip_on_cran()
  skip_if_not_installed("bridgesampling")

  # group 1: published literature with strong selection; group 2: unselected studies
  set.seed(5)
  lit <- sim_mix(K = 30, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(0.05, 0.2, 1), one_sided = TRUE)
  unsel <- sim_mix(K = 30, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                   weights = c(1, 1, 1), one_sided = TRUE)
  y   <- c(lit$y, unsel$y)
  sd  <- c(lit$sds, unsel$sds)
  grp <- rep(1:2, each = 30)
  fit_re   <- suppressWarnings(re_mix(y, sd, M = 1,
                                      chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  fit_flex <- suppressWarnings(sel_flexpb(y, sd, grp, M = 1, steps = c(0.95, 0.975), unselected = 2,
                                          chains = 2, iter = 1500, cores = 2, seed = 1, refresh = 0))
  mu_re   <- rstan::summary(fit_re,   pars = "mu", probs = c(0.025, 0.975))$summary
  mu_flex <- rstan::summary(fit_flex, pars = "mu", probs = c(0.025, 0.975))$summary

  expect_gt(mu_re[, "2.5%"], 0)
  expect_lt(mu_flex[, "2.5%"], 0)
  expect_gt(mu_flex[, "97.5%"], 0)
  omega <- rstan::extract(fit_flex, pars = "omega")$omega
  expect_true(all(omega[, 2, ] == 1))
  expect_lt(mean(omega[, 1, 1]), 0.5)

  p <- model_probs(list(re_mix = fit_re, sel_flexpb = fit_flex))
  expect_gt(p[["sel_flexpb"]], 0.9)
})
