# Tiny fits (20 studies, one short chain) that check the structure of the
# returned stanfit objects, not the quality of the estimates.

set.seed(1)
dat <- sim_mix(K = 20, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
               weights = c(0.5, 0.8, 1), one_sided = TRUE)

test_that("re_mix returns a stanfit with the documented parameters", {
  fit <- suppressWarnings(re_mix(dat$y, dat$sds, M = 2, chains = 1, iter = 300, seed = 1, refresh = 0))
  expect_s4_class(fit, "stanfit")
  draws <- rstan::extract(fit)
  n <- nrow(draws$mu)

  expect_equal(dim(draws$mu), c(n, 2))
  expect_equal(dim(draws$tau), c(n, 2))
  expect_equal(dim(draws$theta), c(n, 2))
  expect_equal(dim(draws$posterior_probs), c(n, 20, 2))
  expect_equal(dim(draws$y_rep), c(n, 20))
  expect_equal(dim(draws$sd_rep), c(n, 20))

  expect_true(all(draws$mu[, 1] <= draws$mu[, 2]))
  expect_true(all(draws$tau > 0))
  expect_equal(rowSums(draws$theta), rep(1, n))
  expect_true(all(abs(apply(draws$posterior_probs, c(1, 2), sum) - 1) < 1e-8))
  expect_equal(draws$log_tau, log(draws$tau))
  expect_true(all(is.finite(draws$y_rep)))
  expect_true(all(draws$sd_rep > 0))
})

test_that("sel_mix returns the selection parameters", {
  fit <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975),
                                  chains = 1, iter = 300, seed = 1, refresh = 0))
  draws <- rstan::extract(fit)
  n <- nrow(draws$omega)

  expect_equal(dim(draws$omega), c(n, 3))
  expect_true(all(abs(draws$omega[, 3] - 1) < 1e-8))
  expect_true(all(draws$omega[, 1] <= draws$omega[, 2]))
  expect_true(all(draws$avg_omega > 0 & draws$avg_omega <= 1))
  expect_equal(rowSums(draws$theta_preselection), rep(1, n))
  expect_equal(dim(draws$y_rep), c(n, 20))
  expect_true(all(is.finite(draws$y_rep)))

  # one step and two-sided selection give two weights
  fit <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = 0.975, one_sided = FALSE,
                                  chains = 1, iter = 300, seed = 1, refresh = 0))
  omega <- rstan::extract(fit, pars = "omega")$omega
  expect_equal(ncol(omega), 2)
  expect_true(all(abs(omega[, 2] - 1) < 1e-8))
})

test_that("sel_flexpb returns one selection model per group", {
  grp <- rep(1:2, each = 10)
  fit <- suppressWarnings(sel_flexpb(dat$y, dat$sds, grp, M = 1, steps = c(0.95, 0.975),
                                     chains = 1, iter = 300, seed = 1, refresh = 0))
  omega <- rstan::extract(fit, pars = "omega")$omega   # draws x G x 3
  expect_equal(dim(omega)[2:3], c(2, 3))
  expect_true(all(abs(omega[, , 3] - 1) < 1e-8))
  expect_false(all(omega[, 2, 1] == 1))

  # with group 2 declared unselected its weights are exactly 1
  fit <- suppressWarnings(sel_flexpb(dat$y, dat$sds, grp, M = 1, steps = c(0.95, 0.975), unselected = 2,
                                     chains = 1, iter = 300, seed = 1, refresh = 0))
  omega <- rstan::extract(fit, pars = "omega")$omega
  expect_true(all(omega[, 2, ] == 1))
  expect_false(all(omega[, 1, 1] == 1))
})
