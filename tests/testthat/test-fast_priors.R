# Prior-only runs (no data enter the model) are almost instant and check that
# the priors implemented in Stan are the ones described in the documentation.
# Data are passed but must be ignored because prior_only = TRUE.

y  <- c(2, 2, 2)
sd <- c(0.1, 0.1, 0.1)

test_that("re_mix samples the default priors when prior_only = TRUE", {
  fit <- suppressWarnings(re_mix(y, sd, M = 2, prior_only = TRUE,
                                 chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  draws <- rstan::extract(fit, pars = c("mu", "tau", "theta"))

  # tau ~ half-normal(0, 0.2), i.e. a normal(0, 0.2) restricted to positive
  # values. The mean of a half-normal with scale s is s * sqrt(2 / pi) (pi is
  # the mathematical constant), here 0.2 * 0.798 = 0.16. With 3000 draws the
  # Monte Carlo error of the mean is about 0.003, so 0.02 is a safe margin.
  expect_lt(abs(mean(draws$tau) - 0.2 * sqrt(2 / pi)), 0.02)

  # The two means have independent normal(0, 1) priors but are ordered, so
  # mu[1] is the smaller and mu[2] the larger of two standard normal draws.
  # The expected minimum of two standard normals is -1 / sqrt(pi) = -0.56 and
  # the expected maximum is +0.56.
  expect_true(all(draws$mu[, 1] <= draws$mu[, 2]))
  expect_lt(abs(mean(draws$mu[, 1]) + 1 / sqrt(pi)), 0.08)
  expect_lt(abs(mean(draws$mu[, 2]) - 1 / sqrt(pi)), 0.08)

  # theta ~ Dirichlet(1, 1): mean 0.5
  expect_lt(abs(mean(draws$theta[, 1]) - 0.5), 0.05)
})

test_that("mu_sd and tau_sd change the priors", {
  fit <- suppressWarnings(re_mix(y, sd, M = 1, mu_sd = 3, tau_sd = 0.5, prior_only = TRUE,
                                 chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  draws <- rstan::extract(fit, pars = c("mu", "tau"))
  expect_lt(abs(mean(draws$mu)), 0.3)
  expect_lt(abs(sd(draws$mu) - 3), 0.3)
  expect_lt(abs(mean(draws$tau) - 0.5 * sqrt(2 / pi)), 0.05)
})

test_that("the gap prior on the means is applied", {
  # gaps ~ lognormal(log(0.2), 0.3). The mean of a lognormal(meanlog, sdlog)
  # is exp(meanlog + sdlog^2 / 2), here 0.2 * exp(0.3^2 / 2) = 0.209.
  fit <- suppressWarnings(re_mix(y, sd, M = 3, mean_prior = "gap", prior_only = TRUE,
                                 chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  draws <- rstan::extract(fit, pars = c("mu1", "mu_gap"))
  expect_lt(abs(mean(draws$mu1)), 0.1)
  expect_lt(abs(mean(draws$mu_gap) - 0.2 * exp(0.3^2 / 2)), 0.01)

  # with gap_min the gaps are floored
  fit <- suppressWarnings(re_mix(y, sd, M = 3, mean_prior = "gap", gap_min = 0.2, prior_only = TRUE,
                                 chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  draws <- rstan::extract(fit, pars = c("mu", "mu_gap"))
  expect_true(all(draws$mu_gap >= 0.2))
  expect_true(all(draws$mu[, 2] - draws$mu[, 1] >= 0.2))
})

test_that("the inverse-gamma prior on tau is applied", {
  # tau ~ inverse-gamma(2.5, 0.15). If tau is inverse-gamma(a, b) then 1 / tau
  # is gamma(a, rate = b), so the median of tau is b / qgamma(0.5, a) = 0.069.
  # The median is used because the inverse-gamma is heavy-tailed, which makes
  # the sample mean unstable.
  fit <- suppressWarnings(re_mix(y, sd, M = 1, tau_prior = "inv_gamma", prior_only = TRUE,
                                 chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  tau <- rstan::extract(fit, pars = "tau")$tau
  expect_lt(abs(median(tau) - 0.15 / qgamma(0.5, 2.5)), 0.01)
  # the density vanishes at 0, so tiny values are practically impossible
  expect_true(all(tau > 0.01))
})

test_that("sel_mix samples the selection-weight prior when prior_only = TRUE", {
  # omega = cumsum(omega_raw) with omega_raw ~ Dirichlet(1, 1, 1), i.e. uniform
  # on the simplex. The first element of a Dirichlet(1, 1, 1) is Beta(1, 2)
  # with mean 1/3, the sum of the first two elements is Beta(2, 1) with mean
  # 2/3, and the sum of all three is exactly 1.
  fit <- suppressWarnings(sel_mix(y, sd, M = 1, steps = c(0.95, 0.975), prior_only = TRUE,
                                  chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  omega <- rstan::extract(fit, pars = "omega")$omega
  expect_equal(ncol(omega), 3)
  expect_lt(abs(mean(omega[, 1]) - 1 / 3), 0.03)
  expect_lt(abs(mean(omega[, 2]) - 2 / 3), 0.03)
  expect_true(all(abs(omega[, 3] - 1) < 1e-8))
  expect_true(all(omega[, 1] <= omega[, 2]))
})

test_that("sel_flexpb fixes the weights of unselected groups at 1", {
  fit <- suppressWarnings(sel_flexpb(y, sd, grp = c(1, 2, 2), M = 1, steps = c(0.95, 0.975),
                                     unselected = 2, prior_only = TRUE,
                                     chains = 1, iter = 4000, warmup = 1000, seed = 1, refresh = 0))
  omega <- rstan::extract(fit, pars = "omega")$omega   # draws x G x 3
  expect_equal(dim(omega)[2:3], c(2, 3))
  expect_lt(abs(mean(omega[, 1, 1]) - 1 / 3), 0.03)
  expect_true(all(omega[, 2, ] == 1))
})
