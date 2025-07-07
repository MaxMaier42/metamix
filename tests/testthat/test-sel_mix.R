test_that("sel_mix recovers mu and tau, and weights (1 component, one-sided)", {
  mus <- c(0)
  taus <- c(0.2)
  weights <- c(0.2, 0.5, 1)

  set.seed(1)
  dat <- sim_mix(K = 500, M = 1, mu = mus[1], tau = taus[1], steps = c(0.95, 0.975), weights = weights, one_sided = TRUE)
  fit <- sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975), chains = 2, cores = 2)
  post_summary <- rstan::summary(fit, pars = c("mu", "tau", "omega"), probs = c(0.001, 0.999))$summary[ , c("mean", "0.1%", "99.9%")]

  post_cis <- post_summary[, c("0.1%", "99.9%")]

  # Check that each true mu is inside the corresponding credible interval
  expect_true(mus[1] >= post_cis[1, 1] && mus[1] <= post_cis[1, 2],
              info = paste0("mu not in 99.8% CI"))

  # Check that each true tau is inside the corresponding credible interval
  expect_true(taus[1] >= post_cis[2, 1] && taus[1] <= post_cis[2, 2],
              info = paste0("tau not in 99.8% CI"))

  for (i in 1:3) {
    expect_true(weights[i] >= post_cis[i + 2, 1] && weights[i] <= post_cis[i + 2, 2],
                info = paste0("omega[", i, "] not in 99.8% CI"))
  }
})


test_that("sel_mix recovers mu and tau, and weights (1 component, two-sided)", {
  mus <- c(0)
  taus <- c(0.2)
  weights <- c(0.2, 0.5, 1)

  set.seed(1)
  dat <- sim_mix(K = 500, M = 1, mu = mus[1], tau = taus[1], steps = c(0.95, 0.975), weights = weights, one_sided = FALSE)
  fit <- sel_mix(dat$y, dat$sds, M = 1, one_sided = FALSE, steps = c(0.95, 0.975), chains = 2, cores = 2)
  post_summary <- rstan::summary(fit, pars = c("mu", "tau", "omega"), probs = c(0.001, 0.999))$summary[ , c("mean", "0.1%", "99.9%")]

  post_cis <- post_summary[, c("0.1%", "99.9%")]

  # Check that each true mu is inside the corresponding credible interval
  expect_true(mus[1] >= post_cis[1, 1] && mus[1] <= post_cis[1, 2],
              info = paste0("mu not in 99.8% CI"))

  # Check that each true tau is inside the corresponding credible interval
  expect_true(taus[1] >= post_cis[2, 1] && taus[1] <= post_cis[2, 2],
              info = paste0("tau not in 99.8% CI"))

  for (i in 1:3) {
    expect_true(weights[i] >= post_cis[i + 2, 1] && weights[i] <= post_cis[i + 2, 2],
                info = paste0("omega[", i, "] not in 99.8% CI"))
  }
})


test_that("sel_mix recovers mu and tau, and weights (3 component, one-sided)", {
  mus <- c(0, 1, 2)
  taus <- c(0.1, 0.2, 0.3)
  weights <- c(0.2, 0.5, 1)

  set.seed(1)
  dat <- sim_mix(K = 500, M = 3, mu = mus, tau = taus, steps = c(0.95, 0.975), weights = weights, one_sided = TRUE)
  fit <- sel_mix(dat$y, dat$sds, M = 3, one_sided = TRUE, steps = c(0.95, 0.975), chains = 2, cores = 2)
  post_summary <- rstan::summary(fit, pars = c("mu", "tau", "omega"), probs = c(0.001, 0.999))$summary[ , c("mean", "0.1%", "99.9%")]

  post_cis <- post_summary[, c("0.1%", "99.9%")]

  for (i in 1:3) {
    expect_true(mus[i] >= post_cis[i, 1] && mus[i] <= post_cis[i, 2],
                info = paste0("mu[", i, "] not in 99.8% CI"))
  }

  # Check that each true tau is inside the corresponding credible interval
  for (i in 1:3) {
    expect_true(taus[i] >= post_cis[i + 3, 1] && taus[i] <= post_cis[i + 3, 2],
                info = paste0("tau[", i, "] not in 99.8% CI"))
  }

  for (i in 1:3) {
    expect_true(weights[i] >= post_cis[i + 6, 1] && weights[i] <= post_cis[i + 6, 2],
                info = paste0("omega[", i, "] not in 99.8% CI"))
  }
})

