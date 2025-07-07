test_that("re_mix recovers mu and tau (1 component)", {
  mus <- c(0)
  taus <- c(0.1)

  set.seed(42)
  dat <- sim_mix(K = 300, M = 1, mu = mus[1], tau = taus[1], steps = c(0.5, 0.95), weights = c(1,1,1), one_sided = T)
  fit <- re_mix(dat$y, dat$sds, M = 1, chains = 2, cores = 2)
  post_summary <- rstan::summary(fit, pars = c("mu", "tau"), probs = c(0.001, 0.999))$summary[ , c("mean", "0.1%", "99.9%")]

  post_cis <- post_summary[, c("0.1%", "99.9%")]

  # Check that each true mu is inside the corresponding credible interval
  expect_true(mus[1] >= post_cis[1, 1] && mus[1] <= post_cis[1, 2],
                info = paste0("mu[", i, "] not in 99.8% CI"))

  # Check that each true tau is inside the corresponding credible interval
    expect_true(taus[1] >= post_cis[2, 1] && taus[1] <= post_cis[2, 2],
                info = paste0("tau[", i, "] not in 99.8% CI"))
})

test_that("re_mix recovers mu and tau (2 components)", {
  mus <- c(0, 1)
  taus <- c(0.1, 0.2)

  set.seed(42)
  dat <- sim_mix(K = 300, M = 2, mu = mus, tau = taus, steps = c(0.5, 0.95), weights = c(1,1,1), one_sided = T)
  fit <- re_mix(dat$y, dat$sds, M = 2, chains = 2, cores = 2)
  post_summary <- rstan::summary(fit, pars = c("mu", "tau"), probs = c(0.001, 0.999))$summary[ , c("mean", "0.1%", "99.9%")]

  post_cis <- post_summary[, c("0.1%", "99.9%")]

  # Check that each true mu is inside the corresponding credible interval
  for (i in 1:2) {
    expect_true(mus[i] >= post_cis[i, 1] && mus[i] <= post_cis[i, 2],
                info = paste0("mu[", i, "] not in 99.8% CI"))
  }

  # Check that each true tau is inside the corresponding credible interval
  for (i in 1:2) {
    expect_true(taus[i] >= post_cis[i + 2, 1] && taus[i] <= post_cis[i + 2, 2],
                info = paste0("tau[", i, "] not in 99.8% CI"))
  }
})

test_that("re_mix recovers mu and tau (3 components)", {
          mus <- c(0, 1, 2)
          taus <- c(0.1, 0.2, 0.3)

          set.seed(42)
          dat <- sim_mix(K = 300, M = 3, mu = mus, tau = taus, steps = c(0.5, 0.95), weights = c(1,1,1), one_sided = T)
          fit <- re_mix(dat$y, dat$sds, M = 3, chains = 2, cores = 2)
          post_summary <- rstan::summary(fit, pars = c("mu", "tau"), probs = c(0.001, 0.999))$summary[ , c("mean", "0.1%", "99.9%")]

          post_cis <- post_summary[, c("0.1%", "99.9%")]

          # Check that each true mu is inside the corresponding credible interval
          for (i in 1:3) {
            expect_true(mus[i] >= post_cis[i, 1] && mus[i] <= post_cis[i, 2],
                        info = paste0("mu[", i, "] not in 99.8% CI"))
          }

          # Check that each true tau is inside the corresponding credible interval
          for (i in 1:3) {
            expect_true(taus[i] >= post_cis[i + 3, 1] && taus[i] <= post_cis[i + 3, 2],
                        info = paste0("tau[", i, "] not in 99.8% CI"))
          }
})
