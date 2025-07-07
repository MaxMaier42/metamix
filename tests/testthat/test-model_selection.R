test_that("no publication bias selected", {
  skip_on_cran()

  mu <- 0
  tau <- 0.2
  set.seed(42)
  dat <- sim_mix(K = 1000, M = 1, mu = mu, tau = tau, steps = c(0.5, 0.95), weights = c(1,1,1), one_sided = T)

  fit1 <- re_mix(dat$y, dat$sds, M = 1, chains = 2, cores = 2)
  fit2 <- suppressWarnings(re_mix(dat$y, dat$sds, M = 2, chains = 2, cores = 2))
  fit1pb <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975), chains = 2, cores = 2))
  fit2pb <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 2, steps = c(0.95, 0.975), chains = 2, cores = 2))

  # Compute bridge sampling estimates
  bridge1 <- bridgesampling::bridge_sampler(fit1, silent = TRUE)
  bridge2 <- bridgesampling::bridge_sampler(fit2, silent = TRUE)
  bridge1pb <- bridgesampling::bridge_sampler(fit1pb, silent = TRUE)
  bridge2pb <- bridgesampling::bridge_sampler(fit2pb, silent = TRUE)

  # Extract log marginal likelihoods
  logmls <- c(
    fit1    = bridge1$logml,
    fit2    = bridge2$logml,
    fit1pb  = bridge1pb$logml,
    fit2pb  = bridge2pb$logml
  )

  # Convert to posterior model probabilities
  logml_shifted <- logmls - max(logmls)  # for numerical stability
  mls <- exp(logml_shifted)
  pmp <- mls / sum(mls)
  #exp(logmls) / sum(exp(logmls))
  expect_true(pmp[1] > 0.5, info = paste0("Correct model not selected"))
})


test_that("publication bias three mixture selected", {
  skip_on_cran()

  mus <- c(0, 1, 2)
  taus <- c(0.1, 0.2, 0.3)
  set.seed(42)
  dat <- sim_mix(K = 1000, M = 3, mu = mus, tau = taus, steps = c(0.95, 0.975), weights = c(0.2,0.5,1), one_sided = T)

  fit2 <- suppressWarnings(re_mix(dat$y, dat$sds, M = 2, chains = 2, cores = 2))
  fit3 <- suppressWarnings(re_mix(dat$y, dat$sds, M = 3, chains = 2, cores = 2))
  fit4 <- suppressWarnings(re_mix(dat$y, dat$sds, M = 4, chains = 2, cores = 2))

  fit2pb <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 2, steps = c(0.95, 0.975), chains = 2, cores = 2))
  fit3pb <- sel_mix(dat$y, dat$sds, M = 3, steps = c(0.95, 0.975), chains = 2, cores = 2)  # warnings not suppressed
  fit4pb <- suppressWarnings(sel_mix(dat$y, dat$sds, M = 4, steps = c(0.95, 0.975), chains = 2, cores = 2))

  # Compute bridge sampling estimates
  bridge2    <- bridgesampling::bridge_sampler(fit2, silent = TRUE)
  bridge3    <- bridgesampling::bridge_sampler(fit3, silent = TRUE)
  bridge4    <- bridgesampling::bridge_sampler(fit4, silent = TRUE)

  bridge2pb  <- bridgesampling::bridge_sampler(fit2pb, silent = TRUE)
  bridge3pb  <- bridgesampling::bridge_sampler(fit3pb, silent = TRUE)
  bridge4pb  <- bridgesampling::bridge_sampler(fit4pb, silent = TRUE)

  # Extract log marginal likelihoods
  logmls <- c(
    fit2    = bridge2$logml,
    fit3    = bridge3$logml,
    fit4    = bridge4$logml,
    fit2pb  = bridge2pb$logml,
    fit3pb  = bridge3pb$logml,
    fit4pb  = bridge4pb$logml
  )

  # Compute posterior model probabilities (assuming equal priors)
  logml_shifted <- logmls - max(logmls)  # numerical stability
  mls <- exp(logml_shifted)
  pmp <- mls / sum(mls)

  #exp(logmls) / sum(exp(logmls))
  expect_true(pmp[5] > 0.5, info = paste0("Correct model not selected"))
})
