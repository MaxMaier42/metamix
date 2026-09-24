test_that("sim_mix returns K studies with columns y and sds", {
  set.seed(1)
  dat <- sim_mix(K = 50, M = 2, mu = c(0, 1), tau = c(0.1, 0.2),
                 steps = c(0.95, 0.975), weights = c(1, 1, 1), one_sided = TRUE)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 50)
  expect_named(dat, c("y", "sds"))
  expect_true(all(is.finite(dat$y)))
  expect_true(all(dat$sds > 0))
})

test_that("sim_mix is reproducible with a seed", {
  set.seed(42)
  a <- sim_mix(K = 30, M = 1, mu = 0.2, tau = 0.1, steps = 0.95, weights = c(0.5, 1), one_sided = TRUE)
  set.seed(42)
  b <- sim_mix(K = 30, M = 1, mu = 0.2, tau = 0.1, steps = 0.95, weights = c(0.5, 1), one_sided = TRUE)
  expect_identical(a, b)
})

test_that("sim_mix only returns studies from intervals with a non-zero publication probability", {
  # Each study falls into a p-value interval and is kept with the publication
  # probability given in `weights` for that interval. Setting all but one
  # weight to 0 therefore forces every returned study into the remaining
  # interval, which checks that the weights are matched to the right intervals
  # (first weight = least significant interval) and that one- and two-sided
  # selection use z and |z| respectively.
  set.seed(1)
  # one-sided, only studies above the second cutoff are published
  dat <- sim_mix(K = 100, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(0, 0, 1), one_sided = TRUE)
  expect_true(all(dat$y / dat$sds > qnorm(0.975)))

  # one-sided, only studies below the first cutoff are published
  dat <- sim_mix(K = 100, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(1, 0, 0), one_sided = TRUE)
  expect_true(all(dat$y / dat$sds < qnorm(0.95)))

  # two-sided, only studies significant in either direction are published
  dat <- sim_mix(K = 100, M = 1, mu = 0, tau = 0.1, steps = 0.975,
                 weights = c(0, 1), one_sided = FALSE)
  expect_true(all(abs(dat$y / dat$sds) > qnorm(0.975)))
})

test_that("without selection the mean of y matches the true mean", {
  set.seed(1)
  dat <- sim_mix(K = 2000, M = 1, mu = 0.5, tau = 0.1, steps = c(0.95, 0.975),
                 weights = c(1, 1, 1), one_sided = TRUE)
  expect_lt(abs(mean(dat$y) - 0.5), 0.02)

  # two components with equal weights: mean of y is the average of the means
  dat <- sim_mix(K = 2000, M = 2, mu = c(0, 1), tau = c(0.1, 0.1), steps = c(0.95, 0.975),
                 weights = c(1, 1, 1), one_sided = TRUE)
  expect_lt(abs(mean(dat$y) - 0.5), 0.05)
})

test_that("sim_mix rejects invalid input", {
  expect_error(sim_mix(10, M = 2, mu = 0, tau = c(0.1, 0.1), steps = 0.95, weights = c(1, 1), one_sided = TRUE),
               "as many means")
  expect_error(sim_mix(10, M = 2, mu = c(0, 1), tau = 0.1, steps = 0.95, weights = c(1, 1), one_sided = TRUE),
               "as many taus")
  expect_error(sim_mix(10, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975), weights = c(1, 1), one_sided = TRUE),
               "one weights for every")
  expect_error(sim_mix(10, M = 1, mu = 0, tau = 0.1, steps = c(0.975, 0.95), weights = c(1, 1, 1), one_sided = TRUE),
               "strictly increasing")
  expect_error(sim_mix(10, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 1), weights = c(1, 1, 1), one_sided = TRUE),
               "between 0 and 1")
  expect_error(sim_mix(10, M = 1, mu = 0, tau = 0.1, steps = 0.95, weights = c(0.5, 2), one_sided = TRUE),
               "between 0 and 1")
  expect_error(sim_mix(10, M = 1, mu = 0, tau = 0.1, steps = 0.95, weights = c(0, 0), one_sided = TRUE),
               "At least one weight")
})
