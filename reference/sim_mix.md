# Simulate Meta-Analytic Data from a Mixture Model with Publication Bias

Simulates `K` published studies from the random-effects mixture model of
[`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md)
combined with the step-function selection model of
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md).
Studies are generated until `K` of them have passed the selection step,
so the returned data set has exactly `K` rows.

## Usage

``` r
sim_mix(
  K,
  M,
  mu,
  tau,
  steps,
  weights,
  one_sided,
  thetas = c(rep(1/M, M)),
  N_low = 25,
  N_high = 500,
  N_shape = 2,
  N_scale = 58
)
```

## Arguments

- K:

  number of (published) studies to return.

- M:

  number of mixture components.

- mu:

  numeric vector of length `M`: the mean true effect of each component.

- tau:

  numeric vector of length `M`: the between-study standard deviation of
  the true effects within each component.

- steps:

  cutoffs of the p-value intervals, in the same format as in
  [`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)
  (cumulative probabilities of the standard normal distribution; e.g.
  `c(0.95, 0.975)` for one-sided p-values of .05 and .025, i.e.
  two-sided .10 and .05 in the expected direction). Strictly increasing
  values in (0, 1).

- weights:

  publication probabilities of the p-value intervals, a vector of length
  `length(steps) + 1` with values in `[0, 1]`. The first element applies
  to the least significant interval, e.g. `c(0.2, 0.5, 1)`.

- one_sided:

  if `TRUE` selection depends on the one-sided p-value (z-statistic
  `y/sd`), if `FALSE` on the two-sided p-value (`abs(y/sd)`).

- thetas:

  mixing weights of the components; defaults to equal weights.

- N_low:

  lower bound on the total sample size of a primary study.

- N_high:

  upper bound on the total sample size of a primary study.

- N_shape:

  shape (`size`) of the negative binomial distribution used to generate
  sample sizes (see Maier et al., 2023).

- N_scale:

  scale of the negative binomial distribution used to generate sample
  sizes; the success probability is `1 / (N_scale + 1)` (see Maier et
  al., 2023).

## Value

A data frame with `K` rows and the columns `y` (observed effect size)
and `sds` (its standard error), ready to be passed to
[`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md),
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)
or
[`sel_flexpb()`](https://maxmaier42.github.io/metamix/reference/sel_flexpb.md).

## Details

See Maier (2026) for details of the model and the simulation design.

## References

Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
mixture modelling. *PsyArXiv*.
[doi:10.31234/osf.io/nkyqm_v2](https://doi.org/10.31234/osf.io/nkyqm_v2)

Maier, M., Bartoš, F., & Wagenmakers, E.-J. (2023). Robust Bayesian
meta-analysis: Addressing publication bias with model-averaging.
*Psychological Methods*, 28(1), 107–122.
[doi:10.1037/met0000405](https://doi.org/10.1037/met0000405)

## See also

[`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md),
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md),
[`sel_flexpb()`](https://maxmaier42.github.io/metamix/reference/sel_flexpb.md)

## Examples

``` r
set.seed(1)
# Two components (null effect and medium effect) without publication bias
dat <- sim_mix(K = 100, M = 2, mu = c(0, 0.5), tau = c(0.05, 0.1),
               steps = c(0.95, 0.975), weights = c(1, 1, 1), one_sided = TRUE)
head(dat)
#>              y       sds
#> 1  0.675906247 0.1638830
#> 2  0.083588527 0.1376880
#> 3  0.286695480 0.1898316
#> 4 -0.095939063 0.1288578
#> 5  0.675788329 0.1910095
#> 6  0.009781839 0.1898770

# Strong one-sided selection: studies that are not significant in the
# expected direction are published with probability .1 only
dat_pb <- sim_mix(K = 100, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
                  weights = c(0.1, 0.5, 1), one_sided = TRUE)
mean(dat_pb$y / dat_pb$sds > qnorm(0.975))
#> [1] 0.27
```
