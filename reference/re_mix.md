# Random-Effects Meta-Analytic Mixture Model

Fits a Bayesian random-effects meta-analysis in which the true study
effects arise from a mixture of `M` normal components, each with its own
mean and between-study heterogeneity (Maier, 2026). With `M = 1` the
model is the ordinary random-effects meta-analysis. No publication bias
adjustment is applied; see
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)
and
[`sel_flexpb()`](https://maxmaier42.github.io/metamix/reference/sel_flexpb.md)
for the corresponding selection models.

## Usage

``` r
re_mix(
  y,
  sd,
  M,
  mu_sd = 1,
  tau_sd = 0.2,
  mean_prior = c("normal", "gap"),
  tau_prior = c("half_normal", "inv_gamma"),
  gap_meanlog = log(0.2),
  gap_sdlog = 0.3,
  gap_min = 0,
  tau_alpha = 2.5,
  tau_beta = 0.15,
  prior_only = FALSE,
  ...
)
```

## Arguments

- y:

  numeric vector of observed effect sizes.

- sd:

  numeric vector of the standard errors of `y` (same length as `y`).

- M:

  number of mixture components (a positive integer).

- mu_sd:

  standard deviation of the normal prior on the component means (used by
  both mean-prior settings: as the prior on every ordered mean when
  `mean_prior = "normal"`, and as the prior on the first mean when
  `mean_prior = "gap"`). Defaults to `1`.

- tau_sd:

  standard deviation of the half-normal prior on the heterogeneity SD
  `tau` (used when `tau_prior = "half_normal"`). Defaults to `0.2`.

- mean_prior:

  prior family for the component means. `"normal"` (default) places an
  independent normal(0, `mu_sd`) prior on each ordered mean. `"gap"`
  uses a lognormal prior on the positive gaps between adjacent means
  (controlled by `gap_meanlog`, `gap_sdlog`, `gap_min`), which
  discourages two components collapsing onto the same location.

- tau_prior:

  prior family for the heterogeneity SD `tau`. `"half_normal"` (default)
  uses a half-normal(0, `tau_sd`) prior. `"inv_gamma"` places an
  inverse-gamma(`tau_alpha`, `tau_beta`) prior directly on `tau`.

- gap_meanlog:

  log-scale location of the repulsive lognormal gap prior (only used
  when `mean_prior = "gap"`). Defaults to `log(0.2)`.

- gap_sdlog:

  log-scale standard deviation of the repulsive lognormal gap prior
  (only used when `mean_prior = "gap"`). Defaults to `0.3`.

- gap_min:

  hard lower floor on the gaps between adjacent means: each gap is
  constrained to be at least `gap_min` and the lognormal gap prior is
  truncated to `[gap_min, Inf)` (only used when `mean_prior = "gap"`).
  Defaults to `0`.

- tau_alpha:

  shape of the inverse-gamma prior on `tau` (only used when
  `tau_prior = "inv_gamma"`). Defaults to `2.5`.

- tau_beta:

  scale of the inverse-gamma prior on `tau` (only used when
  `tau_prior = "inv_gamma"`). Defaults to `0.15`.

- prior_only:

  if `TRUE`, the data are ignored and the prior distribution is sampled
  instead (useful for prior predictive checks).

- ...:

  further arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html),
  e.g. `chains`, `iter`, `cores`, `seed`, `refresh` or `control`.

## Value

An object of class
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) as
returned by
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
It contains posterior draws of

- `mu`:

  vector of length `M`: the ordered component means.

- `tau`:

  vector of length `M`: the between-study standard deviation of each
  component.

- `theta`:

  vector of length `M`: the mixing weights.

- `mu1`, `mu_gap`:

  the underlying parameterisation of the means (the first mean and the
  `M - 1` gaps between adjacent means).

- `log_tau`:

  `log(tau)`, convenient for funnel diagnostics.

- `posterior_probs`:

  `K` by `M` matrix: the posterior probability that study `i` belongs to
  component `m`.

- `y_rep`, `sd_rep`:

  posterior predictive draws of `K` effect sizes and the standard errors
  they were generated with. For each draw a component is sampled from
  `theta`, a study is sampled with probability proportional to its
  responsibility for that component, and a new effect size is generated
  using that study's standard error.

Use e.g. `print(fit, pars = c("mu", "tau", "theta"))`, `summary(fit)` or
[`rstan::extract()`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)
to access them.

## Details

See Maier (2026) for details of the model and the priors.

## References

Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
mixture modelling. *PsyArXiv*.
[doi:10.31234/osf.io/nkyqm_v2](https://doi.org/10.31234/osf.io/nkyqm_v2)

Maier, M., Bartoš, F., & Wagenmakers, E.-J. (2023). Robust Bayesian
meta-analysis: Addressing publication bias with model-averaging.
*Psychological Methods*, 28(1), 107–122.
[doi:10.1037/met0000405](https://doi.org/10.1037/met0000405)

Stan Development Team (2024). RStan: the R interface to Stan.
<https://mc-stan.org/>

## See also

[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)
and
[`sel_flexpb()`](https://maxmaier42.github.io/metamix/reference/sel_flexpb.md)
for the same mixture model combined with step-function selection models;
[`sim_mix()`](https://maxmaier42.github.io/metamix/reference/sim_mix.md)
to simulate data from the model.

Other model fitting functions:
[`sel_flexpb()`](https://maxmaier42.github.io/metamix/reference/sel_flexpb.md),
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)

## Examples

``` r
# \donttest{
# Simulate 40 studies from a two-component mixture without publication bias
set.seed(1)
dat <- sim_mix(K = 40, M = 2, mu = c(0, 0.6), tau = c(0.1, 0.1),
               steps = c(0.95, 0.975), weights = c(1, 1, 1), one_sided = TRUE)

# A single short chain keeps the example fast; use e.g. chains = 4 and
# iter = 2000 (the rstan defaults) for a real analysis.
fit <- re_mix(dat$y, dat$sds, M = 2, chains = 1, iter = 500, refresh = 0)
#> Warning: The largest R-hat is 1.08, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
print(fit, pars = c("mu", "tau", "theta"))
#> Inference for Stan model: random_effects_mix.
#> 1 chains, each with iter=500; warmup=250; thin=1; 
#> post-warmup draws per chain=250, total post-warmup draws=250.
#> 
#>          mean se_mean   sd  2.5%  25%  50%  75% 97.5% n_eff Rhat
#> mu[1]    0.06    0.03 0.15 -0.45 0.01 0.06 0.16  0.33    22 1.08
#> mu[2]    0.52    0.02 0.15  0.27 0.43 0.53 0.59  0.84    63 1.00
#> tau[1]   0.16    0.01 0.10  0.01 0.08 0.14 0.22  0.35    78 1.00
#> tau[2]   0.15    0.02 0.10  0.02 0.07 0.12 0.23  0.36    29 1.02
#> theta[1] 0.49    0.03 0.21  0.02 0.38 0.50 0.60  0.93    43 1.00
#> theta[2] 0.51    0.03 0.21  0.07 0.40 0.50 0.62  0.98    43 1.00
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Sep 24 13:47:53 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).

# Posterior probability that each study belongs to the second component
resp <- rstan::extract(fit, pars = "posterior_probs")$posterior_probs
head(colMeans(resp)[, 2])
#> [1] 0.8760180 0.2089365 0.4848856 0.1406083 0.8657916 0.1791249

# Sample from the prior only (no data enter the model)
prior <- re_mix(dat$y, dat$sds, M = 2, prior_only = TRUE,
                chains = 1, iter = 500, refresh = 0)
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
print(prior, pars = c("mu", "tau"))
#> Inference for Stan model: random_effects_mix.
#> 1 chains, each with iter=500; warmup=250; thin=1; 
#> post-warmup draws per chain=250, total post-warmup draws=250.
#> 
#>         mean se_mean   sd  2.5%   25%   50%  75% 97.5% n_eff Rhat
#> mu[1]  -0.48    0.06 0.77 -1.94 -1.02 -0.44 0.10  1.03   143    1
#> mu[2]   0.60    0.06 0.88 -1.06  0.07  0.58 1.19  2.27   204    1
#> tau[1]  0.16    0.01 0.11  0.00  0.08  0.14 0.23  0.41   289    1
#> tau[2]  0.15    0.01 0.12  0.00  0.06  0.13 0.23  0.39    92    1
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Sep 24 13:47:54 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
