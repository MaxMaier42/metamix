# Meta-Analytic Mixture Model with Group-Specific Selection Models

Fits the meta-analytic mixture model with step-function selection of
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md),
but estimates a separate selection model (vector of selection weights
`omega`) for each group of studies, e.g. for several meta-analyses that
are analysed jointly, or for published versus preregistered studies. The
mixture components (means, heterogeneity and mixing weights) are shared
across groups. Groups that are known to be free of publication bias can
be declared via `unselected`, in which case their selection weights are
fixed at 1.

## Usage

``` r
sel_flexpb(
  y,
  sd,
  grp,
  M,
  steps = c(0.9, 0.95),
  one_sided = TRUE,
  mu_sd = 1,
  tau_sd = 0.2,
  mean_prior = c("normal", "gap"),
  tau_prior = c("half_normal", "inv_gamma"),
  gap_meanlog = log(0.2),
  gap_sdlog = 0.3,
  gap_min = 0,
  tau_alpha = 2.5,
  tau_beta = 0.15,
  unselected = NULL,
  prior_only = FALSE,
  ...
)
```

## Arguments

- y:

  numeric vector of observed effect sizes.

- sd:

  numeric vector of the standard errors of `y` (same length as `y`).

- grp:

  integer vector of the same length as `y` giving the group (e.g.
  meta-analysis or study type) of each study, coded `1, ..., G`.

- M:

  number of mixture components (a positive integer).

- steps:

  cutoffs of the p-value intervals of the selection model, given as
  cumulative probabilities of the standard normal distribution (i.e.
  `qnorm(steps)` are the cutoffs on the z-statistic;
  `steps = c(0.95, 0.975)` corresponds to one-sided p-values of .05 and
  .025, or two-sided p-values of .10 and .05 in the expected direction).
  One or two strictly increasing values in (0, 1) are supported.
  Defaults to `c(0.9, 0.95)`.

- one_sided:

  if `TRUE` (default) selection only operates on p-values in the
  expected direction (z-statistic `y/sd`); to model selection favouring
  negative effects, negate `y` before fitting. If `FALSE` selection is
  two-sided (`abs(y/sd)`), in which case all `steps` must exceed 0.5.

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

- unselected:

  integer vector of group indices assumed NOT subject to publication
  bias: their selection weights `omega` are fixed to 1, so the
  likelihood for their studies reduces to the no-selection model, and
  the selection model is estimated only for the remaining groups (e.g.
  `unselected = 1` if group 1 is an unselected registry sample).
  Defaults to `NULL` (selection estimated for every group).

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

  vector of length `M`: the mixing weights among the published studies.

- `theta_preselection`:

  vector of length `M`: the mixing weights corrected for selection, i.e.
  the estimated proportions of the components among all *conducted*
  studies rather than among the published ones.

- `omega`:

  `G` by `length(steps) + 1` matrix: one row of relative publication
  probabilities of the p-value intervals per group, non-decreasing with
  the last element fixed at 1. Rows of groups listed in `unselected` are
  constant at 1.

- `omega_raw`:

  the underlying simplex of each group,
  `omega[g, ] = cumsum(omega_raw[g, ])`.

- `avg_omega`:

  vector of length `M`: the average selection weight of the studies
  attributed to each component, i.e. how strongly each component is
  affected by selection.

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
  they were generated with.

Use e.g. `print(fit, pars = c("mu", "tau", "omega"))`, `summary(fit)` or
[`rstan::extract()`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)
to access them.

## Details

See Maier (2026) for details of the model and the priors.

## References

Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
mixture modelling. *PsyArXiv*.
[doi:10.31234/osf.io/nkyqm_v2](https://doi.org/10.31234/osf.io/nkyqm_v2)

Vevea, J. L., & Hedges, L. V. (1995). A general linear model for
estimating effect size in the presence of publication bias.
*Psychometrika*, 60(3), 419–435.
[doi:10.1007/BF02294384](https://doi.org/10.1007/BF02294384)

## See also

[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)
for a single selection model,
[`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md)
for no selection,
[`sim_mix()`](https://maxmaier42.github.io/metamix/reference/sim_mix.md)
to simulate data.

Other model fitting functions:
[`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md),
[`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md)

## Examples

``` r
# \donttest{
# Two sets of studies on the same effect: a published literature subject to
# selection (group 1) and preregistered replications that are not (group 2)
set.seed(1)
lit <- sim_mix(K = 15, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
               weights = c(0.2, 0.5, 1), one_sided = TRUE)
rep <- sim_mix(K = 15, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
               weights = c(1, 1, 1), one_sided = TRUE)
y   <- c(lit$y, rep$y)
sd  <- c(lit$sds, rep$sds)
grp <- rep(1:2, each = 15)

# A single short chain keeps the example fast; use e.g. chains = 4 and
# iter = 2000 (the rstan defaults) for a real analysis. Group 2 is known to
# be unselected, so a selection model is estimated for group 1 only.
fit <- sel_flexpb(y, sd, grp, M = 1, steps = c(0.95, 0.975), unselected = 2,
                  chains = 1, iter = 400, refresh = 0)
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
print(fit, pars = c("mu", "tau", "omega"))
#> Inference for Stan model: selection_flexpb.
#> 1 chains, each with iter=400; warmup=200; thin=1; 
#> post-warmup draws per chain=200, total post-warmup draws=200.
#> 
#>            mean se_mean   sd 2.5%  25%  50%  75% 97.5% n_eff Rhat
#> mu[1]      0.25    0.00 0.04 0.18 0.23 0.25 0.28  0.31   105 1.00
#> tau[1]     0.04    0.00 0.03 0.00 0.02 0.03 0.06  0.12   144 1.00
#> omega[1,1] 0.32    0.01 0.18 0.07 0.19 0.28 0.43  0.75   180 1.00
#> omega[1,2] 0.64    0.01 0.21 0.25 0.51 0.66 0.81  0.97   207 1.00
#> omega[1,3] 1.00    0.00 0.00 1.00 1.00 1.00 1.00  1.00   216 0.99
#> omega[2,1] 1.00     NaN 0.00 1.00 1.00 1.00 1.00  1.00   NaN  NaN
#> omega[2,2] 1.00     NaN 0.00 1.00 1.00 1.00 1.00  1.00   NaN  NaN
#> omega[2,3] 1.00     NaN 0.00 1.00 1.00 1.00 1.00  1.00   NaN  NaN
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Sep 24 13:47:54 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
