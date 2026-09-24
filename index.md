# metamixr

metamixr fits Bayesian meta-analytic mixture models: the true effects of
the primary studies are assumed to come from a mixture of M normal
components, each with its own mean and between-study heterogeneity. The
mixture can be combined with a step-function selection model to adjust
for publication bias. Models are estimated with Stan and returned as
stanfit objects, so models can be inspected via rstan and marginal
likelihoods for model comparison can be calculated using bridgesampling.
See Maier (2026) for more information.

## Installation

The package contains Stan models that are compiled when the package is
installed, which takes several minutes and requires a C++ toolchain
(Rtools on Windows, Xcode command line tools on macOS).

``` r

# install.packages("remotes")
remotes::install_github("MaxMaier42/metamix")
```

## Example

``` r

library(metamixr)

# 100 studies from two components (null and medium effect) with one-sided
# publication bias: non-significant studies are published with probability .2,
# marginally significant ones with probability .5
set.seed(1)
dat <- sim_mix(K = 100, M = 2, mu = c(0, 0.5), tau = c(0.1, 0.1),
               steps = c(0.95, 0.975), weights = c(0.2, 0.5, 1), one_sided = TRUE)

# two-component mixture with a selection model
fit <- sel_mix(dat$y, dat$sds, M = 2, steps = c(0.95, 0.975), chains = 4, cores = 4)
print(fit, pars = c("mu", "tau", "theta", "omega"))
```

- [`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md):
  random-effects mixture model without publication bias adjustment.
- [`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md):
  random-effects mixture model with a step-function selection model.
- [`sim_mix()`](https://maxmaier42.github.io/metamix/reference/sim_mix.md):
  simulate data from the mixtures.

## References

Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
mixture modelling. *PsyArXiv*.
<https://doi.org/10.31234/osf.io/nkyqm_v2>

Maier, M., Bartoš, F., & Wagenmakers, E.-J. (2023). Robust Bayesian
meta-analysis: Addressing publication bias with model-averaging.
*Psychological Methods*, 28(1), 107–122.
<https://doi.org/10.1037/met0000405>

Mertens, S., Herberz, M., Hahnel, U. J. J., & Brosch, T. (2022). The
effectiveness of nudging: A meta-analysis of choice architecture
interventions across behavioral domains. *Proceedings of the National
Academy of Sciences*, 119(1), e2107346118.
<https://doi.org/10.1073/pnas.2107346118>
