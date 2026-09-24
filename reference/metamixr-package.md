# metamixr: 'MetaMix': Bayesian Meta-Analytic Mixture Models

Fits Bayesian random-effects meta-analytic mixture models in which the
true study effects are drawn from a mixture of normal components with
component-specific means and between-study heterogeneity. The mixture
model can be combined with step-function selection models (Vevea and
Hedges, 1995,
[doi:10.1007/BF02294384](https://doi.org/10.1007/BF02294384) ) to adjust
for publication bias, either with one selection model for all studies or
with separate selection models for groups of studies. Models are
estimated by Hamiltonian Monte Carlo using 'rstan' and returned as
'stanfit' objects, so posterior summaries, posterior predictive checks
and marginal likelihoods for model comparison (via 'bridgesampling') are
readily available. Also includes a data simulator and the nudging
meta-analysis data of Mertens et al. (2022)
[doi:10.1073/pnas.2107346118](https://doi.org/10.1073/pnas.2107346118) .

## Details

The package implements Bayesian meta-analytic mixture models with
publication bias correction. Three model fitting functions all return a
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
object:

- [`re_mix()`](https://maxmaier42.github.io/metamix/reference/re_mix.md):

  random-effects meta-analytic mixture model without publication bias
  adjustment (the ordinary random-effects meta-analysis when `M = 1`).

- [`sel_mix()`](https://maxmaier42.github.io/metamix/reference/sel_mix.md):

  the same mixture model combined with a step-function selection model
  to adjust for publication bias.

- [`sel_flexpb()`](https://maxmaier42.github.io/metamix/reference/sel_flexpb.md):

  the selection model but different selection processes can be estimated
  for different groups of studies.

[`sim_mix()`](https://maxmaier42.github.io/metamix/reference/sim_mix.md)
simulates data from these models, and
[mertens_nudge](https://maxmaier42.github.io/metamix/reference/mertens_nudge.md)
contains the nudging meta-analysis of Mertens et al. (2022) used in the
vignettes.

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

Useful links:

- <https://github.com/MaxMaier42/metamix>

- Report bugs at <https://github.com/MaxMaier42/metamix/issues>

## Author

**Maintainer**: Maximilian Maier <maximilianmaier0401@gmail.com>
([ORCID](https://orcid.org/0000-0002-9873-6096))

Other contributors:

- Trustees of Columbia University \[copyright holder\]
