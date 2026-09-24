#' @keywords internal
#' @aliases metamixr
#'
#' @details
#' The package implements Bayesian meta-analytic mixture models with publication bias correction.
#' Three model fitting functions all return a
#' [`stanfit`][rstan::stanfit-class] object:
#' \describe{
#'   \item{[re_mix()]}{random-effects meta-analytic mixture model without
#'     publication bias adjustment (the ordinary random-effects meta-analysis
#'     when `M = 1`).}
#'   \item{[sel_mix()]}{the same mixture model combined with a step-function
#'     selection model to adjust for publication bias.}
#'   \item{[sel_flexpb()]}{the selection model but different selection processes
#'   can be estimated for different groups of studies.}
#' }
#' [sim_mix()] simulates data from these models, and [mertens_nudge] contains
#' the nudging meta-analysis of Mertens et al. (2022) used in the vignettes.
#'
#' @references
#' Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
#' mixture modelling. *PsyArXiv*. \doi{10.31234/osf.io/nkyqm_v2}
#'
#' Maier, M., Bartoš, F., & Wagenmakers, E.-J. (2023). Robust Bayesian
#' meta-analysis: Addressing publication bias with model-averaging.
#' *Psychological Methods*, 28(1), 107--122. \doi{10.1037/met0000405}
#'
#' Stan Development Team (2024). RStan: the R interface to Stan.
#' <https://mc-stan.org/>
#'
#' @useDynLib metamixr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats qnorm dnbinom pnbinom rnorm rbinom
"_PACKAGE"
