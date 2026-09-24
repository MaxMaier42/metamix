#' Random-Effects Meta-Analytic Mixture Model
#'
#' Fits a Bayesian random-effects meta-analysis in which the true study effects
#' arise from a mixture of `M` normal components, each with its own mean and
#' between-study heterogeneity (Maier, 2026). With `M = 1` the model is the
#' ordinary random-effects meta-analysis. No publication bias adjustment is
#' applied; see [sel_mix()] and [sel_flexpb()] for the corresponding selection
#' models.
#'
#' See Maier (2026) for details of the model and the priors.
#'
#' @param y numeric vector of observed effect sizes.
#' @param sd numeric vector of the standard errors of `y` (same length as `y`).
#' @param M number of mixture components (a positive integer).
#' @param mu_sd standard deviation of the normal prior on the component means
#'   (used by both mean-prior settings: as the prior on every ordered mean when
#'   `mean_prior = "normal"`, and as the prior on the first mean when
#'   `mean_prior = "gap"`). Defaults to `1`.
#' @param tau_sd standard deviation of the half-normal prior on the heterogeneity
#'   SD `tau` (used when `tau_prior = "half_normal"`). Defaults to `0.2`.
#' @param mean_prior prior family for the component means. `"normal"` (default)
#'   places an independent normal(0, `mu_sd`) prior on each ordered mean.
#'   `"gap"` uses a lognormal prior on the positive gaps between adjacent means
#'   (controlled by `gap_meanlog`, `gap_sdlog`, `gap_min`), which discourages two
#'   components collapsing onto the same location.
#' @param tau_prior prior family for the heterogeneity SD `tau`. `"half_normal"`
#'   (default) uses a half-normal(0, `tau_sd`) prior.
#'   `"inv_gamma"` places an inverse-gamma(`tau_alpha`, `tau_beta`) prior directly
#'   on `tau`.
#' @param gap_meanlog log-scale location of the repulsive lognormal gap prior
#'   (only used when `mean_prior = "gap"`). Defaults to `log(0.2)`.
#' @param gap_sdlog log-scale standard deviation of the repulsive lognormal gap
#'   prior (only used when `mean_prior = "gap"`). Defaults to `0.3`.
#' @param gap_min hard lower floor on the gaps between adjacent means: each gap is
#'   constrained to be at least `gap_min` and the lognormal gap prior is truncated
#'   to `[gap_min, Inf)` (only used when `mean_prior = "gap"`). Defaults to `0`.
#' @param tau_alpha shape of the inverse-gamma prior on `tau` (only used when
#'   `tau_prior = "inv_gamma"`). Defaults to `2.5`.
#' @param tau_beta scale of the inverse-gamma prior on `tau` (only used when
#'   `tau_prior = "inv_gamma"`). Defaults to `0.15`.
#' @param prior_only if `TRUE`, the data are ignored and the prior distribution
#'   is sampled instead (useful for prior predictive checks).
#' @param ... further arguments passed to [rstan::sampling()], e.g. `chains`,
#'   `iter`, `cores`, `seed`, `refresh` or `control`.
#'
#' @return An object of class [`stanfit`][rstan::stanfit-class] as returned by
#' [rstan::sampling()]. It contains posterior draws of
#' \describe{
#'   \item{`mu`}{vector of length `M`: the ordered component means.}
#'   \item{`tau`}{vector of length `M`: the between-study standard deviation of
#'     each component.}
#'   \item{`theta`}{vector of length `M`: the mixing weights.}
#'   \item{`mu1`, `mu_gap`}{the underlying parameterisation of the means (the
#'     first mean and the `M - 1` gaps between adjacent means).}
#'   \item{`log_tau`}{`log(tau)`, convenient for funnel diagnostics.}
#'   \item{`posterior_probs`}{`K` by `M` matrix: the posterior probability that
#'     study `i` belongs to component `m`.}
#'   \item{`y_rep`, `sd_rep`}{posterior predictive draws of `K` effect sizes and
#'     the standard errors they were generated with. For each draw a component
#'     is sampled from `theta`, a study is sampled with probability
#'     proportional to its responsibility for that component, and a new effect
#'     size is generated using that study's standard error.}
#' }
#' Use e.g. `print(fit, pars = c("mu", "tau", "theta"))`, `summary(fit)` or
#' [rstan::extract()] to access them.
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
#' @seealso [sel_mix()] and [sel_flexpb()] for the same mixture model combined
#'   with step-function selection models; [sim_mix()] to simulate data from
#'   the model.
#' @family model fitting functions
#'
#' @examples
#' \donttest{
#' # Simulate 40 studies from a two-component mixture without publication bias
#' set.seed(1)
#' dat <- sim_mix(K = 40, M = 2, mu = c(0, 0.6), tau = c(0.1, 0.1),
#'                steps = c(0.95, 0.975), weights = c(1, 1, 1), one_sided = TRUE)
#'
#' # A single short chain keeps the example fast; use e.g. chains = 4 and
#' # iter = 2000 (the rstan defaults) for a real analysis.
#' fit <- re_mix(dat$y, dat$sds, M = 2, chains = 1, iter = 500, refresh = 0)
#' print(fit, pars = c("mu", "tau", "theta"))
#'
#' # Posterior probability that each study belongs to the second component
#' resp <- rstan::extract(fit, pars = "posterior_probs")$posterior_probs
#' head(colMeans(resp)[, 2])
#'
#' # Sample from the prior only (no data enter the model)
#' prior <- re_mix(dat$y, dat$sds, M = 2, prior_only = TRUE,
#'                 chains = 1, iter = 500, refresh = 0)
#' print(prior, pars = c("mu", "tau"))
#' }
#' @export
re_mix <- function(y, sd, M, mu_sd = 1, tau_sd = 0.2,
                   mean_prior = c("normal", "gap"),
                   tau_prior = c("half_normal", "inv_gamma"),
                   gap_meanlog = log(0.2), gap_sdlog = 0.3, gap_min = 0,
                   tau_alpha = 2.5, tau_beta = 0.15,
                   prior_only = FALSE, ...){
  mean_prior <- match.arg(mean_prior)
  tau_prior  <- match.arg(tau_prior)
  use_gap_prior     <- as.integer(mean_prior == "gap")
  use_inv_gamma_tau <- as.integer(tau_prior == "inv_gamma")

  if(mu_sd <= 0){
    stop("mu_sd must be > 0.")
  }
  if(tau_prior == "half_normal" && tau_sd <= 0){
    stop("tau_sd must be > 0 when tau_prior = 'half_normal'.")
  }
  if(tau_prior == "inv_gamma" && (tau_alpha <= 0 || tau_beta <= 0)){
    stop("tau_alpha and tau_beta must be > 0 when tau_prior = 'inv_gamma'.")
  }
  if(mean_prior == "gap" && gap_sdlog <= 0){
    stop("gap_sdlog must be > 0 when mean_prior = 'gap'.")
  }
  if(gap_min < 0){
    stop("gap_min must be >= 0.")
  }
  if(!is.numeric(M) || length(M) != 1 || is.na(M) || M < 1 || M != round(M)){
    stop("M must be a single positive integer.")
  }
  if(prior_only){
    y  <- numeric(0)
    sd <- numeric(0)
  }
  if(!is.numeric(y) || !is.numeric(sd)){
    stop("y and sd must be numeric vectors.")
  }
  if(length(y) != length(sd)){
    stop("y and sd must have the same length.")
  }
  if(anyNA(y) || anyNA(sd) || any(!is.finite(y)) || any(!is.finite(sd))){
    stop("y and sd must not contain missing or infinite values.")
  }
  if(any(sd <= 0)){
    stop("all standard errors in sd must be > 0.")
  }
  standata <- list(K = length(y),
                   y = y,
                   v = sd^2,
                   M = M,
                   use_gap_prior = use_gap_prior,
                   use_inv_gamma_tau = use_inv_gamma_tau,
                   mu_sd = mu_sd,
                   tau_sd = tau_sd,
                   tau_alpha = tau_alpha,
                   tau_beta = tau_beta,
                   gap_meanlog = gap_meanlog,
                   gap_sdlog = gap_sdlog,
                   gap_min = gap_min)
  out <- rstan::sampling(stanmodels$random_effects_mix, data = standata, ...)
  return(out)
}
