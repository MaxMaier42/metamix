#' Meta-Analytic Mixture Model with a Step-Function Selection Model
#'
#' Fits the random-effects meta-analytic mixture model of [re_mix()] combined
#' with a step-function (weight-function) selection model
#' to adjust for publication bias (Maier, 2026).
#' One common selection model is estimated for all studies.
#' Use [sel_flexpb()] to estimate separate selection models for groups of
#' studies.
#'
#' See Maier (2026) for details of the model and the priors.
#'
#' @inheritParams re_mix
#' @param steps cutoffs of the p-value intervals of the selection model, given as
#'   cumulative probabilities of the standard normal distribution (i.e.
#'   `qnorm(steps)` are the cutoffs on the z-statistic; `steps = c(0.95, 0.975)`
#'   corresponds to one-sided p-values of .05 and .025, or two-sided p-values
#'   of .10 and .05 in the expected direction). One or two strictly increasing
#'   values in (0, 1) are supported. Defaults to `c(0.9, 0.95)`.
#' @param one_sided if `TRUE` (default) selection only operates on p-values in
#'   the expected direction (z-statistic `y/sd`); to model selection favouring
#'   negative effects, negate `y` before fitting. If `FALSE` selection is
#'   two-sided (`abs(y/sd)`), in which case all `steps` must exceed 0.5.
#'
#' @return An object of class [`stanfit`][rstan::stanfit-class] as returned by
#' [rstan::sampling()]. It contains posterior draws of
#' \describe{
#'   \item{`mu`}{vector of length `M`: the ordered component means.}
#'   \item{`tau`}{vector of length `M`: the between-study standard deviation of
#'     each component.}
#'   \item{`theta`}{vector of length `M`: the mixing weights among the
#'     published studies.}
#'   \item{`theta_preselection`}{vector of length `M`: the mixing weights
#'     corrected for selection, i.e. the estimated proportions of the
#'     components among all *conducted* studies rather than among the
#'     published ones.}
#'   \item{`omega`}{vector of length `length(steps) + 1`: the relative
#'     publication probabilities of the p-value intervals, non-decreasing with
#'     the last element fixed at 1.}
#'   \item{`omega_raw`}{the underlying simplex, `omega = cumsum(omega_raw)`.}
#'   \item{`avg_omega`}{vector of length `M`: the average selection weight of
#'     the studies attributed to each component, i.e. how strongly each
#'     component is affected by selection.}
#'   \item{`mu1`, `mu_gap`}{the underlying parameterisation of the means (the
#'     first mean and the `M - 1` gaps between adjacent means).}
#'   \item{`log_tau`}{`log(tau)`, convenient for funnel diagnostics.}
#'   \item{`posterior_probs`}{`K` by `M` matrix: the posterior probability that
#'     study `i` belongs to component `m`.}
#'   \item{`y_rep`, `sd_rep`}{posterior predictive draws of `K` effect sizes and
#'     the standard errors they were generated with.}
#' }
#' Use e.g. `print(fit, pars = c("mu", "tau", "omega"))`, `summary(fit)` or
#' [rstan::extract()] to access them.
#'
#' @references
#' Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
#' mixture modelling. *PsyArXiv*. \doi{10.31234/osf.io/nkyqm_v2}
#'
#' Vevea, J. L., & Hedges, L. V. (1995). A general linear model for estimating
#' effect size in the presence of publication bias. *Psychometrika*, 60(3),
#' 419--435. \doi{10.1007/BF02294384}
#'
#' Maier, M., Bartoš, F., & Wagenmakers, E.-J. (2023). Robust Bayesian
#' meta-analysis: Addressing publication bias with model-averaging.
#' *Psychological Methods*, 28(1), 107--122. \doi{10.1037/met0000405}
#'
#' @seealso [re_mix()] for the model without selection, [sel_flexpb()] for
#'   group-specific selection models, [sim_mix()] to simulate data.
#' @family model fitting functions
#'
#' @examples
#' \donttest{
#' # 30 studies with a true mean of 0.3 under one-sided selection: studies that
#' # are not significant in the expected direction (two-sided p > .10) are
#' # published with probability .2, marginally significant ones with .5
#' set.seed(1)
#' dat <- sim_mix(K = 30, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
#'                weights = c(0.2, 0.5, 1), one_sided = TRUE)
#'
#' # A single short chain keeps the example fast; use e.g. chains = 4 and
#' # iter = 2000 (the rstan defaults) for a real analysis.
#' fit <- sel_mix(dat$y, dat$sds, M = 1, steps = c(0.95, 0.975),
#'                chains = 1, iter = 400, refresh = 0)
#' print(fit, pars = c("mu", "tau", "omega"))
#' }
#' @export
sel_mix <- function(y, sd, M, steps = c(0.9, 0.95), one_sided = TRUE,
                    mu_sd = 1, tau_sd = 0.2,
                    mean_prior = c("normal", "gap"),
                    tau_prior = c("half_normal", "inv_gamma"),
                    gap_meanlog = log(0.2), gap_sdlog = 0.3, gap_min = 0,
                    tau_alpha = 2.5, tau_beta = 0.15,
                    prior_only = FALSE, ...){
  mean_prior <- match.arg(mean_prior)
  tau_prior  <- match.arg(tau_prior)
  use_gap_prior     <- as.integer(mean_prior == "gap")
  use_inv_gamma_tau <- as.integer(tau_prior == "inv_gamma")

  n_step <- length(steps)
  if(!(n_step == 1 | n_step == 2)){
    stop("Please specify one or two steps. The package does not currently support larger step numbers.")
  }
  if(!is.numeric(steps) || anyNA(steps) || any(steps <= 0 | steps >= 1)){
    stop("steps must lie strictly between 0 and 1.")
  }
  if(n_step == 2 && steps[2] <= steps[1]){
    stop("steps must be strictly increasing.")
  }
  if(!is.numeric(M) || length(M) != 1 || is.na(M) || M < 1 || M != round(M)){
    stop("M must be a single positive integer.")
  }

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

  if(steps[1] <= .5 & !one_sided){
    stop("all cutoffs must be larger .5 for two-sided selection.")
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

  ### assign steps
  if(one_sided){
    crit <- y/sd
  } else {
    crit <- abs(y/sd)
  }

  I <- findInterval(crit, qnorm(steps)) + 1

  standata <- list(K = length(y),
                   y = y,
                   v = sd^2,
                   I = as.array(I),
                   M = M,
                   crit_v = as.array(qnorm(steps)),
                   n_step = n_step,
                   one_sided = as.integer(one_sided),
                   use_gap_prior = use_gap_prior,
                   use_inv_gamma_tau = use_inv_gamma_tau,
                   mu_sd = mu_sd,
                   tau_sd = tau_sd,
                   tau_alpha = tau_alpha,
                   tau_beta = tau_beta,
                   gap_meanlog = gap_meanlog,
                   gap_sdlog = gap_sdlog,
                   gap_min = gap_min)
  out <- rstan::sampling(stanmodels$selection_mix, data = standata, ...)
  return(out)
}
