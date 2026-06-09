#' Random Effects Meta-Analytic Mixtures with Selection Models
#'
#' @export
#' @param y vector of primary study effects
#' @param sd standard deviation of primary study effects
#' @param M number of mixture components
#' @param steps p-value cutoffs (one-sided p-values). Currently only supports one or two steps.
#' @param one_sided whether selection is one or two-sided
#' @param mu_sd standard deviation of the normal prior on the component means
#'   (used by both mean-prior settings: as the prior on every ordered mean when
#'   `mean_prior = "normal"`, and as the prior on the first mean when
#'   `mean_prior = "gap"`). Defaults to `1`.
#' @param tau_sd standard deviation of the half-normal prior on the heterogeneity
#'   SD `tau` (used when `tau_prior = "half_normal"`). Defaults to `0.2`.
#' @param mean_prior prior family for the component means. `"normal"` (default)
#'   places an independent normal(0, `mu_sd`) prior on each ordered mean -- the
#'   classic specification. `"gap"` uses the repulsive setup: normal on the first
#'   mean plus a lognormal prior on the positive gaps between adjacent means
#'   (controlled by `gap_meanlog`, `gap_sdlog`, `gap_min`).
#' @param tau_prior prior family for the heterogeneity SD `tau`. `"half_normal"`
#'   (default) uses a half-normal(0, `tau_sd`) prior -- the classic specification.
#'   `"inv_gamma"` places an inverse-gamma(`tau_alpha`, `tau_beta`) prior directly
#'   on `tau`, suppressing near-degenerate (spike) components.
#' @param gap_meanlog log-scale location of the repulsive lognormal gap prior
#'   (only used when `mean_prior = "gap"`). Defaults to `log(0.1)`.
#' @param gap_sdlog log-scale standard deviation of the repulsive lognormal gap
#'   prior (only used when `mean_prior = "gap"`). Defaults to `0.4`.
#' @param gap_min hard lower floor on the gaps between adjacent means: each gap is
#'   constrained to be at least `gap_min` and the lognormal gap prior is truncated
#'   to `[gap_min, Inf)` (only used when `mean_prior = "gap"`). Defaults to `0`.
#' @param tau_alpha shape of the inverse-gamma prior on `tau` (only used when
#'   `tau_prior = "inv_gamma"`). Defaults to `2.5`.
#' @param tau_beta scale of the inverse-gamma prior on `tau` (only used when
#'   `tau_prior = "inv_gamma"`). Defaults to `0.15`.
#' @param prior_only if TRUE, sample from the prior only (no data).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
sel_mix <- function(y, sd, M, steps = c(0.9, 0.95), one_sided = TRUE,
                    mu_sd = 1, tau_sd = 0.2,
                    mean_prior = c("normal", "gap"),
                    tau_prior = c("half_normal", "inv_gamma"),
                    gap_meanlog = log(0.1), gap_sdlog = 0.4, gap_min = 0,
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
