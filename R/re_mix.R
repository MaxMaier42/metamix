#' Random Effects Meta-Analytic Mixtures
#'
#' @export
#' @param y vector of primary study effects
#' @param sd standard deviation of primary study effects
#' @param M number of mixture components
#' @param mu_sd standard deviation of prior on mean
#' @param tau_alpha shape of the inverse-gamma prior placed directly on the
#'   heterogeneity SD `tau` (not on the variance). Larger values concentrate the
#'   prior and lighten its right tail, discouraging large within-component
#'   heterogeneity. Defaults to `2.5`.
#' @param tau_beta scale of the inverse-gamma prior on `tau`. With `tau_alpha`
#'   this sets the location: the prior median is roughly `tau_beta /
#'   qgamma(0.5, tau_alpha)`. Defaults to `0.15` (median tau ~ 0.07, 95%
#'   interval ~ [0.02, 0.36], allowing tau up to ~0.5 only in the extreme tail).
#' @param gap_meanlog log-scale location of the repulsive lognormal prior on the
#'   gaps between adjacent (ordered) component means. Larger values push
#'   components further apart and discourage the overfitted-mixture pathology of
#'   two components collapsing onto the same location. Defaults to `log(0.1)`
#'   (Fisher-z scale); with the default `gap_sdlog` this puts under 5% of prior
#'   mass on separations smaller than Cohen's d = 0.1.
#' @param gap_sdlog log-scale standard deviation of the repulsive lognormal prior
#'   on the gaps between adjacent means. Defaults to `0.4`.
#' @param gap_min hard lower floor on the gaps between adjacent (ordered)
#'   component means: every gap is constrained to be at least `gap_min`, and the
#'   lognormal gap prior is truncated to `[gap_min, Inf)`. Use this to forbid
#'   near-degenerate configurations where two means sit almost on top of each
#'   other. Defaults to `0` (no floor; identical to the previous behaviour).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
re_mix <- function(y, sd, M, mu_sd = 1, tau_alpha = 2.5, tau_beta = 0.15, gap_meanlog = log(0.1), gap_sdlog = 0.4, gap_min = 0, prior_only = FALSE, ...){
  if(mu_sd <= 0){
    stop("mu_sd must be > 0.")
  }
  if(tau_alpha <= 0 | tau_beta <= 0){
    stop("tau_alpha and tau_beta must be > 0.")
  }
  if(gap_sdlog <= 0){
    stop("gap_sdlog must be > 0.")
  }
  if(gap_min < 0){
    stop("gap_min must be >= 0.")
  }
  if(prior_only){
    y  <- numeric(0)
    sd <- numeric(0)
  }
  standata <- list(K = length(y),
                   y = y,
                   v = sd^2,
                   M = M,
                   mu_sd = mu_sd,
                   tau_alpha = tau_alpha,
                   tau_beta = tau_beta,
                   gap_meanlog = gap_meanlog,
                   gap_sdlog = gap_sdlog,
                   gap_min = gap_min)
  out <- rstan::sampling(stanmodels$random_effects_mix, data = standata, ...)
  return(out)
}
