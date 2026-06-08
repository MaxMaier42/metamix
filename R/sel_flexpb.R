#' Random Effects Meta-Analytic Mixtures with Group-Specific Selection Models
#'
#' @export
#' @param y vector of primary study effects
#' @param sd standard deviation of primary study effects
#' @param grp integer vector of group membership for each study (e.g., meta-analysis source), values 1:G
#' @param M number of mixture components
#' @param steps p-value cutoffs (one-sided p-values). Currently only supports one or two steps.
#' @param one_sided whether selection is one or two-sided
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
#' @param prior_only if TRUE, sample from the prior only (no data)
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
sel_flexpb <- function(y, sd, grp, M, steps = c(0.9, 0.95), one_sided = TRUE,
                       mu_sd = 1, tau_alpha = 2.5, tau_beta = 0.15, gap_meanlog = log(0.1),
                       gap_sdlog = 0.4, prior_only = FALSE, ...) {
  n_step <- length(steps)
  if(!(n_step == 1 | n_step == 2)){
    stop("Please specify one or two steps. The package does not currently support larger step numbers.")
  }

  if(mu_sd <= 0){
    stop("mu_sd must be > 0.")
  }

  if(tau_alpha <= 0 | tau_beta <= 0){
    stop("tau_alpha and tau_beta must be > 0.")
  }

  if(gap_sdlog <= 0){
    stop("gap_sdlog must be > 0.")
  }

  if(steps[1] <= .5 & !one_sided){
    stop("all cutoffs must be larger .5 for two-sided selection.")
  }

  grp <- as.integer(grp)
  if(length(grp) != length(y)){
    stop("grp must have the same length as y.")
  }

  G <- max(grp)

  if(prior_only){
    y   <- numeric(0)
    sd  <- numeric(0)
    grp <- integer(0)
  }

  ### assign steps
  if(length(y) > 0){
    if(one_sided){
      crit <- y/sd
    } else {
      crit <- abs(y/sd)
    }
    I <- findInterval(crit, qnorm(steps)) + 1
  } else {
    I <- integer(0)
  }

  standata <- list(K = length(y),
                   y = y,
                   v = sd^2,
                   I = as.array(I),
                   grp = as.array(grp),
                   G = G,
                   M = M,
                   crit_v = as.array(qnorm(steps)),
                   n_step = n_step,
                   one_sided = as.integer(one_sided),
                   mu_sd = mu_sd,
                   tau_alpha = tau_alpha,
                   tau_beta = tau_beta,
                   gap_meanlog = gap_meanlog,
                   gap_sdlog = gap_sdlog)

  out <- rstan::sampling(stanmodels$selection_flexpb, data = standata, ...)
  return(out)
}
