#' Random Effects Meta-Analytic Mixtures
#'
#' @export
#' @param y vector of primary study effects
#' @param sd standard deviation of primary study effects
#' @param M number of mixture components
#' @param mu_sd standard deviation of prior on mean
#' @param tau_sd standard deviation of prior on heterogeneity
#' @param gap_meanlog log-scale location of the repulsive lognormal prior on the
#'   gaps between adjacent (ordered) component means. Larger values push
#'   components further apart and discourage the overfitted-mixture pathology of
#'   two components collapsing onto the same location. Defaults to `log(0.1)`.
#' @param gap_sdlog log-scale standard deviation of the repulsive lognormal prior
#'   on the gaps between adjacent means. Defaults to `0.5`.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
re_mix <- function(y, sd, M, mu_sd = 1, tau_sd = 0.2, gap_meanlog = log(0.1), gap_sdlog = 0.5, prior_only = FALSE, ...){
  if(mu_sd <= 0 | tau_sd <= 0){
    stop("mu_sd and tau_sd must be > 0.")
  }
  if(gap_sdlog <= 0){
    stop("gap_sdlog must be > 0.")
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
                   tau_sd = tau_sd,
                   gap_meanlog = gap_meanlog,
                   gap_sdlog = gap_sdlog)
  out <- rstan::sampling(stanmodels$random_effects_mix, data = standata, ...)
  return(out)
}
