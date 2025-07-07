#' Random Effects Meta-Analytic Mixtures
#'
#' @export
#' @param y vector of primary study effects
#' @param sd standard deviation of primary study effects
#' @param M number of mixture components
#' @param mu_sd standard deviation of prior on mean
#' @param tau_sd standard deviation of prior on heterogeneity
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
re_mix <- function(y, sd, M, mu_sd = 1, tau_sd = 0.2, ...){
  if(mu_sd <= 0 | tau_sd <= 0){
    stop("mu_sd and tau_sd must be > 0.")
  }
  standata <- list(K = length(y),
                   y = y,
                   v = sd^2,
                   M = M,
                   mu_sd = mu_sd,
                   tau_sd = tau_sd)
  out <- rstan::sampling(stanmodels$random_effects_mix, data = standata, ...)
  return(out)
}
