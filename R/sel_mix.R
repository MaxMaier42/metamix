#' Random Effects Meta-Analytic Mixtures with Selection Models
#'
#' @export
#' @param y vector of primary study effects
#' @param sd standard deviation of primary study effects
#' @param M number of mixture components
#' @param steps p-value cutoffs (one-sided p-values). Currently only supports one or two steps.
#' @param one_sided whether selection is one or two-sided
#' @param mu_sd standard deviation of prior on mean
#' @param tau_sd standard deviation of prior on heterogeneity
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
sel_mix <- function(y, sd, M, steps = c(0.5, 0.95), one_sided = TRUE, mu_sd = 1, tau_sd = 0.2, ...){
  n_step <- length(steps)
  if(!(n_step == 1 | n_step == 2)){
    stop("Please specify one or two steps. The package does not currently support larger step numbers.")
  }

  if(mu_sd <= 0 | tau_sd <= 0){
    stop("mu_sd and tau_sd must be > 0.")
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
                   I = I,
                   M = M,
                   crit_v = as.array(qnorm(steps)),
                   n_step = n_step,
                   one_sided = one_sided,
                   mu_sd = mu_sd,
                   tau_sd = tau_sd)
  out <- rstan::sampling(stanmodels$selection_mix, data = standata, ...)
  return(out)
}
