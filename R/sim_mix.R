#' Simulate Meta-Analytic Data from a Mixture Model with Publication Bias
#'
#' Simulates `K` published studies from the random-effects mixture model of
#' [re_mix()] combined with the step-function selection model of [sel_mix()].
#' Studies are generated until `K` of them have passed the selection step, so
#' the returned data set has exactly `K` rows.
#'
#' See Maier (2026) for details of the model and the simulation design.
#'
#' @param K number of (published) studies to return.
#' @param M number of mixture components.
#' @param mu numeric vector of length `M`: the mean true effect of each component.
#' @param tau numeric vector of length `M`: the between-study standard deviation
#'   of the true effects within each component.
#' @param steps cutoffs of the p-value intervals, in the same format as in
#'   [sel_mix()] (cumulative probabilities of the standard normal distribution;
#'   e.g. `c(0.95, 0.975)` for one-sided p-values of .05 and .025, i.e.
#'   two-sided .10 and .05 in the expected direction). Strictly increasing
#'   values in (0, 1).
#' @param weights publication probabilities of the p-value intervals, a vector
#'   of length `length(steps) + 1` with values in `[0, 1]`. The first element
#'   applies to the least significant interval, e.g. `c(0.2, 0.5, 1)`.
#' @param one_sided if `TRUE` selection depends on the one-sided p-value
#'   (z-statistic `y/sd`), if `FALSE` on the two-sided p-value (`abs(y/sd)`).
#' @param thetas mixing weights of the components; defaults to equal weights.
#' @param N_low lower bound on the total sample size of a primary study.
#' @param N_high upper bound on the total sample size of a primary study.
#' @param N_shape shape (`size`) of the negative binomial distribution used to
#'   generate sample sizes (see Maier et al., 2023).
#' @param N_scale scale of the negative binomial distribution used to generate
#'   sample sizes; the success probability is `1 / (N_scale + 1)` (see Maier
#'   et al., 2023).
#'
#' @return A data frame with `K` rows and the columns `y` (observed effect
#'   size) and `sds` (its standard error), ready to be passed to [re_mix()],
#'   [sel_mix()] or [sel_flexpb()].
#'
#' @references
#' Maier, M. (2026). Addressing heterogeneity with Bayesian meta-analytic
#' mixture modelling. *PsyArXiv*. \doi{10.31234/osf.io/nkyqm_v2}
#'
#' Maier, M., Bartoš, F., & Wagenmakers, E.-J. (2023). Robust Bayesian
#' meta-analysis: Addressing publication bias with model-averaging.
#' *Psychological Methods*, 28(1), 107--122. \doi{10.1037/met0000405}
#'
#' @seealso [re_mix()], [sel_mix()], [sel_flexpb()]
#'
#' @examples
#' set.seed(1)
#' # Two components (null effect and medium effect) without publication bias
#' dat <- sim_mix(K = 100, M = 2, mu = c(0, 0.5), tau = c(0.05, 0.1),
#'                steps = c(0.95, 0.975), weights = c(1, 1, 1), one_sided = TRUE)
#' head(dat)
#'
#' # Strong one-sided selection: studies that are not significant in the
#' # expected direction are published with probability .1 only
#' dat_pb <- sim_mix(K = 100, M = 1, mu = 0, tau = 0.1, steps = c(0.95, 0.975),
#'                   weights = c(0.1, 0.5, 1), one_sided = TRUE)
#' mean(dat_pb$y / dat_pb$sds > qnorm(0.975))
#' @export
sim_mix <- function(K, M, mu, tau, steps, weights, one_sided, thetas = c(rep(1/M, M)), N_low = 25, N_high = 500, N_shape = 2, N_scale = 58){

  #simulate sample sizes from negative binomial distribution (see Maier et al., 2023)
  N_seq <- seq(N_low,N_high,1)
  N_den <- dnbinom(N_seq, size = N_shape, prob = 1/(N_scale+1) ) /
    (pnbinom(N_high, size = N_shape, prob = 1/(N_scale+1) ) - pnbinom(N_low - 1, size = N_shape, prob = 1/(N_scale+1) ))

  if(!(length(mu) == M)){
    stop("There must be as many means as mixture components.")
  }

  if(!(length(tau) == M)){
    stop("There must be as many taus as mixture components.")
  }

  if(!(length(weights) == length(steps)+1)){
    stop("There must be one weights for every p-value interval.")
  }
  if(!is.numeric(steps) || anyNA(steps) || any(steps <= 0 | steps >= 1)){
    stop("steps must lie strictly between 0 and 1.")
  }
  if(length(steps) > 1 && any(diff(steps) <= 0)){
    stop("steps must be strictly increasing.")
  }
  if(!is.numeric(weights) || anyNA(weights) || any(weights < 0 | weights > 1)){
    stop("weights must be publication probabilities between 0 and 1.")
  }
  if(all(weights == 0)){
    stop("At least one weight must be > 0, otherwise no study can be published.")
  }

  y <- c()
  sds <- c()
  while(length(y) < K){
    cluster <- sample(1:M, 1, prob = thetas) #select mixture component
    delta_i <- rnorm(1, mu[cluster], tau[cluster]) #simulate true effect sizes

    n_i <- sample(N_seq, 1, TRUE, N_den)/2 #select sample size using negative binomial density
    v_i <- (n_i + n_i)/(n_i*n_i) + delta_i^2/(2*(n_i+n_i)) #Borenstein p.25
    sd_i <- sqrt(v_i)

    y_i <- rnorm(1, delta_i, sd_i) #simulate empirical effect sizes taking sampling variation into account

    ##simulate selection
    if(one_sided){
      crit <- y_i/sd_i
    } else {
      crit <- abs(y_i/sd_i)
    }

    I <- findInterval(crit, qnorm(steps)) + 1

    published <- as.logical(rbinom(1, 1, weights[I]))
    if(published){
      y <- c(y, y_i)
      sds <- c(sds, sd_i)
    }
  }
  return(data.frame(y, sds))
}
