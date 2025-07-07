#' Simulate Data from Random Effects Meta-Analytic Mixtures with Selection Models
#'
#' @export
#' @param K number of studies
#' @param M number of mixture components
#' @param mu vector of mixture means
#' @param tau vector of heterogeneity in each mixture component
#' @param steps p-value cutoffs (one-sided p-values). Currently only supports one or two steps.
#' @param weights relative publication probabilities
#' @param one_sided whether selection is one or two-sided
#' @param thetas mixing weights, defaults to equal weight for each component
#' @param N_low lower bound on primary study sample sizes
#' @param N_high upper bound on primary study sample sizes
#' @param N_shape shape of negative binomial distribution to generate sample sizes (see Maier et al., 2023)
#' @param N_scale scale of negative binomial distribution to generate sample sizes (see Maier et al., 2023)
#' @return A data frame with effect sizes and standard errors
#'
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


    if(I == 1){
      published <- as.logical(rbinom(1,1,weights[1]))
    }
    if(I == 2){
      published <- as.logical(rbinom(1,1,weights[2]))
    }
    if(I == 3){
      published <- as.logical(rbinom(1,1,weights[3]))
    }
    if(published == TRUE){
      y <- c(y, y_i)
      sds <- c(sds, sd_i)
    }
  }
  return(data.frame(y, sds))
}
