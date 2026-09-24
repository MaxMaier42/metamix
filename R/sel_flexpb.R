#' Meta-Analytic Mixture Model with Group-Specific Selection Models
#'
#' Fits the meta-analytic mixture model with step-function selection of
#' [sel_mix()], but estimates a separate selection model (vector of selection
#' weights `omega`) for each group of studies, e.g. for several meta-analyses
#' that are analysed jointly, or for published versus preregistered studies.
#' The mixture components (means, heterogeneity and mixing weights) are shared
#' across groups. Groups that are known to be free of publication bias can be
#' declared via `unselected`, in which case their selection weights are fixed
#' at 1.
#'
#' See Maier (2026) for details of the model and the priors.
#'
#' @inheritParams sel_mix
#' @param grp integer vector of the same length as `y` giving the group
#'   (e.g. meta-analysis or study type) of each study, coded `1, ..., G`.
#' @param unselected integer vector of group indices assumed NOT subject to
#'   publication bias: their selection weights `omega` are fixed to 1, so the
#'   likelihood for their studies reduces to the no-selection model, and the
#'   selection model is estimated only for the remaining groups (e.g.
#'   `unselected = 1` if group 1 is an unselected registry sample). Defaults to
#'   `NULL` (selection estimated for every group).
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
#'   \item{`omega`}{`G` by `length(steps) + 1` matrix: one row of relative
#'     publication probabilities of the p-value intervals per group,
#'     non-decreasing with the last element fixed at 1. Rows of groups listed
#'     in `unselected` are constant at 1.}
#'   \item{`omega_raw`}{the underlying simplex of each group,
#'     `omega[g, ] = cumsum(omega_raw[g, ])`.}
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
#' @seealso [sel_mix()] for a single selection model, [re_mix()] for no
#'   selection, [sim_mix()] to simulate data.
#' @family model fitting functions
#'
#' @examples
#' \donttest{
#' # Two sets of studies on the same effect: a published literature subject to
#' # selection (group 1) and preregistered replications that are not (group 2)
#' set.seed(1)
#' lit <- sim_mix(K = 15, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
#'                weights = c(0.2, 0.5, 1), one_sided = TRUE)
#' rep <- sim_mix(K = 15, M = 1, mu = 0.3, tau = 0.1, steps = c(0.95, 0.975),
#'                weights = c(1, 1, 1), one_sided = TRUE)
#' y   <- c(lit$y, rep$y)
#' sd  <- c(lit$sds, rep$sds)
#' grp <- rep(1:2, each = 15)
#'
#' # A single short chain keeps the example fast; use e.g. chains = 4 and
#' # iter = 2000 (the rstan defaults) for a real analysis. Group 2 is known to
#' # be unselected, so a selection model is estimated for group 1 only.
#' fit <- sel_flexpb(y, sd, grp, M = 1, steps = c(0.95, 0.975), unselected = 2,
#'                   chains = 1, iter = 400, refresh = 0)
#' print(fit, pars = c("mu", "tau", "omega"))
#' }
#' @export
sel_flexpb <- function(y, sd, grp, M, steps = c(0.9, 0.95), one_sided = TRUE,
                       mu_sd = 1, tau_sd = 0.2,
                       mean_prior = c("normal", "gap"),
                       tau_prior = c("half_normal", "inv_gamma"),
                       gap_meanlog = log(0.2), gap_sdlog = 0.3, gap_min = 0,
                       tau_alpha = 2.5, tau_beta = 0.15,
                       unselected = NULL,
                       prior_only = FALSE, ...) {
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
  if(!is.numeric(grp) || anyNA(grp) || any(grp < 1) || any(grp != round(grp))){
    stop("grp must be a vector of positive integers (group indices 1, ..., G).")
  }
  grp <- as.integer(grp)
  if(length(grp) != length(y)){
    stop("grp must have the same length as y.")
  }

  G <- max(grp)

  sel_free <- rep(1L, G)
  if(!is.null(unselected)){
    unselected <- as.integer(unselected)
    if(any(unselected < 1 | unselected > G)){
      stop("unselected must contain group indices between 1 and G.")
    }
    sel_free[unselected] <- 0L
    if(all(sel_free == 0L)){
      stop("All groups are marked unselected; use re_mix() for a model without publication bias.")
    }
  }

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
                   use_gap_prior = use_gap_prior,
                   use_inv_gamma_tau = use_inv_gamma_tau,
                   mu_sd = mu_sd,
                   tau_sd = tau_sd,
                   tau_alpha = tau_alpha,
                   tau_beta = tau_beta,
                   gap_meanlog = gap_meanlog,
                   gap_sdlog = gap_sdlog,
                   gap_min = gap_min,
                   sel_free = as.array(sel_free))

  out <- rstan::sampling(stanmodels$selection_flexpb, data = standata, ...)
  return(out)
}
