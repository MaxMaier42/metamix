data {
  int<lower=0> K; // number of studies
  vector[K] y; // observed effect sizes
  vector<lower=0>[K] v; // variances of observed effects
  int<lower=0> M; //number of mixture components
  int<lower=0, upper=1> use_gap_prior; // 0 = independent normal on each (ordered) mean; 1 = repulsive gap prior
  int<lower=0, upper=1> use_inv_gamma_tau; // 0 = half-normal prior on tau; 1 = inverse-gamma prior on tau
  real<lower=0> mu_sd; //standard deviation of prior on mu
  real<lower=0> tau_sd; // SD of half-normal prior on tau (used when use_inv_gamma_tau = 0)
  real<lower=0> tau_alpha; // shape of inverse-gamma prior on tau (used when use_inv_gamma_tau = 1)
  real<lower=0> tau_beta;  // scale of inverse-gamma prior on tau (used when use_inv_gamma_tau = 1)
  real gap_meanlog; // log-scale location of repulsive lognormal gap prior (used when use_gap_prior = 1)
  real<lower=0> gap_sdlog; // log-scale SD of repulsive lognormal gap prior (used when use_gap_prior = 1)
  real<lower=0> gap_min; // hard lower floor on the gaps between adjacent means (0 = no floor)
}

parameters {
  real mu1; // smallest component mean
  vector<lower=gap_min>[M-1] mu_gap; // gaps between adjacent (ordered) means, floored at gap_min
  vector<lower=0>[M] tau; // between-study heterogeneity SD per component
  simplex[M] theta; //

}

transformed parameters {
  vector[M] mu; // overall mean effect size (ordered via positive gaps)
  mu[1] = mu1;
  for (m in 2:M) mu[m] = mu[m-1] + mu_gap[m-1];
}

model {
  // Priors
  vector[M] log_theta = log(theta);

  // --- prior on the component means ---
  if (use_gap_prior == 1) {
    // Repulsive setup: normal on the first mean + lognormal on the positive gaps,
    // truncated to [gap_min, inf). The truncation correction (constant in the
    // parameters) keeps the prior a proper normalized density for cross-M bridge.
    target += normal_lpdf(mu1 | 0, mu_sd);
    target += lognormal_lpdf(mu_gap | gap_meanlog, gap_sdlog);
    if (gap_min > 0)
      target += -(M - 1) * lognormal_lccdf(gap_min | gap_meanlog, gap_sdlog);
  } else {
    // Classic setup: independent normal prior on each (ordered) component mean.
    // Equivalent to the old `ordered[M] mu; mu ~ normal(0, mu_sd)` (the gap
    // parameterization is a unit-Jacobian reparameterization of the same model).
    target += normal_lpdf(mu | 0, mu_sd);
  }

  // --- prior on the heterogeneity SD tau (lower=0 supplies the log Jacobian) ---
  if (use_inv_gamma_tau == 1) {
    target += inv_gamma_lpdf(tau | tau_alpha, tau_beta);
  } else {
    target += normal_lpdf(tau | 0, tau_sd); // half-normal (tau constrained > 0)
  }

  target += dirichlet_lpdf(theta | rep_vector(1, M));

  // Model for observed data
  for (i in 1:K) {
    vector[M] lps = log_theta;
      for(m in 1:M){
        lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2)); // Likelihood for observed effect sizes
      }
      target += log_sum_exp(lps);
  }
}

generated quantities {
  vector[M] log_tau = log(tau);  // kept for funnel diagnostics (tau on log scale)
  matrix[K, M] posterior_probs;
  vector[K] y_rep;           // posterior predictive effect sizes
  vector[K] sd_rep;          // within-study SDs (sampled from empirical distribution)

  // ===== 1. Compute posterior component probabilities =====
  for (i in 1:K) {
    vector[M] log_weights;

    // Compute log(prior × likelihood) for each component
    for (m in 1:M) {
      log_weights[m] = log(theta[m]) + normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
    }

    // Normalize to get responsibilities
    for (m in 1:M) {
      posterior_probs[i, m] = exp(log_weights[m] - log_sum_exp(log_weights));
    }
  }

  // ===== 2. Generate posterior predictive samples =====
  for (i in 1:K) {
    // Sample component from MIXTURE WEIGHTS (not posterior responsibilities)
    int component = categorical_rng(theta);

    // Sample a study index to borrow its within-study variance structure
    // Weight by how much that study belongs to this component
    vector[K] p_study;
    for (k in 1:K) {
      p_study[k] = posterior_probs[k, component];
    }
    p_study = p_study / sum(p_study);  // normalize
    int study_idx = categorical_rng(p_study);

    // Use the sampled study's within-study variance
    sd_rep[i] = sqrt(v[study_idx]);

    // Generate effect size from the component, using sampled variance structure
    y_rep[i] = normal_rng(mu[component], sqrt(v[study_idx] + square(tau[component])));

  }
}
