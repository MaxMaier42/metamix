data {
  int<lower=0> K; // number of studies
  vector[K] y; // observed effect sizes
  vector<lower=0>[K] v; // variances of observed effects
  int<lower=0> M; //number of mixture components
  real<lower=0> mu_sd; //standard deviation of prior on mu
  real<lower=0> tau_alpha; // shape of inverse-gamma prior on tau (heterogeneity SD)
  real<lower=0> tau_beta;  // scale of inverse-gamma prior on tau (heterogeneity SD)
  real gap_meanlog; // log-scale location of repulsive lognormal prior on gaps between adjacent means
  real<lower=0> gap_sdlog; // log-scale SD of repulsive lognormal prior on gaps between adjacent means
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
  target += normal_lpdf(mu1 | 0, mu_sd);
  // Repulsive lognormal prior on the gaps, truncated to [gap_min, inf). The
  // truncation correction (constant in the parameters) is retained so the prior
  // stays a proper normalized density for the cross-M bridge comparison.
  target += lognormal_lpdf(mu_gap | gap_meanlog, gap_sdlog);
  if (gap_min > 0)
    target += -(M - 1) * lognormal_lccdf(gap_min | gap_meanlog, gap_sdlog);
  // Inverse-gamma prior directly on tau; the lower=0 constraint supplies the
  // log Jacobian automatically (no manual term). target += keeps the prior's
  // normalizing constant, which scales with M and must be retained for the
  // cross-M bridge-sampling marginal-likelihood comparison.
  target += inv_gamma_lpdf(tau | tau_alpha, tau_beta);
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
