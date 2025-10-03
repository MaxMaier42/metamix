data {
  int<lower=0> K; // number of studies
  vector[K] y; // observed effect sizes
  vector<lower=0>[K] v; // variances of observed effects
  int<lower=0> M; //number of mixture components
  real<lower=0> mu_sd; //standard deviation of prior on mu
  real<lower=0> tau_sd; //standard deviation of prior on tau
}

parameters {
  ordered[M] mu; // overall mean effect size
  array[M] real<lower=0> tau; // inverse of between-study variance
  simplex[M] theta; //

}

model {
  // Priors
  vector[M] log_theta = log(theta);
  target += normal_lpdf(mu | 0, mu_sd);
  target += normal_lpdf(tau | 0,  tau_sd);
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
  matrix[K, M] posterior_probs;
  vector[K] y_rep; // posterior predictive samples

  for (i in 1:K) {
    vector[M] log_weights;

    // Compute log(prior × likelihood) for each component
    for (m in 1:M) {
      log_weights[m] = log(theta[m]) + normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
    }

    // Generate posterior predictive sample
    // First, sample which component this observation comes from
    int component = categorical_rng(theta);

    // Then, sample from that component's distribution
    y_rep[i] = normal_rng(mu[component], sqrt(v[i] + tau[component]^2));

    // Normalize using log_sum_exp trick to avoid underflow
    for (m in 1:M) {
      posterior_probs[i, m] = exp(log_weights[m] - log_sum_exp(log_weights));
    }
  }
}
