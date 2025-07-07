data {
  int<lower=1, upper = 2> n_step;
  int<lower=0> K; // number of studies
  int<lower=0> M; //number of mixture components
  vector[K] y; // observed effect sizes
  vector<lower=0>[K] v; // variances of observed effects
  array[n_step] real crit_v;
  array[K] int<lower=1> I; // index for intervals based on p-value
  int<lower=0, upper=1> one_sided;
  real<lower=0> mu_sd; //standard deviation of prior on mu
  real<lower=0> tau_sd; //standard deviation of prior on tau
}

parameters {
  ordered[M] mu; // overall mean effect size
  array[M] real<lower=0> tau; // inverse of between-study variance
  simplex[n_step+1] omega_raw; //
  simplex[M] theta; //
}

transformed parameters {
  vector[n_step + 1] omega;  // (bias-related) publication bias
  omega = cumulative_sum(omega_raw);
}

model {
  vector[M] log_theta = log(theta);
  // Priors
  target += normal_lpdf(mu | 0, mu_sd);
  target += normal_lpdf(tau | 0, tau_sd);
  target += dirichlet_lpdf(omega_raw | rep_vector(1, n_step+1));
  target += dirichlet_lpdf(theta | rep_vector(1, M));

  // Model for observed data
  for (i in 1:K) {
    vector[M] lps = log_theta;
    if(one_sided == 1){
      if(n_step == 1){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps[m] += - log_sum_exp(
            normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[1]),
            normal_lccdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[2])
            );
        }
      }
      if(n_step == 2){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps[m] += - log_sum_exp([
            normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[1]),
            log_diff_exp(normal_lcdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[2]),
            normal_lccdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[3])]
            );
        }
      }
    }
    if(one_sided == 0){
      if(n_step == 1){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps[m] += - log_sum_exp([
            normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[2]),
            log_diff_exp(normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[1]),
            normal_lccdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[2])
            ]);
        }
      }
      if(n_step == 2){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps[m] += - log_sum_exp([
            log_diff_exp(normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[1]),
            log_diff_exp(normal_lcdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[2]),
            log_diff_exp(normal_lccdf(-crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lccdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[2]),
            normal_lccdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[3]),
            normal_lcdf(-crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[3])]
            );
        }
      }
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  matrix[K, M] posterior_probs;

  for (i in 1:K) {
    vector[M] log_weights;

    // Compute log(prior × likelihood) for each component
    for (m in 1:M) {
      log_weights[m] = log(theta[m]) + normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
    }

    // Normalize using log_sum_exp trick to avoid underflow
    for (m in 1:M) {
      posterior_probs[i, m] = exp(log_weights[m] - log_sum_exp(log_weights));
    }
  }
}


