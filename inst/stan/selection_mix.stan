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
    vector[M] lps_normalize;
    
    if(one_sided == 1){
      if(n_step == 1){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps_normalize[m] = log_sum_exp(
            normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[1]),
            normal_lccdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[2])
          );
        }
        lps -= log_sum_exp(log_theta + lps_normalize);
      }
      if(n_step == 2){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps_normalize[m] = log_sum_exp([
            normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[1]),
            log_diff_exp(normal_lcdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[2]),
            normal_lccdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[3])
          ]);
        }
        lps -= log_sum_exp(log_theta + lps_normalize);
      }
    }
    if(one_sided == 0){
      if(n_step == 1){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps_normalize[m] = log_sum_exp([
            normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[2]),
            log_diff_exp(normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[1]),
            normal_lccdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[2])
          ]);
        }
        lps -= log_sum_exp(log_theta + lps_normalize);
      }
      if(n_step == 2){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[I[i]]);
          lps_normalize[m] = log_sum_exp([
            log_diff_exp(normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[1]),
            log_diff_exp(normal_lcdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[2]),
            log_diff_exp(normal_lccdf(-crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lccdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[2]),
            normal_lccdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[3]),
            normal_lcdf(-crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[3])
          ]);
        }
        lps -= log_sum_exp(log_theta + lps_normalize);
      }
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  matrix[K, M] posterior_probs;
  vector[K] y_rep;

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
  
  // Posterior predictive sampling with selection process
  {
    int filled = 0;
    int max_attempts = 50 * K;   // Safety cap; adjust if acceptance low
    int attempts = 0;
    
    while (filled < K && attempts < max_attempts) {
      attempts += 1;

      // 1) Sample component from mixture weights (theta is already pre-selection)
      int component = categorical_rng(theta);
    
      // 2) Sample a study index i conditional on component
      //    (preserves component-specific variance distribution empirically)
      vector[K] p_i;
      for (k in 1:K) p_i[k] = posterior_probs[k, component];
      p_i = p_i / sum(p_i);              // normalize
      int i = categorical_rng(p_i);

      // 3. Draw candidate
      real sigma = sqrt(v[i] + square(tau[component]));
      real y_candidate = normal_rng(mu[component], sigma);

      // 4. Compute z-statistic using only within-study SD sqrt(v[i])
      real z = y_candidate / sqrt(v[i]);
      if (one_sided == 0) {
        z = abs(z);   // two-sided: use absolute value
      }

      // 5. Interval index = 1 + number of cutpoints passed
      int interval_idx = 1;
      for (c in 1:n_step) {
        if (z > crit_v[c])
          interval_idx += 1;
      }

      real p_select = omega[interval_idx];

      if (p_select >= 1) {
        // Always accept (no RNG call)
        filled += 1;
        y_rep[filled] = y_candidate;
      } else {
        if (bernoulli_rng(p_select) == 1) {
          filled += 1;
          y_rep[filled] = y_candidate;
        }
      }
    }

    // Optional: if not all filled (rare if max_attempts high), fill with NaN
    if (filled < K) {
      for (j in (filled + 1):K)
        y_rep[j] = not_a_number();
    }
  }
}

