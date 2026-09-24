data {
  int<lower=1, upper = 2> n_step;
  int<lower=0> K; // number of studies
  int<lower=0> M; //number of mixture components
  int<lower=1> G; // number of groups (e.g., distinct meta-analyses)
  vector[K] y; // observed effect sizes
  vector<lower=0>[K] v; // variances of observed effects
  array[n_step] real crit_v;
  array[K] int<lower=1> I; // index for intervals based on p-value
  array[K] int<lower=1, upper=G> grp; // group membership for each study
  int<lower=0, upper=1> one_sided;
  int<lower=0, upper=1> use_gap_prior; // 0 = independent normal on each (ordered) mean; 1 = repulsive gap prior
  int<lower=0, upper=1> use_inv_gamma_tau; // 0 = half-normal prior on tau; 1 = inverse-gamma prior on tau
  real<lower=0> mu_sd; //standard deviation of prior on mu
  real<lower=0> tau_sd; // SD of half-normal prior on tau (used when use_inv_gamma_tau = 0)
  real<lower=0> tau_alpha; // shape of inverse-gamma prior on tau (used when use_inv_gamma_tau = 1)
  real<lower=0> tau_beta;  // scale of inverse-gamma prior on tau (used when use_inv_gamma_tau = 1)
  real gap_meanlog; // log-scale location of repulsive lognormal gap prior (used when use_gap_prior = 1)
  real<lower=0> gap_sdlog; // log-scale SD of repulsive lognormal gap prior (used when use_gap_prior = 1)
  real<lower=0> gap_min; // hard lower floor on the gaps between adjacent means (0 = no floor)
  array[G] int<lower=0, upper=1> sel_free; // 1 = selection model estimated for this group; 0 = no publication bias (omega fixed to 1)
}

parameters {
  real mu1; // smallest component mean
  vector<lower=gap_min>[M-1] mu_gap; // gaps between adjacent (ordered) means, floored at gap_min
  vector<lower=0>[M] tau; // between-study heterogeneity SD per component
  array[G] simplex[n_step+1] omega_raw; // one selection model per group
  simplex[M] theta; //
}

transformed parameters {
  array[G] vector[n_step + 1] omega;  // (bias-related) publication bias per group
  vector[M] mu; // overall mean effect size (ordered via positive gaps)
  mu[1] = mu1;
  for (m in 2:M) mu[m] = mu[m-1] + mu_gap[m-1];
  // Groups with sel_free = 0 are assumed unselected: omega fixed to 1 (the
  // likelihood then reduces exactly to the no-selection likelihood for them);
  // their omega_raw stays in the model with its Dirichlet prior only.
  for (g in 1:G) {
    if (sel_free[g] == 1)
      omega[g] = cumulative_sum(omega_raw[g]);
    else
      omega[g] = rep_vector(1.0, n_step + 1);
  }
}

model {
  vector[M] log_theta = log(theta);
  // Priors

  // --- prior on the component means ---
  if (use_gap_prior == 1) {
    target += normal_lpdf(mu1 | 0, mu_sd);
    target += lognormal_lpdf(mu_gap | gap_meanlog, gap_sdlog);
    if (gap_min > 0)
      target += -(M - 1) * lognormal_lccdf(gap_min | gap_meanlog, gap_sdlog);
  } else {
    // independent normal on each (ordered) mean; reproduces ordered[M] mu ~ normal
    target += normal_lpdf(mu | 0, mu_sd);
  }

  // --- prior on the heterogeneity SD tau ---
  if (use_inv_gamma_tau == 1) {
    target += inv_gamma_lpdf(tau | tau_alpha, tau_beta);
  } else {
    target += normal_lpdf(tau | 0, tau_sd); // half-normal (tau constrained > 0)
  }

  for (g in 1:G)
    target += dirichlet_lpdf(omega_raw[g] | rep_vector(1, n_step+1));
  target += dirichlet_lpdf(theta | rep_vector(1, M));

  // Model for observed data
  for (i in 1:K) {
    vector[M] lps = log_theta;
    int g = grp[i];

    if(one_sided == 1){
      if(n_step == 1){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[g][I[i]]);
          lps[m] += - log_sum_exp(
            normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][1]),
            normal_lccdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][2])
            );
        }
      }
      if(n_step == 2){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[g][I[i]]);
          lps[m] += - log_sum_exp([
            normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][1]),
            log_diff_exp(normal_lcdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[g][2]),
            normal_lccdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][3])]
            );
        }
      }
    }
    if(one_sided == 0){
      if(n_step == 1){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[g][I[i]]);
          lps[m] += - log_sum_exp([
            normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][2]),
            log_diff_exp(normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[g][1]),
            normal_lccdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][2])
            ]);
        }
      }
      if(n_step == 2){
        for(m in 1:M){
          lps[m] += normal_lpdf(y[i] | mu[m], sqrt(v[i] + tau[m]^2));
          lps[m] += log(omega[g][I[i]]);
          lps[m] += - log_sum_exp([
            log_diff_exp(normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[g][1]),
            log_diff_exp(normal_lcdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lcdf(crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[g][2]),
            log_diff_exp(normal_lccdf(-crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)), normal_lccdf(-crit_v[1]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2))) + log(omega[g][2]),
            normal_lccdf(crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][3]),
            normal_lcdf(-crit_v[2]*sqrt(v[i]) | mu[m], sqrt(v[i] + tau[m]^2)) + log(omega[g][3])]
            );
        }
      }
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  vector[M] log_tau = log(tau);  // kept for funnel diagnostics (tau on log scale)
  matrix[K, M] posterior_probs;
  vector[K] y_rep;
  vector[K] sd_rep;          // within-study SD used for z and as "study SD"
  vector[M] avg_omega;  // Average omega for each component
  simplex[M] theta_preselection;  // Mixture proportions before selection


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

  // Calculate average omega for each component
  // weighted by marginal likelihood of each study under that component
  // Calculate avg_omega[m] using inverse-omega weighting to undo selection
  if (K == 0) {
    // No data: selection cannot be estimated, set to 1 (no correction)
    for (m in 1:M) avg_omega[m] = 1.0;
    theta_preselection = theta;
  } else {
    for (m in 1:M) {
      real num = 0;
      real den = 0;

      for (i in 1:K) {
        real r = posterior_probs[i, m];     // responsibility for component m
        num += r;                           // = r * omega / omega
        den += r / omega[grp[i]][I[i]];    // inverse selection weight
      }

      avg_omega[m] = num / den;
    }

    // Reweight theta to get pre-selection proportions
    vector[M] theta_adjusted;
    for (m in 1:M) {
      theta_adjusted[m] = theta[m] / avg_omega[m];
    }
    theta_preselection = theta_adjusted / sum(theta_adjusted);  // Renormalize to simplex
  }


    int filled = 0;
    int max_attempts = 50 * K;   // Safety cap; adjust if acceptance low
    int attempts = 0;
    while (filled < K && attempts < max_attempts) {
      attempts += 1;

      // 1) Sample component from mixture weights (NOT posterior_probs[i])
      int component = categorical_rng(theta_preselection);

      // 2) Sample a study index i conditional on component
      //    (preserves component-specific variance distribution empirically)
      vector[K] p_i;
      for (k in 1:K) p_i[k] = posterior_probs[k, component];
      p_i = p_i / sum(p_i);              // normalize
      int i = categorical_rng(p_i);

      // 3. Draw candidate
      real sigma = sqrt(v[i] + square(tau[component]));
      real y_candidate = normal_rng(mu[component], sigma);
      real sd_candidate = sqrt(v[i]);

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

      // 6. Apply the selection model for the group of the sampled study
      real p_select = omega[grp[i]][interval_idx];


      if (p_select >= 1) {
        // Always accept (no RNG call)
        filled += 1;
        y_rep[filled] = y_candidate;
        sd_rep[filled] = sd_candidate;
      } else {
        if (bernoulli_rng(p_select) == 1) {
          filled += 1;
          y_rep[filled] = y_candidate;
          sd_rep[filled] = sd_candidate;
        }
      }
    }

    if (filled < K) {
      for (j in (filled + 1):K)
        y_rep[j] = not_a_number();
    }
}
