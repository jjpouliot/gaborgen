data {
  int<lower=0> n;
  int<lower=0> n_observations;
  int<lower=0> n_missing;
  int<lower=0> n_participants;
  int<lower=0> n_trials; 
  int<lower=0> n_phases;
  int<lower=0> n_blocks;
  int<lower=0> n_cues;

  array[n_observations] real amplitude;
  array[n] int participant;
  array[n] int phase;
  array[n] int cue;
  array[n] int trial;
  array[n] int cue_trial_count;
  array[n] int phase_trial_count;
  array[n] int block_trial_count;
  array[n_participants, n_blocks, n_cues] real arousal_pbc_centered;
  array[n_participants, n_blocks, n_cues] real arousal_pbc;
  array[n] int block;
  array[n] int <lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper=n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
}

parameters {
  // raw parameters may be close to zero, hard to sample
  // exp() transform makes them easy to sample
  real sigma_average_raw;
  real sigma_sd_raw;
  real intercept_average;
  real intercept_sd_raw;
  real fatigue_average;
  real fatigue_sd_raw;
  real intercept_fatigue_corr_raw;
  real cue_sd_raw;
  
  array[n_blocks, n_cues] real bcue;
  array[n_participants] real sigma_raw;
  array[n_participants] real intercept;
  array[n_participants] real fatigue;

  array[n_missing] real amplitude_missing;
}

transformed parameters {
  real sigma_average = exp(sigma_average_raw);
  real sigma_sd = exp(sigma_sd_raw);
  array[n_participants] real sigma = exp(sigma_raw);
  real intercept_sd = exp(intercept_sd_raw);
  real fatigue_sd = exp(fatigue_sd_raw);
  real cue_sd = exp(cue_sd_raw);
  
  real intercept_fatigue_corr = -inv_logit(intercept_fatigue_corr_raw);

  cov_matrix[2] Sigma;
  Sigma[1,1] = square(intercept_sd);
  Sigma[2,2] = square(fatigue_sd);
  Sigma[1,2] = intercept_fatigue_corr * intercept_sd * fatigue_sd;
  Sigma[2,1] = Sigma[1,2];
  
  matrix[2,2] Sigma_L = cholesky_decompose(Sigma);
  
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
}

model {
  // Priors
  sigma_average_raw ~ normal(-0.5,1.5);
  sigma_sd_raw ~ normal(-0.5,1.5);
  // it is correct that we are now using sigma_sd
  // can use transformed parameter
  sigma_raw ~ student_t(n_participants-1, sigma_average_raw, sigma_sd);

  intercept_average ~ normal(0, 0.75);
  intercept_sd_raw ~ normal(-0.5,1.5);
  fatigue_average ~ normal(0, 0.01);
  fatigue_sd_raw ~ normal(-0.5,1.5);
  
  intercept_fatigue_corr_raw ~ normal(0, 1.75); // this creates a roughly uniform distribution, if confused visualize in R with hist(boot::inv.logit(rnorm(5000,0,1.75)))
  
  for (p in 1:n_participants) {
    [intercept[p], fatigue[p]]' ~ multi_student_t_cholesky(n_participants - 1,
                                                           [intercept_average, fatigue_average]',
                                                           Sigma_L);
  }
  
  cue_sd_raw ~ normal(-0.5,1.5);

  for (b in 1:n_blocks) {
    for (q in 1:n_cues) {
      bcue[b, q] ~ student_t((n_blocks*n_cues)-1, 0, cue_sd);
    }
  }

  // Likelihood
  for(i in 1:n) {
    real mu = intercept[participant[i]] +
              (fatigue[participant[i]] * trial[i]) +
              bcue[block[i], cue[i]];

    amplitude_all[i] ~ normal(mu, sigma[participant[i]]);

  }
}

generated quantities {
  // These mu_pred and log_lik names cannot change
  
  // For model predictions
  array[n_observations] real mu_pred;
  // This is necessary for cross-validation
  array[n_observations] real log_lik;
  
  for (i in 1:n_observations) {
    mu_pred[i] = intercept[participant[indices_observed[i]]] +
                 (fatigue[participant[indices_observed[i]]] * trial[indices_observed[i]]) +
                 bcue[block[indices_observed[i]], cue[indices_observed[i]]];

    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[indices_observed[i]]]);
  }
}
