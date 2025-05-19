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

transformed data {
   real trial_center = (n_trials + 1) / 2.0;
}

parameters {

  real<lower = 0> sigma_average; //this won't hug zero so it does not need to be exp()
  real sigma_sd_raw;
  real intercept_average;
  real intercept_sd_raw;
  real fatigue_average;
  real fatigue_sd_raw;
  real intercept_fatigue_corr_raw;
  
  real b_arousal_average;       
  real b_arousal_sd_raw;       
  
  array[n_participants] real<lower = 0> sigma; //this won't hug zero so it does not need to be exp()
  array[n_participants] real intercept;
  array[n_participants] real fatigue;
  
  array[n_participants] real b_arousal;       

  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  real sigma_sd = exp(sigma_sd_raw);
  real intercept_sd = exp(intercept_sd_raw);
  real fatigue_sd = exp(fatigue_sd_raw);
  
  real intercept_fatigue_corr = tanh(intercept_fatigue_corr_raw);

  cov_matrix[2] Sigma;
  Sigma[1,1] = square(intercept_sd);
  Sigma[2,2] = square(fatigue_sd);
  Sigma[1,2] = intercept_fatigue_corr * intercept_sd * fatigue_sd;
  Sigma[2,1] = Sigma[1,2];
  
  matrix[2,2] Sigma_L = cholesky_decompose(Sigma);
  
  real b_arousal_sd = exp(b_arousal_sd_raw);       
  
  array[n] real mu;
  
  // Hard-coded number of trials per block (for the transition).
  // Adjust this if your block size differs.
  // 48 trials per block except for block 1.
  real trials_per_block = 48.0;

  for (i in 1:n) {
    real centered_trial = trial[i] - trial_center;
    // Fatigue contribution grows with block_trial_count[i]
    real fatigue_contrib = fatigue[participant[i]] * centered_trial;

    // Compute arousal effect
    if (block[i] == 1) {
      // For the first block, use the block's arousal rating directly
      mu[i] = intercept[participant[i]] 
              + fatigue_contrib
              // + b_arousal[participant[i]] * arousal_pbc[participant[i], block[i], cue[i]];
              + b_arousal[participant[i]] * arousal_pbc_centered[participant[i], block[i], cue[i]];
    } else {
      // For subsequent blocks, smoothly transition from the previous block's rating
      // real previous_arousal = arousal_pbc[participant[i], block[i]-1, cue[i]];
      // real current_arousal  = arousal_pbc[participant[i], block[i], cue[i]];
      real previous_arousal = arousal_pbc_centered[participant[i], block[i]-1, cue[i]];
      real current_arousal  = arousal_pbc_centered[participant[i], block[i], cue[i]];

      // Weight previous block's arousal more at the start of the new block,
      // and gradually shift to the current block's arousal by the block's end.
      real weight_previous = (trials_per_block - block_trial_count[i]) / trials_per_block;
      real weight_current = block_trial_count[i] / trials_per_block;

      real blended_arousal = previous_arousal * weight_previous + current_arousal * weight_current;

      mu[i] = intercept[participant[i]]
              + fatigue_contrib
              + b_arousal[participant[i]] * blended_arousal;
    }
  }
}

model {
  // Priors
  sigma_average ~ normal(1, 0.5);
  sigma_sd_raw ~ normal(-3, 1);

  sigma ~ student_t(n_participants-1, sigma_average, sigma_sd);

  intercept_average ~ normal(0, 0.75);
  intercept_sd_raw ~ normal(-1, 1);
  fatigue_average ~ normal(0, 0.01);
  fatigue_sd_raw ~ normal(-4.5, 0.75);
  
  intercept_fatigue_corr_raw ~ normal(0, 1.75); // this creates a roughly uniform distribution, if confused visualize in R with hist(boot::inv.logit(rnorm(5000,0,1.75)))
  
  for (p in 1:n_participants) {
    [intercept[p], fatigue[p]]' ~ multi_student_t_cholesky(n_participants - 1,
                                                           [intercept_average, fatigue_average]',
                                                           Sigma_L);
  }
  


  b_arousal_average ~ normal(0, 0.01);
  b_arousal_sd_raw ~ normal(-4.5, 0.75);
  // b_arousal_sd_raw ~ normal(-0.5,1.5);
  b_arousal ~ student_t(n_participants-1, b_arousal_average, b_arousal_sd);

  // Likelihood
  for(i in 1:n)
    amplitude_all[i] ~ normal(mu[i], sigma[participant[i]]);
}

generated quantities {
  array[n_observations] real mu_pred = mu[indices_observed];
  array[n_observations] real log_lik;
  
  for (i in 1:n_observations) {
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[indices_observed[i]]]);
  }
}
