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
  array[n] int block;
  array[n] int <lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper=n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
}

parameters {
  real intercept;       // single global intercept
  real fatigue;         // single global fatigue effect
  real<lower=0> sigma;  // single global sigma
  real b_arousal_average;       
  real<lower=0> b_arousal_sd;       
  array[n_participants] real b_arousal;       // single global effect of arousal rating

  // Missing values
  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  array[n] real mu;
  
  // Hard-coded number of trials per block (for the transition).
  // Adjust this if your block size differs.
  // 48 trials per block except for block 1.
  real trials_per_block = 48.0;

  for (i in 1:n) {
    // Fatigue contribution grows with block_trial_count[i]
    real fatigue_contrib = fatigue * trial[i];

    // Compute arousal effect
    if (block[i] == 1) {
      // For the first block, use the block's arousal rating directly
      mu[i] = intercept 
              + fatigue_contrib
              + b_arousal[participant[i]] * arousal_pbc_centered[participant[i], block[i], cue[i]];
    } else {
      // For subsequent blocks, smoothly transition from the previous block's rating
      real previous_arousal = arousal_pbc_centered[participant[i], block[i]-1, cue[i]];
      real current_arousal  = arousal_pbc_centered[participant[i], block[i], cue[i]];

      // Weight previous block's arousal more at the start of the new block,
      // and gradually shift to the current block's arousal by the block's end.
      real weight_previous = (trials_per_block - block_trial_count[i]) / trials_per_block;
      real weight_current = block_trial_count[i] / trials_per_block;

      real blended_arousal = previous_arousal * weight_previous + current_arousal * weight_current;

      mu[i] = intercept
              + fatigue_contrib
              + b_arousal[participant[i]] * blended_arousal;
    }
  }
}

model {
  // Priors
  intercept ~ normal(0, 1);
  fatigue ~ normal(0, 0.1);
  sigma ~ normal(1, 1);
  b_arousal_average ~ normal(0, 0.1);
  b_arousal_sd ~ normal(0, 0.1);
  
  b_arousal ~ normal(b_arousal_average, b_arousal_sd);

  // Likelihood
  amplitude_all ~ normal(mu, sigma);
}

generated quantities {
  array[n_observations] real mu_pred = mu[indices_observed];
  array[n_observations] real log_lik;
  
  for (i in 1:n_observations) {
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma);
  }
}
