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
  real b_arousal;       // single global effect of arousal rating
  real<lower=0> sigma;  // single global sigma

  // Missing values
  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  array[n] real mu;
  for (i in 1:n) {
    // Simple linear predictor:
    // intercept + (fatigue * trial[i]) + (b_arousal * arousal rating)
    mu[i] = intercept 
            + fatigue * trial[i]
            + b_arousal * arousal_pbc_centered[participant[i], block[i], cue[i]];
  }
}

model {
  // Priors
  intercept ~ normal(0, 1);
  fatigue ~ normal(0, 0.1);
  b_arousal ~ normal(0, 1);
  sigma ~ normal(1, 1);

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
