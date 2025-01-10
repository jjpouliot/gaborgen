data {
  int<lower=0> n;
  // int<lower=0> n_observations;
  // int<lower=0> n_missing;
  // int<lower=0> n_channels;
  int<lower=0> n_participants;
  int<lower=0> n_trials; // should be the number of possible trials per participant, 176 for gaborgen24 day 1
  int<lower=0> n_phases;
  // int<lower=0> n_blocks;
  int<lower=0> n_cues;

  array[n] real amplitude;
  // array[n] real <lower=0> amplitude;
  array[n] int participant;
  // array[n_trials] int channel;
  array[n] int phase;
  array[n] int cue;
  array[n] int trial;
  // array[n] int block;
  // array[n] real time;
  array[n] int <lower=0, upper=1> paired;
  // array[n_observations] int<lower=1, upper= n_observations + n_missing> indices_observed;
  // array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
}

parameters {

  real<lower=0> sigma_average;
  real<lower=0> sigma_sd;
  array[n_participants] real<lower=0> sigma;

  real intercept_average;
  real<lower=0> intercept_sd;
  array[n_participants, n_cues] real intercept;
  real fatigue_average;
  real<lower=0> fatigue_sd;
  array[n_participants] real fatigue;
  real<lower=0> cue_sd;
  array[n_phases - 1, n_cues] real bcue;

  // missing values get listed as parameters, where as observed values (amplitude) was specified as data
  // array[n_missing] real <lower=0> amplitude_missing;
}


model {
  // Priors
  intercept_average ~ normal(0,.5);
  intercept_sd ~ normal(0,.25);
  for (p in 1:n_participants) {
    for (q in 1:n_cues) {
      intercept[p,q] ~ student_t(n_participants*n_cues, intercept_average, intercept_sd);
    }
  }
  
  fatigue_average ~ normal(0,.025);
  fatigue_sd ~ normal(0,.025);
  fatigue ~ student_t(n_participants, fatigue_average, fatigue_sd);
  
  cue_sd ~ normal(0,.25);
  for (ph in 1:(n_phases - 1)) {
    for (q in 1:n_cues) {
      bcue[ph, q] ~ student_t((n_phases - 1)*n_cues, 0, cue_sd);
    }
  }

  sigma_average ~ normal(0,1);
  sigma_sd ~ normal(0,1);
  sigma ~ student_t(n_participants, sigma_average, sigma_sd); // Multilevel prior for each participant

  // Likelihood
  for (i in 1:n) {
    real mu;
    if (phase[i] == 1) {
      mu = intercept[participant[i], cue[i]] +
              (fatigue[participant[i]] * trial[i]);
    } else {
      mu = intercept[participant[i], cue[i]] +
           (fatigue[participant[i]] * trial[i]) +
           bcue[phase[i] - 1, cue[i]];
    }
    

    amplitude[i] ~ normal(mu, sigma[participant[i]]);

  }
}

generated quantities {
  // These mu_pred and log_lik names cannot change
  
  // For model predictions
  array[n] real mu_pred;
  // This is necessary for cross-validation
  array[n] real log_lik;
  
  for (i in 1:n) {
    if (phase[i] == 1) {
      mu_pred[i] = intercept[participant[i], cue[i]] +
              (fatigue[participant[i]] * trial[i]);
    } else {
      mu_pred[i] = intercept[participant[i], cue[i]] +
           (fatigue[participant[i]] * trial[i]) +
           bcue[phase[i] - 1, cue[i]];
    }
    
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[i]]);
  }
}
