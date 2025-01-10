data {
  int<lower=0> n;
  int<lower=0> n_observations;
  int<lower=0> n_missing;
  // int<lower=0> n_channels;
  int<lower=0> n_participants;
  int<lower=0> n_trials; // should be the number of possible trials per participant, 176 for gaborgen24 day 1
  int<lower=0> n_phases;
  int<lower=0> n_blocks;
  int<lower=0> n_cues;

  array[n_observations] real amplitude;
  // array[n] real <lower=0> amplitude;
  array[n] int participant;
  // array[n_trials] int channel;
  array[n] int phase;
  array[n] int cue;
  array[n] int trial;
  array[n] int cue_trial_count;
  array[n] int phase_trial_count;
  array[n] int block_trial_count;
  array[n_participants, n_blocks, n_cues] real arousal_pbc_centered;
  array[n] int block;
  // array[n] real time;
  array[n] int <lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper= n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
}

parameters {

  real<lower=0> sigma_average;
  real<lower=0> sigma_sd;
  array[n_participants] real<lower=0> sigma;

  real intercept_average;
  real<lower=0> intercept_sd;
  real intercept_par_sd;
  array[n_participants] real<lower=0> intercept_par;
  array[n_participants, n_cues] real<lower=0> intercept_par_cue;
  real fatigue_average;
  real<lower=0> fatigue_sd;
  array[n_participants] real fatigue;
  real cue_average;
  real<lower=0> cue_sd;
  array[n_participants] real bcue;

  // missing values get listed as parameters, where as observed values (amplitude) was specified as data
  array[n_missing] real <lower=0> amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  array[n] real mu;
  real current_fatigue;
  real current_cue_effect;
  real previous_cue_effect;
  real combined_effect;
  for(i in 1:n) {
    if (trial[i] == 1){
      current_fatigue = 0.0;
      current_cue_effect = 0.0;
    }
    if (block[i] == 1){
      current_fatigue = current_fatigue + fatigue[participant[i]];
      current_cue_effect = bcue[participant[i]] * 
                           arousal_pbc_centered[participant[i], 
                                                block[i], 
                                                cue[i]];
      combined_effect = current_fatigue + current_cue_effect;
      
      mu[i] = intercept_par[participant[i]] +
              combined_effect;
    } else {
      current_fatigue = current_fatigue + fatigue[participant[i]];
      
      previous_cue_effect = bcue[participant[i]] * 
                            arousal_pbc_centered[participant[i], 
                                                 (block[i] - 1), 
                                                 cue[i]];
                                                
      current_cue_effect = bcue[participant[i]] * 
                           arousal_pbc_centered[participant[i], 
                                                block[i], 
                                                cue[i]];
      // smooth transition from previous block rating
      // hard coding trial count of 48 per block after block 1
      combined_effect = current_fatigue + 
                        ((previous_cue_effect * ((48.0 - block_trial_count[i])/48.0)) +
                        (current_cue_effect * (block_trial_count[i]/48.0)));
      
      mu[i] = intercept_par[participant[i]] +
              combined_effect;
    }
  }
}

model {
  // Priors
  intercept_average ~ normal(0,.5);
  intercept_sd ~ normal(0,.25);
  intercept_par ~ student_t(n_participants-1, intercept_average, intercept_sd); // Multilevel prior for each participant
  // intercept_par_sd ~ normal(0,.25);
  // for (p in 1:n_participants){
  //   for (q in 1:n_cues){
  //     intercept_par_cue[p, q] ~ student_t((n_participants * n_cues)-1, intercept_par[p], intercept_par_sd);
  //   }
  // }

  fatigue_average ~ normal(0,.025);
  fatigue_sd ~ normal(0,.025);
  fatigue ~ student_t(n_participants-1, fatigue_average, fatigue_sd); // Multilevel prior for each participant

  cue_sd ~ normal(0,.5);
  cue_average ~ normal(0, .25);
  bcue ~ student_t(n_participants-1, cue_average, cue_sd);

  sigma_average ~ normal(1,1);
  sigma_sd ~ normal(0,.5);
  sigma ~ student_t(n_participants-1, sigma_average, sigma_sd); // Multilevel prior for each participant

  // Likelihood, mu already defined in tranformed parameters
  for(i in 1:n) {
    amplitude_all[i] ~ normal(mu[i], sigma[participant[i]]);
  }
}

generated quantities {
  // These mu_pred and log_lik names cannot change
  
  // For model predictions
  array[n_observations] real mu_pred;
  mu_pred = mu[indices_observed];
  // This is necessary for cross-validation
  array[n_observations] real log_lik;
  
  for (i in 1:n_observations) {
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[indices_observed[i]]]);
  }
}
