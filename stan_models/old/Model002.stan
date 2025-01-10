data {
  int<lower=0> n;
  int<lower=0> n_observations;
  int<lower=0> n_missing;
  // int<lower=0> n_channels;
  int<lower=0> n_participants;
  int<lower=0> n_trials; // should be the number of possible trials per participant, 176 for gaborgen24 day 1
  int<lower=0> n_phases;
  // int<lower=0> n_blocks;
  int<lower=0> n_cues;

  array[n_observations] real amplitude;
  // array[n] real <lower=0> amplitude;
  array[n] int participant;
  // array[n_trials] int channel;
  array[n] int phase;
  array[n] int cue;
  array[n] int trial;
  // array[n] int block;
  // array[n] real time;
  array[n] int <lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper= n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
}

parameters {

  real<lower=0> sigma_average;
  real<lower=0> sigma_sd;
  real<lower=0> intercept_average;
  real<lower=0> intercept_sd;
  real<upper=0> fatigue_average;
  real<lower=0> fatigue_sd;
  real<lower=0, upper=1> learning_rate_average;
  real<lower=0> learning_rate_sd;
  real<lower=0> scaling_sd;
  
  array[n_participants] real<lower=0> sigma;
  array[n_participants] real<lower=0> intercept;
  array[n_participants] real<upper=0> fatigue;
  array[n_participants] real<lower=0, upper=1> learning_rate;
  array[n_cues] real scaling;
  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  array[n] real<lower=-0.000001,upper=1> CSP_associative_strength;
  array[n-1] real<lower=-1, upper=1> CSP_change;

  // Estimate CSP_associative_strength, Rescola Wagner type equation
  for(i in 1:(n-1)) {
    if (phase[i] == 1){ // Habituation phase
      CSP_change[i] = 0;
      CSP_associative_strength[i] = 0;
        CSP_associative_strength[i+1] = fmin(fmax(CSP_associative_strength[i+1], 0.0), 1.0);
    } else if (cue[i] == 1 && phase[i] != 1) {
        CSP_change[i] = learning_rate[participant[i]] * 
                        (paired[i] - CSP_associative_strength[i]);
        CSP_associative_strength[i+1] = CSP_associative_strength[i] + CSP_change[i];

        // Ensure CSP_associative_strength stays within the bounds
        CSP_associative_strength[i+1] = fmin(fmax(CSP_associative_strength[i+1], 0.0), 1.0);
    } else {
       CSP_change[i] = 0;
       CSP_associative_strength[i+1] = CSP_associative_strength[i]; // maybe should scale learning_rate by some parameter for other cues
    }
  }  
}

model {
  // Priors
  intercept_average ~ normal(0,.5);
  intercept_sd ~ normal(0,.25);
  intercept ~ student_t(n_participants, intercept_average, intercept_sd); 

  fatigue_average ~ normal(0,.025);
  fatigue_sd ~ normal(0,.025);
  fatigue ~ student_t(n_participants, fatigue_average, fatigue_sd); 

  scaling_sd ~ normal(0,.25);
  for (q in 1:n_cues) {
    scaling[q] ~ student_t(n_cues, 0, scaling_sd); //What if no multilevel prior here?
  }

  sigma_average ~ normal(1,1);
  sigma_sd ~ normal(0,.5);
  sigma ~ student_t(n_participants, sigma_average, sigma_sd); 

  learning_rate_average ~ normal(0, 1) T[0, 1];
  learning_rate_sd ~ normal(0, 1);
  learning_rate ~ student_t(n_participants, learning_rate_average, learning_rate_sd) T[0,1]; 


  // Likelihood
  for(i in 1:n) {
  real mu = intercept[participant[i]] +
            (fatigue[participant[i]] * trial[i]) +
            (scaling[cue[i]] * CSP_associative_strength[i]);

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
                 (scaling[cue[indices_observed[i]]] * CSP_associative_strength[indices_observed[i]]);
                 
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[indices_observed[i]]]);
  }
}
