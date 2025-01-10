data {
  int<lower=0> n;
  int<lower=0> n_observations;
  int<lower=0> n_missing;
  int<lower=0> n_participants;
  int<lower=0> n_trials; // should be the number of possible trials per participant, e.g., 176 for gaborgen24 day 1
  int<lower=0> n_phases;
  int<lower=0> n_blocks;
  int<lower=0> n_cues;
  int<lower=0> n_knots; // number of knots for B-splines
  int<lower=0> spline_degree; // degree of the B-splines
  
  array[n_observations] real amplitude;
  array[n] int participant;
  array[n] int phase;
  array[n] int cue;
  array[n] int trial;
  array[n] int block;
  array[n] real time;
  array[n] int<lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper=n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
  
  matrix[n_trials, n_knots] spline_basis; // Spline basis matrix for trials
}

parameters {
  real<lower=0> sigma_average;
  real<lower=0> sigma_sd;
  array[n_participants] real<lower=0> sigma;

  real intercept_average;
  real<lower=0> intercept_sd;
  array[n_participants] real intercept;
  real<upper=0> fatigue_average;
  real<lower=0> fatigue_sd;
  array[n_participants] real<upper=0> fatigue;

  matrix[n_cues, n_knots] spline_coefs; // Coefficients for B-splines per cue

  // missing values get listed as parameters, whereas observed values (amplitude) were specified as data
  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
}

model {
  // Priors
  intercept_average ~ normal(0, 1);
  intercept_sd ~ normal(0, 0.5);
  intercept ~ student_t(n_participants, intercept_average, intercept_sd); // Multilevel prior for each participant

  fatigue_average ~ normal(0, 0.25);
  fatigue_sd ~ normal(0, 0.25);
  fatigue ~ student_t(n_participants, fatigue_average, fatigue_sd) T[, 0]; // Multilevel prior for each participant

  for (q in 1:n_cues) {
    spline_coefs[q] ~ student_t(3, 0, 1); // Multilevel prior for regularization of spline coefficients
  } // Regularizing spline coefficients

  sigma_average ~ normal(0, 1);
  sigma_sd ~ normal(0, 0.25);
  sigma ~ student_t(n_participants, sigma_average, sigma_sd); // Multilevel prior for each participant

  // Likelihood
  for (i in 1:n) {
    // Compute the spline effect for the current cue and trial
    real spline_effect = dot_product(spline_basis[trial[i]], to_array_1d(spline_coefs[cue[i]]));

    real mu = intercept[participant[i]] +
              (fatigue[participant[i]] * trial[i]) +
              spline_effect;

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
    int current_index = indices_observed[i];
    real spline_effect = dot_product(spline_basis[trial[current_index]], to_array_1d(spline_coefs[cue[current_index]]));
    mu_pred[i] = intercept[participant[current_index]] +
                 (fatigue[participant[current_index]] * trial[current_index]) +
                 spline_effect;
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[current_index]]);
  }
}
