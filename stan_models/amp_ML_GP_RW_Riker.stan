data {
  int<lower=0> n_trials;
  int<lower=0> n_observations;
  int<lower=0> n_missing;
  int<lower=0> n_channels;
  int<lower=0> n_participants;
  int<lower=0> n_cues;

  array[n_observations] real amplitude;
  array[n_trials] int participant;
  array[n_trials] int channel;
  array[n_trials] int phase;
  array[n_trials] int cue;
  array[n_trials] int trial;
  array[n_trials] int <lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper= n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
  // array[n_observations] real k; //cue distance from CS+
}

transformed data {
  // This is used to fill the covariance matrices later on by figuring out where each correlation needs to go
  int<lower=1> num_off_diagonals = n_channels * (n_channels - 1) %/% 2;
  // Precompute indices for the off-diagonal elements
  array[num_off_diagonals] int row_idx;
  array[num_off_diagonals] int col_idx;
  {
    int idx = 1;
    for (i in 1:(n_channels - 1)) {
      for (j in (i + 1):n_channels) {
        row_idx[idx] = i;
        col_idx[idx] = j;
        idx += 1;
      }
    }
  }
}

parameters {
  real<lower=0> sigma_average;
  real<lower=0> sigma_sd;
  array[n_participants] real<lower=0> sigma;

  real<lower=0,upper=1> learning_rate_average;
  real<lower=0,upper=1> learning_rate_sd;
  array[n_participants] real<lower=0,upper=1> learning_rate;

  vector[n_channels] intercepts_average;
  array[n_channels] real<lower=0> channel_intercept_sd;
  matrix[n_participants, n_channels] intercept;
  real<lower=-1,upper=1> intercept_Rho_average;
  real<lower=0> intercept_Rho_sd;
  vector[num_off_diagonals] intercept_channel_Rhos;
  
  real<lower=0> cue_sd_average;
  real<lower=0> cue_sd_tau;
  vector[n_channels] cue_sd;
  real average_of_cues_all_chan;
  real<lower=0> average_of_cues_all_chan_sd;
  vector[n_channels] average_of_cues_per_channel;
  matrix[n_cues, n_channels] cue_scaling_average;
  array[n_channels] real<lower=0> channel_scaling_sd;
  array[n_participants, n_cues] vector[n_channels] microvolt_scaling;
  real<lower=-1,upper=1> scaling_Rho_average;
  real<lower=0> scaling_Rho_sd;
  vector[num_off_diagonals] scaling_channel_Rhos;

  // missing values get listed as parameters, where as observed values (amplitude) was specified as data
  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n_trials] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  array[n_trials] real<lower=0,upper=1> CSP_associative_strength;
  array[n_trials-1] real<lower=-1, upper=1> CSP_change;

  // Estimate CSP_associative_strength, Rescola Wagner type equation
  for(i in 1:n_trials-1) {
    if (phase[i] == 1){ // Habituation phase
      CSP_change[i] = 0;
      CSP_associative_strength[i] = 0;
    } else if (cue[i] == 1) {
        CSP_change[i] = learning_rate[participant[i]] * 
                        (paired[i] - CSP_associative_strength[i]);
    
        CSP_associative_strength[i+1] = CSP_associative_strength[i] + CSP_change[i];
    } else {
       CSP_change[i] = 0;
       CSP_associative_strength[i+1] = CSP_associative_strength[i]; // maybe should scale learning_rate by some parameter for other cues
    }
  }  
    
  // correlation, covariance, and cholesky matrix construction
  matrix[n_channels, n_channels] intercept_correlation_matrix;
  matrix[n_channels, n_channels] intercept_cov_mat;
  cholesky_factor_cov[n_channels, n_channels] L_intercept_cov_mat;
  
  // Initialize the correlation matrix
  for (i in 1:n_channels) {
    for (j in 1:n_channels) {
      intercept_correlation_matrix[i, j] = 0;
    }
    intercept_correlation_matrix[i, i] = 1;
  }

  // Fill in the off-diagonal elements using precomputed indices
  for (n in 1:num_off_diagonals) {
    int i = row_idx[n];
    int j = col_idx[n];
    intercept_correlation_matrix[i, j] = intercept_channel_Rhos[n];
    intercept_correlation_matrix[j, i] = intercept_channel_Rhos[n];  // Symmetry
  }

  // Construct the covariance matrix using the standard deviations
  for (i in 1:n_channels) {
    for (j in 1:n_channels) {
      intercept_cov_mat[i, j] = intercept_correlation_matrix[i, j] * channel_intercept_sd[i] * channel_intercept_sd[j];
    }
  }
  // Cholesky decomposition for efficiency in sampling
  L_intercept_cov_mat = cholesky_decompose(intercept_cov_mat);
  
  
  // repeat for scaling matrices
  matrix[n_channels, n_channels] scaling_correlation_matrix;
  matrix[n_channels, n_channels] scaling_cov_mat;
  cholesky_factor_cov[n_channels, n_channels] L_scaling_cov_mat;
  
  // Initialize the correlation matrix
  for (i in 1:n_channels) {
    for (j in 1:n_channels) {
      scaling_correlation_matrix[i, j] = 0;
    }
    scaling_correlation_matrix[i, i] = 1;
  }

  // Fill in the off-diagonal elements using precomputed indices
  for (n in 1:num_off_diagonals) {
    int i = row_idx[n];
    int j = col_idx[n];
    scaling_correlation_matrix[i, j] = scaling_channel_Rhos[n];
    scaling_correlation_matrix[j, i] = scaling_channel_Rhos[n];  // Symmetry
  }

  // Construct the covariance matrix using the standard deviations
  for (i in 1:n_channels) {
    for (j in 1:n_channels) {
      scaling_cov_mat[i, j] = scaling_correlation_matrix[i, j] * channel_scaling_sd[i] * channel_scaling_sd[j];
    }
  }
  // Cholesky decomposition for efficiency in sampling
  L_scaling_cov_mat = cholesky_decompose(scaling_cov_mat);

}

model {
  // Simple priors that don't involve channel
  sigma_average ~ normal(0,2);
  sigma_sd ~ normal(0,2);
  sigma ~ normal(sigma_average, sigma_sd); // Multilevel prior for each participant
  learning_rate_average ~ beta(1.1, 1.1);
  learning_rate_sd ~ beta(1.1, 1.1);
  learning_rate ~ normal(learning_rate_average, learning_rate_sd) T[0,1]; // Multilevel prior for each participant, truncated normal distribution

  // Priors that include channel take into account how channels are correlated
  // Set uninformative priors on MVNormal mean of intercept for all channels (A) and variance (B), perhaps these should also be multilevel.
  // The correlation between channels per trial is estimated from a multilevel prior, 
  // that in turn estimates the average per-trial correlation between all the channels (C).
  // Some of these priors are truncated normal distributions to keep them within the bounds
  // of correlations parameters (-1 to 1). The number of unique correlations to estimate is 
  // large (n_channels *(n_channels -1) / 2), so regularization is needed. Their indices are 
  // precomputed in the transformed data section. In the transformed parameters section, the 
  // average channel variance (channel_intercept_sd) and correlations (intercept_channel_Rhos) 
  // are used to estimate a covarinance matrix and its noncentered Cholesky transformation 
  // (L_intercept_cov_mat) that is easier to fit. This culminates in specifying that the 
  // intercept per participant and channel is sampled from a multivariate normal distribution (D).
  // A
  intercepts_average ~ normal(0, 1); //vector n_channels long
  // B
  channel_intercept_sd ~ normal(0, 2); //vector n_channels long
  // C Priors for the off-diagonal elements (lower triangular) correlation matrix
  intercept_Rho_average ~ normal(0,1) T[-1,1];
  intercept_Rho_sd ~ normal(0, 1);
  intercept_channel_Rhos ~ normal(intercept_Rho_average, intercept_Rho_sd) T[-1,1];
  // D
  for (p in 1:n_participants) {
    intercept[p] ~ multi_normal_cholesky(intercepts_average, L_intercept_cov_mat);
  }

  // The goal here was to a parameter (microvolt_scaling) per participant, cue, 
  // and channel that scales CS+ (CSP) associative strength to the brain amplitude 
  // of interest. The construction is similar to what was just done for intercept,
  // but with the added dimension of cue.

  // A
  cue_sd_average ~ normal(0,2);
  cue_sd_tau ~ normal(0,2);
  cue_sd ~ normal(cue_sd_average, cue_sd_tau);
  average_of_cues_all_chan ~ normal(0,1);
  average_of_cues_all_chan_sd ~ normal(0,2);
  average_of_cues_per_channel ~ normal(average_of_cues_all_chan, average_of_cues_all_chan_sd);
  for (q in 1:n_cues) {
    for (c in 1:n_channels) {
      cue_scaling_average[cue[q], channel[c]] ~ student_t(n_cues, average_of_cues_per_channel[channel[c]], cue_sd[channel[c]]);
    }
  }
  // B
  channel_scaling_sd ~ normal(0, 2);
  // C Priors for the off-diagonal elements (lower triangular) correlation matrix
  scaling_Rho_average ~ normal(0, 1) T[-1,1];
  scaling_Rho_sd ~ normal(0, 1);
  scaling_channel_Rhos ~ normal(scaling_Rho_average, scaling_Rho_sd) T[-1,1];
  
  // D
  // Prior for microvolt_scaling (deviations from intercept)
  for (p in 1:n_participants) {
    for (q in 1:n_cues) {
      microvolt_scaling[p, q] ~ multi_normal_cholesky(cue_scaling_average[q], L_scaling_cov_mat);
    }
  }
  // Likelihood
  for(i in 1:n_trials) {
    real mu = intercept[participant[i], channel[i]] + 
              (microvolt_scaling[participant[i], cue[i], channel[i]] * 
              CSP_associative_strength[i]);

    amplitude_all[i] ~ normal(mu, sigma[participant[i]]);

    // real wavelet = beta2[participant[i], channel[i]] * 
    //                (exp(-0.5 * square(k[n_trials] / sigma_wavelet)) * cos(frequency_wavelet * k[n_trials]));
    // real mu = beta1[participant[i], channel[i]] + wavelet;
    
    // real ricker = beta2[participant[i], channel[i]] * 
    //               ((1 - (2 * pow(pi(), 2) * pow(sigma_ricker[channel_id_vec[i]], 2) * pow(k[i], 2))) * 
    //               exp(-(2 * pow(pi(), 2) * pow(sigma_ricker[channel_id_vec[i]], 2) * pow(k[i], 2))));
    // real mu = beta1[participant[i], channel[i]] + ricker;
    
    // amplitude_all[i] ~ normal(mu, sigma[participant[i]])
  }
}

generated quantities {
  // These mu_pred and log_lik names cannot change
  
  // For model predictions
  array[n_observations] real mu_pred;
  // This is necessary for cross-validation
  array[n_observations] real log_lik;
  
  // for (i in 1:n_observations) {
  //   mu_pred[i] = (1 - (2 * pow(pi(), 2) * pow(sigma[channel_id_vec[i]], 2) * pow(k[i], 2))) * 
  //                exp(-(2 * pow(pi(), 2) * pow(sigma[channel_id_vec[i]], 2) * pow(k[i], 2)));
  //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma_error[channel_id_vec[i]]);
  // }
}
