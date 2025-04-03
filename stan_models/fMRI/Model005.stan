functions {
  // This function computes the log-likelihood for a slice of observations.
  // "amplitude_slice" and "f_array" are arrays corresponding to a slice of the data.
  // "eps" is a tiny noise value.
  real log_lik_reduce(array[] real amplitude_slice,
                      int start, int end,
                      array[] real f_array,
                      real eps) {
    real lp = 0;
    for (i in 1:size(amplitude_slice)) {
      int global_index = start + i - 1;
      lp += normal_lpdf(amplitude_slice[i] | f_array[global_index], eps);
    }
    return lp;
  }
}

data {
  int<lower=0> n_amplitude;
  int<lower=0> n_censor;
  int<lower=0> n_DM_cols;
  matrix[1070, n_DM_cols] design_matrix;
  vector[n_amplitude] amplitude_no_censor;
  array[n_amplitude - n_censor] int censor;
}

transformed data {
  // Apply censoring.
  matrix[1070 - n_censor, n_DM_cols] DM_censored = design_matrix[censor,];
  // Convert the censored amplitude vector to an array.
  array[n_amplitude - n_censor] real amplitude = to_array_1d(amplitude_no_censor[censor]);
}

parameters {
  real<lower=0> sigma;
  real delta_raw;
  real rho_raw;
  vector[n_DM_cols] Betas;
  // Latent standard normals for the correlated process.
  vector[1070 - n_censor] eta;
}

transformed parameters {
  real delta = inv_logit(delta_raw);
  real rho   = inv_logit(rho_raw);
  
  // Build full AR(1)-style correlation matrix for 1070 time points.
  matrix[1070, 1070] Cor;
  matrix[1070, 1070] L_Cor;
  for (i in 1:1070) {
    for (j in 1:1070) {
      if (i == j)
        Cor[i,j] = 1;
      else
        Cor[i,j] = delta * pow(rho, abs(i - j));
    }
  }
  L_Cor = cholesky_decompose(Cor);
  
  // Extract the censored portion.
  matrix[1070 - n_censor, 1070 - n_censor] L_Cor_censor = L_Cor[censor, censor];
  
  // Fixed effects: linear predictor.
  vector[1070 - n_censor] Mu = DM_censored * Betas;
  
  // Here is the key non-centered parameterization:
  // f has marginal distribution MVN(Mu, sigma^2 * Cor_censor).
  vector[1070 - n_censor] f = Mu + sigma * (L_Cor_censor * eta);
  
  // Convert f to an array for reduce_sum.
  // (If you have Stan version 2.26 or later, to_array_1d is available.)
  array[1070 - n_censor] real f_array = to_array_1d(f);
}

model {
  // Priors.
  sigma ~ normal(0, 0.5);
  delta_raw ~ normal(0, 1.75);
  rho_raw ~ normal(0, 1.75);
  Betas ~ normal(0, 1);
  eta ~ std_normal();
  
  // Likelihood via reduce_sum.
  // Use a very small eps (e.g., 1e-6) so that the observation is essentially equal to f.
  target += reduce_sum(log_lik_reduce,
                       amplitude,      // the sliced amplitude array
                       1,              // grainsize (can be tuned)
                       f_array,        // shared latent function values
                       1e-6);         // nearly zero observation noise
}
