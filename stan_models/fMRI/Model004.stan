functions {
  // This function computes the log-likelihood contribution for observations
  // from index "start" to "end". The first argument is the slice of spike.
  real log_lik_reduce(array[] real amplitude,
                      int start, int end,
                      array[] real f,
                      real sigma) {
    real lp = 0;
    // Loop over the slice; size(spike) equals (end - start + 1)
    for (i in 1:size(amplitude)) {
      int global_index = start + i - 1;
      lp += normal_lpdf(amplitude[i] | f[global_index], sigma);
    }
    return lp;
  }
}

data {
  int<lower=0> n_amplitude;
  int<lower=0> n_censor;
  int<lower=0> n_DM_cols;

  matrix[1070, n_DM_cols] design_matrix;
  array[n_amplitude] real amplitude_no_censor;
  array[n_amplitude - n_censor] int censor;



}

transformed data {
  matrix[1070 - n_censor, n_DM_cols] DM_censored = design_matrix[censor,];

  array[n_amplitude - n_censor] real amplitude = amplitude_no_censor[censor];
}

parameters {
  real <lower=0> sigma;
  real delta_raw;
  real rho_raw;
  vector[n_DM_cols] Betas; //intercept 
  vector[1070 - n_censor] eta;
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;

  real delta = inv_logit(delta_raw);
  real rho = inv_logit(rho_raw);



}

model {
  // Priors
  // for(p in 1:n_participants) {
  //   for(r in 1:n_ROIs) {
  //     sigma[p,r] ~ normal(0,.5);
  //     delta_raw[p,r] ~ normal(0,1.75);
  //     rho_raw[p,r] ~ normal(0,1.75);
  //   }
  // }
  
  sigma ~ normal(0,.5);
  delta_raw ~ normal(0,1.75);
  rho_raw ~ normal(0,1.75);
  // Betas[1] ~ normal(100,.5); //intercept 
  Betas ~ normal(0,1);



  matrix[1070, 1070] Cor;        // Correlation matrix
  matrix[1070, 1070] L_Cor;      // Cholesky factor of the correlation matrix
  // matrix[1070, 1070] Cov;
  // matrix[1070, 1070] L_Cov;
  
  for (i in 1:1070) {
    for (j in 1:1070) {
      // This is the AR(1)-style correlation:  Cor(i,j) = rho^|i-j|
      if (i == j) {
        Cor[i,j] = 1;
      } else {
        Cor[i,j] = delta * pow(rho, abs(i-j));
      }
    }
  }
  
  L_Cor = cholesky_decompose(Cor);


  matrix[1070 - n_censor, 1070 - n_censor] L_Cor_censor = L_Cor[censor, censor];

  eta ~ std_normal();
  
  vector[rows(DM_censored)] Mu = DM_censored * Betas;

  array[1070-n_censor] real f = to_array_1d(Mu + (L_Cor_censor * eta));



  // Likelihood
  target += reduce_sum(log_lik_reduce,
                       amplitude,               // This array will be sliced.
                       1,                   // Grainsize.
                       f,  // Shared data (not sliced automatically).
                       sigma);

}

// generated quantities {
//   // array[n_observations] real mu_pred = mu[indices_observed];
//   // array[n_observations] real log_lik;
  
//   // for (i in 1:n_observations) {
//   //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
//   // }
// }
