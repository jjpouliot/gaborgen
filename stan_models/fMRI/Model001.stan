data {
  int<lower=0> n_amplitude;
  int<lower=0> n_censor;
  int<lower=0> n_DM_cols;

  matrix[1070, n_DM_cols] design_matrix;
  vector[n_amplitude] amplitude_no_censor;
  array[n_amplitude - n_censor] int censor;



}

transformed data {
  matrix[1070 - n_censor, n_DM_cols] DM_censored = design_matrix[censor,];

  vector[n_amplitude - n_censor] amplitude = amplitude_no_censor[censor];
}

parameters {
  real <lower=0> sigma;
  real delta_raw;
  real rho_raw;
  vector[n_DM_cols] Betas; //intercept 
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;

  real delta = inv_logit(delta_raw);
  real rho = inv_logit(rho_raw);


  // matrix[1070, 1070] Cor;        // Correlation matrix
  // matrix[1070, 1070] L_Cor;      // Cholesky factor of the correlation matrix
  // matrix[1070, 1070] Cov;
  // matrix[1070, 1070] L_Cov;
  
  // for (i in 1:1070) {
  //   for (j in 1:1070) {
  //     // This is the AR(1)-style correlation:  Cor(i,j) = rho^|i-j|
  //     if (i == j) {
  //       Cor[i,j] = 1;
  //     } else {
  //       Cor[i,j] = delta * pow(rho, abs(i-j));
  //     }
  //   }
  // }
  
  // L_Cor = cholesky_decompose(Cor);

  // // Build the diagonal scaling matrix with σ (standard deviation).
  // // Since σ is already the standard deviation, we don't square it.
  // matrix[1070,1070] D = diag_matrix(rep_vector(sigma, 1070));
  
  // // Form the covariance matrix: Cov = D * Cor * D.
  // // Cov = D * Cor * D;
  
  // // Obtain the Cholesky factor of the covariance matrix:
  // // L_Cov = D * L_Cor satisfies L_Cov * L_Cov' = Cov.
  // L_Cov = D * L_Cor;

  // matrix[1070 - n_censor, 1070 - n_censor] L_Cov_censor = L_Cov[censor, censor];

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
  Betas[1] ~ normal(100,.5); //intercept 
  Betas[2:n_DM_cols] ~ normal(0,.5);



  matrix[1070, 1070] Cor;        // Correlation matrix
  matrix[1070, 1070] L_Cor;      // Cholesky factor of the correlation matrix
  // matrix[1070, 1070] Cov;
  matrix[1070, 1070] L_Cov;
  
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

  // Build the diagonal scaling matrix with σ (standard deviation).
  // Since σ is already the standard deviation, we don't square it.
  matrix[1070,1070] D = diag_matrix(rep_vector(sigma, 1070));
  
  // Form the covariance matrix: Cov = D * Cor * D.
  // Cov = D * Cor * D;
  
  // Obtain the Cholesky factor of the covariance matrix:
  // L_Cov = D * L_Cor satisfies L_Cov * L_Cov' = Cov.
  L_Cov = D * L_Cor;

  matrix[1070 - n_censor, 1070 - n_censor] L_Cov_censor = L_Cov[censor, censor];


  // Likelihood
  // for(p in 1:n_participants) {
  //   for(r in 1:n_ROIs) {
  vector[rows(DM_censored)] Mu = DM_censored * Betas;
      // amplitude_all[participant[p], ROI[i],] ~ multi_normal_cholesky(Mu[i], Sigma_L[participant[i]]);
      amplitude ~ multi_normal_cholesky(Mu, L_Cov_censor);
    // }
  // }



}

// generated quantities {
//   // array[n_observations] real mu_pred = mu[indices_observed];
//   // array[n_observations] real log_lik;
  
//   // for (i in 1:n_observations) {
//   //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
//   // }
// }
