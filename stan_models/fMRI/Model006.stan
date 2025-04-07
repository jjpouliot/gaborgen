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
  
  matrix[1070,1070] rho_power_matrix;

  for (i in 1:1070) {
    for (j in 1:1070) {
      rho_power_matrix[i,j] = abs(i-j);
    }
  }

  matrix[1070,1070] zero_diag_matrix;

  for (i in 1:1070) {
    for (j in 1:1070) {
      if (i == j) {
        zero_diag_matrix[i,j] = 0;
      } else {
        zero_diag_matrix[i,j] = 1;
      }
    }
  }
}


parameters {
  real <lower =0> sigma_raw; // this still needs to be half normal
  real delta_raw;
  real rho_raw;
  vector[n_DM_cols] Betas; //intercept 
  vector[1070-n_censor] alpha;
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;
  
  real sigma = sigma_raw *.5;
  real delta = inv_logit(delta_raw * 1.75);
  real rho = inv_logit(rho_raw * 1.75);

  matrix[1070,1070] Cor;
  
  Cor = add_diag((delta .* zero_diag_matrix),1) .* (rho .^ rho_power_matrix);
  
  matrix[1070 - n_censor, 1070 - n_censor] Cor_censor = Cor[censor, censor];
  
  matrix[1070 - n_censor, 1070 - n_censor] L_Cor_censor = cholesky_decompose(Cor_censor);
  
  

  vector[rows(DM_censored)] Mu = (DM_censored * Betas) + // matrix[n,p] times vector[p] = vector [p]
                           rep_vector(sigma, 1070 - n_censor) .* (L_Cor_censor * alpha); 

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
  
  sigma_raw ~ std_normal();
  delta_raw ~ std_normal();
  rho_raw ~ std_normal();
  Betas ~ std_normal();
  alpha ~ std_normal();
  
  amplitude ~ normal(Mu, 1);

  // sigma ~ normal(0,.5);
  // delta_raw ~ normal(0,1.75);
  // rho_raw ~ normal(0,1.75);
  // Betas[1] ~ normal(100,.5); //intercept 

    // matrix[1070, 1070] Cor;        // Correlation matrix
  // matrix[1070, 1070] L_Cor;      // Cholesky factor of the correlation matrix
  // matrix[1070, 1070] Cov;
  // // matrix[1070, 1070] L_Cov;
  // 
  // for (i in 1:1070) {
  //   for (j in 1:1070) {
  //     if (i == j) {
  //       Cov[i,j] = pow(sigma,2);
  //     } else {
  //       Cov[i,j] = delta * pow(rho, abs(i-j)) * pow(sigma, 2);
  //     }
  //   }
  // }
  
  // matrix[1070 - n_censor, 1070 - n_censor] Cov_censor = Cov[censor, censor];

  // Likelihood
  // for(p in 1:n_participants) {
  //   for(r in 1:n_ROIs) {
  // vector[rows(DM_censored)] Mu = DM_censored * Betas;
  //     // amplitude_all[participant[p], ROI[i],] ~ multi_normal_cholesky(Mu[i], Sigma_L[participant[i]]);
  //     amplitude ~ multi_normal(Mu, Cov_censor);
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
