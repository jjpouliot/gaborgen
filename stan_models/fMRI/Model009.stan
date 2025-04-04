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
  real <lower =0> mu_sigma_raw; // this still needs to be half normal
  real mu_delta_raw;
  real mu_rho_time_raw;
  array[n_participants, vector[n_DM_cols] Betas; //intercept 
  // vector[1070-n_censor] alpha;
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;
  
  real mu_sigma = mu_sigma_raw *.5;
  real mu_delta = inv_logit(mu_delta_raw * 1.75);
  real mu_rho_time = inv_logit(mu_rho_time_raw * 1.75);



  // map_rect here?
  array[n_par, n_roi] vector[rows(DM_censored)] Mu; 
  array [n_par, n_roi] matrix[1070, 1070] Cov_time;
  array [n_par, n_roi] [1070, 1070] L_time; // have to pass in censor at some point
  
  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      Mu[p, r] = (DM_censored[p, r] * Betas[p, r]); 
      
      Cov_time[p, r] = pow(sigma[p, r], 2) .* add_diag((delta[p, r] .* zero_diag_matrix),1) .* (rho_time[p, r] .^ rho_time_power_matrix);
      
      L_time[p, r] = cholesky_decompose(Cov_time[p, r]);
    }
  }

}

model {
  // Priors
  //vectors per roi?
  mu_sigma_raw ~ std_normal();
  mu_delta_raw ~ std_normal();
  mu_rho_time_raw ~ std_normal();
  mu_rho_d ~ std_normal();
  mu_B ~ std_normal();
  
  
  // multilevel priors
  tau
  for(p in 1:n_participants) {
    // capturing correlations between rois per participant
    // should I be taking raw variants here?
    sigma[p, ] ~ multi_normal_cholesky(mu_sigma, L_sigma);
    delta[p, ] ~ multi_normal_cholesky(mu_delta, L_delta);
    rho_time[p, ] ~ multi_normal_cholesky(mu_rho_time, L_rho_time);
    for (b in 1:n_betas)
      Betas[p, ][b] ~ ~ multi_normal_cholesky(mu_B[b], L_B[b]);
  }
  
  // likelihood
  for(p in 1:n_participants) {
    for(r in 1:n_roi) {
      // map_rect with shards of p and r?
      // or just reduce sum this if it remains multi_normal
      matrix [1070-n_censor[p], 1070-n_censor[p]] L = L_time[censor[p], censor[p]];
      bold[p, r, ] ~ multi_normal_cholesky(Mu[p, r], L);
    }
  }


}

// generated quantities {
//   // array[n_observations] real mu_pred = mu[indices_observed];
//   // array[n_observations] real log_lik;
  
//   // for (i in 1:n_observations) {
//   //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
//   // }
// }
