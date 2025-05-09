
functions {
  // This function 
  vector participant_roi_ll(vector theta, // reused for every shared
                            vector phi, // shard specific parameters
                            array[] real data_shard, // data (real numbers) used for this shard
                            array[] int int_data_shard) { // data (integers) used for this shard
    // Expected structure of theta:
    //   theta[1:K]         : beta coefficients (K = n_DM_cols)
    //   theta[K+1]         : sigma (standard deviation)
    //   theta[K+2]         : delta
    //   theta[K+3]         : rho_time
    int n_bold = int_data_shard[1];       // number of BOLD observations for this shard
    int K = int_data_shard[2];            // number of DM columns

    // Reconstruct the design matrix DM for this shard.
    // We assume the first (n_bold*K) entries of data_shard are the design matrix (row-major order)
    matrix[n_bold, K] DM;
    {
      int pos = 1;
      for (i in 1:n_bold)
        for (j in 1:K) {
          DM[i, j] = data_shard[pos];
          pos += 1;
        }
    }
    
    // Next, the remaining n_bold entries of data_shard are the BOLD time series.
    vector[n_bold] bold;
    {
      int offset = n_bold * K;
      for (i in 1:n_bold)
        bold[i] = data_shard[offset + i];
    }
    
    // Compute the mean vector (Mu) for the BOLD data
    vector[n_bold] mu = DM * theta[1:K];
    
    // Extract shard-specific parameters
    real sigma    = theta[K+1];
    real delta    = theta[K+2];
    real rho_time = theta[K+3];
    
    // Construct the covariance matrix Cov_time.
    // For illustration, we build a simple structure using the supplied parameters.
    matrix[n_bold, n_bold] Cov_time;
    for (i in 1:n_bold) {
      for (j in 1:n_bold) {
        // For example, let the off-diagonals be modulated by delta and rho_time:
        if (i == j) {
          Cov_time[i, j] = square(sigma);
        } else {
          Cov_time[i, j] = square(sigma) * delta * pow(rho_time, fabs(i - j));
        }
      }
    }
    
    // Compute the Cholesky factor of Cov_time.
    matrix[n_bold, n_bold] L_time = cholesky_decompose(Cov_time);
    
    // Compute and return the log likelihood for this shard.
    // Note: multi_normal_cholesky_lpdf returns a real, but we need to return a vector.
    return [ multi_normal_cholesky_lpdf(bold | mu, L_time) ]';
  }
}


data {
  int<lower=0> n_bold;
  int<lower=0> n_censor;
  int<lower=0> n_par;
  int<lower=0> n_roi;
  int<lower=0> n_DM_cols; 

  matrix[1070, n_DM_cols] design_matrix;
  array[n_par, n_roi] vector[n_bold] bold;
  array[n_par, n_roi, n_bold] int censor;
}

transformed data {

  matrix[1070 - n_censor, n_DM_cols] DM_censored = design_matrix[censor[],];

  vector[n_bold - n_censor] amplitude = amplitude_no_censor[censor];
  
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
  array[n_roi] vector[n_DM_cols] mu_beta;
  
  array[n_par, n_roi] real <lower =0> sigma_raw; // this still needs to be half normal
  array[n_par, n_roi] real delta_raw;
  array[n_par, n_roi] real rho_time_raw;
  array[n_par, n_roi] vector[n_DM_cols] Betas;
  
  
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;
  
  real mu_sigma = mu_sigma_raw *.5;
  real mu_delta = inv_logit(mu_delta_raw * 1.75);
  real mu_rho_time = inv_logit(mu_rho_time_raw * 1.75);

  matrix[n_roi, n_roi] Cov_sigma;
  matrix[n_roi, n_roi] Cov_delta;
  matrix[n_roi, n_roi] Cov_rho_time;
  array[n_betas] matrix[n_roi, n_roi] Cov_betas;
  
  for (r_1 in 1:n_roi){
    for (r_2 in 1:n_roi) {
      if (r_1 == r_2) {
        Cov_sigma[r_1, r_2] = sigma_sigma[r_1] * sigma_sigma[r_2];
        Cov_delta[r_1, r_2] = delta_sigma[r_1] * delta_sigma[r_2];
        Cov_rho_time[r_1, r_2] = rho_time_sigma[r_1] * rho_time_sigma[r_2];
        for (b in 1:n_betas) {
          Cov_betas[r_1, r_2][b] = beta_sigma[r_1, b] * beta_sigma[r_2, b];
        }
      } else {
        Cov_sigma[r_1, r_2] = sigma_sigma[r_1] * 
                              sigma_sigma[r_2] *
                              sigma_rho;
        Cov_delta[r_1, r_2] = delta_sigma[r_1] * 
                              delta_sigma[r_2] *
                              delta_rho;
        Cov_rho_time[r_1, r_2] = rho_time_sigma[r_1] * 
                                 rho_time_sigma[r_2] *
                                 rho_time_rho;
        for (b in 1:n_betas) {
          Cov_betas[r_1, r_2][b] = beta_sigma[r_1, b] * 
                                   beta_sigma[r_2, b] *
                                   beta_rho[b];
        }
      }
    }
  }
  
  matrix[n_roi, n_roi] L_sigma = cholesky_decompose(Cov_sigma);
  matrix[n_roi, n_roi] L_delta = cholesky_decompose(Cov_delta);
  matrix[n_roi, n_roi] L_rho_time = cholesky_decompose(Cov_rho_time);
  array[n_betas] matrix[n_roi, n_roi] L_betas;
  for (b in 1:n_betas) {
    L_betas[1070,1070][b] = cholesky_decompose(Cov_betas);
  }

  // map_rect here?
  array[n_par, n_roi] vector[rows(DM_censored)] Mu; 
  array [n_par, n_roi] matrix[1070, 1070] Cov_time;
  array [n_par, n_roi] matrix [1070, 1070] L_time; // have to pass in censor at some point
  
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
  mu_betas ~ std_normal();
  
  
  // multilevel priors
  for(p in 1:n_participants) {
    // capturing correlations between rois per participant
    // should I be taking raw variants here?
    sigma_raw[p, ] ~ multi_normal_cholesky(mu_sigma_raw, L_sigma);
    delta_raw[p, ] ~ multi_normal_cholesky(mu_delta_raw, L_delta);
    rho_time_raw[p, ] ~ multi_normal_cholesky(mu_rho_time_raw, L_rho_time);
    for (b in 1:n_betas)
      vector[n_roi] Beta_vec = to_vector(Betas[p, ][b]);
      Beta_vec ~ multi_normal_cholesky(mu_betas[b], L_betas[b]);
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
  // or
  // Assume:
  //   n_shards = n_par * n_roi
  //   theta_dim = n_DM_cols + 3; // K beta coefficients + sigma + delta + rho_time

  // You would typically build these arrays in transformed data or via helper functions.
  // For illustration, assume we have pre-built:
  //   vector[theta_dim] theta_shards[n_shards];       // parameters for each shard
  //   real data_shards[n_shards][data_length];          // data for each shard
  //   int int_data_shards[n_shards][int_data_length];     // integer meta-data for each shard

  // Set priors for global parameters here (or on the individual theta_shards if modeled hierarchically).

  // Parallelize the likelihood across shards using map_rect.
  target += sum(map_rect(participant_roi_ll, theta_shards, data_shards, int_data_shards));


}

// generated quantities {
//   // array[n_observations] real mu_pred = mu[indices_observed];
//   // array[n_observations] real log_lik;
  
//   // for (i in 1:n_observations) {
//   //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
//   // }
// }
