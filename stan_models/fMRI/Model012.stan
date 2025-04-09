functions {
  // This function calculates the necessary log-likelihoods for each ROI
  // This function is put in map_rect so that each "shard" is a participant
  // This allows for parallel compututing such that multiple CPU threads can be used per chain
  // Matrix and vector operations are used so that GPU accelleration can be used with OpenCL
  vector participant_ll(vector phi, // reused for every shared
                        vector theta, // shard specific parameters
                        array[] real data_shard_real, // data (real numbers) used for this shard
                        array[] int data_shard_int) { // data (integers) used for this shard

  // unpack data_shard_int
  int n_bold = data_shard_int[1];
  int n_roi = data_shard_int[2];
  int n_censor = data_shard_int[3];
  int n_beta = data_shard_int[4];
  array[n_bold] int usable_bold_indices = data_shard_int[5:(n_bold+4)];

  array[n_bold - n_censor] int uncensored_indices;
  int keep_index = 1;
  for (i in 1:n_bold) {
    if (usable_bold_indices[i] == 1) {
      uncensored_indices[keep_index] = i;
      keep_index = keep_index + 1;
    }
  }

  // upack data_shard_real
  matrix[n_bold, n_roi] bold = to_matrix(data_shard_real[1:(n_bold * n_roi)], n_bold, n_roi);
  matrix[n_bold, n_beta] design_matrix = to_matrix(data_shard_real[(n_bold*n_roi + 1):(n_bold*n_roi + n_bold*n_beta)], n_bold, n_beta);


  // unpack parameters
  int start_index = 1;
  int stop_index = n_roi;
  array[n_roi] real sigma = to_array_1d(theta[start_index:stop_index]);

  start_index = stop_index + 1;
  stop_index = stop_index + n_roi;
  array[n_roi] real delta = to_array_1d(theta[start_index:stop_index]);

  start_index = stop_index + 1;
  stop_index = stop_index + n_roi;
  array[n_roi] real rho_time = to_array_1d(theta[start_index:stop_index]);
  
  array[n_roi] vector[n_beta] betas;
  for (r in 1:n_roi) {
    start_index = stop_index + 1;
    stop_index = stop_index + n_beta;
    betas[r] = theta[start_index:stop_index];
  }

  // create operation matrices  
  matrix[n_bold, n_bold] zero_diag_matrix =
    rep_matrix(1, n_bold, n_bold) - diag_matrix(rep_vector(1, n_bold));

  vector[n_bold] idx = linspaced_vector(n_bold, 1, n_bold);
  
  // Create a matrix where each row is the index i.
  // v * rep_row_vector(1, n_bold) produces an n_bold x n_bold matrix
  // where row i consists of the value idx[i] repeated n_bold times.
  matrix[n_bold, n_bold] row_idx = idx * rep_row_vector(1, n_bold);
  
  // Similarly, rep_column_vector(1, n_bold) * idx'
  // produces an n_bold x n_bold matrix where column j consists of the value idx[j] repeated.
  matrix[n_bold, n_bold] col_idx = rep_vector(1, n_bold) * idx';
  
  // Now compute the matrix of absolute differences:
  matrix[n_bold, n_bold] rho_power_matrix = abs(row_idx - col_idx);
  

  // Create mu predictions and covariance matrices per roi, then calc likelihood
  matrix[n_bold - n_censor, n_beta] design_matrix_censor = design_matrix[uncensored_indices,];
  array[n_roi] vector[n_bold - n_censor] Mu; 
  array[n_roi] matrix[n_bold, n_bold] Cov;
  array[n_roi] matrix[n_bold - n_censor, n_bold - n_censor] Cov_censor;
  array[n_roi] matrix[n_bold - n_censor, n_bold - n_censor] L_Cov_censor;

  matrix[n_bold - n_censor, n_roi] bold_censored = bold[uncensored_indices,1:n_roi];

  vector[1] log_lik;
  
  for (r in 1:n_roi) {
    Mu[r]= (design_matrix_censor * betas[r]);

    Cov[r] = pow(sigma[r], 2) .* add_diag((delta[r] .* zero_diag_matrix),1) .* (rho_time[r] .^ rho_power_matrix);
  
    Cov_censor[r] = Cov[r][uncensored_indices, uncensored_indices];

    L_Cov_censor[r] = cholesky_decompose(Cov_censor[r]);

    log_lik += multi_normal_cholesky_lpdf(bold_censored[,r] | Mu[r], L_Cov_censor[r]);

  }

  return log_lik;
}  
}



data {
  int<lower=0> n;
  int<lower=0> n_par;
  int<lower=0> n_roi;
  int<lower=0> n_bold;
  int<lower=0> n_beta; 

  array[n] int par;
  array[n] int roi;
  array[n_par, n_bold*n_roi] real bold;
  array[n_par, n_bold] int usable_bold_indices;
  array[n_par] int n_censor;
  array[n_par, n_bold, n_beta] real design_array;
}

transformed data {
  // these matrices are used for more efficient operations
  // matrix[n_bold, n_bold] zero_diag_matrix;

  // for (i in 1:n_bold) {
  //   for (j in 1:n_bold) {
  //     if (i == j) {
  //       zero_diag_matrix[i,j] = 0;
  //     } else {
  //       zero_diag_matrix[i,j] = 1;
  //     }
  //   }
  // }
  
  // matrix[n_bold,n_bold] rho_power_matrix;

  // for (i in 1:n_bold) {
  //   for (j in 1:n_bold) {
  //     rho_power_matrix[i,j] = abs(i-j);
  //   }
  // }



  // For map_rect: 
  // create data_shard_real for each participant
  // array[n_par, n_bold * n_roi + n_bold * n_beta]
  array[n_par, n_bold * n_roi + n_bold * n_beta] real data_shard_real;
  for (p in 1:n_par) {
    data_shard_real[p, 1:(n_bold * n_roi)] = bold[p, ];
    // I think this is row-major order 
    // design_array[p,1,1]
    // design_array[p,1,2]
    // design_array[p,2,1]
    // design_array[p,2,2]
    data_shard_real[p, (n_bold*n_roi + 1):(n_bold*n_roi + n_bold*n_beta)] = to_array_1d(design_array[p,,]);
  }



  // row size: n_bold rows
  // column size: n_roi + n_beta + rows(zero_diag_matrix) + rows(rho_power_matrix)
  // bold per roi
  // design matrix
  // operation matrices
  // array[n_par] matrix[n_bold,n_roi+n_beta+rows(zero_diag_matrix)+rows(rho_power_matrix)] data_shard_real;
  // for (p in 1:n_par) {
  //   for (r in 1:n_roi) {
  //     data_shard_real[p][,r] = to_vector(bold[p, r, ]);
  //   }
  //   data_shard_real[p][, (n_roi + 1):(n_roi + n_beta)] = to_matrix(design_array[p, , ]);

  //   data_shard_real[p][, (n_roi + n_beta + 1):(n_roi + n_beta + rows(zero_diag_matrix))] = 
  //     zero_diag_matrix;
      
  //   data_shard_real[p][, (n_roi+n_beta+rows(zero_diag_matrix)+1):(n_roi+n_beta+rows(zero_diag_matrix)+rows(rho_power_matrix))] = 
  //     rho_power_matrix;

  // }

  // create data_shard_int for each participant
  // array[n_bold + 3]
  // n_bold
  // n_roi
  // n_censor
  // censor
  array[n_par, n_bold + 4] int data_shard_int;
  for(p in 1:n_par){
    data_shard_int[p, 1] = n_bold;
    data_shard_int[p, 2] = n_roi;
    data_shard_int[p, 3] = n_censor[p];
    data_shard_int[p, 4] = n_beta;
    data_shard_int[p, 5:(n_bold + 4)] = usable_bold_indices[p,];
  }

}


parameters {
  // real <lower =0> mu_sigma_raw; // this still needs to be half normal
  // real mu_delta_raw;
  // real mu_rho_time_raw;
  // array[n_roi] vector[n_DM_cols] mu_beta;

  array[n_par, n_roi] real <lower =0> sigma_raw; // this still needs to be half normal
  array[n_par, n_roi] real delta_raw;
  array[n_par, n_roi] real rho_time_raw;
  array[n_par, n_roi] vector[n_beta] betas;

  // array[n_par, n_roi] real <lower =0> sigma_raw; // this still needs to be half normal
  // array[n_par, n_roi] real delta_raw;
  // array[n_par, n_roi] real rho_time_raw;
  // array[n_par, n_roi] vector[n_DM_cols] Betas;
  
  
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;

  array[n_par, n_roi] real sigma;
  array[n_par, n_roi] real delta;
  array[n_par, n_roi] real rho_time;

  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma[p, r] = sigma_raw[p, r] *.5;
      delta[p, r] = inv_logit(delta_raw[p, r] * 1.75);
      rho_time[p, r] = inv_logit(rho_time_raw[p, r] * 1.75);
    }
  }

  // real mu_sigma = mu_sigma_raw *.5;
  // real mu_delta = inv_logit(mu_delta_raw * 1.75);
  // real mu_rho_time = inv_logit(mu_rho_time_raw * 1.75);

  // matrix[n_roi, n_roi] Cov_sigma;
  // matrix[n_roi, n_roi] Cov_delta;
  // matrix[n_roi, n_roi] Cov_rho_time;
  // array[n_betas] matrix[n_roi, n_roi] Cov_betas;
  
  // for (r_1 in 1:n_roi){
  //   for (r_2 in 1:n_roi) {
  //     if (r_1 == r_2) {
  //       Cov_sigma[r_1, r_2] = sigma_sigma[r_1] * sigma_sigma[r_2];
  //       Cov_delta[r_1, r_2] = delta_sigma[r_1] * delta_sigma[r_2];
  //       Cov_rho_time[r_1, r_2] = rho_time_sigma[r_1] * rho_time_sigma[r_2];
  //       for (b in 1:n_betas) {
  //         Cov_betas[r_1, r_2][b] = beta_sigma[r_1, b] * beta_sigma[r_2, b];
  //       }
  //     } else {
  //       Cov_sigma[r_1, r_2] = sigma_sigma[r_1] * 
  //                             sigma_sigma[r_2] *
  //                             sigma_rho;
  //       Cov_delta[r_1, r_2] = delta_sigma[r_1] * 
  //                             delta_sigma[r_2] *
  //                             delta_rho;
  //       Cov_rho_time[r_1, r_2] = rho_time_sigma[r_1] * 
  //                                rho_time_sigma[r_2] *
  //                                rho_time_rho;
  //       for (b in 1:n_betas) {
  //         Cov_betas[r_1, r_2][b] = beta_sigma[r_1, b] * 
  //                                  beta_sigma[r_2, b] *
  //                                  beta_rho[b];
  //       }
  //     }
  //   }
  // }
  
  // matrix[n_roi, n_roi] L_sigma = cholesky_decompose(Cov_sigma);
  // matrix[n_roi, n_roi] L_delta = cholesky_decompose(Cov_delta);
  // matrix[n_roi, n_roi] L_rho_time = cholesky_decompose(Cov_rho_time);
  // array[n_betas] matrix[n_roi, n_roi] L_betas;
  // for (b in 1:n_betas) {
  //   L_betas[1070,1070][b] = cholesky_decompose(Cov_betas);
  // }

  // // map_rect here?
  // array[n_par, n_roi] vector[rows(DM_censored)] Mu; 
  // array [n_par, n_roi] matrix[1070, 1070] Cov_time;
  // array [n_par, n_roi] matrix [1070, 1070] L_time; // have to pass in censor at some point
  
  // for (p in 1:n_par) {
  //   for (r in 1:n_roi) {
  //     Mu[p, r] = (DM_censored[p, r] * Betas[p, r]); 
      
  //     Cov_time[p, r] = pow(sigma[p, r], 2) .* add_diag((delta[p, r] .* zero_diag_matrix),1) .* (rho_time[p, r] .^ rho_time_power_matrix);
      
  //     L_time[p, r] = cholesky_decompose(Cov_time[p, r]);
  //   }
  // }

}

model {
  // Priors
  //vectors per roi?
  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma_raw[p, r] ~ std_normal();
      delta_raw[p, r] ~ std_normal();
      rho_time_raw[p, r] ~ std_normal();
      betas[p, r] ~ std_normal();
    }
  }

  // local data shard parameters theta for each participant
  // vector[n_roi * 3 + n_roi * n_beta] because there is a parameter per roi
  // sigma 
  // delta
  // rho_time
  // betas
  array[n_par] vector[(n_roi * 3) + (n_roi * n_beta)] theta;

  for (p in 1:n_par) {
    int start_index = 1;
    int stop_index = n_roi;
    theta[p][start_index:stop_index] = to_vector(sigma[p,]);

    start_index = stop_index + 1;
    stop_index = stop_index + n_roi;
    theta[p][start_index:stop_index] = to_vector(delta[p,]);

    start_index = stop_index + 1;
    stop_index = stop_index + n_roi;
    theta[p][start_index:stop_index] = to_vector(rho_time[p,]);

    for (r in 1:n_roi) {
      start_index = stop_index + 1;
      stop_index = stop_index + n_beta;
      theta[p][start_index:stop_index] = to_vector(betas[p,r]);
    }
  }

vector[0] phi;

target += sum(map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int));

  // rho_d ~ std_normal();
  // mu_sigma_raw ~ std_normal();
  // mu_delta_raw ~ std_normal();
  // mu_rho_time_raw ~ std_normal();
  // mu_betas ~ std_normal();
  // mu_rho_d ~ std_normal();
  

  
  
  // multilevel priors
  // for(p in 1:n_participants) {
  //   // capturing correlations between rois per participant
  //   // should I be taking raw variants here?
  //   sigma_raw[p, ] ~ multi_normal_cholesky(mu_sigma_raw, L_sigma);
  //   delta_raw[p, ] ~ multi_normal_cholesky(mu_delta_raw, L_delta);
  //   rho_time_raw[p, ] ~ multi_normal_cholesky(mu_rho_time_raw, L_rho_time);
  //   for (b in 1:n_betas)
  //     vector[n_roi] Beta_vec = to_vector(Betas[p, ][b]);
  //     Beta_vec ~ multi_normal_cholesky(mu_betas[b], L_betas[b]);
  // }
  
  // // likelihood
  // for(p in 1:n_participants) {
  //   for(r in 1:n_roi) {
  //     // map_rect with shards of p and r?
  //     // or just reduce sum this if it remains multi_normal
  //     matrix [1070-n_censor[p], 1070-n_censor[p]] L = L_time[censor[p], censor[p]];
  //     bold[p, r, ] ~ multi_normal_cholesky(Mu[p, r], L);
  //   }
  // }
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
  // target += sum(map_rect(participant_roi_ll, phi, theta_shards, data_shards, int_data_shards));


}

// generated quantities {
//   // array[n_observations] real mu_pred = mu[indices_observed];
//   // array[n_observations] real log_lik;
  
//   // for (i in 1:n_observations) {
//   //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
//   // }
// }
