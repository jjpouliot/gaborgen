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
  array[n_bold] int uncensored_indices_padded = data_shard_int[5:(n_bold+4)];
  array[n_bold - n_censor] int uncensored_indices = uncensored_indices_padded[1:(n_bold - n_censor)];


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

  // real lp_sigma = 0;
  real lp_mvn = 0;
  
  for (r in 1:n_roi) {

    Mu[r]= (design_matrix_censor * betas[r]);

    Cov[r] = pow(sigma[r], 2) .* add_diag((delta[r] .* zero_diag_matrix),1) .* (rho_time[r] .^ rho_power_matrix);
  
    Cov_censor[r] = Cov[r][uncensored_indices, uncensored_indices];

    L_Cov_censor[r] = cholesky_decompose(Cov_censor[r]);

    lp_mvn += multi_normal_cholesky_lpdf(bold_censored[,r] | Mu[r], L_Cov_censor[r]);

  }

  return [lp_mvn]';
}  
}



data {
  int<lower=0> n;
  int<lower=0> n_par;
  int<lower=0> n_roi;
  int<lower=0> n_bold;
  int<lower=0> n_beta; 
  // int<lower=0> n_beta_stim; // assuming that the first n_beta_stim coefficients in the desing matrix need a multilevel prior 

  array[n] int par;
  array[n] int roi;
  array[n_par, n_bold*n_roi] real bold;
  array[n_par, n_bold] int usable_bold_indices_one_is_true;
  array[n_par] int n_censor;
  array[n_par, n_bold, n_beta] real design_array;
}

transformed data {

  // map_rect needs inputs to be same size, so zeros are padded at end per participant
  array[n_par, n_bold] int uncensored_indices_padded;
  for (p in 1:n_par) {
    int keep_index = 1;
    for (i in 1:n_bold) {
      if (usable_bold_indices_one_is_true[p,i] == 1) {
        uncensored_indices_padded[p, keep_index] = i;
        keep_index = keep_index + 1;
      }
    }

    if(n_censor[p] > 0) {
      array [n_censor[p]] int pad = rep_array(0, n_censor[p]);
      uncensored_indices_padded[p, keep_index:n_bold] = pad;
    }
  }


  // For map_rect: 
  // create data_shard_real for each participant
  // array[n_par, n_bold * n_roi + n_bold * n_beta]
  array[n_par, n_bold * n_roi + n_bold * n_beta] real data_shard_real;
  for (p in 1:n_par) {
    data_shard_real[p, 1:(n_bold * n_roi)] = bold[p, ];
    // pack from column-major order 
    // data_shard[p, 1] = to_array_1d(to_matrix(design_array[p,1,1]))
    // data_shard[p, 2] = to_array_1d(to_matrix(design_array[p,2,1]))
    // data_shard[p, 3] = to_array_1d(to_matrix(design_array[p,3,1]))
    // ...
    // data_shard[p, n_bold + 1] = to_array_1d(to_matrix(design_array[p,1,2]))
    // data_shard[p, n_bold + 2] = to_array_1d(to_matrix(design_array[p,2,2]))
    // data_shard[p, n_bold + 3] = to_array_1d(to_matrix(design_array[p,3,2]))
    data_shard_real[p, (n_bold*n_roi + 1):(n_bold*n_roi + n_bold*n_beta)] = to_array_1d(to_matrix(design_array[p,,]));
  }

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
    data_shard_int[p, 5:(n_bold + 4)] = uncensored_indices_padded[p,];
  }

}


parameters {
  array[n_roi] real <lower =0> mu_sigma_raw; // this still needs to be half normal
  array[n_roi] real mu_delta_raw;
  array[n_roi] real mu_rho_time_raw;
  array[n_roi] vector[n_beta] mu_betas;
  // array[n_roi] vector[n_beta_stim] mu_betas;

  // think about what the priors should be fore these
  array[n_roi] real <lower =0> tau_sigma_raw;
  array[n_roi] real <lower =0>  tau_delta_raw;
  array[n_roi] real <lower =0> tau_rho_time_raw;
  array[n_roi, n_beta] real <lower =0> tau_betas;
  // array[n_roi] vector[n_beta_stim] tau_betas;

  array[n_par, n_roi] real sigma_z;
  array[n_par, n_roi] real delta_z;
  array[n_par, n_roi] real rho_time_z;
  array[n_par, n_roi] vector[n_beta] betas_z;
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;

  vector[0] phi; // parameters shared between all data shards


  array[n_roi] real <lower =0>  mu_sigma;
  array[n_roi] real <lower =0>  mu_delta;
  array[n_roi] real <lower =0>  mu_rho_time;

  for (r in 1:n_roi) {
    mu_sigma[r] = mu_sigma_raw[r] *.5;
    mu_delta[r] = inv_logit(mu_delta_raw[r] * 1.75);
    mu_rho_time[r] = inv_logit(mu_rho_time_raw[r] * 1.75);
  }


  array[n_par, n_roi] real <lower =0>  sigma;
  array[n_par, n_roi] real <lower =0>  delta;
  array[n_par, n_roi] real <lower =0>  rho_time;
    array[n_par, n_roi] vector[n_beta] betas;


  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma[p, r] = mu_sigma[r] + tau_sigma_raw[r] * sigma_z[p, r];
      delta[p, r] = inv_logit((mu_delta_raw[r] * 1.75) + tau_delta_raw[r] * delta_z[p, r]);
      rho_time[p, r] = inv_logit((mu_rho_time_raw[r] * 1.75) + tau_rho_time_raw[r] * rho_time_z[p, r]);
      for (b in 1:n_beta){
        betas[p,r][b] = mu_betas[r][b] + tau_betas[r,b] * betas_z[p,r][b];
      }
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

  // save the log_lik per shard for elpd_loo later
  // this is bad slow way of doing this, makes it a parameter that gets a lot of tracking
  // vector[n_par] log_lik_shards = map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int);

}

model {
  // Priors
  for (r in 1:n_roi) {
    mu_sigma_raw[r] ~ std_normal();
    mu_delta_raw[r] ~ std_normal();
    mu_rho_time_raw[r] ~ std_normal();
    tau_sigma_raw[r] ~ std_normal();
    tau_delta_raw[r] ~ std_normal();
    tau_rho_time_raw[r] ~ std_normal();
    for (b in 1:n_beta){
      mu_betas[r][b] ~ std_normal();
      tau_betas[r][b] ~ std_normal();
    }
  }

  // Multilevel priors (or z-score priors for the multilevel)
  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma_z[p, r] ~ std_normal();
      delta_z[p, r] ~ std_normal();
      rho_time_z[p, r] ~ std_normal();
      betas_z[p,r] ~ std_normal();
      
      // for (b in 1:n_beta){
      // sigma_raw[p, r] ~ normal(mu_sigma_raw[r], tau_sigma_raw[r]);
      // delta_raw[p, r] ~ normal(mu_delta_raw[r], tau_delta_raw[r]);
      // rho_time_raw[p, r] ~ normal(mu_rho_time_raw[r], tau_rho_time_raw[r]);
      
      // for (b in 1:n_beta){
      //   betas[p,r][b] ~ normal(mu_betas[r][b], tau_betas[r, b]);
      // only regression coefficents for stim not movement get multilevel prior
      // for (b in 1:n_beta){
      //   if (b < n_beta_stim) {
      //     betas[p, r][b] ~ std_normal();
      //   } else {
      //     betas[p, r][b] ~ std_normal();
      //   }
    }
  }


target += map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int);

// target += sum(map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int));

}

generated quantities {
  // array[n_observations] real mu_pred = mu[indices_observed];
  vector[n_par] log_lik;
  
  log_lik = map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int);


  // for (i in 1:n_observations) {
  //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
  // }
}
