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
  int n_censor = data_shard_int[2];
  int n_beta = data_shard_int[3];
  int p = data_shard_int[4];
  int r = data_shard_int[5];
  array[n_bold] int uncensored_indices_padded = 
    data_shard_int[6:(n_bold+5)];

  array[n_bold - n_censor] int uncensored_indices = 
    uncensored_indices_padded[1:(n_bold - n_censor)];

  // upack data_shard_real
  vector[n_bold] bold = to_vector(data_shard_real[1:n_bold]);
  matrix[n_bold, n_beta] design_matrix = 
    to_matrix(data_shard_real[(1 + n_bold):
                              (n_bold + n_bold * n_beta)],
              n_bold, n_beta);

  matrix[n_bold, n_bold] rho_power_matrix = 
    to_matrix(data_shard_real[(1 + n_bold + n_bold * n_beta):
                             (n_bold + n_bold * n_beta + n_bold * n_bold)], 
              n_bold, n_bold);


  // unpack parameters


  real sigma_z = theta[1];
  real delta_z = theta[2];
  real rho_time_z = theta[3];
  vector[n_beta] betas_z = theta[4:(3 + n_beta)];

  // log-likelihoods of priors


  // Create mu predictions and covariance matrices per roi, then calc likelihood
  matrix[n_bold - n_censor, n_beta] design_matrix_censor = design_matrix[uncensored_indices,];
  vector[n_bold - n_censor] Mu = (design_matrix_censor * betas_z);
  matrix[n_bold, n_bold] delta_mat = 
  // 1 on diagonal, delta off
    add_diag(
      add_diag(
        rep_matrix(delta_z, n_bold, n_bold), 
        - delta_z),
      1);

  matrix [n_bold, n_bold] rho_mat = rho_time_z .^ rho_power_matrix;
  matrix[n_bold, n_bold] Cov = pow(sigma_z, 2) .*  delta_mat .* rho_mat;
  matrix[n_bold - n_censor, n_bold - n_censor] Cov_censor = Cov[uncensored_indices, uncensored_indices];
  matrix[n_bold - n_censor, n_bold - n_censor] L_Cov_censor = cholesky_decompose(Cov_censor);

  vector[n_bold - n_censor] bold_censored = bold[uncensored_indices];

  real lp = multi_normal_cholesky_lpdf(bold_censored | Mu, L_Cov_censor);


  return [lp]';
}  

  vector multi_normal_elementwise_log_lik(vector phi, // reused for every shared
                                          vector theta, // shard specific parameters
                                          array[] real data_shard_real, // data (real numbers) used for this shard
                                          array[] int data_shard_int) { // data (integers) used for this shard


  // unpack data_shard_int
  int n_bold = data_shard_int[1];
  int n_censor = data_shard_int[2];
  int n_beta = data_shard_int[3];
  int p = data_shard_int[4];
  int r = data_shard_int[5];
  array[n_bold] int uncensored_indices_padded = 
    data_shard_int[6:(n_bold+5)];

  array[n_bold - n_censor] int uncensored_indices = 
    uncensored_indices_padded[1:(n_bold - n_censor)];

  // upack data_shard_real
  vector[n_bold] bold = to_vector(data_shard_real[1:n_bold]);
  matrix[n_bold, n_beta] design_matrix = 
    to_matrix(data_shard_real[(1 + n_bold):
                              (n_bold + n_bold * n_beta)],
              n_bold, n_beta);

  matrix[n_bold, n_bold] rho_power_matrix = 
    to_matrix(data_shard_real[(1 + n_bold + n_bold * n_beta):
                             (n_bold + n_bold * n_beta + n_bold * n_bold)], 
              n_bold, n_bold);


  // unpack parameters
  real sigma = theta[1];
  real delta = theta[2];
  real rho_time = theta[3];
  vector[n_beta] betas = theta[4:(3 + n_beta)];
  // vector[n_bold] bold_z_padded = theta[(4 + n_beta):
  //                                      (3 + n_beta + n_bold)];
  // vector[n_bold - n_censor] bold_z = bold_z_padded[1:(n_bold-n_censor)];

  // Create mu predictions and covariance matrices per roi, then calc likelihood
  matrix[n_bold - n_censor, n_beta] design_matrix_censor = design_matrix[uncensored_indices,];
  vector[n_bold - n_censor] Mu = (design_matrix_censor * betas);
  matrix[n_bold, n_bold] delta_mat = 
  // 1 on diagonal, delta off
    add_diag(
      add_diag(
        rep_matrix(delta, n_bold, n_bold), 
        - delta),
      1);
  matrix [n_bold, n_bold] rho_mat = rho_time .^ rho_power_matrix;
  matrix[n_bold, n_bold] Cov = pow(sigma, 2) .*  delta_mat .* rho_mat;
  matrix[n_bold - n_censor, n_bold - n_censor] Cov_censor = Cov[uncensored_indices, uncensored_indices];
  matrix[n_bold - n_censor, n_bold - n_censor] L_Cov_censor = cholesky_decompose(Cov_censor);

  vector[n_bold - n_censor] bold_censored = bold[uncensored_indices];

  vector[n_bold - n_censor] z = mdivide_left_tri_low(L_Cov_censor, bold_censored - Mu);

  vector[n_bold - n_censor] comps;

  for (k in 1:(n_bold - n_censor)) {
      // normal_lpdf(z[k] | 0, 1) == -0.5*log(2*pi()) - 0.5*z[k]^2
      comps[k] = normal_lpdf(z[k] | 0, 1) - log(L_Cov_censor[k,k]);
  }
  return comps;
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
  array[n_par, n_roi, n_bold] real bold;
  array[n_par, n_bold] int usable_bold_indices_one_is_true;
  array[n_par] int n_censor;
  array[n_par, n_bold, n_beta] real design_array;
}

transformed data {
  int n_shards = n_par * n_roi;

  // Create a matrix where each row is the index i.
  // v * rep_row_vector(1, n_bold) produces an n_bold x n_bold matrix
  // where row i consists of the value idx[i] repeated n_bold times.
  vector[n_bold] idx = linspaced_vector(n_bold, 1, n_bold);
  matrix[n_bold, n_bold] row_idx = idx * rep_row_vector(1, n_bold);
  
  // Similarly, rep_column_vector(1, n_bold) * idx'
  // produces an n_bold x n_bold matrix where column j consists of the value idx[j] repeated.
  matrix[n_bold, n_bold] col_idx = rep_vector(1, n_bold) * idx';
  
  // Now compute the matrix of absolute differences:
  matrix[n_bold, n_bold] rho_power_matrix = abs(row_idx - col_idx);
  


  // map_rect needs inputs to be same size, so zeros are padded at end per participant
  array[n_par, n_bold] int uncensored_indices_padded;
  int total_usable_TRs = 0;
  for (p in 1:n_par) {
    int keep_index = 1;
    for (i in 1:n_bold) {
      if (usable_bold_indices_one_is_true[p,i] == 1) {
        uncensored_indices_padded[p, keep_index] = i;
        keep_index = keep_index + 1;
        total_usable_TRs += 1;
      }
    }

    if(n_censor[p] > 0) {
      array [n_censor[p]] int pad = rep_array(0, n_censor[p]);
      uncensored_indices_padded[p, keep_index:n_bold] = pad;
    }
  }

  int total_usable_TRs_all = total_usable_TRs * n_roi;


  // For map_rect: 
  // create data_shard_int for each roi by participant
  // n_bold; 1
  // n_censor; 2
  // n_beta; 3
  // p (participant); 4 
  // r (roi); 5
  // uncensored_indices_padded; 5 + n_bold
    array[n_shards, 5+n_bold] int data_shard_int;
//   array[n_shards, 3+n_bold] int data_shard_int;
  int shard_index_1stD = 1;
  for(p in 1:n_par) {
    for(r in 1:n_roi) {
      // pack scalars
      data_shard_int[shard_index_1stD, 1] = n_bold;
      data_shard_int[shard_index_1stD, 2] = n_censor[p];
      data_shard_int[shard_index_1stD, 3] = n_beta;
      data_shard_int[shard_index_1stD, 4] = p;
      data_shard_int[shard_index_1stD, 5] = r;
      // pack uncensored indices padded vector
      data_shard_int[shard_index_1stD, 
                     6:(5 + n_bold)] = 
                     uncensored_indices_padded[p,];


      shard_index_1stD = shard_index_1stD + 1;
    }
  }

  // For map_rect: 
  // create data_shard_real for each roi by participant
  // bold; n_bold
  // design_matrix; n_bold + n_bold * n_beta
  // rho_power_matrix; n_bold + n_bold * n_beta + n_bold * n_bold
  array[n_shards, n_bold + n_bold * n_beta + n_bold * n_bold] real data_shard_real;
  shard_index_1stD = 1;
  for (p in 1:n_par) {
    int bold_start = 1;
    int bold_stop = n_bold;
    for (r in 1:n_roi) {
      // pack bold
      int data_shard_2ndD_start = 1;
      int data_shard_2ndD_stop = n_bold;
      data_shard_real[shard_index_1stD, 
                      data_shard_2ndD_start:data_shard_2ndD_stop] = 
                      bold[p, r, 1:n_bold];

      // pack design matrix from column-major order 
      data_shard_2ndD_start = data_shard_2ndD_start + n_bold;
      data_shard_2ndD_stop = data_shard_2ndD_stop + n_bold * n_beta;
      data_shard_real[shard_index_1stD, 
                      data_shard_2ndD_start:data_shard_2ndD_stop] = 
                      to_array_1d(to_matrix(design_array[p,,]));

      // pack rho_power_matrix for Toeplitz
      data_shard_2ndD_start = data_shard_2ndD_start + n_bold * n_beta;
      data_shard_2ndD_stop = data_shard_2ndD_stop + n_bold * n_bold;
      data_shard_real[shard_index_1stD, 
                     data_shard_2ndD_start:data_shard_2ndD_stop] = 
                     to_array_1d(rho_power_matrix);
      
      shard_index_1stD = shard_index_1stD + 1;
    }
  }
}

parameters {

  // parameters shared between all data shards
  // total number of parameters in phi
  // raw versions get transformed later in map_rect
  // mu_sigma_raw; 1  * n_roi
  // tau_sigma_raw; 2  * n_roi
  // mu_delta_raw; 3  * n_roi
  // tau_delta_raw; 4  * n_roi
  // mu_rho_time_raw; 5  * n_roi
  // tau_rho_time; 6 * n_roi
  // mu_betas; 6 * n_roi + 1 * n_beta * n_roi 
  // tau_betas_raw; 6 * n_roi + 2 * n_beta * n_roi 
  // vector[6 * n_roi + 2 * n_beta * n_roi] phi; 
  // jk phi is empty again, but maybe there should be another function for duplicating this for each parameter
  vector[0] phi; // parameters shared between all data shards; sometimes empty
  

  // for theta
  array[n_par] vector[n_roi] sigma_z;
  array[n_par] vector[n_roi] delta_z_raw;
  array[n_par] vector[n_roi] rho_time_z_raw;
  array[n_par, n_beta] vector[n_roi] betas_z;
  
  vector <lower=0> [n_roi] mu_sigma_raw; // this still needs to be half normal
  // array[n_roi] real <lower=0> mu_sigma_raw; // this still needs to be half normal
  vector <lower=0> [n_roi] tau_sigma_raw;
  vector[n_roi] mu_delta_raw;
  vector <lower=0> [n_roi]  tau_delta_raw;
  vector[n_roi] mu_rho_time_raw;
  vector <lower=0> [n_roi] tau_rho_time_raw;
  array[n_beta] vector[n_roi] mu_betas;
  array[n_beta] vector <lower=0> [n_roi] tau_betas;

  // old used to be packed into theta, now done in map rect
  // array[n_par, n_roi] real sigma_z;
  // array[n_par, n_roi] real delta_z;
  // array[n_par, n_roi] real rho_time_z;
  // array[n_par, n_roi] vector[n_beta] betas_z;
  
  // these are the pairwise correlations between the ROIs
  // they get used to make covariance matrices later on
  // rho_z because it goes through tanh to become a real correlation
  // now I don't think it makes sense to find these since there is but one group mu and tau per roi, so finding the correlation here doesn't make sense
//   vector[(n_roi*(n_roi-1))%/%2] rho_z_mu_sigma_raw;
//   vector[(n_roi*(n_roi-1))%/%2] rho_z_tau_sigma_raw;
//   vector[(n_roi*(n_roi-1))%/%2] rho_z_mu_delta_raw;
//   vector[(n_roi*(n_roi-1))%/%2] rho_z_tau_delta_raw;
//   vector[(n_roi*(n_roi-1))%/%2] rho_z_mu_rho_time_raw;
//   vector[(n_roi*(n_roi-1))%/%2] rho_z_tau_rho_time;
//   array[n_beta] vector[(n_roi*(n_roi-1))%/%2] rho_z_mu_betas;
//   array[n_beta] vector[(n_roi*(n_roi-1))%/%2] rho_z_tau_betas_raw;

  vector[(n_roi*(n_roi-1))%/%2] rho_z_sigma_z_raw;
  vector[(n_roi*(n_roi-1))%/%2] rho_z_delta_z_raw;
  vector[(n_roi*(n_roi-1))%/%2] rho_z_rho_time_z_raw; //rho_z_rho is not a typo, it is how rho over time is correlated between ROIs
  array[n_beta] vector[(n_roi*(n_roi-1))%/%2] rho_z_betas_z;

}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;






  array[n_par, n_roi] real <lower=0>  sigma;
  array[n_par, n_roi] real <lower=0>  delta;
  array[n_par, n_roi] real <lower=0>  rho_time;
  array[n_par, n_roi] vector[n_beta] betas;

  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma[p, r] = mu_sigma_raw[r] * 0.5 + tau_sigma_raw[r] * sigma_z[p][r];
      delta[p, r] = inv_logit((mu_delta_raw[r] * 1.75) + tau_delta_raw[r] * delta_z_raw[p][r]);
      rho_time[p, r] = inv_logit((mu_rho_time_raw[r] * 1.75) + tau_rho_time_raw[r] * rho_time_z_raw[p][r]);
      for (b in 1:n_beta){
        //notice ordering changes from betas_z to betas
        betas[p,r][b] = mu_betas[b][r] + tau_betas[b][r] * betas_z[p,b][r];
      }
    }
  }

  
  // local shard parameters theta for participant by roi
  // sigma_z; 1 
  // delta_z; 2
  // rho_time_z; 3
  // betas_z; 3 + n_beta
  array[n_shards] vector[(3 + n_beta)] theta;


  for (i in 1:1) {
    // has to be in loop to make this integer
    int Pshard_index_1stD = 1;
    for (p in 1:n_par) {
      for (r in 1:n_roi) {
        theta[Pshard_index_1stD][1] = sigma[p,r];
        theta[Pshard_index_1stD][2] = delta[p,r];
        theta[Pshard_index_1stD][3] = rho_time[p,r];
        theta[Pshard_index_1stD][4:(3 + n_beta)] = to_vector(betas[p,r]);
      
        // int current_theta_index = 4 + n_beta;
        // int useable_bold_index = 1;
        // for (b in 1:n_bold) {
        //   if (usable_bold_indices_one_is_true[p,b] == 1) {
        //     theta[Pshard_index_1stD][current_theta_index] = bold_z[useable_bold_index];
        //     current_theta_index += 1;
        //     useable_bold_index += 1;
        //   }
        // }
        // if (n_censor[p] > 0) {
        //   theta[Pshard_index_1stD][current_theta_index:current_theta_index + n_censor[p] - 1] = 
        //     rep_vector(0, n_censor[p]);
        // }
        Pshard_index_1stD = Pshard_index_1stD + 1;
      }
    }
  }


//   theta[1] = to_vector(sigma_z);

//   int theta_start_idx = 1;
//   int theta_stop_idx = n_roi;
  
  

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
    mu_betas[r] ~ std_normal();
    tau_betas[r] ~ std_normal();
  }

  // now found through mvn
  // Multilevel priors (or z-score priors for the multilevel)
  // for (p in 1:n_par) {
  //   for (r in 1:n_roi) {
  //     sigma_z[p, r] ~ std_normal();
  //     delta_z[p, r] ~ std_normal();
  //     rho_time_z[p, r] ~ std_normal();
  //     betas_z[p,r] ~ std_normal();
  //   }
  // }

  // uninformative priors for pairwise-ROI correlations for each parameter
//   rho_z_mu_sigma_raw ~ std_normal();
//   rho_z_tau_sigma_raw ~ std_normal();
//   rho_z_mu_delta_raw ~ std_normal();
//   rho_z_tau_delta_raw ~ std_normal();
//   rho_z_mu_rho_time_raw ~ std_normal();
//   rho_z_tau_rho_time ~ std_normal();

  
  rho_z_sigma_z_raw ~ std_normal();
  rho_z_delta_z_raw ~ std_normal();
  rho_z_rho_time_z_raw ~ std_normal();

  for (b in 1:n_beta) {
    // rho_z_mu_betas[b] ~ std_normal();
    // rho_z_tau_betas_raw[b] ~ std_normal();

    rho_z_betas_z[b] ~ std_normal();
  }    

  // Use rhos for correlation matrices
//   matrix[n_roi, n_roi] Corr_mu_sigma_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_tau_sigma_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_mu_delta_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_tau_delta_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_mu_rho_time_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_tau_rho_time_raw = rep_matrix(1, n_roi, n_roi);
  
  matrix[n_roi, n_roi] Corr_sigma_z_raw = rep_matrix(1, n_roi, n_roi);
  matrix[n_roi, n_roi] Corr_delta_z_raw = rep_matrix(1, n_roi, n_roi);
  matrix[n_roi, n_roi] Corr_rho_time_z_raw = rep_matrix(1, n_roi, n_roi);

//   array[n_beta] matrix[n_roi, n_roi] Corr_mu_betas = rep_array(rep_matrix(1, n_roi, n_roi), n_beta);
//   array[n_beta] matrix[n_roi, n_roi] Corr_tau_betas_raw = rep_array(rep_matrix(1, n_roi, n_roi), n_beta);
  
  array[n_beta] matrix[n_roi, n_roi] Corr_betas_z = rep_array(rep_matrix(1, n_roi, n_roi), n_beta);
  
  int rho_idx = 1;
  for (i in 1:(n_roi-1)) {
    for (j in (i + 1):n_roi) {
      // group
    //   Corr_mu_sigma_raw[i,j] = tanh(rho_z_mu_sigma_raw[rho_idx]);
    //   Corr_mu_sigma_raw[j,i] = tanh(rho_z_mu_sigma_raw[rho_idx]);

    //   Corr_tau_sigma_raw[i,j] = tanh(rho_z_tau_sigma_raw[rho_idx]);
    //   Corr_tau_sigma_raw[j,i] = tanh(rho_z_tau_sigma_raw[rho_idx]);
      
    //   Corr_mu_delta_raw[i,j] = tanh(rho_z_mu_delta_raw[rho_idx]);
    //   Corr_mu_delta_raw[j,i] = tanh(rho_z_mu_delta_raw[rho_idx]);

    //   Corr_tau_delta_raw[i,j] = tanh(rho_z_tau_delta_raw[rho_idx]);
    //   Corr_tau_delta_raw[j,i] = tanh(rho_z_tau_delta_raw[rho_idx]);

    //   Corr_mu_rho_time_raw[i,j] = tanh(rho_z_mu_rho_time_raw[rho_idx]);
    //   Corr_mu_rho_time_raw[j,i] = tanh(rho_z_mu_rho_time_raw[rho_idx]);
      
    //   Corr_tau_rho_time_raw[i,j] = tanh(rho_z_tau_delta_raw[rho_idx]);
    //   Corr_tau_rho_time_raw[j,i] = tanh(rho_z_tau_delta_raw[rho_idx]);
      
      // local
      Corr_sigma_z_raw[i,j] = tanh(rho_z_sigma_z_raw[rho_idx]);
      Corr_sigma_z_raw[j,i] = tanh(rho_z_sigma_z_raw[rho_idx]);

      Corr_delta_z_raw[i,j] = tanh(rho_z_delta_z_raw[rho_idx]);
      Corr_delta_z_raw[j,i] = tanh(rho_z_delta_z_raw[rho_idx]);

      Corr_rho_time_z_raw[i,j] = tanh(rho_z_rho_time_z_raw[rho_idx]);
      Corr_rho_time_z_raw[j,i] = tanh(rho_z_rho_time_z_raw[rho_idx]);
      for (b in 1:n_beta) {
        // group
        // Corr_mu_betas[b][i,j] = tanh(rho_z_mu_betas[b][rho_idx]);
        // Corr_mu_betas[b][j,i] = tanh(rho_z_mu_betas[b][rho_idx]);
        
        // Corr_tau_betas_raw[b][i,j] = tanh(rho_z_tau_betas_raw[b][rho_idx]);
        // Corr_tau_betas_raw[b][j,i] = tanh(rho_z_tau_betas_raw[b][rho_idx]);
        //local
        Corr_betas_z[b][i,j] = tanh(rho_z_betas_z[b][rho_idx]);
        Corr_betas_z[b][j,i] = tanh(rho_z_betas_z[b][rho_idx]);
      }
      rho_idx += 1;
    }
  }

//   matrix[n_roi, n_roi] L_Corr_mu_sigma = cholesky_decompose(Corr_mu_sigma_raw);
//   matrix[n_roi, n_roi] L_Corr_tau_sigma = cholesky_decompose(Corr_tau_sigma_raw);
//   matrix[n_roi, n_roi] L_Corr_mu_delta_raw = cholesky_decompose(Corr_mu_delta_raw);
//   matrix[n_roi, n_roi] L_Corr_tau_delta_raw = cholesky_decompose(Corr_tau_delta_raw);
//   matrix[n_roi, n_roi] L_Corr_mu_rho_time_raw = cholesky_decompose(Corr_mu_rho_time_raw);
//   matrix[n_roi, n_roi] L_Corr_tau_rho_time_raw = cholesky_decompose(Corr_tau_rho_time_raw);
  
  matrix[n_roi, n_roi] L_Corr_sigma_z = cholesky_decompose(Corr_sigma_z_raw);
  matrix[n_roi, n_roi] L_Corr_delta_z_raw = cholesky_decompose(Corr_delta_z_raw);
  matrix[n_roi, n_roi] L_Corr_rho_time_z_raw = cholesky_decompose(Corr_rho_time_z_raw);
  
//   array[n_beta] matrix[n_roi, n_roi] L_Corr_mu_betas;
//   array[n_beta] matrix[n_roi, n_roi] L_Corr_tau_betas_raw;
  
  array[n_beta] matrix[n_roi, n_roi] L_Corr_betas_z;

  for (b in 1:n_beta) {
    // L_Corr_mu_betas[b] = cholesky_decompose(Corr_mu_betas[b]);
    // L_Corr_tau_betas_raw[b] = cholesky_decompose(Corr_tau_betas_raw[b]);
    
    L_Corr_betas_z[b] = cholesky_decompose(Corr_betas_z[b]);
  }    

  // The ROI correlational structure is used for better priors
  // shared group parameters
//   mu_sigma_raw ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_mu_sigma);
//   tau_sigma_raw ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_tau_sigma);
//   mu_delta_raw ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_mu_delta_raw);
//   tau_delta_raw ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_tau_delta_raw);
//   mu_rho_time_raw ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_mu_rho_time_raw);
//   tau_rho_time_raw ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_tau_rho_time_raw);
//   for (b in 1:n_beta) {
//     mu_betas[b] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_mu_betas[b]);
//     tau_betas[b] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_tau_betas_raw[b]);
//   }

  // Theta local participant and ROI parameters
  // The noncentered parameterization is used, thus z-scores get sampled, but this is still multilevel

  for(p in 1:n_par){
    sigma_z[p] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_sigma_z);
    delta_z_raw[p] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_delta_z_raw);
    rho_time_z_raw[p] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_rho_time_z_raw);
    for (b in 1:n_beta) {
      betas_z[p,b] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_betas_z[b]);
    }
  }




  target += map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int);

}

generated quantities {
  // array[n_observations] real mu_pred = mu[indices_observed];

  // these need to be done from phi / theta now
  // array[n_roi] real <lower =0>  mu_sigma;
  // array[n_roi] real <lower =0>  mu_delta;
  // array[n_roi] real <lower =0>  mu_rho_time;

  // for (r in 1:n_roi) {
  //   mu_sigma[r] = mu_sigma_raw[r] *.5;
  //   mu_delta[r] = inv_logit(mu_delta_raw[r] * 1.75);
  //   mu_rho_time[r] = inv_logit(mu_rho_time_raw[r] * 1.75);
  // }

  vector[total_usable_TRs_all] log_lik;
  int shard_index = 1;
  int log_lik_start = 1;
  int log_lik_stop = n_bold - n_censor[1];
  for(s in 1:n_shards) {
      if (s > 1) {
        int current_par = data_shard_int[s  ][4];
        int last_par    = data_shard_int[s-1][4];
        int current_n_censor = n_censor[current_par];
        int last_n_censor = n_censor[last_par];
        log_lik_start += n_bold - last_n_censor; 
        log_lik_stop += n_bold - current_n_censor;
      }
      log_lik[log_lik_start:log_lik_stop] = 
        multi_normal_elementwise_log_lik(phi, 
                                         theta[s], 
                                         data_shard_real[s], 
                                         data_shard_int[s]);
    }

  // log_lik = map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int);


  // for (i in 1:n_observations) {
  //   log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], L_Cov[participant[i]]);
  // }
}
