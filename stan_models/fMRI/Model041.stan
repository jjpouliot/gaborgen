functions {

  // --- Helper: solve for MA(1) parameter theta from (rho, delta) ---
  real arma11_theta_from_rho_delta(real rho, real delta) {
    // Handles edge cases; returns invertible MA root (|theta| < 1) when possible.
    real theta;
    if (delta > 0.9999) {
      // delta -> 1: pure AR(1) case
      theta = 0;
    } else {
      real denom = rho * (1 - delta);
      if (denom < 1e-9) denom = 1e-9;
      real B = (1 + square(rho) * (1 - 2 * delta)) / denom;
      real disc = B * B - 4;
      if (disc < 0) disc = 0;
      real t1 = (-B + sqrt(disc)) / 2;
      real t2 = (-B - sqrt(disc)) / 2;

      // choose root with |theta| < 1 using squares to avoid fabs()
      if (square(t1) < 1) {
        theta = t1;
      } else {
        theta = t2;
      }
      // final guard if numerical issues
      if (square(theta) >= 1) {
        if (theta > 0) theta = 0.999;
        else theta = -0.999;
      }
    }
    return theta;
  }

  // --- Helper: map (rho, delta, sigma) -> innovation sd sigma_e ---
  real arma11_sigma_e_from_rho_delta_sigma(real rho, real delta, real sigma) {
    real theta = arma11_theta_from_rho_delta(rho, delta);
    real denom = 1 + square(theta) + 2 * rho * theta;
    if (denom < 1e-12) denom = 1e-12;
    real num = 1 - square(rho);
    if (num < 1e-12) num = 1e-12;
    return sigma * sqrt(num / denom);
  }

  // --- Core: exact ARMA(1,1) innovations log-likelihood (handles censoring) ---
  real arma11_innovations_loglik(vector y, vector mu,
                             real rho, real delta, real sigma,
                             array[] int is_obs) {
    int N = num_elements(y);

    // ARMA(1,1) parameters that reproduce your Toeplitz Sigma
    real Phi     = rho;
    real theta   = arma11_theta_from_rho_delta(rho, delta);
    real sigma_e = arma11_sigma_e_from_rho_delta_sigma(rho, delta, sigma);

    // State is alpha_t = [ y_t - mu_t , epsilon_t ]'
    // Transition T = [[Phi, theta],[0,0]]
    // Process cov Q = sigma_e^2 * [[1,1],[1,1]]
    // Observation: y_t - mu_t = [1, 0] * alpha_t, measurement noise = 0

    // State mean init
    real a1 = 0;
    real a2 = 0;

    // Stationary covariance init P solves P = T P T' + Q (closed form for this T, Q)
    real s2 = square(sigma_e);
    real p11 = s2 * (1 + 2 * Phi * theta + square(theta)) / (1 - square(Phi));
    real p12 = s2;
    real p22 = s2;

    real lp = 0;
    real eps = 1e-12;

    for (t in 1:N) {
      // Predict step
      real a1_pred = Phi * a1 + theta * a2;
      real a2_pred = 0;

      // P_pred = T P T' + Q  (worked out elementwise for this T, Q)
      real p11_pred = square(Phi) * p11 + 2 * Phi * theta * p12 + square(theta) * p22 + s2;
      real p12_pred = s2;     // equals p21_pred
      real p22_pred = s2;

      if (is_obs[t] == 1) {
        // Innovation and its variance
        real v = (y[t] - mu[t]) - a1_pred;
        real F = p11_pred + eps;  // H=0

        lp += -0.5 * (log(2 * pi()) + log(F) + v * v / F);

        // Kalman gain (2x1)
        real k1 = p11_pred / F;
        real k2 = p12_pred / F;

        // Update state mean
        a1 = a1_pred + k1 * v;
        a2 = a2_pred + k2 * v;

        // Joseph form update for covariance to preserve symmetry, H=0, Z=[1,0]
        // P_new = (I - KZ) P_pred (I - KZ)'
        real one_minus_k1 = 1 - k1;
        real p11_new = square(one_minus_k1) * p11_pred;
        real p12_new = (one_minus_k1) * (p12_pred - k2 * p11_pred);
        real p22_new = square(k2) * p11_pred - 2 * k2 * p12_pred + p22_pred;

        p11 = p11_new;
        p12 = p12_new;
        p22 = p22_new;
      } else {
        // Censored/missing: prediction only
        a1 = a1_pred;
        a2 = a2_pred;
        p11 = p11_pred;
        p12 = p12_pred;
        p22 = p22_pred;
      }
    }
    return lp;
  }

  // ---------- MAP_RECT WORKER: now calls innovations lp (no Cholesky) ----------
  vector participant_ll(vector phi, vector theta,
                        array[] real data_shard_real,
                        array[] int data_shard_int) {

    // unpack data_shard_int
    int n_bold  = data_shard_int[1];
    int n_censor = data_shard_int[2];
    int n_beta  = data_shard_int[3];
    int p       = data_shard_int[4];
    int r       = data_shard_int[5];
    array[n_bold] int uncensored_indices_padded =
      data_shard_int[6:(n_bold + 5)];
    array[n_bold - n_censor] int uncensored_indices =
      uncensored_indices_padded[1:(n_bold - n_censor)];

    // upack data_shard_real
    vector[n_bold] bold = to_vector(data_shard_real[1:n_bold]);
    matrix[n_bold, n_beta] design_matrix =
      to_matrix(data_shard_real[(1 + n_bold):(n_bold + n_bold * n_beta)],
                n_bold, n_beta);

    // unpack parameters
    real sigma     = theta[1];
    real delta     = theta[2];
    real rho_time  = theta[3];
    vector[n_beta] betas = theta[4:(3 + n_beta)];

    // Build full-length mean
    vector[n_bold] Mu_full = design_matrix * betas;

    // Build observation mask (1=keep, 0=censored)
    array[n_bold] int is_obs;
    {
    //   int i;
      for (i in 1:n_bold) is_obs[i] = 0;
      for (i in 1:(n_bold - n_censor)) is_obs[uncensored_indices[i]] = 1;
    }

    // Exact MVN log-lik via ARMA(1,1) innovations (O(N), no Cholesky)
    real lp = arma11_innovations_loglik(bold, Mu_full, rho_time, delta, sigma, is_obs);

    return [lp]';
  }

  // ---------- Per-observation log-lik components (for loo etc.) ----------
  vector multi_normal_elementwise_log_lik(vector phi, vector theta,
                                          array[] real data_shard_real,
                                          array[] int data_shard_int) {

    // unpack data_shard_int
    int n_bold  = data_shard_int[1];
    int n_censor = data_shard_int[2];
    int n_beta  = data_shard_int[3];
    int p       = data_shard_int[4];
    int r       = data_shard_int[5];
    array[n_bold] int uncensored_indices_padded =
      data_shard_int[6:(n_bold + 5)];
    array[n_bold - n_censor] int uncensored_indices =
      uncensored_indices_padded[1:(n_bold - n_censor)];

    // upack data_shard_real
    vector[n_bold] bold = to_vector(data_shard_real[1:n_bold]);
    matrix[n_bold, n_beta] design_matrix =
      to_matrix(data_shard_real[(1 + n_bold):(n_bold + n_bold * n_beta)],
                n_bold, n_beta);

    // unpack parameters
    real sigma     = theta[1];
    real delta     = theta[2];
    real rho_time  = theta[3];
    vector[n_beta] betas = theta[4:(3 + n_beta)];

    // Full-length mean and obs mask
    vector[n_bold] Mu_full = design_matrix * betas;
    array[n_bold] int is_obs;
    {
    //   int i;
      for (i in 1:n_bold) is_obs[i] = 0;
      for (i in 1:(n_bold - n_censor)) is_obs[uncensored_indices[i]] = 1;
    }

    // Re-run the innovations loop but store per-observation contributions
    int N_keep = n_bold - n_censor;
    vector[N_keep] comps;

    // ARMA params
    real Phi     = rho_time;
    real theta_ma = arma11_theta_from_rho_delta(rho_time, delta);
    real sigma_e = arma11_sigma_e_from_rho_delta_sigma(rho_time, delta, sigma);
    real s2 = square(sigma_e);

    // State init
    real a1 = 0;
    real a2 = 0;
    real p11 = s2 * (1 + 2 * Phi * theta_ma + square(theta_ma)) / (1 - square(Phi));
    real p12 = s2;
    real p22 = s2;

    int k = 1; // index into comps
    real eps = 1e-12;

    for (t in 1:n_bold) {
      // Predict
      real a1_pred = Phi * a1 + theta_ma * a2;
      real a2_pred = 0;

      real p11_pred = square(Phi) * p11 + 2 * Phi * theta_ma * p12 + square(theta_ma) * p22 + s2;
      real p12_pred = s2;
      real p22_pred = s2;

      if (is_obs[t] == 1) {
        real v = (bold[t] - Mu_full[t]) - a1_pred;
        real F = p11_pred + eps;

        comps[k] = -0.5 * (log(2 * pi()) + log(F) + v * v / F);
        k = k + 1;

        // Gain
        real k1 = p11_pred / F;
        real k2 = p12_pred / F;

        // Update state
        a1 = a1_pred + k1 * v;
        a2 = a2_pred + k2 * v;

        // Joseph covariance update
        real one_minus_k1 = 1 - k1;
        real p11_new = square(one_minus_k1) * p11_pred;
        real p12_new = (one_minus_k1) * (p12_pred - k2 * p11_pred);
        real p22_new = square(k2) * p11_pred - 2 * k2 * p12_pred + p22_pred;

        p11 = p11_new;
        p12 = p12_new;
        p22 = p22_new;
      } else {
        // Missing
        a1 = a1_pred;
        a2 = a2_pred;
        p11 = p11_pred;
        p12 = p12_pred;
        p22 = p22_pred;
      }
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
  // array[n_par, n_bold*n_roi] real bold;
   array[n_par, n_roi, n_bold] real bold;
  array[n_par, n_bold] int usable_bold_indices_one_is_true;
  array[n_par] int n_censor;
  array[n_par, n_bold, n_beta] real design_array;
}

transformed data {
  int n_beta_motion = 6;
  int n_beta_stim = n_beta - n_beta_motion;
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


  // old from 18
  // For map_rect: 
  // create data_shard_int for each roi by participant
  // n_bold 1
  // n_roi 2
  // n_censor 3
  // n_beta 4
  // censor 4 + n_bold
  // rho_power_matrix 4 + n_bold + n_bold * n_beta
  // array[n_par * n_roi,
  //       4 + n_bold + n_bold * n_beta] int data_shard_int;
  // int shard_index_1stD = 1;
  // for(p in 1:n_par) {
  //   for(r in 1:n_roi) {
  //     // pack scalars
  //     data_shard_int[shard_index_1stD, 1] = n_bold;
  //     data_shard_int[shard_index_1stD, 2] = n_roi;
  //     data_shard_int[shard_index_1stD, 3] = n_censor[p];
  //     data_shard_int[shard_index_1stD, 4] = n_beta;
  //     // pack uncensored indices padded vector
  //     data_shard_int[shard_index_1stD, 
  //                    5:(4 + n_bold)] = 
  //                    uncensored_indices_padded[p,];


  //     shard_index_1stD = shard_index_1stD + 1;
  //   }
  // }
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



 // old from 18
  // For map_rect: 
  // create data_shard_real for each roi by participant
  // bold n_bold
  // design_matrix n_bold + n_bold * n_beta
  // rho_power_matrix n_bold + n_bold * n_beta + n_bold * n_bold
  // array[n_par * n_roi, n_bold + n_bold * n_beta + n_bold * n_bold] real data_shard_real;
  // shard_index_1stD = 1;
  // for (p in 1:n_par) {
  //   int bold_start = 1;
  //   int bold_stop = n_bold;
  //   for (r in 1:n_roi) {
  //     // pack bold
  //     int data_shard_2ndD_start = 1;
  //     int data_shard_2ndD_stop = n_bold;
  //     data_shard_real[shard_index_1stD, 
  //                     data_shard_2ndD_start:data_shard_2ndD_stop] = 
  //                     bold[p, bold_start:bold_stop];

  //     // pack design matrix from column-major order 
  //     data_shard_2ndD_start = data_shard_2ndD_start + n_bold;
  //     data_shard_2ndD_stop = data_shard_2ndD_stop + n_bold * n_beta;
  //     data_shard_real[shard_index_1stD, 
  //                     data_shard_2ndD_start:data_shard_2ndD_stop] = 
  //                     to_array_1d(to_matrix(design_array[p,,]));

  //     // pack rho_power_matrix for Toeplitz
  //     data_shard_2ndD_start = data_shard_2ndD_start + n_bold * n_beta;
  //     data_shard_2ndD_stop = data_shard_2ndD_stop + n_bold * n_bold;
  //     data_shard_real[shard_index_1stD, 
  //                    data_shard_2ndD_start:data_shard_2ndD_stop] = 
  //                    to_array_1d(rho_power_matrix);
      
  //     shard_index_1stD = shard_index_1stD + 1;
  //   }
  // }

  // new from 38
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
  array[n_roi] real <lower =0> mu_sigma_raw; // this still needs to be half normal
  array[n_roi] real mu_delta_raw;
  array[n_roi] real mu_rho_time_raw;
  array[n_roi] vector[n_beta_stim] mu_betas;
  // array[n_roi] vector[n_beta_stim] mu_betas;

  // think about what the priors should be fore these
  array[n_roi] real <lower =0> tau_sigma_raw;
  array[n_roi] real <lower =0>  tau_delta_raw;
  array[n_roi] real <lower =0> tau_rho_time_raw;
  array[n_roi, n_beta_stim] real <lower =0> tau_betas;
  // array[n_roi] vector[n_beta_stim] tau_betas;

  array[n_par, n_roi] real sigma_z;
  array[n_par, n_roi] real delta_z;
  array[n_par, n_roi] real rho_time_z;
  array[n_par, n_roi] vector[n_beta_stim] betas_z;
  array[n_par, n_roi] vector[n_beta_motion] betas_motion;
  // non centered mvn
  // vector[total_usable_TRs] bold_z;

  vector[0] phi; // parameters shared between all data shards; sometimes empty
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;

  array[n_par, n_roi] real <lower =0>  sigma;
  array[n_par, n_roi] real <lower =0>  delta;
  array[n_par, n_roi] real <lower =0>  rho_time;
  array[n_par, n_roi] vector[n_beta] betas;

  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma[p, r] = mu_sigma_raw[r] * 0.5 + tau_sigma_raw[r] * sigma_z[p, r];
      delta[p, r] = inv_logit((mu_delta_raw[r] * 1.75) + tau_delta_raw[r] * delta_z[p, r]);
      rho_time[p, r] = inv_logit((mu_rho_time_raw[r] * 1.75) + tau_rho_time_raw[r] * rho_time_z[p, r]);
      int b_motion_index = 1;
      for (b in 1:n_beta){ 
        if (b <= n_beta_stim) {
          betas[p,r][b] = mu_betas[r][b] + tau_betas[r,b] * betas_z[p,r][b];
        } else {
          betas[p,r][b] = betas_motion[p,r][b_motion_index];
          b_motion_index += 1;
        }
      }
    }
  }


  // local data shard parameters theta for each roi per participant
  // sigma 1 
  // delta 2
  // rho_time 3
  // betas 3 + n_beta
  array[n_par * n_roi] vector[(3 + n_beta)] theta;
  
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

}

model {
  // Priors
  // eventually move these into phi
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

  // Multilevel priors (or z-score priors for the multilevel)
  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma_z[p, r] ~ std_normal();
      delta_z[p, r] ~ std_normal();
      rho_time_z[p, r] ~ std_normal();
      betas_z[p,r] ~ std_normal();
      betas_motion[p,r] ~ normal(0, 3);
    }
  }

  target += map_rect(participant_ll, phi, theta, data_shard_real, data_shard_int);

}

generated quantities {
  // array[n_observations] real mu_pred = mu[indices_observed];
  // old model 18
  // array[n_roi] real <lower =0>  mu_sigma;
  // array[n_roi] real <lower =0>  mu_delta;
  // array[n_roi] real <lower =0>  mu_rho_time;

  // for (r in 1:n_roi) {
  //   mu_sigma[r] = mu_sigma_raw[r] *.5;
  //   mu_delta[r] = inv_logit(mu_delta_raw[r] * 1.75);
  //   mu_rho_time[r] = inv_logit(mu_rho_time_raw[r] * 1.75);
  // }

  // vector[total_usable_TRs] log_lik;
  // int shard_index = 1;
  // int log_lik_start = 1;
  // int log_lik_stop = n_bold - n_censor[1];
  // for (p in 1:n_par){
  //   for (r in 1:n_roi){
  //     if (shard_index > 1) {
  //       log_lik_start += n_bold - n_censor[p-1]; // this doesn't make sense once more rois are added
  //       log_lik_stop += n_bold - n_censor[p];
  //     }
  //     log_lik[log_lik_start:log_lik_stop] = 
  //       multi_normal_elementwise_log_lik(phi, 
  //                                        theta[shard_index], 
  //                                        data_shard_real[shard_index], 
  //                                        data_shard_int[shard_index]);
  //     shard_index += 1;
  //   }
  // }

  // from model 38
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
