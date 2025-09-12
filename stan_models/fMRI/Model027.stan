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
  array[n_par, n_roi, n_bold] real bold;
  array[n_par, n_bold] int usable_bold_indices_one_is_true;
  array[n_par] int n_censor;
  array[n_par, n_bold, n_beta] real design_array;
}

transformed data {
  int n_cue = 4;
  int n_shards = n_par * n_roi;
  int n_motion_beta = 6;
  int n_shock_per_run = 15;
  int n_trials_per_cue = 8 + 12 + 12 + 12;
  array[n_shock_per_run] int shock_linespace = linspaced_int_array(n_shock_per_run,1,n_shock_per_run);
  array[n_trials_per_cue] int cue_linespace = linspaced_int_array(n_trials_per_cue,1,n_trials_per_cue);

  array[n_trials_per_cue] int csp_beta_indices;
  array[n_trials_per_cue] int gs1_beta_indices;
  array[n_trials_per_cue] int gs2_beta_indices;
  array[n_trials_per_cue] int gs3_beta_indices;
  array[n_shock_per_run] int shock_beta_indices;

  //habituation 8 trials
  int current_start_index = 1;
  int current_stop_index = 8;
  csp_beta_indices[1:8] = linspaced_int_array(8, current_start_index, current_stop_index);
  current_start_index = current_start_index + 8;
  current_stop_index = current_stop_index + 8;
  gs1_beta_indices[1:8] = linspaced_int_array(8, current_start_index, current_stop_index);
  current_start_index = current_start_index + 8;
  current_stop_index = current_stop_index + 8;
  gs2_beta_indices[1:8] = linspaced_int_array(8, current_start_index, current_stop_index);
  current_start_index = current_start_index + 8;
  current_stop_index = current_stop_index + 8;
  gs3_beta_indices[1:8] = linspaced_int_array(8, current_start_index, current_stop_index);
  //acquisition #1 12 trials
  current_start_index = current_start_index + 8;
  current_stop_index = current_stop_index + 12;
  csp_beta_indices[1+8:8+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs1_beta_indices[1+8:8+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs2_beta_indices[1+8:8+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs3_beta_indices[1+8:8+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  //acquisition #2 12 trials
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  csp_beta_indices[1+8+12:8+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs1_beta_indices[1+8+12:8+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs2_beta_indices[1+8+12:8+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs3_beta_indices[1+8+12:8+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  //extinction 12 trials
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  csp_beta_indices[1+8+12+12:8+12+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs1_beta_indices[1+8+12+12:8+12+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs2_beta_indices[1+8+12+12:8+12+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 12;
  gs3_beta_indices[1+8+12+12:8+12+12+12] = linspaced_int_array(12, current_start_index, current_stop_index);
  //shock 15 times
  current_start_index = current_start_index + 12;
  current_stop_index = current_stop_index + 15;
  shock_beta_indices = linspaced_int_array(15, current_start_index, current_stop_index);

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
  // array[n_par, n_beta] vector[n_roi] betas_z;
  array[n_par, n_motion_beta] vector[n_roi] betas_motion_z;
  // now a shared group parameter
  vector[n_beta - n_motion_beta] betas_z;
  
  vector <lower=0> [n_roi] mu_sigma_raw; // this still needs to be half normal
  // array[n_roi] real <lower=0> mu_sigma_raw; // this still needs to be half normal
  vector <lower=0> [n_roi] tau_sigma_raw;
  vector[n_roi] mu_delta_raw;
  vector <lower=0> [n_roi]  tau_delta_raw;
  vector[n_roi] mu_rho_time_raw;
  vector <lower=0> [n_roi] tau_rho_time_raw;
  // array[n_beta] vector[n_roi] mu_betas;
  // array[n_beta] vector <lower=0> [n_roi] tau_betas;
  array[n_motion_beta] vector[n_roi] mu_betas_motion;
  array[n_motion_beta] vector <lower=0> [n_roi] tau_betas_motion;

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

  // not going to worry about correlated rois
  // vector[(n_roi*(n_roi-1))%/%2] rho_z_sigma_z_raw;
  // vector[(n_roi*(n_roi-1))%/%2] rho_z_delta_z_raw;
  // vector[(n_roi*(n_roi-1))%/%2] rho_z_rho_time_z_raw; //rho_z_rho is not a typo, it is how rho over time is correlated between ROIs
  // array[n_beta] vector[(n_roi*(n_roi-1))%/%2] rho_z_betas_z;

  real<lower=0> csp_rho;
  real<lower=0> gs1_rho;
  real<lower=0> gs2_rho;
  real<lower=0> gs3_rho;
  real<lower=0> shock_rho;
  real<lower=0> csp_alpha;
  real<lower=0> gs1_alpha;
  real<lower=0> gs2_alpha;
  real<lower=0> gs3_alpha;
  real<lower=0> shock_alpha;
  real<lower=0> csp_sigma;
  real<lower=0> gs1_sigma;
  real<lower=0> gs2_sigma;
  real<lower=0> gs3_sigma;
  real<lower=0> shock_sigma;
}

transformed parameters {
  // array[n] real amplitude_all;
  // amplitude_all[indices_observed] = amplitude;
  // amplitude_all[indices_missing] = amplitude_missing;

  array[n_par, n_roi] real <lower=0>  sigma;
  array[n_par, n_roi] real <lower=0>  delta;
  array[n_par, n_roi] real <lower=0>  rho_time;
  array[n_par, n_roi] vector[6] betas_motion;

  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma[p, r] = mu_sigma_raw[r] * 0.5 + tau_sigma_raw[r] * sigma_z[p][r];
      delta[p, r] = inv_logit((mu_delta_raw[r] * 1.75) + tau_delta_raw[r] * delta_z_raw[p][r]);
      rho_time[p, r] = inv_logit((mu_rho_time_raw[r] * 1.75) + tau_rho_time_raw[r] * rho_time_z_raw[p][r]);
      for (b in 1:n_motion_beta){
        //notice ordering changes from betas_z to betas
        betas_motion[p,r][b] = mu_betas_motion[b][r] + tau_betas_motion[b][r] * betas_motion_z[p,b][r];
      }
    }
  }

  


  matrix[n_trials_per_cue, n_trials_per_cue] csp_L_K;
  matrix[n_trials_per_cue, n_trials_per_cue] gs1_L_K;
  matrix[n_trials_per_cue, n_trials_per_cue] gs2_L_K;
  matrix[n_trials_per_cue, n_trials_per_cue] gs3_L_K;
  matrix[n_shock_per_run, n_shock_per_run] shock_L_K;

  matrix[n_trials_per_cue, n_trials_per_cue] csp_K = gp_exp_quad_cov(cue_linespace, csp_alpha, csp_rho);
  matrix[n_trials_per_cue, n_trials_per_cue] gs1_K = gp_exp_quad_cov(cue_linespace, gs1_alpha, gs1_rho);
  matrix[n_trials_per_cue, n_trials_per_cue] gs2_K = gp_exp_quad_cov(cue_linespace, gs2_alpha, gs2_rho);
  matrix[n_trials_per_cue, n_trials_per_cue] gs3_K = gp_exp_quad_cov(cue_linespace, gs3_alpha, gs3_rho);
  matrix[n_shock_per_run, n_shock_per_run] shock_K = gp_exp_quad_cov(shock_linespace, shock_alpha, shock_rho);
  
  // assuming 1 roi lazily
  vector[n_beta - n_motion_beta] betas;
  for (t in 1:n_trials_per_cue) {
    csp_K[t, t] = csp_K[t, t] + square(csp_sigma);
    gs1_K[t, t] = gs1_K[t, t] + square(gs1_sigma);
    gs2_K[t, t] = gs2_K[t, t] + square(gs2_sigma);
    gs3_K[t, t] = gs3_K[t, t] + square(gs3_sigma);
  }
  
  for (t in 1:n_shock_per_run) {
    shock_K[t, t] = shock_K[t, t] + square(shock_sigma);
  }

  // csp_L_K = cholesky_decompose(csp_K);
  // gs1_L_K = cholesky_decompose(gs1_K);
  // gs2_L_K = cholesky_decompose(gs2_K);
  // gs3_L_K = cholesky_decompose(gs3_K);
  // shock_L_K = cholesky_decompose(shock_K);
  //helps matrix be invertible 
  csp_L_K  = cholesky_decompose( add_diag(csp_K, 1e-12) );
  gs1_L_K  = cholesky_decompose( add_diag(gs1_K, 1e-12) );
  gs2_L_K  = cholesky_decompose( add_diag(gs2_K, 1e-12) );
  gs3_L_K  = cholesky_decompose( add_diag(gs3_K, 1e-12) );
  shock_L_K= cholesky_decompose( add_diag(shock_K, 1e-12) );


  // just assuming 1 roi for now lazily
  betas[csp_beta_indices] = csp_L_K * betas_z[csp_beta_indices];
  betas[gs1_beta_indices] = gs1_L_K * betas_z[gs1_beta_indices];
  betas[gs2_beta_indices] = gs2_L_K * betas_z[gs2_beta_indices];
  betas[gs3_beta_indices] = gs3_L_K * betas_z[gs3_beta_indices];
  betas[shock_beta_indices] = shock_L_K * betas_z[shock_beta_indices];



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
        theta[Pshard_index_1stD][4:(3 + n_beta)] = append_row(betas, to_vector(betas_motion[p,r]));
      
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
  for (r in 1:n_roi) {
    mu_sigma_raw[r] ~ std_normal();
    mu_delta_raw[r] ~ std_normal();
    mu_rho_time_raw[r] ~ std_normal();
    tau_sigma_raw[r] ~ std_normal();
    tau_delta_raw[r] ~ std_normal();
    tau_rho_time_raw[r] ~ std_normal();
  }
  for(b in 1:n_motion_beta) {
    mu_betas_motion[b] ~ std_normal();
    tau_betas_motion[b] ~ std_normal();
  }
  
  // Gaussian process priors
  csp_rho ~ lognormal(log(8), .5);
  gs1_rho ~ lognormal(log(8), .5);
  gs2_rho ~ lognormal(log(8), .5);
  gs3_rho ~ lognormal(log(8), .5);
  shock_rho ~ lognormal(log(8), .5);
  csp_alpha ~ std_normal();
  gs1_alpha ~ std_normal();
  gs2_alpha ~ std_normal();
  gs3_alpha ~ std_normal();
  shock_alpha ~ std_normal();
  csp_sigma ~ std_normal();
  gs1_sigma ~ std_normal();
  gs2_sigma ~ std_normal();
  gs3_sigma ~ std_normal();
  shock_sigma ~ std_normal();

  betas_z[csp_beta_indices] ~ std_normal();
  betas_z[gs1_beta_indices] ~ std_normal();
  betas_z[gs2_beta_indices] ~ std_normal();
  betas_z[gs3_beta_indices] ~ std_normal();
  betas_z[shock_beta_indices] ~ std_normal();




  // Multilevel priors (or z-score priors for the multilevel)
  for (p in 1:n_par) {
    for (r in 1:n_roi) {
      sigma_z[p, r] ~ std_normal();
      delta_z_raw[p, r] ~ std_normal();
      rho_time_z_raw[p, r] ~ std_normal();
      for (b in 1:n_motion_beta){
        betas_motion_z[p,b][r] ~ std_normal();
      }
    }
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

  // rho_z_sigma_z_raw ~ std_normal();
  // rho_z_delta_z_raw ~ std_normal();
  // rho_z_rho_time_z_raw ~ std_normal();

  // for (b in 1:n_beta) {
  //   // rho_z_mu_betas[b] ~ std_normal();
  //   // rho_z_tau_betas_raw[b] ~ std_normal();

  //   rho_z_betas_z[b] ~ std_normal();
  // }    

  // Use rhos for correlation matrices
//   matrix[n_roi, n_roi] Corr_mu_sigma_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_tau_sigma_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_mu_delta_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_tau_delta_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_mu_rho_time_raw = rep_matrix(1, n_roi, n_roi);
//   matrix[n_roi, n_roi] Corr_tau_rho_time_raw = rep_matrix(1, n_roi, n_roi);
  
  // matrix[n_roi, n_roi] Corr_sigma_z_raw = rep_matrix(1, n_roi, n_roi);
  // matrix[n_roi, n_roi] Corr_delta_z_raw = rep_matrix(1, n_roi, n_roi);
  // matrix[n_roi, n_roi] Corr_rho_time_z_raw = rep_matrix(1, n_roi, n_roi);

//   array[n_beta] matrix[n_roi, n_roi] Corr_mu_betas = rep_array(rep_matrix(1, n_roi, n_roi), n_beta);
//   array[n_beta] matrix[n_roi, n_roi] Corr_tau_betas_raw = rep_array(rep_matrix(1, n_roi, n_roi), n_beta);
  
  // array[n_beta] matrix[n_roi, n_roi] Corr_betas_z = rep_array(rep_matrix(1, n_roi, n_roi), n_beta);
  
  // int rho_idx = 1;
  // for (i in 1:(n_roi-1)) {
  //   for (j in (i + 1):n_roi) {
  //     // group
  //   //   Corr_mu_sigma_raw[i,j] = tanh(rho_z_mu_sigma_raw[rho_idx]);
  //   //   Corr_mu_sigma_raw[j,i] = tanh(rho_z_mu_sigma_raw[rho_idx]);

  //   //   Corr_tau_sigma_raw[i,j] = tanh(rho_z_tau_sigma_raw[rho_idx]);
  //   //   Corr_tau_sigma_raw[j,i] = tanh(rho_z_tau_sigma_raw[rho_idx]);
      
  //   //   Corr_mu_delta_raw[i,j] = tanh(rho_z_mu_delta_raw[rho_idx]);
  //   //   Corr_mu_delta_raw[j,i] = tanh(rho_z_mu_delta_raw[rho_idx]);

  //   //   Corr_tau_delta_raw[i,j] = tanh(rho_z_tau_delta_raw[rho_idx]);
  //   //   Corr_tau_delta_raw[j,i] = tanh(rho_z_tau_delta_raw[rho_idx]);

  //   //   Corr_mu_rho_time_raw[i,j] = tanh(rho_z_mu_rho_time_raw[rho_idx]);
  //   //   Corr_mu_rho_time_raw[j,i] = tanh(rho_z_mu_rho_time_raw[rho_idx]);
      
  //   //   Corr_tau_rho_time_raw[i,j] = tanh(rho_z_tau_delta_raw[rho_idx]);
  //   //   Corr_tau_rho_time_raw[j,i] = tanh(rho_z_tau_delta_raw[rho_idx]);
      
  //     // local
  //     Corr_sigma_z_raw[i,j] = tanh(rho_z_sigma_z_raw[rho_idx]);
  //     Corr_sigma_z_raw[j,i] = tanh(rho_z_sigma_z_raw[rho_idx]);

  //     Corr_delta_z_raw[i,j] = tanh(rho_z_delta_z_raw[rho_idx]);
  //     Corr_delta_z_raw[j,i] = tanh(rho_z_delta_z_raw[rho_idx]);

  //     Corr_rho_time_z_raw[i,j] = tanh(rho_z_rho_time_z_raw[rho_idx]);
  //     Corr_rho_time_z_raw[j,i] = tanh(rho_z_rho_time_z_raw[rho_idx]);
  //     for (b in 1:n_beta) {
  //       // group
  //       // Corr_mu_betas[b][i,j] = tanh(rho_z_mu_betas[b][rho_idx]);
  //       // Corr_mu_betas[b][j,i] = tanh(rho_z_mu_betas[b][rho_idx]);
        
  //       // Corr_tau_betas_raw[b][i,j] = tanh(rho_z_tau_betas_raw[b][rho_idx]);
  //       // Corr_tau_betas_raw[b][j,i] = tanh(rho_z_tau_betas_raw[b][rho_idx]);
  //       //local
  //       Corr_betas_z[b][i,j] = tanh(rho_z_betas_z[b][rho_idx]);
  //       Corr_betas_z[b][j,i] = tanh(rho_z_betas_z[b][rho_idx]);
  //     }
  //     rho_idx += 1;
  //   }
  // }

//   matrix[n_roi, n_roi] L_Corr_mu_sigma = cholesky_decompose(Corr_mu_sigma_raw);
//   matrix[n_roi, n_roi] L_Corr_tau_sigma = cholesky_decompose(Corr_tau_sigma_raw);
//   matrix[n_roi, n_roi] L_Corr_mu_delta_raw = cholesky_decompose(Corr_mu_delta_raw);
//   matrix[n_roi, n_roi] L_Corr_tau_delta_raw = cholesky_decompose(Corr_tau_delta_raw);
//   matrix[n_roi, n_roi] L_Corr_mu_rho_time_raw = cholesky_decompose(Corr_mu_rho_time_raw);
//   matrix[n_roi, n_roi] L_Corr_tau_rho_time_raw = cholesky_decompose(Corr_tau_rho_time_raw);
  
  // matrix[n_roi, n_roi] L_Corr_sigma_z = cholesky_decompose(Corr_sigma_z_raw);
  // matrix[n_roi, n_roi] L_Corr_delta_z_raw = cholesky_decompose(Corr_delta_z_raw);
  // matrix[n_roi, n_roi] L_Corr_rho_time_z_raw = cholesky_decompose(Corr_rho_time_z_raw);
  
//   array[n_beta] matrix[n_roi, n_roi] L_Corr_mu_betas;
//   array[n_beta] matrix[n_roi, n_roi] L_Corr_tau_betas_raw;
  
  // array[n_beta] matrix[n_roi, n_roi] L_Corr_betas_z;

  // for (b in 1:n_beta) {
  //   // L_Corr_mu_betas[b] = cholesky_decompose(Corr_mu_betas[b]);
  //   // L_Corr_tau_betas_raw[b] = cholesky_decompose(Corr_tau_betas_raw[b]);
    
  //   L_Corr_betas_z[b] = cholesky_decompose(Corr_betas_z[b]);
  // }    

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

  // for(p in 1:n_par){
  //   sigma_z[p] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_sigma_z);
  //   delta_z_raw[p] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_delta_z_raw);
  //   rho_time_z_raw[p] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_rho_time_z_raw);
  //   for (b in 1:n_beta) {
  //     betas_z[p,b] ~ multi_normal_cholesky(rep_vector(0,n_roi), L_Corr_betas_z[b]);
  //   }
  // }




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
