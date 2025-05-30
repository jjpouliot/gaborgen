data {
  int<lower=0> n;
  int<lower=0> n_observations;
  int<lower=0> n_missing;
  int<lower=0> n_participants;
  int<lower=0> n_trials; 
  int<lower=0> n_phases;
  int<lower=0> n_blocks;
  int<lower=0> n_cues;

  array[n_observations] real amplitude;
  array[n] int participant;
  array[n] int phase;
  array[n] int cue;
  array[n] int trial;
  array[n] int cue_trial_count;
  array[n] int phase_trial_count;
  array[n] int block_trial_count;
  array[n_participants, n_blocks, n_cues] real arousal_pbc_centered;
  array[n_participants, n_blocks, n_cues] real arousal_pbc;
  array[n] int block;
  array[n] int <lower=0, upper=1> paired;
  array[n_observations] int<lower=1, upper=n_observations + n_missing> indices_observed;
  array[n_missing] int<lower=1, upper=n_observations + n_missing> indices_missing;
}

parameters {

  real<lower = 0> sigma_average; //this won't hug zero so it does not need to be exp() /1
  real sigma_sd_raw; //2
  real intercept_average; //3
  real intercept_sd_raw; //4
  real fatigue_average; //5
  real fatigue_sd_raw; //6
  real intercept_fatigue_corr_raw; //7
  real<lower=0, upper=1> learning_paired_average; //8
  real<lower=0, upper=1> learning_unpaired_average; //9
  real learning_paired_sd_raw; //10
  real learning_unpaired_sd_raw;  //11
  real scaling_mean;  //12
  real scaling_sd_raw;  //13
  
  array[n_participants] real<lower = 0> sigma; //this won't hug zero so it does not need to be exp() //37 if n_participants is 24
  array[n_participants] real intercept; //61
  array[n_participants] real fatigue; //85
  array[n_participants] real<lower=0, upper=1> learning_paired; //109
  array[n_participants] real<lower=0, upper=1> learning_unpaired; //133
  array[n_cues] real scaling; //137
  array[n_missing] real amplitude_missing;
}

transformed parameters {
  array[n] real amplitude_all;
  amplitude_all[indices_observed] = amplitude;
  amplitude_all[indices_missing] = amplitude_missing;
  
  real sigma_sd = exp(sigma_sd_raw);
  real intercept_sd = exp(intercept_sd_raw);
  real fatigue_sd = exp(fatigue_sd_raw);
  
  real intercept_fatigue_corr = -inv_logit(intercept_fatigue_corr_raw);

  cov_matrix[2] Sigma;
  Sigma[1,1] = square(intercept_sd);
  Sigma[2,2] = square(fatigue_sd);
  Sigma[1,2] = intercept_fatigue_corr * intercept_sd * fatigue_sd;
  Sigma[2,1] = Sigma[1,2];
  
  matrix[2,2] Sigma_L = cholesky_decompose(Sigma);
  
  real learning_paired_sd = exp(learning_paired_sd_raw);
  real learning_unpaired_sd = exp(learning_unpaired_sd_raw);
  real scaling_sd = exp(scaling_sd_raw);
  
  
  array[n] real<lower=-0.000001,upper=1> CSP_associative_strength;
  array[n-1] real<lower=-1, upper=1> CSP_change;

  // Estimate CSP_associative_strength, Rescola Wagner type equation
  for(i in 1:(n-1)) {
    if (phase[i] == 1){ // Habituation phase
      CSP_change[i] = 0;
      CSP_associative_strength[i] = 0;
        CSP_associative_strength[i+1] = fmin(fmax(CSP_associative_strength[i+1], 0.0), 1.0);
    } else if (cue[i] == 1 && phase[i] != 1) {
        if(paired[i] == 1){
          CSP_change[i] = learning_paired[participant[i]] * 
                          (paired[i] - CSP_associative_strength[i]);
        } else {
          CSP_change[i] = learning_unpaired[participant[i]] * 
                          (paired[i] - CSP_associative_strength[i]);
        }
        CSP_associative_strength[i+1] = CSP_associative_strength[i] + CSP_change[i];

        // Ensure CSP_associative_strength stays within the bounds
        CSP_associative_strength[i+1] = fmin(fmax(CSP_associative_strength[i+1], 0.0), 1.0);
    } else {
       CSP_change[i] = 0;
       CSP_associative_strength[i+1] = CSP_associative_strength[i];
    }
  }  
}

model {
  // Priors
  sigma_average ~ normal(1, 0.5);
  sigma_sd_raw ~ normal(-3, 1);

  sigma ~ student_t(n_participants-1, sigma_average, sigma_sd);

  intercept_average ~ normal(0, 0.75);
  intercept_sd_raw ~ normal(-1, 1);
  fatigue_average ~ normal(0, 0.01);
  fatigue_sd_raw ~ normal(-4.5, 0.75);
  
  intercept_fatigue_corr_raw ~ normal(0, 1.75); // this creates a roughly uniform distribution, if confused visualize in R with hist(boot::inv.logit(rnorm(5000,0,1.75)))
  
  for (p in 1:n_participants) {
    [intercept[p], fatigue[p]]' ~ multi_student_t_cholesky(n_participants - 1,
                                                           [intercept_average, fatigue_average]',
                                                           Sigma_L);
  }

  scaling_mean ~ normal(0,1);
  scaling_sd_raw ~ normal(-0.5,1.5);
  // scaling_sd_raw ~ normal(-0.5,1);
  for (q in 1:n_cues) {
    scaling[q] ~ student_t(n_cues-1, scaling_mean, scaling_sd); 
  }

  learning_paired_average ~ beta(1.1, 1.1);
  learning_paired_sd_raw ~ normal(-0.5,1.5);
  // learning_paired_sd_raw ~ normal(-0.5,1);
  learning_paired ~ student_t(n_participants-1, 
                                  learning_paired_average, 
                                  learning_paired_sd); 

  learning_unpaired_average ~ beta(1.1, 1.1);
  learning_unpaired_sd_raw ~ normal(-0.5,1.5);
  // learning_unpaired_sd_raw ~ normal(-0.5,1);
  learning_unpaired ~ student_t(n_participants-1, 
                                    learning_unpaired_average, 
                                    learning_unpaired_sd); 


  // Likelihood
  for(i in 1:n) {
    real mu = intercept[participant[i]] +
              (fatigue[participant[i]] * trial[i]) +
              (scaling[cue[i]] * CSP_associative_strength[i]);

    amplitude_all[i] ~ normal(mu, sigma[participant[i]]);
  }
}

generated quantities {
  // These mu_pred and log_lik names cannot change
  
  // For model predictions
  array[n_observations] real mu_pred;
  // This is necessary for cross-validation
  array[n_observations] real log_lik;
  
  for (i in 1:n_observations) {
    mu_pred[i] = intercept[participant[indices_observed[i]]] +
                 (fatigue[participant[indices_observed[i]]] * trial[indices_observed[i]]) +
                 (scaling[cue[indices_observed[i]]] * CSP_associative_strength[indices_observed[i]]);
                 
    log_lik[i] = normal_lpdf(amplitude[i] | mu_pred[i], sigma[participant[indices_observed[i]]]);
  }
}
