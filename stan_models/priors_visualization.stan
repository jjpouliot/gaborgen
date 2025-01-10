data {
real intecept_prior_sd;
real fatigue_prior_sd;
}

parameters {
  real intercept;
  real fatigue;
  // real intercept_sd_raw;
  // real fatigue_sd_raw;
  
  real intercept_fatigue_corr_raw;
  
  // real intercept;
  // real fatigue;
  
}

transformed parameters {
  // Convert raw SD parameters to positive scale
  // real intercept_sd = exp(intercept_sd_raw);
  // real fatigue_sd   = exp(fatigue_sd_raw);
  real intercept_fatigue_corr = -inv_logit(intercept_fatigue_corr_raw);
  
  cov_matrix[2] Sigma;
  Sigma[1,1] = square(intecept_prior_sd);
  Sigma[2,2] = square(fatigue_prior_sd);
  Sigma[1,2] = intercept_fatigue_corr * intecept_prior_sd * fatigue_prior_sd;
  Sigma[2,1] = Sigma[1,2];
  
  matrix[2,2] Sigma_L = cholesky_decompose(Sigma);
}

model {
  
  intercept_fatigue_corr_raw ~ normal(0, 1.75); 
  
  [intercept, fatigue]' ~ multi_normal_cholesky([0, 0]', Sigma_L);
  
}

