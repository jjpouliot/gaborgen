data {

}

parameters {
  real intercept_average;
  real fatigue_average;
  real intercept;
  real fatigue;
  real intercept_sd_raw;
  real fatigue_sd_raw;
  
  real intercept_fatigue_corr_raw;
  
  // real intercept;
  // real fatigue;
  
}

transformed parameters {
  // Convert raw SD parameters to positive scale
  real intercept_sd = exp(intercept_sd_raw);
  real fatigue_sd   = exp(fatigue_sd_raw);
  real intercept_fatigue_corr = -inv_logit(intercept_fatigue_corr_raw);
  
  cov_matrix[2] Sigma;
  Sigma[1,1] = square(intercept_sd);
  Sigma[2,2] = square(fatigue_sd);
  Sigma[1,2] = intercept_fatigue_corr * intercept_sd * fatigue_sd;
  Sigma[2,1] = Sigma[1,2];
  
  matrix[2,2] Sigma_L = cholesky_decompose(Sigma);
}

model {
  
  intercept_average ~ normal(0,0.2);
  intercept_sd_raw ~ normal(-1,1);
  fatigue_average ~ normal(0,0.01);
  fatigue_sd_raw ~ normal(-4.5,.75);
  intercept_fatigue_corr_raw ~ normal(0, 1.75); 

[intercept, fatigue]' ~ multi_student_t_cholesky(23,
                                                 [intercept_average, fatigue_average]', 
                                                 Sigma_L);

}

