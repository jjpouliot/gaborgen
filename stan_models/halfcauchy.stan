data {
  real<lower=0> mu;  // Scale parameter for Half-Cauchy
  real<lower=0> scale;  // Scale parameter for Half-Cauchy
  // real<lower=0> alpha;  // Scale parameter for Half-Cauchy
  // real<lower=0> beta;  // Scale parameter for Half-Cauchy
}

parameters {
  real<lower=0> y;  // Samples 
}

model {
  y ~ cauchy(mu, scale);  // Half-Cauchy distribution
  // y ~ gamma(alpha, beta);  // gamma distribution
}

