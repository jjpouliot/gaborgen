library(tidyverse)
library(cmdstanr)
library(patchwork)

# we recommend running this in a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

model013_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_1.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_2.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_3.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_4.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_5.csv"
  )
)

model013_fit_mpi_summary <- model013_fit_mpi$summary()


model013_mpi_draws <- model013_fit_mpi$draws(format = "df")

model013_loo <- model013_fit_mpi$loo()


model016_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/model016_chain_65887620_1.csv"
  )
)

model016_draws <- model016_fit_mpi$draws(format = "df")

model018_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/model018_chain_65957320_1.csv"
  )
)

model018_draws <- model018_fit_mpi$draws(format = "df")

names(model016_draws)
names(model018_draws)

model016_draws %>%
  select(starts_with('betas[')) %>%
  select(ends_with(c(",1,5]", ",1,13]"))) %>%
  # pivot_longer(cols = ends_with(",1,5]")) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model018_draws %>%
  select(starts_with('betas[')) %>%
  select(ends_with(c(",1,5]", ",1,13]"))) %>%
  # pivot_longer(cols = ends_with(",1,5]")) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model016_draws %>%
  select(starts_with(c(
    'betas[8,1,5]',
    'betas[13,1,5]',
    'betas[14,1,5]',
    'betas[15,1,5]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,5]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model018_draws %>%
  select(starts_with(c(
    'betas[8,1,5]',
    'betas[13,1,5]',
    'betas[14,1,5]',
    'betas[15,1,5]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,5]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model016_draws %>%
  select(starts_with(c(
    'betas[8,1,13]',
    'betas[13,1,13]',
    'betas[14,1,13]',
    'betas[15,1,13]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model018_draws %>%
  select(starts_with(c(
    'betas[8,1,13]',
    'betas[13,1,13]',
    'betas[14,1,13]',
    'betas[15,1,13]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))


# can I get the per element log_lik

model018_draws$`log_lik[1]`[1]

fmri_stan_list$usable_bold_indices_one_is_true[1, ]

cov_mat <- matrix(nrow = 1070, ncol = 1070)

for (i in 1:1070) {
  for (j in 1:1070) {
    if (i == j) {
      cov_mat[i, j] = (model018_draws$`sigma[1,1]`[1])^2
    } else {
      cov_mat[i, j] = model018_draws$`delta[1,1]`[1] *
        ((model018_draws$`rho_time[1,1]`[1])^abs(i - j)) *
        ((model018_draws$`sigma[1,1]`[1])^2)
    }
  }
}

l_cov_mat <- t(chol(cov_mat[
  fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1,
  fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1
]))

dm <- fmri_stan_list$design_array[1, , ] %>% matrix(nrow = 1070)

beta_vector <- model018_draws %>%
  select(starts_with("betas[1,")) %>%
  slice(1) %>%
  unlist() %>%
  matrix(ncol = 1)

# visualize
(dm %*% beta_vector) %>% plot(type = "l")
fmri_stan_list$bold[1, ] %>% plot(type = "l")

mu <- dm %*% beta_vector
bold <- fmri_stan_list$bold[1, ]

mu_censor <- mu[fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1]
bold_censor <- bold[fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1]

mu_censor %>% plot(type = "l")
bold_censor %>% plot(type = "l")

data.frame(
  time = seq(0, 1069 * 2, by = 2),
  usable = fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1,
  mu,
  bold
) %>%
  # filter(usable ==T) %>%
  ggplot() +
  geom_line(aes(x = time, y = bold)) +
  geom_line(aes(x = time, y = mu), color = "red") +
  theme_bw()

bold_censor_z <- solve(l_cov_mat, (bold_censor - mu_censor))

bold_censor_z %>% plot(type = "l")

log_liks <- vector(length = length(bold_censor))

for (i in 1:length(bold_censor)) {
  log_liks[i] <- dnorm(bold_censor_z[i], log = T) - log(l_cov_mat[i, i])
}

# matches perfectly
sum(log_liks)
model018_draws$`log_lik[1]`[1]

# vectorized cov_mat for fun
# pull out your parameters
sigma_val <- model018_draws$`sigma[1,1]`[1]
delta_val <- model018_draws$`delta[1,1]`[1]
rho_val <- model018_draws$`rho_time[1,1]`[1]
s2 <- sigma_val^2

# index vector
n <- 1070
idx <- seq_len(n)

# 1) build the |i-j| distance matrix
dist_mat <- abs(outer(idx, idx, FUN = "-"))

# 2) fill in covariances: delta*rho^|i-j|*sigma^2 everywhere
cov_mat <- s2 * delta_val * rho_val^dist_mat

# 3) override the diagonal to be sigma^2
diag(cov_mat) <- s2
