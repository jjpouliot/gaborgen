library(tidyverse)
library(cmdstanr)
library(patchwork)


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
