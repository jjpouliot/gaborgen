library(tidyverse)
library(cmdstanr)
library(R.matlab)


# Load data ####
parent_folder <- "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/"
single_trial_mat_folder <- paste0(parent_folder,"single_trial_timeseries_FFTs")

all_mat_filepaths <- list.files(path = single_trial_mat_folder,full.names = T)

habituation_raw_fft_paths <- all_mat_filepaths[
  str_detect(all_mat_filepaths, pattern = "conditionS1[1-4]_") & str_detect(all_mat_filepaths, pattern = "ChanFrequency")
]



hold_data <- 
  readMat()
