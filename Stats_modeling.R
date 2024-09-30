library(tidyverse)
library(cmdstanr)
library(R.matlab)


# Load data ####
parent_folder <- "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI"
git_repository <- "/home/andrewf/Repositories/gaborgen"

participants_without_complete_trial_list_in_dat_file <- c("110")

channel_information <- readMat(paste0(git_repository, "/chanlocs_ssv4att_MRI.mat")) 

single_trial_mat_folder <- paste0(parent_folder,"/single_trial_timeseries_FFTs")

all_mat_filepaths <- list.files(path = single_trial_mat_folder,full.names = T)
all_mat_filepaths <- all_mat_filepaths[!grepl(x = all_mat_filepaths,
                                              pattern = paste0("participant", 
                                                                participants_without_complete_trial_list_in_dat_file))]

raw_fft_paths <- all_mat_filepaths[str_detect(all_mat_filepaths, 
                                              pattern = "rawChanFrequency")]

sliding_window_fft_paths <- all_mat_filepaths[str_detect(all_mat_filepaths, 
                                                         pattern = "slidingWindow_ChanFrequency")]




channel_labels_in_order <- channel_information$chanlocs[1,,] %>% unlist()
ssVEP_frequencyHz <- 15

raw_fft_df <- data.frame("participant" = character(), 
                         "channel" = factor(levels = channel_labels_in_order), 
                         "trial" = numeric(), 
                         "phase" = factor(levels = c("habituation","acquisition","extinction")), 
                         "cue" = factor(levels = c("CSP", "GS1", "GS2", "GS3")), 
                         "paired" = factor(levels = c("shock", "no_shock")), 
                         "amplitude" = numeric())

for (path_index in 1:length(raw_fft_paths)) {
  mat_data <- readMat(raw_fft_paths[path_index])
  
  ssVEP_frequency_index <- which.min(abs(mat_data$rawFreqs - ssVEP_frequencyHz))

  amplitude <- mat_data$currentRawFFT[,ssVEP_frequency_index]

  participant <- stringr::str_extract(raw_fft_paths[path_index], "(?<=participant)[0-9]+")
  
  trial <- stringr::str_extract(raw_fft_paths[path_index], 
                                "(?<=trial)[0-9]+") %>% 
           as.numeric()

  cue_marker <- stringr::str_extract(raw_fft_paths[path_index], 
                                "(?<=condition)S[0-9]+")
  
  if(cue_marker == "S121"){
    paired <- rep("shock", length(amplitude)) %>% factor(levels = c("shock", "no_shock"))

  } else {
    paired <- rep("no_shock", length(amplitude)) %>% factor(levels = c("shock", "no_shock"))
  }
  
  if(cue_marker %in% c("S11", "S12", "S13", "S14")) {
    phase <- factor(rep("habituation", length(amplitude), levels = c("habituation","acquisition","extinction")))
  } else if (cue_marker %in% c("S121", "S21", "S22", "S23", "S24")){
    phase <- factor(rep("acquisition", length(amplitude), levels = c("habituation","acquisition","extinction")))
  } else if (cue_marker %in% c("S31", "S32", "S33", "S34")){
    phase <- factor(rep("extinction", length(amplitude), levels = c("habituation","acquisition","extinction")))

  }

  cue <- case_when(cue_marker %in% c("S11", "S21", "S31", "S121") ~ "CSP",
                   cue_marker %in% c("S12", "S22", "S32")         ~ "GS1",
                   cue_marker %in% c("S13", "S23", "S33")         ~ "GS2",
                   cue_marker %in% c("S14", "S24", "S34")         ~ "GS3")
  
  cue <- factor(cue, levels = c("CSP", "GS1", "GS2", "GS3"))


  current_raw_fft_df <- data.frame("participant" = rep(participant, length(amplitude)),
                                   "channel" = factor(channel_labels_in_order, 
                                                      levels = channel_labels_in_order), 
                                   "trial" = rep(trial, length(amplitude)), 
                                   "phase" = phase, 
                                   "cue" = cue, 
                                   "paired" = paired, 
                                   "amplitude" = amplitude)

  raw_fft_df <- rbind.data.frame(raw_fft_df, current_raw_fft_df)
  
}

sliding_window_fft_df <- data.frame("participant" = character(), 
                                    "channel" = factor(levels = channel_labels_in_order), 
                                    "trial" = numeric(), 
                                    "phase" = factor(levels = c("habituation","acquisition","extinction")), 
                                    "cue" = factor(levels = c("CSP", "GS1", "GS2", "GS3")), 
                                    "paired" = factor(levels = c("shock", "no_shock")), 
                                    "amplitude" = numeric())
                                    
for (path_index in 1:length(sliding_window_fft_paths)) {
  mat_data <- readMat(sliding_window_fft_paths[path_index])
  
  ssVEP_frequency_index <- which.min(abs(mat_data$slidingWindowFreqs - ssVEP_frequencyHz))

  amplitude <- mat_data$currentSlidingWindowFFT[,ssVEP_frequency_index]

  participant <- stringr::str_extract(sliding_window_fft_paths[path_index], "(?<=participant)[0-9]+")
  
  trial <- stringr::str_extract(sliding_window_fft_paths[path_index], 
                                "(?<=trial)[0-9]+") %>% 
           as.numeric()

  cue_marker <- stringr::str_extract(sliding_window_fft_paths[path_index], 
                                "(?<=condition)S[0-9]+")
  
  if(cue_marker == "S121"){
    paired <- rep("shock", length(amplitude)) %>% factor(levels = c("shock", "no_shock"))

  } else {
    paired <- rep("no_shock", length(amplitude)) %>% factor(levels = c("shock", "no_shock"))
  }
  
 if(cue_marker %in% c("S11", "S12", "S13", "S14")) {
    phase <- factor(rep("habituation", length(amplitude), levels = c("habituation","acquisition","extinction")))
  } else if (cue_marker %in% c("S121", "S21", "S22", "S23", "S24")){
    phase <- factor(rep("acquisition", length(amplitude), levels = c("habituation","acquisition","extinction")))
  } else if (cue_marker %in% c("S31", "S32", "S33", "S34")){
    phase <- factor(rep("extinction", length(amplitude), levels = c("habituation","acquisition","extinction")))

  }

  cue <- case_when(cue_marker %in% c("S11", "S21", "S31", "S121") ~ "CSP",
                   cue_marker %in% c("S12", "S22", "S32")         ~ "GS1",
                   cue_marker %in% c("S13", "S23", "S33")         ~ "GS2",
                   cue_marker %in% c("S14", "S24", "S34")         ~ "GS3")
  
  cue <- factor(cue, levels = c("CSP", "GS1", "GS2", "GS3"))

  current_sliding_window_fft_df <- data.frame("participant" = rep(participant, length(amplitude)),
                                              "channel" = factor(channel_labels_in_order, 
                                                      levels = channel_labels_in_order), 
                                              "trial" = rep(trial, length(amplitude)), 
                                              "phase" = phase, 
                                              "cue" = cue, 
                                              "paired" = paired, 
                                              "amplitude" = amplitude)

  sliding_window_fft_df <- rbind.data.frame(sliding_window_fft_df, current_sliding_window_fft_df)
  
}

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant, phase, cue, paired) %>% 
  summarise(mean_amp = mean(amplitude),
            se_amplitude = plotrix::std.error(amplitude),
            n()) %>% 
  print(n = 9999)

sliding_window_fft_df %>% 
  filter(channel %in% c("Oz")) %>% 
  group_by(phase, cue, paired) %>% 
  summarise(mean_amp = mean(amplitude),
            se_amplitude = plotrix::std.error(amplitude))

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  ggplot() +
  geom_point(aes(x = trial, y = amplitude, color = cue))

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(trial) %>% 
  summarise_all(mean) %>% 
  ggplot() +
  geom_line(aes(x = trial, y = amplitude)) +
  ggtitle("Average Oz ssVEP per Trial") +
  theme_bw() +
  theme(text = element_text(size = 22)) 

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(trial,participant) %>% 
  summarise_all(mean) %>% 
  ggplot() +
  geom_line(aes(x = trial, y = amplitude, color = participant))

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude))

raw_fft_df %>%
  filter(channel %in% c("Oz", "O1", "O2", "POz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  theme_bw()

sliding_window_fft_df %>%
  filter(channel %in% c("Oz", "O1", "O2", "POz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  theme_bw()

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  facet_grid(. ~ phase) +
  ggtitle("Average ssVEP per cue by phase") +
  theme_bw() +
  theme(text = element_text(size = 22)) 

sliding_window_fft_df %>%
  filter(channel %in% c("Oz", "O1", "O2", "POz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  facet_grid(. ~ phase) +
  theme_bw()

sliding_window_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  theme_bw()

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  theme_bw()




raw_fft_df %>%
  filter(channel %in% c("Oz"),
         paired == "no_shock") %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  ungroup() %>% 
  select(-c("channel", "phase", "paired", "participant")) %>% 
  group_by(trial) %>% 
  summarise_all(mean) %>% 
  ggplot() +
  geom_line(aes(x = trial, y = z_scored_amplitude))

sliding_window_fft_df %>%
  filter(channel %in% c("Oz"),
         paired == "no_shock") %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  ungroup() %>% 
  select(-c("channel", "phase", "paired", "participant")) %>% 
  group_by(trial) %>% 
  summarise_all(mean) %>% 
  ggplot() +
  geom_line(aes(x = trial, y = z_scored_amplitude))

sliding_window_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude))



# habituation_raw_fft_paths <- all_mat_filepaths[
#   str_detect(all_mat_filepaths, 
#     pattern = "conditionS1[1-4]_") & str_detect(all_mat_filepaths, 
#                                                 pattern = "rawChanFrequency")
# ]

# acquisition_raw_fft_paths <- all_mat_filepaths[
#   str_detect(all_mat_filepaths, 
#     pattern = "2[1-4]_samp") & str_detect(all_mat_filepaths, 
#                                                 pattern = "rawChanFrequency")
# ]

# Modeling ####

# Prepare data for modeling ####
# must sort df some of these models
# going to have to load dat files to know what is missing

raw_stan_list <- list()

raw_fft_df_Oz <- raw_fft_df %>% 
  filter(channel %in% c("Oz")) %>% 
  arrange(participant, trial)


raw_stan_list$n_observations <- nrow(raw_fft_df_Oz)

raw_stan_list$n_participants <- unique(raw_fft_df_Oz$participant) %>% length()

raw_stan_list$n_trials <- unique(raw_fft_df_Oz$trial) %>% length()

participant_by_n <- raw_fft_df_Oz$participant %>% unique() %>% rep(each = raw_stan_list$n_trials)

raw_stan_list$participant <- participant_by_n %>% as.factor() %>% as.integer()

raw_stan_list$n <- raw_stan_list$n_trial * raw_stan_list$n_participants
raw_stan_list$n_missing <- raw_stan_list$n - raw_stan_list$n_observations

## Add in missing trials
indices_observed <- integer()
indices_missing <- integer()
phase <- integer()
block <- integer()
cue <- integer()
trial <- integer()
paired <- integer()
time <- double()
current_trial_number <- 1
current_observation <- 1

for (i in 1:raw_stan_list$n) {

  current_participant <- participant_by_n[i]
  current_participant_log_filepaths <- list.files(path = paste0(parent_folder, "/raw_data/"),
                                              pattern = paste0(current_participant, "_logfile.dat$"), recursive = T, full.names = T)

  current_day1_log_filepath <- current_participant_log_filepaths[grepl(pattern = "Day1", 
                                                                       x = current_participant_log_filepaths)]
  
  if (length(current_day1_log_filepath) != 1) {
    stop(paste0("more or less than one day1 log file for participant ", current_participant))
  }

  current_day1_log_file <- read.table(current_day1_log_filepath, header = T, sep = ",")

  if(nrow(current_day1_log_file) != raw_stan_list$n_trials){
    stop(paste0("Participant ", current_participant, " does not have all trials in log file. Add this person as a useable participant at top of script."))

  }

  if(current_trial_number == raw_fft_df_Oz$trial[current_observation]){
    indices_observed <- c(indices_observed, i)

    phase <- c(phase, current_day1_log_file$phase[current_trial_number])
    block <- c(block, current_day1_log_file$block[current_trial_number])
    cue <- c(cue, current_day1_log_file$stim[current_trial_number])
    trial <- c(trial, current_day1_log_file$trial[current_trial_number])
    paired <- c(paired, current_day1_log_file$paired[current_trial_number])
    time <- c(time, current_day1_log_file$timeSinceFirstTR[current_trial_number])

    current_observation <- current_observation + 1
    current_trial_number <- current_trial_number + 1
    if(current_trial_number > raw_stan_list$n_trial){
      current_trial_number <- 1
    }
  } else {
    indices_missing <- c(indices_missing, i)

    phase <- c(phase, current_day1_log_file$phase[current_trial_number])
    block <- c(block, current_day1_log_file$block[current_trial_number])
    cue <- c(cue, current_day1_log_file$stim[current_trial_number])
    trial <- c(trial, current_day1_log_file$trial[current_trial_number])
    paired <- c(paired, current_day1_log_file$paired[current_trial_number])
    time <- c(time, current_day1_log_file$timeSinceFirstTR[current_trial_number])

    current_trial_number <- current_trial_number + 1
    if(current_trial_number > raw_stan_list$n_trial){
      current_trial_number <- 1
    }
  }
  if(is.na(phase[i])){
    stop()
  }
}

raw_stan_list$indices_observed <- indices_observed %>% as.integer()
raw_stan_list$indices_missing <- indices_missing %>% as.integer()

raw_stan_list$phase <- phase %>% as.integer()
raw_stan_list$block <- block %>% as.integer()
raw_stan_list$cue <- cue %>% as.integer()
raw_stan_list$trial <- trial %>% as.integer()
raw_stan_list$paired <- paired %>% as.integer()
raw_stan_list$time <- time

raw_stan_list$n_phases <- unique(raw_stan_list$phase) %>% length()
raw_stan_list$n_blocks <- unique(raw_stan_list$block) %>% length()
raw_stan_list$n_cues <- unique(raw_stan_list$cue) %>% length()

raw_stan_list$amplitude <- raw_fft_df_Oz$amplitude


## Model settings
number_of_chains <- 4
number_of_parallel_chains  <- 4
posterior_samples_per_chain <- 1000
warmup_samples_per_chain <- 1000
where_to_save_chains <- paste0(parent_folder,"/stan_chains")

## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_fatigue_single_channel.stan")

#force_recompile = T is sometimes helpful
cue_by_phase_1chan_fatigue <- cmdstanr::cmdstan_model(model_path, 
                                                      force_recompile = T)

#Model source code
# cue_by_phase_1chan_fatigue$print()

# Clear previous chains
list.files(pattern = "amp_ML_fatigue_single_channel",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

cue_by_phase_1chan_fatigue_fit <- cue_by_phase_1chan_fatigue$sample(data = raw_stan_list,
                                                                    refresh = 250,
                                                                    seed = 2,
                                                                    iter_warmup = warmup_samples_per_chain, 
                                                                    iter_sampling = posterior_samples_per_chain, 
                                                                    save_warmup = F, 
                                                                    show_messages = T,
                                                                    output_dir = where_to_save_chains,
                                                                    chains = number_of_chains,
                                                                    parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
cue_by_phase_1chan_fatigue_fit_meta_data <- cue_by_phase_1chan_fatigue_fit$metadata()

cue_by_phase_1chan_fatigue_fit_relevant_parameters <- cue_by_phase_1chan_fatigue_fit_meta_data$model_params[
  !str_detect(cue_by_phase_1chan_fatigue_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

cue_by_phase_1chan_fatigue_summary <- cue_by_phase_1chan_fatigue_fit$summary(variables = cue_by_phase_1chan_fatigue_fit_relevant_parameters)

# Get the actual posterior samples


cue_by_phase_1chan_fatigue_samples <- posterior::as_draws_df(
  cue_by_phase_1chan_fatigue_fit$draws(variables = cue_by_phase_1chan_fatigue_fit_relevant_parameters))

# rm(mcmc_long, cue_by_phase_1chan_fatigue_samples_df)
# cue_by_phase_1chan_fatigue_samples_df <- as.data.frame(cue_by_phase_1chan_fatigue_samples)


# mcmc_long <- cue_by_phase_1chan_fatigue_samples %>%
#   select(.draw, .chain, .iteration, cue_by_phase_1chan_fatigue_fit_relevant_parameters) %>%
#   tidyr::pivot_longer(cols = all_of(cue_by_phase_1chan_fatigue_fit_relevant_parameters),
#                       names_to = "parameter", values_to = "value")


cue_by_phase_1chan_fatigue_samples %>% 
  ggplot() +
  geom_line(aes(x = .iteration, y = sigma_average, color = as.factor(.chain), linewidth = as.factor(.chain))) +
  theme_classic()

cue_by_phase_1chan_fatigue_samples %>% 
  ggplot() +
  geom_line(aes(x = .iteration, y = intercept_average, color = as.factor(.chain), linewidth = as.factor(.chain))) +
  theme_classic()


library(bayesplot)
library(RColorBrewer)

# Use colorRampPalette to generate 16 colors from the "Set3" palette
custom_colors <- colorRampPalette(brewer.pal(12, "Set3"))(16)

bayesplot::color_scheme_set(custom_colors)

bayesplot::mcmc_trace(cue_by_phase_1chan_fatigue_samples[,,"sigma_average"]) 

bayesplot::mcmc_trace(cue_by_phase_1chan_fatigue_samples[,,"intercept_average"])

bayesplot::mcmc_trace(cue_by_phase_1chan_fatigue_samples[,,"fatigue_average"])

bayesplot::mcmc_trace(cue_by_phase_1chan_fatigue_samples[,,66:70])



raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(trial) %>% 
  summarise_all(mean) %>% 
  ggplot() +
  geom_line(aes(x = trial, y = amplitude))

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(trial, participant) %>% 
  summarise_all(mean) %>% 
  ggplot() +
  geom_line(aes(x = trial, y = amplitude, color = participant))

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>% 
  group_by(cue, phase, paired) %>% 
  summarise(mean_z_amp = mean(z_scored_amplitude),
            se_amp = plotrix::std.error(z_scored_amplitude)) %>% 
  ggplot() +
  geom_pointrange(aes(color = paired, x = cue, y = mean_z_amp,
                      ymax = mean_z_amp + se_amp, ymin = mean_z_amp - se_amp)) +
  facet_grid(.~ phase)
  theme_bw()

raw_fft_df %>%
  filter(channel %in% c("Oz")) %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude))) %>%
  ggplot() +
  geom_line(aes(x = trial, y = z_scored_amplitude, color = participant), linewidth = .2) +
  theme_classic()

## Z_score within participant because of excess variation

raw_fft_df_Oz<- raw_fft_df_Oz %>% 
  group_by(participant) %>%
  mutate(z_scored_amplitude = as.vector(scale(amplitude)))

raw_stan_list$amplitude <- raw_fft_df_Oz %>% 
  pull(z_scored_amplitude)

## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_fatigue_single_channel_normalized_participant.stan")

#force_recompile = T is sometimes helpful
cue_by_phase_1chan_fatigue_norm <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
# cue_by_phase_1chan_fatigue$print()

# Clear previous chains
list.files(pattern = "amp_ML_fatigue_single_channel_normalized_participant",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

cue_by_phase_1chan_fatigue_norm_fit <- cue_by_phase_1chan_fatigue_norm$sample(data = raw_stan_list,
                                                                              refresh = 500,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
cue_by_phase_1chan_fatigue_norm_fit_meta_data <- cue_by_phase_1chan_fatigue_norm_fit$metadata()

cue_by_phase_1chan_fatigue_norm_fit_relevant_parameters <- cue_by_phase_1chan_fatigue_norm_fit_meta_data$model_params[
  !str_detect(cue_by_phase_1chan_fatigue_norm_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

cue_by_phase_1chan_fatigue_norm_fit_summary <- 
  cue_by_phase_1chan_fatigue_norm_fit$summary(variables = cue_by_phase_1chan_fatigue_norm_fit_relevant_parameters)


cue_by_phase_1chan_fatigue_norm_fit_loo <- cue_by_phase_1chan_fatigue_norm_fit$loo()



# Get the actual posterior samples


cue_by_phase_1chan_fatigue_norm_samples <- 
  cue_by_phase_1chan_fatigue_norm_fit$draws(variables = cue_by_phase_1chan_fatigue_norm_fit_relevant_parameters,
  format = "df")

cue_by_phase_1chan_fatigue_norm_long <- cue_by_phase_1chan_fatigue_norm_samples %>% 
  pivot_longer(cols = cue_by_phase_1chan_fatigue_norm_fit_relevant_parameters)


## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_fatigue_time_single_channel_normalized_participant.stan")

#force_recompile = T is sometimes helpful
cue_by_phase_1chan_fatigue_norm_time <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
# cue_by_phase_1chan_fatigue$print()

# Clear previous chains
list.files(pattern = "amp_ML_fatigue_time_single_channel_normalized_participant",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

cue_by_phase_1chan_fatigue_norm_time_fit <- cue_by_phase_1chan_fatigue_norm_time$sample(data = raw_stan_list,
                                                                              refresh = 500,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
cue_by_phase_1chan_fatigue_norm_time_fit_meta_data <- cue_by_phase_1chan_fatigue_norm_time_fit$metadata()

cue_by_phase_1chan_fatigue_norm_time_fit_relevant_parameters <- cue_by_phase_1chan_fatigue_norm_time_fit_meta_data$model_params[
  !str_detect(cue_by_phase_1chan_fatigue_norm_time_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

cue_by_phase_1chan_fatigue_norm_time_fit_summary <- 
  cue_by_phase_1chan_fatigue_norm_time_fit$summary(variables = cue_by_phase_1chan_fatigue_norm_time_fit_relevant_parameters)


cue_by_phase_1chan_fatigue_norm_time_fit_loo <- cue_by_phase_1chan_fatigue_norm_time_fit$loo()

loo::loo_compare(cue_by_phase_1chan_fatigue_norm_time_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_fit_loo)

cue_by_phase_1chan_fatigue_norm_time_fit_samples <- 
  cue_by_phase_1chan_fatigue_norm_time_fit$draws(variables = cue_by_phase_1chan_fatigue_norm_time_fit_relevant_parameters,
  format = "df")

cue_by_phase_1chan_fatigue_norm_time_fit_long <- cue_by_phase_1chan_fatigue_norm_time_fit_samples %>% 
  pivot_longer(cols = cue_by_phase_1chan_fatigue_norm_time_fit_relevant_parameters)





loo::loo_model_weights(list(cue_by_phase_1chan_fatigue_norm_time_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_RW_fit_loo))

loo::loo_model_weights(list(cue_by_phase_1chan_fatigue_norm_fit_loo,
                            cue_by_phase_1chan_fatigue_norm_RW_fit_loo))

cue_by_phase_1chan_fatigue_norm_RW_fit_samples <- 
  cue_by_phase_1chan_fatigue_norm_RW_fit$draws(variables = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters,
  format = "df")

cue_by_phase_1chan_fatigue_norm_RW_fit_long <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>% 
  pivot_longer(cols = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters)



## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_RW_norm_single_channel.stan")

#force_recompile = T is sometimes helpful
cue_by_phase_1chan_fatigue_norm_RW <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
# cue_by_phase_1chan_fatigue$print()

# Clear previous chains
list.files(pattern = "amp_ML_RW_norm_single_channel",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

cue_by_phase_1chan_fatigue_norm_RW_fit <- cue_by_phase_1chan_fatigue_norm_RW$sample(data = raw_stan_list,
                                                                              refresh = 50,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
cue_by_phase_1chan_fatigue_norm_RW_fit_meta_data <- cue_by_phase_1chan_fatigue_norm_RW_fit$metadata()

cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters <- cue_by_phase_1chan_fatigue_norm_RW_fit_meta_data$model_params[
  !str_detect(cue_by_phase_1chan_fatigue_norm_RW_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

cue_by_phase_1chan_fatigue_norm_RW_fit_summary <- 
  cue_by_phase_1chan_fatigue_norm_RW_fit$summary(variables = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters)


cue_by_phase_1chan_fatigue_norm_RW_fit_loo <- cue_by_phase_1chan_fatigue_norm_RW_fit$loo()

loo::loo_compare(cue_by_phase_1chan_fatigue_norm_time_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_RW_fit_loo)

loo::loo_model_weights(list(cue_by_phase_1chan_fatigue_norm_time_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_fit_loo,
                 cue_by_phase_1chan_fatigue_norm_RW_fit_loo))

loo::loo_model_weights(list(cue_by_phase_1chan_fatigue_norm_fit_loo,
                            cue_by_phase_1chan_fatigue_norm_RW_fit_loo))

cue_by_phase_1chan_fatigue_norm_RW_fit_samples <- 
  cue_by_phase_1chan_fatigue_norm_RW_fit$draws(variables = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters,
  format = "df")

cue_by_phase_1chan_fatigue_norm_RW_fit_long <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>% 
  pivot_longer(cols = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters)



## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_scale_by_par_norm_single_channel.stan")

#force_recompile = T is sometimes helpful
amp_ML_scale_by_par_norm_single_channel <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
# cue_by_phase_1chan_fatigue$print()

# Clear previous chains
list.files(pattern = "amp_ML_scale_by_par_norm_single_channel",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

amp_ML_scale_by_par_norm_single_channel_fit <- amp_ML_scale_by_par_norm_single_channel$sample(data = raw_stan_list,
                                                                              refresh = 50,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
amp_ML_scale_by_par_norm_single_channel_fit_meta_data <- amp_ML_scale_by_par_norm_single_channel_fit$metadata()

amp_ML_scale_by_par_norm_single_channel_fit_relevant_parameters <- amp_ML_scale_by_par_norm_single_channel_fit_meta_data$model_params[
  !str_detect(amp_ML_scale_by_par_norm_single_channel_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

amp_ML_scale_by_par_norm_single_channel_fit_summary <- 
  amp_ML_scale_by_par_norm_single_channel_fit$summary(variables = amp_ML_scale_by_par_norm_single_channel_fit_relevant_parameters)


amp_ML_scale_by_par_norm_single_channel_fit_loo <- amp_ML_scale_by_par_norm_single_channel_fit$loo()

loo::loo_compare(cue_by_phase_1chan_fatigue_norm_RW_fit_loo,
                 amp_ML_scale_by_par_norm_single_channel_fit_loo)

## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_RW_norm_single_channel_larger_LR.stan")

#force_recompile = T is sometimes helpful
amp_ML_RW_norm_single_channel_larger_LR <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
amp_ML_RW_norm_single_channel_larger_LR$print()

# Clear previous chains
list.files(pattern = "amp_ML_RW_norm_single_channel_larger_LR",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

amp_ML_RW_norm_single_channel_larger_LR_fit <- amp_ML_RW_norm_single_channel_larger_LR$sample(data = raw_stan_list,
                                                                              refresh = 50,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
amp_ML_RW_norm_single_channel_larger_LR_fit_meta_data <- amp_ML_RW_norm_single_channel_larger_LR_fit$metadata()

amp_ML_RW_norm_single_channel_larger_LR_fit_relevant_parameters <- amp_ML_RW_norm_single_channel_larger_LR_fit_meta_data$model_params[
  !str_detect(amp_ML_scale_by_par_norm_single_channel_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

amp_ML_RW_norm_single_channel_larger_LR_fit_summary <- 
  amp_ML_RW_norm_single_channel_larger_LR_fit$summary(variables = amp_ML_RW_norm_single_channel_larger_LR_fit_relevant_parameters)


amp_ML_RW_norm_single_channel_larger_LR_fit_loo <- amp_ML_RW_norm_single_channel_larger_LR_fit$loo()

loo::loo_compare(cue_by_phase_1chan_fatigue_norm_RW_fit_loo,
                 amp_ML_scale_by_par_norm_single_channel_fit_loo,
                 amp_ML_RW_norm_single_channel_larger_LR_fit_loo)


## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_RW_norm_single_channel_scale_phase.stan")

#force_recompile = T is sometimes helpful
amp_ML_RW_norm_single_channel_scale_phase <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
amp_ML_RW_norm_single_channel_scale_phase$print()

# Clear previous chains
list.files(pattern = "amp_ML_RW_norm_single_channel_scale_phase",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

amp_ML_RW_norm_single_channel_scale_phase_fit <- amp_ML_RW_norm_single_channel_scale_phase$sample(data = raw_stan_list,
                                                                              refresh = 50,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
amp_ML_RW_norm_single_channel_scale_phase_fit_meta_data <- amp_ML_RW_norm_single_channel_scale_phase_fit$metadata()

amp_ML_RW_norm_single_channel_scale_phase_fit_relevant_parameters <- amp_ML_RW_norm_single_channel_scale_phase_fit_meta_data$model_params[
  !str_detect(amp_ML_RW_norm_single_channel_scale_phase_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

amp_ML_RW_norm_single_channel_scale_phase_fit_summary <- 
  amp_ML_RW_norm_single_channel_scale_phase_fit$summary(variables = amp_ML_RW_norm_single_channel_scale_phase_fit_relevant_parameters)


amp_ML_RW_norm_single_channel_scale_phase_fit_loo <- amp_ML_RW_norm_single_channel_scale_phase_fit$loo()

loo::loo_compare(cue_by_phase_1chan_fatigue_norm_RW_fit_loo,
                 amp_ML_scale_by_par_norm_single_channel_fit_loo,
                 amp_ML_RW_norm_single_channel_larger_LR_fit_loo,
                 amp_ML_RW_norm_single_channel_scale_phase_fit_loo)

## cue by phase with fatigue one channel ####
model_path <- paste0(git_repository, "/stan_models/amp_ML_norm_single_channel_par_cue_phase_mat.stan")

#force_recompile = T is sometimes helpful
amp_ML_norm_single_channel_par_cue_phase_mat <- cmdstanr::cmdstan_model(model_path, 
                                                           force_recompile = T)

#Model source code
amp_ML_norm_single_channel_par_cue_phase_mat$print()

# Clear previous chains
list.files(pattern = "amp_ML_norm_single_channel_par_cue_phase_mat",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

amp_ML_norm_single_channel_par_cue_phase_mat_fit <- amp_ML_norm_single_channel_par_cue_phase_mat$sample(data = raw_stan_list,
                                                                              refresh = 50,
                                                                              seed = 2,
                                                                              iter_warmup = warmup_samples_per_chain, 
                                                                              iter_sampling = posterior_samples_per_chain, 
                                                                              save_warmup = F, 
                                                                              show_messages = T,
                                                                              output_dir = where_to_save_chains,
                                                                              chains = number_of_chains,
                                                                              parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
amp_ML_norm_single_channel_par_cue_phase_mat_fit_meta_data <- amp_ML_norm_single_channel_par_cue_phase_mat_fit$metadata()

amp_ML_norm_single_channel_par_cue_phase_mat_fit_relevant_parameters <- amp_ML_norm_single_channel_par_cue_phase_mat_fit_meta_data$model_params[
  !str_detect(amp_ML_norm_single_channel_par_cue_phase_mat_fit_meta_data$model_params, "log_lik|mu_pred|amplitude")]

amp_ML_norm_single_channel_par_cue_phase_mat_fit_summary <- 
  amp_ML_norm_single_channel_par_cue_phase_mat_fit$summary(variables = amp_ML_norm_single_channel_par_cue_phase_mat_fit_relevant_parameters)


amp_ML_norm_single_channel_par_cue_phase_mat_fit_loo <- amp_ML_norm_single_channel_par_cue_phase_mat_fit$loo()

loo::loo_compare(cue_by_phase_1chan_fatigue_norm_RW_fit_loo,
                 amp_ML_scale_by_par_norm_single_channel_fit_loo,
                 amp_ML_RW_norm_single_channel_larger_LR_fit_loo,
                 amp_ML_RW_norm_single_channel_scale_phase_fit_loo,
                 amp_ML_norm_single_channel_par_cue_phase_mat_fit_loo)

cue_by_phase_1chan_fatigue_norm_RW_fit_samples <- 
  cue_by_phase_1chan_fatigue_norm_RW_fit$draws(variables = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters,
  format = "df")

cue_by_phase_1chan_fatigue_norm_RW_fit_long <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>% 
  pivot_longer(cols = cue_by_phase_1chan_fatigue_norm_RW_fit_relevant_parameters)



library(ggridges)
library(patchwork)

raw_fft_df_Oz %>%
  filter(participant == "115") %>% 
  ggplot() +
  geom_line(aes(x = trial, y= z_scored_amplitude))

cue_by_phase_1chan_fatigue_norm_long %>% 
  filter(grepl(pattern = "sigma", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # scale_x_continuous(limits = c(0, 3)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_long %>% 
  filter(grepl(pattern = "sigma\\[", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # scale_x_continuous(limits = c(0, 2)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_long %>% 
  filter(grepl(pattern = "intercept\\[", .data$name)) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name),, fill = "black",
   alpha = .4) +
  scale_x_continuous(name = "Z-Scored ssVEP") +
  ggtitle("ssVEP Intercept at Trial 0") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank()) +

cue_by_phase_1chan_fatigue_norm_long %>% 
  filter(grepl(pattern = "fatigue\\[", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),, fill = "black",
   alpha = .4) +
  scale_x_continuous(name = "Z-Scored ssVEP") +
  ggtitle("Slope with Trials") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())

cue_by_phase_1chan_fatigue_norm_long %>% 
  filter(grepl(pattern = "bcue", .data$name)) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name),, fill = "black",
   alpha = .4) +
  scale_x_continuous(name = "Z-Scored ssVEP", breaks = seq(-.75, .75, by = .25)) +
  scale_y_discrete(labels = c("Hab, CS+",
                              "Hab, GS1",
                              "Hab, GS2",
                              "Hab, GS3",
                              "Acq, CS+",
                              "Acq, GS1",
                              "Acq, GS2",
                              "Acq, GS3",
                              "Ext, CS+",
                              "Ext, GS1",
                              "Ext, GS2",
                              "Ext, GS3")) +
  ggtitle("Cue by phase") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())

cue_by_phase_1chan_fatigue_norm_long %>% 
  filter(grepl(pattern = "bcue", .data$name)) %>% 
  group_by(name) %>% 
  summarize((sum(value - 0 < 0)/ 12000)*100)



cue_by_phase_1chan_fatigue_norm_time_fit_long %>% 
  filter(grepl(pattern = "sigma", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # scale_x_continuous(limits = c(0, 3)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_time_fit_long %>% 
  filter(grepl(pattern = "sigma\\[", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # scale_x_continuous(limits = c(0, 2)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_time_fit_long %>% 
  filter(grepl(pattern = "bcue", .data$name)) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name),, fill = "black",
   alpha = .4) +
  # scale_x_continuous(limits = c(-1.5, 1.5)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_time_fit_long %>% 
  filter(grepl(pattern = "bcue", .data$name)) %>% 
  group_by(name) %>% 
  summarize((sum(value - 0 < 0)/ 12000)*100)


cue_by_phase_1chan_fatigue_norm_time_fit_long %>% 
  filter(grepl(pattern = "intercept", .data$name)) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name),, fill = "black",
   alpha = .4) +
  # scale_x_continuous(limits = c(-1.5, 1.5)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_time_fit_long %>% 
  filter(grepl(pattern = "fatigue\\[", .data$name)) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name),, fill = "black",
   alpha = .4) +
  # scale_x_continuous(limits = c(-1.5, 1.5)) +
  theme_bw()



cue_by_phase_1chan_fatigue_norm_RW_fit_long  %>% 
  filter(grepl(pattern = "sigma", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # scale_x_continuous(limits = c(0, 3)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_RW_fit_long  %>% 
  filter(grepl(pattern = "sigma\\[", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # scale_x_continuous(limits = c(0, 2)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_RW_fit_long  %>% 
  filter(grepl(pattern = "intercept", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # coord_cartesian(xlim = c(0, .1)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_RW_fit_long  %>% 
  filter(grepl(pattern = "fatigue", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  # coord_cartesian(xlim = c(0, .1)) +
  theme_bw()

cue_by_phase_1chan_fatigue_norm_RW_fit_long  %>% 
  filter(grepl(pattern = "scaling", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  scale_y_discrete(labels = c("Average of Cues", "SD of Cues", "CS+", "GS1", "GS2", "GS3")) +
  scale_x_continuous(name = "Z-Scored ssVEP") +
  ggtitle("Cue Scaling with CS+ Associative Strength") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())

cue_by_phase_1chan_fatigue_norm_RW_fit_long  %>% 
  filter(grepl(pattern = "learning_rate", .data$name)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  scale_x_continuous(name = "Bound between 0 and 1") +
  coord_cartesian(xlim = c(0, .075)) +
  ggtitle("Learning Rate by Participant") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())




CSP_scaling <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>%
    pull(`scaling[1]`)

GS1_scaling <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>%
    pull(`scaling[2]`)

GS2_scaling <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>%
    pull(`scaling[3]`)

GS3_scaling <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>%
    pull(`scaling[4]`)

scaling_plot_df <- data.frame("participant" = double(),
                              "iteration" = double(),
                              "trial"= double(), 
                              "name" = character(), 
                              "value" = double())


start_trial_index <- 1
stop_trial_index <- 176

for(p in 1:4){
  
  for(i in start_trial_index:stop_trial_index){
    
    current_CSP_variable_name <- paste0("CSP_associative_strength",'[', i,']') 
    
    current_associative_strength <- cue_by_phase_1chan_fatigue_norm_RW_fit_samples %>%
    pull(current_CSP_variable_name)
    
    n_samples <- length(current_associative_strength)
    
    CSP_scaling_mod <- current_associative_strength * CSP_scaling
    GS1_scaling_mod <- current_associative_strength * GS1_scaling
    GS2_scaling_mod <- current_associative_strength * GS2_scaling
    GS3_scaling_mod <- current_associative_strength * GS3_scaling
    
    current_df_to_rbind <- rbind(data.frame(participant = p,
                                            iteration = 1:n_samples,
                                            trial = rep(i, n_samples),
                                            name = "CSP_scaling_mod",
                                            value = CSP_scaling_mod),
                                 data.frame(participant = p,
                                            iteration = 1:n_samples,
                                            trial = rep(i, n_samples),
                                            name = "GS1_scaling_mod",
                                            value = GS1_scaling_mod),
                                 data.frame(participant = p,
                                            iteration = 1:n_samples,
                                            trial = rep(i, n_samples),
                                            name = "GS2_scaling_mod",
                                            value = GS2_scaling_mod),
                                 data.frame(participant = p,
                                            iteration = 1:n_samples,
                                            trial = rep(i, n_samples),
                                            name = "GS3_scaling_mod",
                                            value = GS3_scaling_mod))
            
    scaling_plot_df <- rbind(scaling_plot_df, current_df_to_rbind)
            
  }
  start_trial_index <- start_trial_index + 176
  stop_trial_index <- stop_trial_index + 176
}

max_samples <- scaling_plot_df$iteration %>% max()
n_samples_to_use <- 200
iterations_to_use <- sample(x = 1:max_samples, 
                            size = n_samples_to_use)


scaling_plot_df %>% 
  filter(iteration %in% iterations_to_use) %>% 
  ggplot() +
  geom_line(aes(x = trial, 
                y = value, 
                color = name, 
                group = interaction(iteration, name)),
            alpha = .3,
            linewidth = .5) +
  scale_y_continuous(name = "Z-Scored ssVEP") +
  scale_x_continuous(name = "Trial", breaks = seq(0, 176*5, by = 50)) +
  scale_color_discrete(name = "Condition",
                       labels = c("CS+", "GS1", "GS2", "GS3")) +
  ggtitle(paste0(n_samples_to_use, " posterior draws for CS+ associative strength on ssVEP",
                 "\nFirst 4 participants")) + 
  theme_classic() +
  theme(text = element_text(size = 22)) +
  guides(color = guide_legend(override.aes = list(linewidth = 6)))
 
# save.image()
# getwd()
# load('.RData')
