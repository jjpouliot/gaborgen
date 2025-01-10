library(tidyverse)
library(cmdstanr)
library(R.matlab)
library(ggbeeswarm)
library(patchwork)
library(ggridges)

# Load data ####

parent_folder <- "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI"
git_repository <- "/home/andrewf/Repositories/gaborgen"

participants_without_complete_trial_list_in_dat_file <- c("110")

channel_information <- readMat(paste0(git_repository, "/chanlocs_ssv4att_MRI.mat")) 

single_trial_mat_folder <- paste0(parent_folder,"/single_trial_timeseries_FFTs_CSD") # old is CSD transformed, new is average reference

all_mat_filepaths <- list.files(path = single_trial_mat_folder,full.names = T)
all_mat_filepaths <- all_mat_filepaths[!grepl(x = all_mat_filepaths,
                                              pattern = paste0("participant", 
                                                               participants_without_complete_trial_list_in_dat_file))]

raw_fft_paths <- all_mat_filepaths[str_detect(all_mat_filepaths, 
                                              pattern = "rawChanFrequency")]

sliding_window_fft_paths <- all_mat_filepaths[str_detect(all_mat_filepaths, 
                                                         pattern = "slidingWindow_ChanFrequency")]

sliding_window_fft_paths <- sliding_window_fft_paths[str_detect(sliding_window_fft_paths, 
                                                                pattern = "600Hz.mat")]




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
  
  amplitude <- mat_data$rawAmp[,ssVEP_frequency_index]
  
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
  
  amplitude <- mat_data$slidingWindowAmp[,ssVEP_frequency_index]
  
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


raw_fft_df <- raw_fft_df %>%
  rename(amplitude_15Hz_fft = amplitude)



sliding_window_fft_df <- sliding_window_fft_df %>%
  rename(amplitude_15Hz_sliding_window_fft_upsampled_600Hz = amplitude)

fft_df <- merge(x = raw_fft_df, y = sliding_window_fft_df,
                by = c("participant", "channel", "trial", "phase", "cue", "paired"))


#csPerm1 is when 15 degrees is CSP (most vertical), 35, 55, 75. csPerm2 is opposite order
csPerm1 <- c(101,110,111,114,116,117,119,125,129,131,133,134,135,137,138,139,140,141,143)
csPerm2 <- c(102,112,113,115,118,120,121,122,123,124,126,127,128,130,132,136,141)


fft_df <- fft_df %>% 
  mutate(csPerm = if_else(participant %in% as.character(csPerm1), 1,2))


Day1_ratings_paths <- list.files(path = parent_folder, 
                                 pattern = "Day1.*ratings.dat$",
                                 recursive = T,
                                 full.names = T)



for(i in 1:length(Day1_ratings_paths)) {
  current_ratings <- read.csv(Day1_ratings_paths[i])
  if(i == 1){
    ratings_df <- current_ratings
  } else {
    ratings_df <- rbind.data.frame(ratings_df, current_ratings)
  }
}

ratings_df$partInd <- as.character(ratings_df$partInd)

# experiment mislabels a rating
names(ratings_df)[9] <- "ar_gs2"

ratings_df <-ratings_df %>% 
  select(!contains("Dur")) %>% 
  pivot_longer(
    cols = c(val_csp, val_gs1, val_gs2, val_gs3,
             ar_csp, ar_gs1, ar_gs2, ar_gs3,
             exp_csp, exp_gs1, exp_gs2, exp_gs3),
    names_to = c(".value", "condition"),
    names_sep = "_")

ratings_df <-ratings_df %>% 
  mutate(condition =  case_when(condition == "csp" ~ "CSP",
                                 condition == "gs1" ~ "GS1",
                                 condition == "gs2" ~ "GS2",
                                 condition == "gs3" ~ "GS3"))

fft_df <- fft_df %>% 
  mutate(.after = phase, 
         block = case_when(trial <= 32 ~ 1,
                           trial > 32 & trial <= (32+48) ~ 2,
                           trial > (32+48) & trial <= (32+48+48) ~ 3,
                           trial > (32+48+48) & trial <= (32+48+48+48) ~ 4))

fft_df <- merge(x = fft_df, y = ratings_df, 
              by.x = c("participant","block", "cue"), by.y = c("partInd", "ratInd", "condition"),
              all.x = T)

# Extract unique participant IDs
unique_participants <- unique(fft_df$participant)
# Ask before showing new page (new plot)
devAskNewPage(TRUE)

for (p in unique_participants) {
  # Create the plot
  # p_plot <-   fft_df %>% 
  #   filter(channel %in% c("Oz"), participant == p) %>%
  #   ggplot() +
  #   geom_quasirandom(aes(x = cue, y = amplitude_15Hz_fft, color = cue))+ 
  #   facet_wrap(~phase)
  p_plot <-   fft_df %>%
    filter(channel %in% c("Oz"), participant == p) %>%
    group_by(participant, channel, cue, phase) %>%
    summarise(mean_amp = mean(amplitude_15Hz_fft),
              se_amp = plotrix::std.error(amplitude_15Hz_fft)) %>%
    ggplot() +
    geom_line(aes(x = cue, y = mean_amp, group = channel, color = channel)) +
    geom_errorbar(aes(x = cue,
                      ymin = mean_amp - se_amp,
                      ymax = mean_amp + se_amp,
                      group = participant, color = channel),width = 0) +
    facet_wrap(~phase)
  
  # Print the plot
  print(p_plot)
  
  fft_df %>% 
    filter(channel %in% c("Oz"), participant == p) %>%
    group_by(participant,channel, csPerm,phase,cue) %>% 
    summarise(mean_amp = mean(amplitude_15Hz_fft),
              se_amp = plotrix::std.error(amplitude_15Hz_fft),
              mean_ar =mean(ar),
              mean_val =mean(val),
              mean_exp =mean(exp),
              n()) %>% 
    print(n = 999)
  
  user_input <- readline(prompt = "Press [enter] for the next plot or 'q' to quit: ")
  
  # Check if the user wants to quit early
  if (tolower(user_input) == "q") {
    message("Exiting early as requested.")
    break
  }
}
devAskNewPage(FALSE)

Oz_fft_df <- fft_df %>% 
  filter(channel == "Oz") %>% 
  filter(!participant %in% c("123","136"))

## temp data visualization ####



Oz_fft_df %>% 
  group_by(csPerm, block, paired, cue) %>% 
  summarise(n()) %>% 
  print(n = 99)


Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = scale(amplitude_15Hz_fft),
         zcsp_ar = scale(ar_csp)) %>% 
  filter(block == 4,cue == "CSP") %>% 
  ggplot() +
  geom_quasirandom(aes(x = zcsp_ar, y = zamp)) +
  geom_smooth(aes(x = zcsp_ar, y = zamp),method = "lm")

Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = scale(amplitude_15Hz_fft)) %>%
  ggplot() +
  geom_smooth(aes(x = trial, y = zamp, color = cue))


Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
  group_by(phase, cue, paired) %>% 
  summarise(mean_zamp = mean(zamp),
            se_zamp = plotrix::std.error(zamp)) %>% 
  ggplot() +
  geom_pointrange(aes(x = cue, y = mean_zamp, ymax = mean_zamp + se_zamp,
                      ymin = mean_zamp - se_zamp, color = paired)) +
  facet_wrap(~ phase)


Oz_fft_df %>% 
  ggplot() +
  geom_smooth(aes(x = trial, y = amplitude_15Hz_fft, color = cue)) +
  facet_wrap(~csPerm)


Oz_fft_df %>% 
  group_by(phase, cue, paired) %>% 
  summarise(mean_amp = mean(amplitude_15Hz_fft),
            se_amp = plotrix::std.error(amplitude_15Hz_fft)) %>% 
  ggplot() +
  geom_pointrange(aes(x = cue, y = mean_amp, ymax = mean_amp + se_amp,
                      ymin = mean_amp - se_amp, color = paired)) +
  facet_wrap(~ phase)

Oz_fft_df %>% 
  # group_by(participant) %>% 
  # mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
  group_by(block, cue, paired, csPerm) %>% 
  summarise(mean_amp = mean(amplitude_15Hz_fft),
            se_amp = plotrix::std.error(amplitude_15Hz_fft)) %>% 
  ggplot() +
  geom_pointrange(aes(x = cue, y = mean_amp, ymax = mean_amp + se_amp,
                      ymin = mean_amp - se_amp, color = paired)) +
  facet_wrap(~ csPerm * block, ncol = 4)

Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
  group_by(block, cue, paired, csPerm) %>% 
  summarise(mean_zamp = mean(zamp),
            se_zamp = plotrix::std.error(zamp)) %>% 
  ggplot() +
  geom_pointrange(aes(x = cue, y = mean_zamp, ymax = mean_zamp + se_zamp,
                      ymin = mean_zamp - se_zamp, color = paired)) +
  facet_wrap(~ csPerm * block, ncol = 4)

Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
  # group_by(phase, cue, paired, csPerm) %>% 
  # summarise(mean_zamp = mean(zamp),
  #           se_zamp = plotrix::std.error(zamp)) %>% 
  ggplot() +
  geom_quasirandom(aes(x = cue, y = zamp, color = paired)) +
  # geom_pointrange(aes(x = cue, y = mean_zamp, ymax = mean_zamp + se_zamp,
  #                     ymin = mean_zamp - se_zamp, color = paired)) +
  facet_wrap(~ csPerm * phase)


# Create stan list ####
gaborgen_stan_list <- list()

gaborgen_stan_list$n <- nrow(Oz_fft_df)

gaborgen_stan_list$trial <- Oz_fft_df$trial %>% as.integer()
gaborgen_stan_list$n_trials <- 176

gaborgen_stan_list$participant <- Oz_fft_df$participant %>% as.factor() %>% as.integer()
gaborgen_stan_list$n_participants <- gaborgen_stan_list$participant %>% unique() %>% length()

gaborgen_stan_list$cue <- Oz_fft_df$cue %>% as.integer()
gaborgen_stan_list$n_cues <- gaborgen_stan_list$cue %>% unique() %>% length()

gaborgen_stan_list$phase <- Oz_fft_df$phase %>% as.integer()
gaborgen_stan_list$n_phases <- gaborgen_stan_list$phase %>% unique() %>% length()

gaborgen_stan_list$paired <- ifelse(Oz_fft_df$paired == "shock", 1, 0) %>% as.integer()

gaborgen_stan_list$amplitude <- Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(amplitude = scale(amplitude_15Hz_fft)) %>% 
  pull(amplitude) %>% 
  as.vector()




# Sampling settings ####
number_of_chains <- 4
warmup_samples_per_chain <- 2000
posterior_samples_per_chain <- 2000
where_to_save_chains <- paste0(parent_folder,"/stan_chains")

# Fit Models ####

## Model 001 - Descriptive model: Fatigue slope, cue by phase ####


model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model001.stan'

#force_recompile = T is sometimes helpful
model001 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
model001$print()

# Clear previous chains
list.files(pattern = "Model001",
           path = where_to_save_chains, 
           full.names = T) %>% 
  file.remove()

model001_fit <- model001$sample(data = gaborgen_stan_list,
                                refresh = 1000,
                                seed = 2,
                                iter_warmup = warmup_samples_per_chain, 
                                iter_sampling = posterior_samples_per_chain, 
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                chains = number_of_chains,
                                parallel_chains = number_of_chains)

# Check summary to see if everything was estimated as expected
model001_fit_meta_data <- model001_fit$metadata()

model001_fit_relevant_parameters <- model001_fit_meta_data$model_params[
  !str_detect(model001_fit_meta_data$model_params, "log_lik|mu_pred")]

model001_fit_summary <- 
  model001_fit$summary(variables = model001_fit_relevant_parameters)

model001_fit_loo <- model001_fit$loo()

model001_df <- posterior::as_draws_df(
  model001_fit$draws(variables = model001_fit_relevant_parameters))

model001_df %>%
  select(starts_with("sigma")) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model001_df %>%
  select(starts_with("intercept")) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model001_df %>%
  select(starts_with("fatigue")) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))
  
model001_df %>%
  select(starts_with("bcue")) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name), alpha = .5) +
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
  
## Model 002 - Rescola-Wagner model ####



## Model 003 - Ratings model: Slope per cue + slope per pairing ####
