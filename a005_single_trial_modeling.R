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

single_trial_mat_folder <- paste0(parent_folder,"/single_trial_timeseries_FFTs_CSD") 
# single_trial_mat_folder <- paste0(parent_folder,"/single_trial_timeseries_FFTs_CSD") #0 to 2000

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
    facet_wrap(~phase) +
  ggtitle(paste0(p))
  # p_plot <-   fft_df %>%
  #   filter(channel %in% c("Oz"), participant == p) %>%
  #   # group_by(participant) %>%
  #   mutate(zamp = scale(amplitude_15Hz_fft) %>% as.vector()) %>%
  #   ggplot() +
  #   geom_point(aes(x = trial, y = amplitude_15Hz_fft, color = cue)) +
  #   scale_x_continuous(limits = c(0,176)) +
  #   ggtitle(paste0(p))
  # geom_errorbar(aes(x = cue,
  #                   ymin = mean_amp - se_amp,
  #                   ymax = mean_amp + se_amp,
  #                   group = participant, color = channel),width = 0) +
  # facet_wrap(~phase)

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


## temp data visualization ####

# Oz_fft_df %>% 
#   select(csPerm,stan_par_id) %>% 
#   group_by(stan_par_id,csPerm) %>% 
#   reframe(n()) %>% 
#   print(n = 9999)
# 
# 
# Oz_fft_df %>%
#   filter(stan_par_id %in% c(1)) %>%
#   # filter(stan_par_id %in% c(5)) %>%
#   # filter(participant %in% c("118")) %>% 
#   ggplot() +
#   geom_beeswarm(aes(x = cue, y = amplitude_15Hz_fft),width = .1) +
#   # geom_quasirandom(aes(x = cue, y = amplitude_15Hz_fft),width = .1) +
#   facet_wrap(~block, nrow = 1)
# 
# (Oz_fft_df %>% 
#     group_by(participant) %>% 
#     mutate(zamp = scale(amplitude_15Hz_fft) %>% as.vector()) %>% 
#     ungroup() %>% 
#     filter(stan_par_id %in% c(1)) %>%
#     # filter(stan_par_id %in% c(5)) %>%
#     # filter(participant %in% c("118")) %>% 
#     ggplot() +
#     # geom_line(aes(x = trial, y = zamp)) +
#     geom_point(aes(x = trial, y = zamp, color = cue)) +
#     geom_smooth(aes(x = trial, y = zamp, color = cue), se = 0) +
#     # geom_quasirandom(aes(x = cue, y = zamp),width = .1) +
#     scale_color_manual(values = cue_color) +
#     theme_classic()) / 
#   (Oz_fft_df %>% 
#      group_by(participant) %>% 
#      mutate(zamp = scale(amplitude_15Hz_fft) %>% as.vector()) %>% 
#      ungroup() %>% 
#      filter(stan_par_id %in% c(5)) %>%
#      # filter(stan_par_id %in% c(5)) %>%
#      # filter(participant %in% c("118")) %>% 
#      ggplot() +
#      # geom_line(aes(x = trial, y = zamp)) +
#      geom_point(aes(x = trial, y = zamp, color = cue)) +
#      geom_smooth(aes(x = trial, y = zamp, color = cue), se = 0) +
#      # geom_quasirandom(aes(x = cue, y = zamp),width = .1) +
#      scale_color_manual(values = cue_color) +
#      theme_classic()) /
#   (Oz_fft_df %>% 
#      group_by(participant) %>% 
#      mutate(zamp = scale(amplitude_15Hz_fft) %>% as.vector()) %>% 
#      ungroup() %>% 
#      filter(stan_par_id %in% c(21)) %>%
#      # filter(stan_par_id %in% c(5)) %>%
#      # filter(participant %in% c("118")) %>% 
#      ggplot() +
#      # geom_line(aes(x = trial, y = zamp)) +
#      geom_point(aes(x = trial, y = zamp, color = cue)) +
#      geom_smooth(aes(x = trial, y = zamp, color = cue), se = 0) +
#      # geom_quasirandom(aes(x = cue, y = zamp),width = .1) +
#      scale_color_manual(values = cue_color) +
#      scale_x_continuous(breaks = seq(0,180,by = 20)) +
#      theme_classic()) 

# 
# Oz_fft_df %>% 
#   group_by(csPerm, block, paired, cue) %>% 
#   summarise(n()) %>% 
#   print(n = 99)
# 
# 
# Oz_fft_df %>% 
#   group_by(participant) %>% 
#   mutate(zamp = scale(amplitude_15Hz_fft),
#          zcsp_ar = scale(ar_csp)) %>% 
#   filter(block == 4,cue == "CSP") %>% 
#   ggplot() +
#   geom_quasirandom(aes(x = zcsp_ar, y = zamp)) +
#   geom_smooth(aes(x = zcsp_ar, y = zamp),method = "lm")
# 
# Oz_fft_df %>% 
#   group_by(participant) %>% 
#   mutate(zamp = scale(amplitude_15Hz_fft)) %>%
#   ggplot() +
#   geom_smooth(aes(x = trial, y = zamp, color = cue))
# 
# 
# Oz_fft_df %>% 
#   group_by(participant) %>% 
#   mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
#   group_by(phase, cue, paired) %>% 
#   summarise(mean_zamp = mean(zamp),
#             se_zamp = plotrix::std.error(zamp)) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = cue, y = mean_zamp, ymax = mean_zamp + se_zamp,
#                       ymin = mean_zamp - se_zamp, color = paired)) +
#   facet_wrap(~ phase)
# 
# 
# Oz_fft_df %>% 
#   ggplot() +
#   geom_smooth(aes(x = trial, y = amplitude_15Hz_fft, color = cue)) +
#   facet_wrap(~csPerm)
# 
# 
# Oz_fft_df %>% 
#   group_by(phase, cue, paired) %>% 
#   summarise(mean_amp = mean(amplitude_15Hz_fft),
#             se_amp = plotrix::std.error(amplitude_15Hz_fft)) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = cue, y = mean_amp, ymax = mean_amp + se_amp,
#                       ymin = mean_amp - se_amp, color = paired)) +
#   facet_wrap(~ phase)
# 
# Oz_fft_df %>% 
#   # group_by(participant) %>% 
#   # mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
#   group_by(block, cue, paired, csPerm) %>% 
#   summarise(mean_amp = mean(amplitude_15Hz_fft),
#             se_amp = plotrix::std.error(amplitude_15Hz_fft)) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = cue, y = mean_amp, ymax = mean_amp + se_amp,
#                       ymin = mean_amp - se_amp, color = paired)) +
#   facet_wrap(~ csPerm * block, ncol = 4)
# 
# Oz_fft_df %>% 
#   group_by(participant) %>% 
#   mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
#   group_by(block, cue, paired, csPerm) %>% 
#   summarise(mean_zamp = mean(zamp),
#             se_zamp = plotrix::std.error(zamp)) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = cue, y = mean_zamp, ymax = mean_zamp + se_zamp,
#                       ymin = mean_zamp - se_zamp, color = paired)) +
#   facet_wrap(~ csPerm * block, ncol = 4)
# 
# Oz_fft_df %>% 
#   group_by(participant) %>% 
#   mutate(zamp = scale(amplitude_15Hz_fft)) %>% 
#   # group_by(phase, cue, paired, csPerm) %>% 
#   # summarise(mean_zamp = mean(zamp),
#   #           se_zamp = plotrix::std.error(zamp)) %>% 
#   ggplot() +
#   geom_quasirandom(aes(x = cue, y = zamp, color = paired)) +
#   # geom_pointrange(aes(x = cue, y = mean_zamp, ymax = mean_zamp + se_zamp,
#   #                     ymin = mean_zamp - se_zamp, color = paired)) +
#   facet_wrap(~ csPerm * phase)


# Create stan list ####
#111 horrible quality
#112 horrible quality
#117 couldn't see
#123 weird at other sensor, but good at Oz, fmri not recording for first 10 mins according to master sheet?
#134 fell asleep
#136 missing too much
Oz_fft_df <- fft_df %>% 
  filter(channel == "Oz") %>% 
  filter(!participant %in% c("111", "112","117","123","134","136")) 

# need to sort df, then find trials missing
Oz_fft_df <- Oz_fft_df %>% 
  arrange(participant,trial)

learned_relationship_at_end_of_study <- c("101",
                                          "102",
                                          "103",
                                          "107",
                                          "111", #kind of
                                          "113",
                                          "114",
                                          "115",
                                          "116",
                                          "117",
                                          "118",
                                          "119",
                                          "120",
                                          "121",
                                          "122",
                                          "123",
                                          "124",
                                          "127",
                                          "129",
                                          "131",
                                          "132",
                                          "135",
                                          "136",
                                          "137",
                                          "138",
                                          "139",
                                          "140",
                                          "142",# thought a shock came with wrong pattern
                                          "144",
                                          "145")

Oz_fft_df <- Oz_fft_df %>% 
  mutate(learned_at_end = if_else(participant %in% learned_relationship_at_end_of_study,T, F))
Oz_fft_df <- Oz_fft_df %>% 
  mutate(stan_par_id = participant %>% 
           as.factor() %>% 
           as.integer())

gaborgen_stan_list <- list()


gaborgen_stan_list$participant <- Oz_fft_df$participant %>% 
                                    unique() %>% 
                                    as.factor() %>% 
                                    as.integer() %>%
                                    rep(each = 176)

gaborgen_stan_list$n_participants <- gaborgen_stan_list$participant %>% unique() %>% length()

gaborgen_stan_list$n_trials <- 176
# gaborgen_stan_list$trial <- rep(1:176, gaborgen_stan_list$n_participants)

gaborgen_stan_list$n <- (gaborgen_stan_list$n_trials * 
                           gaborgen_stan_list$n_participants)

gaborgen_stan_list$n_observations <- nrow(Oz_fft_df)
gaborgen_stan_list$n_missing <- gaborgen_stan_list$n - 
                                   gaborgen_stan_list$n_observations




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
participant_by_n <- Oz_fft_df$participant %>% unique() %>% rep(each = gaborgen_stan_list$n_trials)
observed_trials <- Oz_fft_df$trial

for (i in 1:(gaborgen_stan_list$n)) {
  
  if(!exists("current_participant")){
    
    current_participant <- participant_by_n[i]
    
    current_participant_log_filepaths <- list.files(path = paste0(parent_folder, "/raw_data/"),
                                                    pattern = paste0(current_participant, "_logfile.dat$"), recursive = T, full.names = T)
    
    current_day1_log_filepath <- current_participant_log_filepaths[grepl(pattern = "Day1", 
                                                                         x = current_participant_log_filepaths)]
    
    if (length(current_day1_log_filepath) != 1) {
      stop(paste0("more or less than one day1 log file for participant ", current_participant))
    }
    
    current_day1_log_file <- read.table(current_day1_log_filepath, header = T, sep = ",")
    
    if(nrow(current_day1_log_file) != gaborgen_stan_list$n_trials){
      stop(paste0("Participant ", current_participant, " does not have all trials in log file. Add this person as a useable participant at top of script."))
    }
    
  } else if (current_participant != participant_by_n[i]) {
    current_participant <- participant_by_n[i]
    
    current_participant_log_filepaths <- list.files(path = paste0(parent_folder, "/raw_data/"),
                                                    pattern = paste0(current_participant, "_logfile.dat$"), recursive = T, full.names = T)
    
    current_day1_log_filepath <- current_participant_log_filepaths[grepl(pattern = "Day1", 
                                                                         x = current_participant_log_filepaths)]
    
    if (length(current_day1_log_filepath) != 1) {
      stop(paste0("more or less than one day1 log file for participant ", current_participant))
    }
    
    current_day1_log_file <- read.table(current_day1_log_filepath, header = T, sep = ",")
    
    if(nrow(current_day1_log_file) != gaborgen_stan_list$n_trials){
      stop(paste0("Participant ", current_participant, " does not have all trials in log file. Add this person as a useable participant at top of script."))
      
    }
  }
  
  #only necessary if the very last trial is missing
  if (length(observed_trials) < current_observation){
    indices_missing <- c(indices_missing, i)
    
    phase <- c(phase, current_day1_log_file$phase[current_trial_number])
    block <- c(block, current_day1_log_file$block[current_trial_number])
    cue <- c(cue, current_day1_log_file$stim[current_trial_number])
    trial <- c(trial, current_day1_log_file$trial[current_trial_number])
    paired <- c(paired, current_day1_log_file$paired[current_trial_number])
    time <- c(time, current_day1_log_file$timeSinceFirstTR[current_trial_number])
    
    current_trial_number <- current_trial_number + 1
    if(current_trial_number > gaborgen_stan_list$n_trial){
      current_trial_number <- 1
    }
  }  else if(current_trial_number == observed_trials[current_observation]){
    indices_observed <- c(indices_observed, i)
    
    phase <- c(phase, current_day1_log_file$phase[current_trial_number])
    block <- c(block, current_day1_log_file$block[current_trial_number])
    cue <- c(cue, current_day1_log_file$stim[current_trial_number])
    trial <- c(trial, current_day1_log_file$trial[current_trial_number])
    paired <- c(paired, current_day1_log_file$paired[current_trial_number])
    time <- c(time, current_day1_log_file$timeSinceFirstTR[current_trial_number])
    
    current_observation <- current_observation + 1
    current_trial_number <- current_trial_number + 1
    if(current_trial_number > gaborgen_stan_list$n_trial){
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
    if(current_trial_number > gaborgen_stan_list$n_trial){
      current_trial_number <- 1
    }
  }
  if(is.na(phase[i])){
    stop()
  }
}


gaborgen_stan_list$indices_observed <- indices_observed %>% as.integer()
gaborgen_stan_list$indices_missing <- indices_missing %>% as.integer()

gaborgen_stan_list$phase <- phase %>% as.integer()
block4 <- numeric()
for (i in 1:length(phase)){
  if(phase[i] == 1){
    block4[i] <- 1
  }
  if(phase[i] == 2 & block[i] == 1){
    block4[i] <- 2
  }
  if(phase[i] == 2 & block[i] == 2){
    block4[i] <- 3
  }
  if(phase[i] == 3){
    block4[i] <- 4
  }
}

gaborgen_stan_list$block <- block4 %>% as.integer()
gaborgen_stan_list$cue <- cue %>% as.integer()
gaborgen_stan_list$trial <- trial %>% as.integer()
gaborgen_stan_list$paired <- paired %>% as.integer()
gaborgen_stan_list$time <- time

gaborgen_stan_list$n_phases <- unique(gaborgen_stan_list$phase) %>% length()
gaborgen_stan_list$n_blocks <- unique(gaborgen_stan_list$block) %>% length()
gaborgen_stan_list$n_cues <- unique(gaborgen_stan_list$cue) %>% length()

# cue_trial_count
cue_trial_count <- c()
for (i in 1:gaborgen_stan_list$n) {
  if (gaborgen_stan_list$trial[i] == 1){
    cue_count_vector <- rep(0,gaborgen_stan_list$n_cues)
  }
  cue_count_vector[gaborgen_stan_list$cue[i]] <- 
    cue_count_vector[gaborgen_stan_list$cue[i]] + 1
  
  cue_trial_count[i] <- cue_count_vector[gaborgen_stan_list$cue[i]]
}

gaborgen_stan_list$cue_trial_count <- cue_trial_count %>% as.integer()

phase_trial_count <- c()
for (i in 1:gaborgen_stan_list$n) {
  if (gaborgen_stan_list$trial[i] == 1){
    phase_count_vector <- rep(0,gaborgen_stan_list$n_cues)
  }
  phase_count_vector[gaborgen_stan_list$phase[i]] <- 
    phase_count_vector[gaborgen_stan_list$phase[i]] + 1
  
  phase_trial_count[i] <- phase_count_vector[gaborgen_stan_list$phase[i]]
}

gaborgen_stan_list$phase_trial_count <- phase_trial_count %>% as.integer()

block_trial_count <- c()
for (i in 1:gaborgen_stan_list$n) {
  if (gaborgen_stan_list$trial[i] == 1){
    block_count_vector <- rep(0,gaborgen_stan_list$n_cues)
  }
  block_count_vector[gaborgen_stan_list$block[i]] <- 
    block_count_vector[gaborgen_stan_list$block[i]] + 1
  
  block_trial_count[i] <- block_count_vector[gaborgen_stan_list$block[i]]
}

gaborgen_stan_list$block_trial_count <- block_trial_count %>% as.integer()

# 
# 
# 
# gaborgen_stan_list$cue <- Oz_fft_df$cue %>% as.integer()
# gaborgen_stan_list$n_cues <- gaborgen_stan_list$cue %>% unique() %>% length()
# 
# gaborgen_stan_list$phase <- Oz_fft_df$phase %>% as.integer()
# gaborgen_stan_list$n_phases <- gaborgen_stan_list$phase %>% unique() %>% length()
# 
# gaborgen_stan_list$paired <- ifelse(Oz_fft_df$paired == "shock", 1, 0) %>% as.integer()

gaborgen_stan_list$amplitude <- Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(amplitude = scale(amplitude_15Hz_fft)) %>% 
  pull(amplitude) %>% 
  as.vector()

gaborgen_stan_list$arousal_pbc <- ratings_df %>% 
  filter(partInd %in% unique(Oz_fft_df$participant)) %>% 
  arrange(condition,ratInd,partInd) %>% 
  pull(ar) %>% 
  array(dim = c(gaborgen_stan_list$n_participants,
                gaborgen_stan_list$n_blocks,
                gaborgen_stan_list$n_cues))

gaborgen_stan_list$arousal_pbc_centered <- ratings_df %>% 
  filter(partInd %in% unique(Oz_fft_df$participant)) %>% 
  mutate(arousal_centered = ar - mean(ar)) %>% 
  arrange(condition,ratInd,partInd) %>% 
  pull(arousal_centered) %>% 
  array(dim = c(gaborgen_stan_list$n_participants,
                gaborgen_stan_list$n_blocks,
                gaborgen_stan_list$n_cues))


gaborgen_stan_list$expect_pbc <- ratings_df %>% 
  filter(partInd %in% unique(Oz_fft_df$participant)) %>% 
  arrange(condition,ratInd,partInd) %>% 
  mutate(exp_div_10 = exp /10) %>% 
  pull(exp_div_10) %>% 
  array(dim = c(gaborgen_stan_list$n_participants,
                gaborgen_stan_list$n_blocks,
                gaborgen_stan_list$n_cues))

gaborgen_stan_list$expect_pbc_centered <- ratings_df %>% 
  filter(partInd %in% unique(Oz_fft_df$participant)) %>% 
  mutate(exp_div_10 = exp /10) %>%
  mutate(expect_centered = exp_div_10 - mean(exp_div_10)) %>% 
  arrange(condition,ratInd,partInd) %>% 
  pull(expect_centered) %>% 
  array(dim = c(gaborgen_stan_list$n_participants,
                gaborgen_stan_list$n_blocks,
                gaborgen_stan_list$n_cues))

  

# Sampling settings ####
number_of_chains <- 8
number_of_parallel_chains <- 4
warmup_samples_per_chain <- 5000
posterior_samples_per_chain <- 5000
where_to_save_chains <- paste0(parent_folder,"/stan_chains")
# used for saved .Rdata
# number_of_chains <- 8
# number_of_parallel_chains <- 4
# warmup_samples_per_chain <- 5000
# posterior_samples_per_chain <- 5000
# where_to_save_chains <- paste0(parent_folder,"/stan_chains")

# Fit Models ####

## Model 001 - Descriptive model: Fatigue slope, cue by phase ####

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model001.stan'

#force_recompile = T is sometimes helpful
model001 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
model001$print()

# Clear previous chains
# list.files(pattern = "Model001",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model001_fit <- model001$sample(data = gaborgen_stan_list,
                                refresh = 500,
                                seed = 2,
                                iter_warmup = warmup_samples_per_chain, 
                                iter_sampling = posterior_samples_per_chain, 
                                save_warmup = F, 
                                show_messages =T,
                                output_dir = where_to_save_chains,
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
model001_fit_meta_data <- model001_fit$metadata()

model001_fit_relevant_parameters <- model001_fit_meta_data$model_params[
  !str_detect(model001_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model001_fit_summary <- 
  model001_fit$summary(variables = model001_fit_relevant_parameters)

model001_fit_loo <- model001_fit$loo()

loo::loo_compare(model001_fit_loo)


model001_df <- posterior::as_draws_df(
  model001_fit$draws(variables = model001_fit_relevant_parameters))

model001_df %>%
  select(starts_with("sigma[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  scale_x_continuous(breaks = seq(0,1.2,by = .2))

model001_df %>%
  select(starts_with("intercept")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model001_df %>%
  select(starts_with("intercept[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))
  
model001_df %>%
  select(starts_with("fatigue")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))
  
model001_df %>%
  select(starts_with("fatigue[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))
  
model001_df %>%
  select(starts_with("bcue")) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, y = name), alpha = .5) +
  scale_x_continuous(name = "Z-Scored ssVEP", breaks = seq(-.75, .75, by = .25)) +
  coord_cartesian(xlim = c(-.5,.5)) +
  # scale_y_discrete(labels = c("Hab, CS+",
  #                             "Hab, GS1",
  #                             "Hab, GS2",
  #                             "Hab, GS3",
  #                             "Acq, CS+",
  #                             "Acq, GS1",
  #                             "Acq, GS2",
  #                             "Acq, GS3",
  #                             "Ext, CS+",
  #                             "Ext, GS1",
  #                             "Ext, GS2",
  #                             "Ext, GS3")) +
  ggtitle("Cue by phase") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())


  
## Model 002 - Rescola-Wagner model: 1 learning rate per participant ####

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model002.stan'

#force_recompile = T is sometimes helpful
model002 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
# model002$print()

# Clear previous chains
# list.files(pattern = "Model002",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model002_fit <- model002$sample(data = gaborgen_stan_list,
                                refresh = 100,
                                seed = 2,
                                iter_warmup = warmup_samples_per_chain, 
                                iter_sampling = posterior_samples_per_chain, 
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)

# Check summary to see if everything was estimated as expected
model002_fit_meta_data <- model002_fit$metadata()

model002_fit_relevant_parameters <- model002_fit_meta_data$model_params[
  !str_detect(model002_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|CSP_|S|raw")]

model002_fit_summary <- 
  model002_fit$summary(variables = model002_fit_relevant_parameters)

model002_fit_loo <- model002_fit$loo()


loo::loo_compare(model001_fit_loo, 
                 model002_fit_loo)

#high k-s don't seem to matter
model001_fit_loo$pointwise[,1][model002_fit_loo$pointwise[,5] > .7]
model002_fit_loo$pointwise[,1][model002_fit_loo$pointwise[,5] > .7]


loo::loo_model_weights(list(model001_fit_loo, 
                            model002_fit_loo))

model002_df <- posterior::as_draws_df(
  model002_fit$draws(variables = model002_fit_relevant_parameters))

model002_df %>%
  select(starts_with("sigma[")) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  scale_x_continuous(breaks = seq(0,1.1, by = .1))

model002_df %>%
  select(starts_with("intercept[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model002_df %>%
  select(starts_with("fatigue[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model002_df %>%
  select(starts_with("learning_rate[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),
                      stat = "binline", bins = 400) +
  coord_cartesian(xlim = c(0,.3))
                  
model002_df %>%
  select(starts_with("learning_rate[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),
                      stat = "binline", bins = 400) +
  coord_cartesian(xlim = c(0,1))


model002_df %>%
  select(starts_with("scaling")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(-2, 5))


## Model 003 - Rescola-Wagner: different learning rates for paired or unpaired per participant ####


model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model003.stan'

#force_recompile = T is sometimes helpful
model003 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
# model003$print()

# Clear previous chains
# list.files(pattern = "Model003",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model003_fit <- model003$sample(data = gaborgen_stan_list,
                                refresh = 100,
                                seed = 2,
                                iter_warmup = warmup_samples_per_chain, 
                                iter_sampling = posterior_samples_per_chain, 
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)


model003_fit_meta_data <- model003_fit$metadata()

model003_fit_relevant_parameters <- model003_fit_meta_data$model_params[
  # !str_detect(model003_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|CSP_|S|raw")]
  !str_detect(model003_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|CSP_|S")]

model003_fit_summary <- 
  model003_fit$summary(variables = model003_fit_relevant_parameters)


model003_fit_loo <- model003_fit$loo()

#high k-s don't seem to matter
model001_fit_loo$pointwise[,1][model003_fit_loo$pointwise[,5] > .7]
model002_fit_loo$pointwise[,1][model003_fit_loo$pointwise[,5] > .7]
model003_fit_loo$pointwise[,1][model003_fit_loo$pointwise[,5] > .7]

loo::loo_compare(model001_fit_loo, 
                 # model002_fit_loo,
                 model003_fit_loo)

loo::loo_model_weights(list(model001_fit_loo, 
                            # model002_fit_loo,
                            model003_fit_loo))

loo::loo_compare(model002_fit_loo,
                 model003_fit_loo)

loo::loo_model_weights(list(model002_fit_loo,
                            model003_fit_loo))

model003_fit_relevant_parameters <- model003_fit_meta_data$model_params[
  !str_detect(model003_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|CSP_|S|raw")]

model003_df <- posterior::as_draws_df(
  model003_fit$draws(variables = model003_fit_relevant_parameters))

model003_df %>%
  select(starts_with("sigma[")) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model003_df %>%
  select(starts_with("intercept[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model003_df %>%
  select(starts_with("fatigue[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model003_df %>%
  select(starts_with("learning_paired_average")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density(aes(x = value)) +
  coord_cartesian(xlim = c(0, 1)) |
  
  model003_df %>%
  select(starts_with("learning_unpaired_average")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density(aes(x = value)) +
  coord_cartesian(xlim = c(0, 1))

model003_df %>%
  select(starts_with("learning_paired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(0, 1)) |
  
  model003_df %>%
  select(starts_with("learning_unpaired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(0, 1))

model003_df %>%
  select(starts_with("learning_paired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),
                      stat = "binline", bins = 100) +
  coord_cartesian(xlim = c(0, 1)) |
  
  model003_df %>%
  select(starts_with("learning_unpaired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),
                      stat = "binline", bins = 100) +
  coord_cartesian(xlim = c(0, 1))

model003_df %>%
  select(starts_with("scaling")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(-1, 1.5))

## Model 004 - Rescola-Wagner: model003 without inv_logit of learning rates ####


model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model004.stan'

#force_recompile = T is sometimes helpful
model004 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
# model004$print()

# Clear previous chains
# list.files(pattern = "Model004",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model004_fit <- model004$sample(data = gaborgen_stan_list,
                                refresh = 100,
                                seed = 2,
                                iter_warmup = warmup_samples_per_chain, 
                                iter_sampling = posterior_samples_per_chain, 
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)


model004_fit_meta_data <- model004_fit$metadata()

model004_fit_relevant_parameters <- model004_fit_meta_data$model_params[
  !str_detect(model004_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|raw|CSP_|S")]

model004_fit_summary <- 
  model004_fit$summary(variables = model004_fit_relevant_parameters)


model004_fit_loo <- model004_fit$loo()

loo::loo_compare(model001_fit_loo, 
                 model002_fit_loo,
                 model003_fit_loo,
                 model004_fit_loo)

loo::loo_model_weights(list(model001_fit_loo, 
                            model002_fit_loo,
                            model003_fit_loo,
                            model004_fit_loo))

model004_df <- posterior::as_draws_df(
  model004_fit$draws(variables = model004_fit_relevant_parameters))

model004_df %>%
  select(starts_with("sigma")) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model004_df %>%
  select(starts_with("sigma[")) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model004_df %>%
  select(starts_with("intercept")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model004_df %>%
  select(starts_with("intercept[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model004_df %>%
  select(starts_with("fatigue")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model004_df %>%
  select(starts_with("fatigue[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model004_df %>%
  select(starts_with("learning_paired_average")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density(aes(x = value)) +
  coord_cartesian(xlim = c(0, 1)) |
  
  model004_df %>%
  select(starts_with("learning_unpaired_average")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density(aes(x = value)) +
  coord_cartesian(xlim = c(0, 1))

model004_df %>%
  select(starts_with("learning_paired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(0, 1)) |
  
  model004_df %>%
  select(starts_with("learning_unpaired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(0, 1))

model004_df %>%
  select(starts_with("learning_paired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),
                      stat = "binline", bins = 50) +
  coord_cartesian(xlim = c(0, 1)) |
  
  model004_df %>%
  select(starts_with("learning_unpaired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name),
                      stat = "binline", bins = 50) +
  coord_cartesian(xlim = c(0, 1))

model004_df %>%
  select(starts_with("scaling")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(-1, 2.5))


## Model 005 - Arousal rating model ####

# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model005.stan'
#arousal centered version
model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model006.stan'

#force_recompile = T is sometimes helpful
model005 <- cmdstanr::cmdstan_model(model_path,
                                    force_recompile = T)


#Model source code
# model005$print()

# Clear previous chains
# list.files(pattern = "Model005",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model005_fit <- model005$sample(data = gaborgen_stan_list,
                                refresh = 50,
                                seed = 3,
                                iter_warmup = warmup_samples_per_chain, 
                                iter_sampling = posterior_samples_per_chain, 
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)

# model006_fit <- model006$sample(data = gaborgen_stan_list,
#                                 refresh = 50,
#                                 seed = 3,
#                                 iter_warmup = warmup_samples_per_chain, 
#                                 iter_sampling = posterior_samples_per_chain, 
#                                 save_warmup = F, 
#                                 show_messages = T,
#                                 output_dir = where_to_save_chains,
#                                 chains = number_of_chains,
#                                 parallel_chains = number_of_parallel_chains)


model005_fit_meta_data <- model005_fit$metadata()

model005_fit_relevant_parameters <- model005_fit_meta_data$model_params[
  !str_detect(model005_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model005_fit_summary <- 
  model005_fit$summary(variables = model005_fit_relevant_parameters)

model005_fit_loo <- model005_fit$loo()

model006_fit_meta_data <- model006_fit$metadata()

model006_fit_relevant_parameters <- model006_fit_meta_data$model_params[
  !str_detect(model006_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model006_fit_summary <- 
  model006_fit$summary(variables = model006_fit_relevant_parameters)

model006_fit_loo <- model006_fit$loo()

loo::loo_compare(model005_fit_loo,
                 model006_fit_loo)

loo::loo_compare(model001_fit_loo, 
                 model002_fit_loo,
                 model003_fit_loo,
                 model004_fit_loo,
                 model005_fit_loo,
                 model006_fit_loo)

loo::loo_model_weights(list(model001_fit_loo,
                            model002_fit_loo,
                            model003_fit_loo,
                            model004_fit_loo,
                            model005_fit_loo))

loo::loo_model_weights(list(model001_fit_loo,
                            model002_fit_loo))

model005_df <- model005_fit$draws(variables = model005_fit_relevant_parameters,
                                  format = "df")
model006_df <- model006_fit$draws(variables = model006_fit_relevant_parameters,
                                  format = "df")


model005_df %>%
  select(starts_with("sigma[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model005_df %>%
  select(starts_with("intercept")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model005_df %>%
  select(starts_with("fatigue")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model005_df %>%
  select(starts_with("b_arousal")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

(model005_df$b_arousal_average > 0) %>% 
  sum()/ (number_of_chains*posterior_samples_per_chain)

model006_df %>%
  select(starts_with("b_arousal")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

(model005_df$b_arousal_average > 0) %>% 
  sum()/ (number_of_chains*posterior_samples_per_chain)

(model006_df$b_arousal_average > 0) %>% 
  sum()/ (number_of_chains*posterior_samples_per_chain)


## Model 007 - Arousal and expectancy rating model ####
model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model007.stan'

#force_recompile = T is sometimes helpful
model007 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
# model007$print()

# Clear previous chains
# list.files(pattern = "Model007",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model007_fit <- model007$sample(data = gaborgen_stan_list,
                                refresh = 50,
                                seed = 3,
                                # iter_warmup = 1000, 
                                # iter_sampling = 1000, 
                                iter_warmup = warmup_samples_per_chain,
                                iter_sampling = posterior_samples_per_chain,
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                # chains = 4,
                                # parallel_chains = 4)
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)


model007_fit_meta_data <- model007_fit$metadata()

model007_fit_relevant_parameters <- model007_fit_meta_data$model_params[
  !str_detect(model007_fit_meta_data$model_params, "log_lik|mu|amplitude|S")]
  # !str_detect(model007_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model007_fit_summary <- 
  model007_fit$summary(variables = model007_fit_relevant_parameters)

model007_fit_loo <- model007_fit$loo()

loo::loo_compare(model001_fit_loo,
                 model003_fit_loo,
                 model005_fit_loo,
                 model006_fit_loo,
                 model007_fit_loo)

loo::loo_model_weights(list(model001_fit_loo,
                            model007_fit_loo))

model007_df <- model007_fit$draws(variables = model007_fit_relevant_parameters,
                                  format = "df")


model007_df %>%
  select(starts_with("sigma[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model007_df %>%
  select(starts_with("intercept[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model007_df %>%
  select(starts_with("fatigue[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

(model005_df %>%
  select(starts_with("b_arousal_average")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(-.02, .06))) /
  
  model007_df %>%
  select(starts_with("b_arousal_average")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +
  coord_cartesian(xlim = c(-.02, .06))

model007_df %>%
  select(starts_with("b_expect_average")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model007_df %>%
  select(starts_with("b_expect[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model007_df %>%
  select(starts_with("b_arousal_expect[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

model007_df %>%
  select(starts_with("b_arousal[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +

model007_df %>%
  select(starts_with("b_expect[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name)) +

model007_df %>%
  select(starts_with("b_arousal_expect[")) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name,
                       levels = rev(unique(name)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, y = name))

(model007_df$b_arousal_average > 0) %>% 
  sum()/ (number_of_chains*posterior_samples_per_chain)

(model007_df$b_expect_average > 0) %>% 
  sum()/ (number_of_chains*posterior_samples_per_chain)


## Model 008 - just adaptation slop ####
model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/Model008.stan'

#force_recompile = T is sometimes helpful
model008 <- cmdstanr::cmdstan_model(model_path, 
                                    force_recompile = T)

#Model source code
# model008$print()

# Clear previous chains
# list.files(pattern = "Model008",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

model008_fit <- model008$sample(data = gaborgen_stan_list,
                                refresh = 50,
                                seed = 3,
                                # iter_warmup = 1000, 
                                # iter_sampling = 1000, 
                                iter_warmup = warmup_samples_per_chain,
                                iter_sampling = posterior_samples_per_chain,
                                save_warmup = F, 
                                show_messages = T,
                                output_dir = where_to_save_chains,
                                # chains = 4,
                                # parallel_chains = 4)
                                chains = number_of_chains,
                                parallel_chains = number_of_parallel_chains)


model008_fit_meta_data <- model008_fit$metadata()

model008_fit_relevant_parameters <- model008_fit_meta_data$model_params[
  !str_detect(model008_fit_meta_data$model_params, "log_lik|mu|amplitude|S")]
# !str_detect(model007_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model008_fit_summary <- 
  model008_fit$summary(variables = model008_fit_relevant_parameters)

model008_fit_loo <- model008_fit$loo()

loo::loo_compare(model001_fit_loo,
                 model008_fit_loo)


# Save data and fits ####
save(fft_df,
     Oz_fft_df,
     gaborgen_stan_list,
     model001_fit,
     model001_fit_summary,
     model001_fit_loo,
     model002_fit,
     model002_fit_summary,
     model002_fit_loo,
     model003_fit,
     model003_fit_summary,
     model003_fit_loo,
     model004_fit,
     model004_fit_summary,
     model004_fit_loo,
     model005_fit,
     model005_fit_summary,
     model005_fit_loo,
     model006_fit,
     model006_fit_summary,
     model006_fit_loo,
     model007_fit,
     model007_fit_summary,
     model007_fit_loo,
     model008_fit,
     model008_fit_summary,
     model008_fit_loo,
     file = paste0(parent_folder,"/gaborgen_eeg_manuscript_data_models.RData"))

# Load to create manuscript stats and figures ####
library(tidyverse)
library(cmdstanr)
library(R.matlab)
library(ggbeeswarm)
library(patchwork)
library(ggridges)

parent_folder <- "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI"
git_repository <- "/home/andrewf/Repositories/gaborgen"

load(file = paste0(parent_folder,"/gaborgen_eeg_manuscript_data_models.RData"))

## Load posterior sample data frames for figures ####
model001_fit_meta_data <- model001_fit$metadata()

model001_fit_relevant_parameters <- model001_fit_meta_data$model_params[
  !str_detect(model001_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model001_df <- posterior::as_draws_df(
  model001_fit$draws(variables = model001_fit_relevant_parameters))

model003_fit_meta_data <- model003_fit$metadata()

model003_fit_relevant_parameters <- model003_fit_meta_data$model_params[
  !str_detect(model003_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|CSP_|S|raw")]

model003_df <- posterior::as_draws_df(
  model003_fit$draws(variables = model003_fit_relevant_parameters))

model008_fit_meta_data <- model008_fit$metadata()

model008_fit_relevant_parameters <- model008_fit_meta_data$model_params[
  !str_detect(model008_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model008_df <- model008_fit$draws(variables = model008_fit_relevant_parameters,
                                  format = "df")

# Sample demographics ####
master_excel <- readxl::read_excel(
  paste0(parent_folder,"/GABORGEN24_MRI_mastersheet.xlsx"),
  sheet = 2,
  skip = 2,
  col_names = c("participant", "blank", "GUID", "date", "time", 
                "experimenters", "sex", "gender", "race", "ethnicity", 
                "age_years", "age_months", "handedness", "seizures",
                "vision", "glasses_contacts", "head_size_cm", 
                "shock_level_mA", "learned_response", "day1_notes",
                "day2_date","day2_time","day1_to_2_changes", "day2_notes",
                "blank"))

sample_demographics <- master_excel %>% 
  select(!starts_with("blank")) %>% 
  filter(participant %in% (Oz_fft_df$participant %>% unique())) %>% 
  mutate(gender = case_when(participant %in% c(139,145) ~ "Male",
                            participant %in% c(140:144) ~"Female",
                            .default = gender)) %>% 
  mutate(race = case_when(participant %in% c(139:142,144) ~ "White",
                          participant %in% c(143) ~"Asian",
                          participant %in% c(121,131,145) ~"Multiracial",
                          .default = race)) %>% 
  mutate(ethnicity = case_when(participant %in% c(139,140,142,143,144,145) ~ "Non-Hispanic",
                               participant %in% c(141) ~"Hispanic",
                               .default = ethnicity)) %>% 
  mutate(age_years = case_when(participant %in% c(139,140,143,144) ~ 18,
                               participant %in% c(141) ~ 20,
                               participant %in% c(142) ~ 21,
                               participant %in% c(145) ~ 28,
                               .default = age_years)) %>% 
  mutate(age_months = age_years *12) %>% 
  mutate(across(where(is.character), ~ na_if(.x, "N/A")))

  
sample_demographics %>% 
  reframe(avg_age = mean(age_years),
          sd_age = sd(age_years))

sample_demographics %>% 
  group_by(sex) %>% 
  reframe(n())

sample_demographics %>%
  group_by(race) %>% 
  reframe(n())

sample_demographics %>%
  group_by(ethnicity) %>% 
  reframe(n())

sample_demographics %>% 
  mutate(shock_level_mA = as.numeric(shock_level_mA)) %>% 
  reframe(avg_shock = mean(shock_level_mA, na.rm = T),
          sd_shock = sd(shock_level_mA, na.rm = T))

## general plot settings ####
cue_color <- c("red1","green1", "purple1", "blue1")
text_font <- "Arial"
text_size <- 15
axis_line_thickness <- .75

## Overall means and standard errors
annotation_text_size <- 9
fig1_dot_size <- .8
range_linewidth <- 1.5

Oz_fft_df %>% 
  reframe(n(),
          n()/gaborgen_stan_list$n_participants,
          plotrix::std.error(n()/gaborgen_stan_list$n_participants),
          gaborgen_stan_list$n_participants * gaborgen_stan_list$n_trials)

Oz_fft_df %>% 
  group_by(participant) %>% 
  reframe(trial_count = n()) %>% 
  reframe(avg_trial_count = mean(trial_count),
          se_trial_count = sd(trial_count))

Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = amplitude_15Hz_fft %>% scale() %>% as.vector()) %>% 
  group_by(cue, block) %>%  
  reframe(mean_zamp = mean(zamp),
          se_zamp = plotrix::std.error(zamp)) %>% 
  ggplot() +
  geom_pointrange(aes(x = cue, y = mean_zamp,
                      ymin = mean_zamp - se_zamp,
                      ymax = mean_zamp + se_zamp,
                      group = cue),
                  size = fig1_dot_size +.2,
                  linewidth = range_linewidth +.75) +
  geom_pointrange(aes(x = cue, y = mean_zamp,
                      ymin = mean_zamp - se_zamp,
                      ymax = mean_zamp + se_zamp,
                      color = cue),
                  size = fig1_dot_size,
                  linewidth = range_linewidth) +
  scale_color_manual(values = cue_color) +
  scale_y_continuous(name = "Z-Scored ssVEP",
                     expand = c(0,0),
                     breaks = seq(-.5,.5,by = .25)) +
  coord_cartesian(ylim = c(-.6,.6)) +
  facet_wrap(~ block, nrow = 1, 
             labeller = labeller(
               block = c("1" = "Habituation",
                         "2" = "Acquisition #1",
                         "3" = "Acquisition #2",
                         "4" = "Extinction"))) +
  geom_text(data = data.frame(cue = c(2.5),
                              mean_zamp = c(.5,
                                            .5-.125,
                                            .5-.125-.125,
                                            .5-.125-.125-.125),
                              block = c(4,4,4,4),
                              anno_label = c("CS+","GS1","GS2","GS3")),
              # aes(x = 2.5, y = .5, 
              #   block = "4", label = "CS+")
            aes(x = cue, 
                y = mean_zamp, 
                label = anno_label),
            color = cue_color,
            family = "Arial",
            size = annotation_text_size) +
  # annotate(
  #   "text",
  #   x = 2.5, y = .5, block = 4,
  #   label = "CS+",
  #   size = annotation_text_size,
  #   color = cue_color[1]) +
  theme_classic() +
  theme(text = element_text(family = text_font,
                            size = text_size,
                            color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position = "none",
        axis.line = element_line(linewidth = axis_line_thickness,
                                 lineend = "square"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
    )
  

ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure1.png",
       dpi = 300, units = "in", 
       height = 1, width = 2,
       scale = 3.5)


## Loo stats ####


loo::loo_compare(model008_fit_loo, #model 1 in paper
                 model001_fit_loo, #model 2
                 model003_fit_loo) #model 3

loo::loo_model_weights(list(model008_fit_loo,
                            model001_fit_loo,
                            model003_fit_loo))

loo::loo_compare(model008_fit_loo,
                 model001_fit_loo)

loo::loo_model_weights(list(model008_fit_loo,
                            model001_fit_loo))

loo::loo_compare(model008_fit_loo,
                 model001_fit_loo,
                 model003_fit_loo,
                 model002_fit_loo, # one learning rate
                 model004_fit_loo, # 2 learning rates no invlogit
                 model005_fit_loo, # block arousal slope
                 # model006_fit_loo, # block centered-arousal slope 
                 model007_fit_loo) # block arousal, expectancy, and arousal-expectancy slope

loo::loo_model_weights(list(model008_fit_loo,
                            model001_fit_loo,
                            model003_fit_loo,
                            model002_fit_loo, # one learning rate
                            model004_fit_loo, # 2 learning rates no invlogit
                            model005_fit_loo, # block arousal slope
                            model006_fit_loo, # block centered-arousal slope 
                            model007_fit_loo)) # block arousal, expectancy, and arousal-expectancy slope


loo::loo_compare(model001_fit_loo,
                 model008_fit_loo)

loo::loo_model_weights(list(model001_fit_loo,
                            model008_fit_loo))

Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(participant) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  mutate(exp_mean = exp(-1*mean_elpd_1_minus_3)) %>% 
  print(n = 999)



## Loo contributions figure ####
loo_dot_size <- .75
loo_dot_alpha <- 1
loo_par_dot_size <- .8
loo_par_range_linewidth <- 1.5
loo_figure_text <-20


CV_by_cue_by_block_plot <-
  Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(block, cue) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  mutate(sum_better_for_mod3 = sum_elpd_1_minus_3 < 0) %>% 
  ggplot() +
  geom_vline(aes(xintercept = cue),
             color = "lightgray",
             linetype = "dotted") +
  geom_hline(aes(yintercept = 0), 
             linewidth = axis_line_thickness) +
  geom_pointrange(aes(x = cue,
                      y = sum_elpd_1_minus_3,
                      ymin = sum_elpd_1_minus_3 - sum_se,
                      ymax = sum_elpd_1_minus_3 + sum_se,
                      color = sum_better_for_mod3),
                  size = loo_par_dot_size,
                  linewidth = loo_par_range_linewidth) +
  facet_wrap(~block, nrow = 1, 
             labeller = labeller(
               block = c("1" = "Habituation",
                         "2" = "Acquisition #1",
                         "3" = "Acquisition #2",
                         "4" = "Extinction"))) +
    geom_text(data = data.frame(cue = c(2.5),
                                elpd_value = c(4.5),
                                block = c(1,2,3,4),
                                anno_label = c("Habituation","Acquisition #1","Acquisition #2","Extinction")),
              # aes(x = 2.5, y = .5,
              #   block = "4", label = "CS+")
              aes(x = cue,
                  y = elpd_value,
                  label = anno_label),
              color = "black",
              family = "Arial",
              size = 7.5) +
  scale_color_manual(values = c("red1","blue1")) +
  scale_y_continuous(name = "ELPD") +
  coord_cartesian(ylim = c(-19, 4)) +
  scale_x_discrete(name = "Cue",labels = c("CS+", "GS1", "GS2", "GS3")) +
  ggtitle("Cross-Validation Accuracy by Cue and Block") +
  theme_classic() +
  theme(text = element_text(family = text_font,
                            size = loo_figure_text,
                            color = "black"),
        axis.line = element_line(linewidth = axis_line_thickness,
                                 lineend = "square"),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )

CV_by_participant_plot <-
  Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(participant) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  mutate(sum_better_for_mod3 = sum_elpd_1_minus_3 < 0) %>% 
  ggplot() +
  geom_vline(aes(xintercept = participant),
             color = "lightgray",
             linetype = "dotted") +
  geom_hline(aes(yintercept = 0), 
             linewidth = axis_line_thickness) +
  geom_pointrange(aes(x = participant,
                      y = sum_elpd_1_minus_3,
                      ymin = sum_elpd_1_minus_3 - sum_se,
                      ymax = sum_elpd_1_minus_3 + sum_se,
                      color = sum_better_for_mod3),
                  size = loo_par_dot_size,
                  linewidth = loo_par_range_linewidth) +
    annotate(geom = "text", 
             x = 19, 
             y = -13,
             label = "Cue by Block\nModel Better",
             family = "Arial",
             color = "red1",
             lineheight = 0.8, # Adjust this value to decrease spacing
             size = 12)+
    annotate(geom = "text", 
             x = 19, 
             y = -17,
             label = "Learning\nModel Better",
             family = "Arial",
             color = "blue1",
             lineheight = 0.8, # Adjust this value to decrease spacing
             size = 12)+
  scale_color_manual(values = c("red1","blue1")) +
  scale_y_continuous(name = "ELPD") +
  coord_cartesian(ylim = c(-19, 4)) +
  scale_x_discrete(name = "Participant") +
  ggtitle("Cross-Validation Accuracy By Participant") +
  theme_classic() +
  theme(text = element_text(family = text_font,
                            size = loo_figure_text,
                            color = "black"),
        axis.line = element_line(linewidth = axis_line_thickness,
                                 lineend = "square"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )

layout_grid <- c('
AAAABBBB
')

CV_by_cue_by_block_plot +
  CV_by_participant_plot +
  plot_annotation(title = "Comparing Cross-Validation Contributions",
                  theme = theme(
                    plot.title = element_text(family = "Arial",
                                              size = 27,
                                              color = "black",
                                              hjust = 0.5,
                                              face = "bold")
                  )) +
  plot_layout(design = layout_grid)



ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure2.png",
       dpi = 300, units = "in", 
       height = 1.5, width = 3.75,
       scale = 5)


## Adaptation ####
# plot specifics
number_of_samples_to_plot <- 4000
number_of_chains_prior <- 4
prior_samples_per_chain <- number_of_samples_to_plot / number_of_chains_prior
line_width_adapt <- .1
line_alpha_adapt <- .1

## Median participant adaptation

# median_model001_adaptation_df <- model001_df %>% 
#   select(starts_with("intercept["), starts_with("fatigue[")) %>% 
#   pivot_longer(
#     cols = everything(),
#     names_to = c("parameter_type", "index"),
#     names_pattern = "(intercept|fatigue)\\[(.*)\\]",
#     values_to = "value") %>%
#   group_by(parameter_type, index) %>% 
#   reframe(median_value = median(value)) %>% 
#   pivot_wider(
#     names_from = parameter_type,
#     values_from = median_value) %>% 
#   mutate(index = as.numeric(index)) %>% 
#   arrange(index) %>% 
#   mutate(index = factor(index, 
#                         levels = 1:gaborgen_stan_list$n_participants))
# 
# median_model003_adaptation_df <- model003_df %>% 
#   select(starts_with("intercept["), starts_with("fatigue[")) %>% 
#   pivot_longer(
#     cols = everything(),
#     names_to = c("parameter_type", "index"),
#     names_pattern = "(intercept|fatigue)\\[(.*)\\]",
#     values_to = "value") %>%
#   group_by(parameter_type, index) %>% 
#   reframe(median_value = median(value)) %>% 
#   pivot_wider(
#     names_from = parameter_type,
#     values_from = median_value) %>% 
#   mutate(index = as.numeric(index)) %>% 
#   arrange(index) %>% 
#   mutate(index = factor(index, 
#                         levels = 1:gaborgen_stan_list$n_participants))
# 
# median_model008_adaptation_df <- model008_df %>% 
#   select(starts_with("intercept["), starts_with("fatigue[")) %>% 
#   pivot_longer(
#     cols = everything(),
#     names_to = c("parameter_type", "index"),
#     names_pattern = "(intercept|fatigue)\\[(.*)\\]",
#     values_to = "value") %>%
#   group_by(parameter_type, index) %>% 
#   reframe(median_value = median(value)) %>% 
#   pivot_wider(
#     names_from = parameter_type,
#     values_from = median_value) %>% 
#   mutate(index = as.numeric(index)) %>% 
#   arrange(index) %>% 
#   mutate(index = factor(index, 
#                         levels = 1:gaborgen_stan_list$n_participants))



# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/priors_visualization_2.stan'
# 
# model_priors <- cmdstanr::cmdstan_model(model_path, 
#                                         force_recompile = T)
# 
# #Model source code
# # model_priors$print()
# 
# model_priors_fit <- model_priors$sample(refresh = 1000,
#                                         seed = 4,
#                                         iter_warmup = prior_samples_per_chain, 
#                                         iter_sampling = prior_samples_per_chain, 
#                                         save_warmup = F, 
#                                         show_messages = T,
#                                         output_dir = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains",
#                                         chains = number_of_chains_prior,
#                                         parallel_chains = 4)
# 
# model_priors_fit_df <- model_priors_fit$draws(format = "df")

number_of_samples_to_plot <- 4000

set.seed(0)
samples_to_plot <- sample(1:nrow(model001_df),
                          size = number_of_samples_to_plot,
                          replace = F)

## Adapt plot settings####
line_width_adapt <- .1
line_alpha_adapt <- .1
text_size <- 15
median_line_width_adapt <-.7
trial_breaks <- c(1, seq(20,180, by = 20))

### Prior plot####
# adaptation_prior_plot <-
#   model_priors_fit_df %>% 
#   ggplot() +
#   geom_abline(aes(intercept = intercept_average - (fatigue_average * (gaborgen_stan_list$n_trials/2)), 
#                   slope = fatigue_average),
#               linewidth = line_width_adapt, 
#               alpha = line_alpha_adapt) +
#   scale_y_continuous(limits =c(-2,2), 
#                      breaks = seq(-4, 4, by = 1),
#                      name = "Z-scored ssVEP") +
#   scale_x_continuous(limits =c(1,176), 
#                      breaks = trial_breaks,
#                      expand = c(0,0),
#                      name = "Trial") +
#   theme_classic() +
#   theme(text = element_text(family = "Arial", 
#                             size = text_size),
#         axis.ticks.y = element_blank())

adaptation_model1_avg_plot <- model008_df[samples_to_plot,] %>%
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 36, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48 + 48,
             color = "black", linetype = "dashed") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = trial_breaks,
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma") +
  ggtitle("Average Adaptation") +
  theme_classic() +
  theme(text = element_text(family = "Arial", 
                            size = text_size),
        plot.title = element_text(family = "Arial",
                             size = text_size + 5,
                             hjust = .5),
        axis.ticks.y = element_blank())


cue_block_average_adapt_plot <- model001_df[samples_to_plot,] %>%
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 36, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48 + 48,
             color = "black", linetype = "dashed") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = trial_breaks,
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma") +
  ggtitle("Average Adaptation") +
  theme_classic() +
  theme(text = element_text(family = "Arial", 
                            size = text_size),
        plot.title = element_text(family = "Arial",
                                  size = text_size + 5,
                                  hjust = .5),
        axis.ticks.y = element_blank())

learning_average_adapt_plot <- model003_df[samples_to_plot,] %>%
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  geom_vline(xintercept = 36, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48 + 48,
             color = "black", linetype = "dashed") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = trial_breaks,
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma") +
  ggtitle("Average Adaptation") +
  theme_classic() +
  theme(text = element_text(family = "Arial", 
                            size = text_size),
        plot.title = element_text(family = "Arial",
                                  size = text_size + 5,
                                  hjust = .5),
        axis.ticks.y = element_blank())

### Median participant adaptation####

median_model001_adaptation_df <- model001_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))

median_model003_adaptation_df <- model003_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))

median_model008_adaptation_df <- model008_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))


by_cue_adapt_par_median_plot <- median_model001_adaptation_df %>% 
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 36, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48 + 48, 
             color = "black", linetype = "dashed") +
  geom_abline(aes(intercept = intercept, 
                  slope = fatigue,
                  color = index),
              linewidth = median_line_width_adapt, 
              alpha = 1) +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = trial_breaks,
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma", name = "Participant") +
  ggtitle("Participant Posterior Median Adaptation") +
  theme_classic() +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(
                                linewidth = median_line_width_adapt + 1))) +
  theme(text = element_text(family = "Arial", 
                            size = text_size),
        plot.title = element_text(family = "Arial",
                                  size = text_size + 5,
                                  hjust = .5),
        legend.position = "right",
        legend.direction = "vertical",
        axis.ticks.y = element_blank())

  
  


learning_adapt_median_par_plot <- median_model003_adaptation_df %>% 
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 36, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48 + 48, 
             color = "black", linetype = "dashed") +
  geom_abline(aes(intercept = intercept, 
                  slope = fatigue,
                  color = index),
              linewidth = median_line_width_adapt, 
              alpha = 1) +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = trial_breaks,
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma", name = "Participant") +
  ggtitle("Participant Posterior Median Adaptation") +
  theme_classic() +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(
                                linewidth = median_line_width_adapt + 1))) +
  theme(text = element_text(family = "Arial", 
                            size = text_size),
        plot.title = element_text(family = "Arial",
                                  size = text_size + 5,
                                  hjust = .5),
        legend.position = "right",
        legend.direction = "vertical",
        axis.ticks.y = element_blank())


null_adapt_median_par_plot <- median_model008_adaptation_df %>% 
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_vline(xintercept = 36, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48, 
             color = "black", linetype = "dashed") +
  geom_vline(xintercept = 36 + 48 + 48, 
             color = "black", linetype = "dashed") +
  geom_abline(aes(intercept = intercept, 
                  slope = fatigue,
                  color = index),
              linewidth = median_line_width_adapt, 
              alpha = 1) +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = trial_breaks,
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma", name = "Participant") +
  ggtitle("Participant Posterior Median") +
  theme_classic() +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(
                                linewidth = median_line_width_adapt + 1))) +
  theme(text = element_text(family = "Arial", 
                            size = text_size),
        plot.title = element_text(family = "Arial",
                                  size = text_size + 5,
                                  hjust = .5),
        legend.position = "right",
        legend.direction = "vertical",
        axis.ticks.y = element_blank())


title_size <- 20
modelavg_title <- ggplot() + 
  ggtitle("Average Adaptation") +
  theme_void() +
  theme(plot.title = element_text(family = "Arial", 
                                  size = title_size, 
                                  face = "bold", 
                                  hjust = 0.5))
modelpar_title <- ggplot() + 
  ggtitle("Participant Median") +
  theme_void() +
  theme(plot.title = element_text(family = "Arial", 
                                  size = title_size, 
                                  face = "bold", 
                                  hjust = 0.5))

model1_title <- ggplot() + 
  ggtitle("Model 1: Only Adaptation") +
  theme_void() +
  theme(plot.title = element_text(family = "Arial", 
                                  size = title_size,
                                  hjust = 0.5))

model2_title <- ggplot() + 
  ggtitle("Model 2: Cue by Block") +
  theme_void() +
  theme(plot.title = element_text(family = "Arial", 
                                  size = title_size,
                                  hjust = 0.5))

model3_title <- ggplot() + 
  ggtitle("Model 3: Learning Model") +
  theme_void() +
  theme(plot.title = element_text(family = "Arial", 
                                  size = title_size,
                                  hjust = 0.5))


layout_grid <- c('
aaabbb
cccccc
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
dddiii
eeeeee
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
fffjjj
gggggg
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
hhhkkk
')


# free(adaptation_prior_plot) +
modelavg_title +
  modelpar_title +
model1_title +
  free(adaptation_model1_avg_plot +
         theme(axis.title.x = element_blank(),
               plot.title = element_blank())) +
  model2_title +
  free(cue_block_average_adapt_plot +
     theme(plot.title = element_blank(),
           axis.title.x = element_blank())) +
  model3_title +
  free(learning_average_adapt_plot +
         theme(plot.title = element_blank())) +
  free(null_adapt_median_par_plot +
         theme(axis.title = element_blank(),
               plot.title = element_blank())) +
  free(by_cue_adapt_par_median_plot+
         theme(plot.title = element_blank(),
               axis.title = element_blank())) +
  free(learning_adapt_median_par_plot +
         theme(plot.title = element_blank(),
               axis.title.y = element_blank())) + 
  plot_annotation(title = "Models Estimate Adaptation Over Trials Similarly",
                  theme = theme(
                    plot.title = element_text(family = "Arial",
                                              size = 25,
                                              color = "black",
                                              hjust = 0.5,
                                              face = "bold")
                  ))   +
  plot_layout(
    design = layout_grid, 
    guides = 'collect'
  ) & 
  theme(
    legend.position = 'right',
    legend.text = element_text(hjust = .5),
    legend.title = element_text(hjust = .5),
    legend.margin = margin(t = .5, r = .5, b = .5, l = .5),  # Reduce margin around the legend
    legend.box.spacing = unit(0.01, "cm")
  )


ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/FigureS1.png",
       dpi = 300, units = "in", 
       height = 3, width = 3,
       scale = 3)


## Learning CS+ fig####
# Plot settings
par_axis_line_thickness <- .75
text_size <- 15
line_width_csp <- .2
line_alpha_csp <- .2
line_width_csp_average <- 1.5

learning_paired_plot <-
  model003_df %>%
  select(starts_with("learning_paired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot(aes(x = value, 
             y = factor(name, 
                        levels = rev(levels(name)))
  )) +
  geom_hline(yintercept = c(1:24),
             linewidth = par_axis_line_thickness,
             linetype = "dashed",
             color = "gray")+
  # geom_density_ridges(aes(height = after_stat(scaled)),
  #                     fill = "black",
  #                     scale = .9,
  #                     stat = "density",
  #                     # rel_min_height = 0.01,
  #                     # from = 0,
  #                     # to = 1,
  #                     panel_scaling = FALSE) +
  stat_density_ridges(fill = "black",
                      color = "black",
                      scale = 1,
                      rel_min_height = 0.01,
                      from = 0,
                      panel_scaling = FALSE,
                      to = 1) +
  scale_y_discrete(labels = paste0("Participant ", 
                                   gaborgen_stan_list$n_participants:1)) +
  scale_x_continuous(breaks = seq(0,1, by = .5),
                     labels = c("0", ".5", "1")) +
  coord_cartesian(#xlim = c(-0.01,1.01), 
                  expand = 0,
                  ylim = c(1,25))+
  ggtitle("Learn Paired") +
  theme_classic() +
  theme(text = element_text(family = "Arial",
                            size = text_size),
        axis.title = element_blank(),
        axis.text.y = element_text(vjust = -.5),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5,vjust = 0))
  


# model003_df %>%
#   select(starts_with("learning_paired[1]")) %>% 
#   pivot_longer(cols = everything()) %>% 
#   mutate(name = factor(name,
#                        levels = unique(name))) %>% 
#   ggplot() +
#   geom_histogram(aes(x = value),bins = 100)

# learning_paired_plot + loo_r2_plot

learning_unpaired_plot <- 
  model003_df %>%
  select(starts_with("learning_unpaired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot(aes(x = value, y = factor(name, levels = rev(levels(name))))) +
  geom_hline(yintercept = c(1:24),
             linewidth = par_axis_line_thickness,
             linetype = "dashed",
             color = "gray")+
  # geom_density_ridges(aes(height = after_stat(scaled)),
  #                     fill = "black",
  #                     scale = .9,
  #                     stat = "density",
  #                     # rel_min_height = 0.01,
  #                     # from = 0,
  #                     # to = 1,
  #                     panel_scaling = FALSE) +
  stat_density_ridges(fill = "black",
                      color = "black",
                      scale = 1,
                      rel_min_height = 0.01,
                      panel_scaling = FALSE,
                      from = 0,
                      to = 1) +
  scale_y_discrete(labels = paste0("Participant ", 
                                   gaborgen_stan_list$n_participants:1)) +
  scale_x_continuous(breaks = seq(0,1, by = .5),
                     labels = c("0", ".5", "1")) +
  coord_cartesian(#xlim = c(-0.01,1.01), 
                  expand = 0,
                  ylim = c(1,25))+
  ggtitle("Learn Unpaired") +
  theme_classic() +
  theme(text = element_text(family = "Arial",
                              size = text_size),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5,vjust = 0))

# learning_paired_plot + learning_unpaired_plot
# learning_unpaired_plot + loo_r2_plot

### CS+ asso####

model003_fit_meta_data <- model003_fit$metadata()

CSP_parameters <- model003_fit_meta_data$model_params[
  str_detect(model003_fit_meta_data$model_params, "CSP_ass|scaling\\[")] 

CSP_df <- posterior::as_draws_df(
  model003_fit$draws(variables = CSP_parameters))

number_of_samples_to_plot <- 4000

set.seed(0)
samples_to_plot <- sample(1:nrow(CSP_df),
                          size = number_of_samples_to_plot,
                          replace = F)



CSP_plot_df <- CSP_df %>% 
  pivot_longer(starts_with("CSP_ass")) %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  mutate(participant = rep(gaborgen_stan_list$participant, 
                           nrow(CSP_df)),
         trial = rep(gaborgen_stan_list$trial, 
                     nrow(CSP_df)),
         cue = rep(gaborgen_stan_list$cue, 
                   nrow(CSP_df))) %>% 
  group_by(participant, trial) %>% 
  mutate(mean_value = mean(value),
         median_value = median(value),.before = 1) %>% 
  ungroup() %>% 
  filter(.draw %in% samples_to_plot)

CSP_plot <- CSP_plot_df %>%
  filter(#participant %in% c(1,2,3),
    trial >= 30) %>% 
  ggplot() +
  geom_line(aes(x = trial, 
                y = value, 
                group = .draw),
            linewidth = line_width_csp,
            alpha = line_alpha_csp) +
  geom_line(aes(x = trial, 
                y = median_value),
            linewidth = line_width_csp_average,
            # alpha = line_alpha_csp,
            color = "red"
  ) +
  scale_x_continuous(breaks = seq(40,180, by = 20),
                     name = "Trial") +
  coord_cartesian(ylim = c(-0.05,1.05),expand = F) +
  theme_classic() +
  facet_grid(participant ~., 
             # scales = "free_y",
             space = "free") +
  ggtitle("Associate Value Posterior Draws and Median") +
  theme(text = element_text(family = "Arial",
                            size = text_size),
        strip.background = element_blank(),
        strip.text = element_blank(),
        # strip.text.y = element_blank(),
        # panel.spacing = unit(0, "lines"),
        panel.spacing = unit(0.15, "lines"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5,vjust = 0.1))



# layout_grid <- c('
# ABCCCC
# ')
# 
# learning_paired_plot +
#   learning_unpaired_plot +
#   CSP_plot +
# plot_annotation(title = expression("Model 3 Estimated Learning Rates, Associate Value, and Cross-Validation R"^2),
#                 theme = theme(
#                   plot.title = element_text(family = "Arial",
#                                             size = 25,
#                                             color = "black",
#                                             hjust = 0.5,
#                                             face = "bold")
#                 ))   +
#   plot_layout(
#     design = layout_grid)
# 
# 
# ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure4.png",
#        dpi = 300, units = "in", 
#        height = 5, width = 5,
#        scale = 2.75)


### R-squared####

loo_R2_cmdstanr <- function(
    observed_data,     # numeric vector of length N with observed outcomes
    predicted_means,   # numeric matrix (draws x N) of posterior E[y] for each obs
    log_lik_matrix     # numeric matrix (draws x N) of log-likelihood for each obs
) {
  # Make sure your inputs align:
  #   - predicted_means and log_lik_matrix each have the same number of rows = total MCMC draws
  #   - both have the same number of columns = N (length(observed_data))
  #   - observed_data is length N
  # browser()
  if (!is.vector(observed_data)) {
    stop("`observed_data` must be a numeric vector.")
  }
  if (!(is.matrix(predicted_means) && is.matrix(log_lik_matrix))) {
    stop("`predicted_means` and `log_lik_matrix` must be numeric matrices.")
  }
  if (nrow(predicted_means) != nrow(log_lik_matrix)) {
    stop("Row counts of `predicted_means` and `log_lik_matrix` differ.")
  }
  if (ncol(predicted_means) != length(observed_data)) {
    stop("Number of columns in `predicted_means` != length(observed_data).")
  }
  if (any(dim(predicted_means) != dim(log_lik_matrix))) {
    stop("`predicted_means` and `log_lik_matrix` must have the same dimensions.")
  }
  
  # Number of draws (M) and observations (N)
  M <- nrow(predicted_means)
  N <- length(observed_data)
  
  # Convert log_lik_matrix to an exponentiated form for relative_eff()
  # You also need a vector of chain IDs. If your total draws M is from, say, 4 chains,
  # each with (M/4) post-warmup draws, you can do something like:
  chains <- 8  # <-- replace as needed
  draws_per_chain <- M / chains
  chain_id <- rep(seq_len(chains), each = draws_per_chain)
  
  # r_eff requires the 'loo' package
  r_eff <- loo::relative_eff(
    x        = exp(log_lik_matrix),
    chain_id = chain_id
  )
  
  # Build a PSIS object
  psis_object <- loo::psis(
    log_ratios = -log_lik_matrix,
    r_eff      = r_eff
  )
  
  # E_loo() gives the PSIS-LOO expectation of ypred
  #   i.e. the "leave-one-out prediction" for each observation.
  ypred_loo <- loo::E_loo(
    x           = predicted_means,
    psis_object = psis_object,
    log_ratios  = -log_lik_matrix
  )$value
  
  # Compute the LOO residuals
  eloo <- ypred_loo - observed_data
  
  # Now we do the Dirichlet resampling approach to get a distribution for R²
  # For details, see: https://avehtari.github.io/bayes_R2/bayes_R2.html
  # (You need `rdirichlet` or `rudirichlet` from the 'MCMCpack' or 'extraDistr' packages.)
  
  # Number of resamples:
  nsim <- 4000
  
  # Sample Dirichlet weights for N items
  # e.g. from library("extraDistr") => `rdiri(n, alpha)`
  # or library("MCMCpack") => `rdirichlet(n, alpha)`
  # For uniform alpha=1 => each draw is Dirichlet(1,...,1)
  
  if (!requireNamespace("extraDistr", quietly = TRUE)) {
    stop("Package 'extraDistr' must be installed for Dirichlet draws (or use another method).")
  }
  
  rd <- extraDistr::rdirichlet(nsim, rep(1, N))  # shape: (nsim x N)
  
  # We compute the variance of the observed data, y:
  #   var(y) ~ sum(weights * y^2) - [sum(weights * y)]^2
  # multiplied by (N/(N-1)) correction
  # rowSums() is for each row of rd
  vary <- (
    rowSums(sweep(rd, 2, observed_data^2, FUN = "*")) -
      rowSums(sweep(rd, 2, observed_data, FUN = "*"))^2
  ) * (N / (N - 1))
  
  # The variance of the LOO residuals, e_loo
  vareloo <- (
    rowSums(sweep(rd, 2, eloo^2, FUN = "*")) -
      rowSums(sweep(rd, 2, eloo, FUN = "*"))^2
  ) * (N / (N - 1))
  
  # The LOO-based R²
  looR2 <- 1 - (vareloo / vary)
  
  # Clip to [-1, 1]
  looR2[looR2 < -1] <- -1
  looR2[looR2 > 1]  <- 1
  
  return(looR2)
}

model001_log_lik_array <- model001_fit$draws("log_lik", format = "draws_matrix") 
model003_log_lik_array <- model003_fit$draws("log_lik", format = "draws_matrix") 
model008_log_lik_array <- model008_fit$draws("log_lik", format = "draws_matrix") 

model001_predicted_means <- model001_fit$draws("mu_pred", format = "draws_matrix")
model003_predicted_means <- model003_fit$draws("mu_pred", format = "draws_matrix")
model008_predicted_means <- model008_fit$draws("mu_pred", format = "draws_matrix")

loo_r2_df <- data.frame("participant" = integer(),
                        "model" = character(),
                        "loo_r2" = numeric())

for(current_par in 1:gaborgen_stan_list$n_participants){
  
  # 1) Extract your observed_data (y) from your original data or from the same data list
  current_observed_data <- 
    gaborgen_stan_list$amplitude[
      gaborgen_stan_list$participant[
        gaborgen_stan_list$indices_observed] == current_par]
  
  # 2) Extract the log_lik draws:
  #    For this, your Stan model must have a 'log_lik[...]' in the parameters or generated quantities
  # shape: (#draws) x (#observations) 
  # or possibly (#draws) x something else. Check dimension.
  current_model001_log_lik_array <- model001_log_lik_array[,
                                                           gaborgen_stan_list$participant[
                                                             gaborgen_stan_list$indices_observed] == current_par]
  
  current_model003_log_lik_array <- model003_log_lik_array[,
                                                           gaborgen_stan_list$participant[
                                                             gaborgen_stan_list$indices_observed] == current_par]

  current_model008_log_lik_array <- model008_log_lik_array[,
                                                           gaborgen_stan_list$participant[
                                                             gaborgen_stan_list$indices_observed] == current_par]
  
  # 3) Extract the posterior prediction draws:
  #    Suppose your Stan model has a parameter or generated quantity called 'y_hat' or 'mu_pred'
  current_model001_predicted_means <- model001_predicted_means[,
                                                               gaborgen_stan_list$participant[
                                                                 gaborgen_stan_list$indices_observed] == current_par]
  
  current_model003_predicted_means <- model003_predicted_means[,
                                                               gaborgen_stan_list$participant[
                                                                 gaborgen_stan_list$indices_observed] == current_par]

  current_model008_predicted_means <- model008_predicted_means[,
                                                               gaborgen_stan_list$participant[
                                                                 gaborgen_stan_list$indices_observed] == current_par]
  
  
  
  # shape: (#draws) x (#observations)
  
  # 4) Call your LOO-based R² function
  current_model001_loo_r2 <- loo_R2_cmdstanr(
    observed_data   = current_observed_data,
    predicted_means = current_model001_predicted_means,
    log_lik_matrix  = current_model001_log_lik_array
  )
  
  current_model003_loo_r2 <- loo_R2_cmdstanr(
    observed_data   = current_observed_data,
    predicted_means = current_model003_predicted_means,
    log_lik_matrix  = current_model003_log_lik_array
  )

  current_model008_loo_r2 <- loo_R2_cmdstanr(
    observed_data   = current_observed_data,
    predicted_means = current_model008_predicted_means,
    log_lik_matrix  = current_model008_log_lik_array
  )
  
  current_loo_df_entry <- rbind(
    data.frame(participant = current_par,
               model = "adapt",
               loo_r2 = current_model008_loo_r2),
    data.frame(participant = current_par,
               model = "block_by_cue",
               loo_r2 = current_model001_loo_r2),
    data.frame(participant = current_par,
               model = "learning",
               loo_r2 = current_model003_loo_r2))
  
  loo_r2_df <- rbind(loo_r2_df,
                     current_loo_df_entry)
}

### Plot loo R^2####
loo_r2_plot <-
loo_r2_df %>% 
  ggplot(aes(x = loo_r2, 
             y = factor(participant, levels = rev(unique(participant))),
             color = model,
             fill = model)) +
  geom_vline(xintercept = 0,
             linewidth = par_axis_line_thickness,
             linetype = "dashed",) +
  geom_hline(yintercept = c(1:24),
             linewidth = par_axis_line_thickness,
             linetype = "dashed",
             color = "gray")+
  geom_density_ridges(aes(height = after_stat(scaled)),
                      # stat_density_ridges(aes(height = after_stat(scaled)),
                      alpha =.2,
                      scale = .9,
                      stat = "density",
                      rel_min_height = 0.01,
                      panel_scaling = FALSE) +
  annotate(geom = "text", 
           x = .4, 
           y = 3.5,
           label = "Model 1",
           family = "Arial",
           color = "green1",
           lineheight = 0.8, # Adjust this value to decrease spacing
           size = 8)+
  annotate(geom = "text", 
           x = .4, 
           y = 2.5,
           label = "Model 2",
           family = "Arial",
           color = "red1",
           lineheight = 0.8, # Adjust this value to decrease spacing
           size = 8)+
  annotate(geom = "text", 
           x = .4, 
           y = 1.5,
           label = "Model 3",
           family = "Arial",
           color = "blue1",
           lineheight = 0.8, # Adjust this value to decrease spacing
           size = 8)+
  scale_y_discrete(labels = paste0("Participant ", 
                                   gaborgen_stan_list$n_participants:1)) +
  scale_color_manual(values = c("green1","red1", "blue1")) +
  scale_fill_manual(values = c("green1","red1", "blue1")) +
  # scale_x_continuous(breaks = seq(0,1, by = .5),
  #                    labels = c("0", ".5", "1")) +
  coord_cartesian(#xlim = c(-0.01,1.01), 
                  expand = 0,
                  ylim = c(1,25))+
  ggtitle(expression("LOO Adjusted R"^2)) +
  theme_classic() +
  theme(text = element_text(family = "Arial",
                            size = text_size),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = .5,vjust = 0))




layout_grid <- c('
AABBCCCCCCCDDD
')

# common_theme <- theme(
#   plot.title.position = "plot",        # Title is positioned relative to the entire plot, not just the panel
#   plot.margin = margin(t = 5, r = 5, b = 5, l = 5)  # Same margin on all sides
# )


learning_paired_plot +
  learning_unpaired_plot + 
  CSP_plot +
  loo_r2_plot +
# free(learning_paired_plot) +
#   free(learning_unpaired_plot) + 
#   free(CSP_plot) +
#   free(loo_r2_plot) +
  plot_annotation(title = expression("Model 3 Estimated Learning Rates, Associate Value, and Cross-Validation R"^2),
                  theme = theme(
                    plot.title = element_text(family = "Arial",
                                              size = 25,
                                              color = "black",
                                              hjust = 0.5,
                                              vjust = -0.1,
                                              face = "bold")
                  ))   +
  plot_layout(
    design = layout_grid)


ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure3.png",
       dpi = 300, units = "in", 
       height = 4.5, width = 6,
       scale = 2.2)

## Effect of cue ####
fig_cue_scaling_text <- 22


bcue_df <- model001_df %>%
  select(starts_with("bcue")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(block = case_when(grepl(pattern = "\\[1", x = name) ~ "Habituation",
                           grepl(pattern = "\\[2", x = name) ~ "Acquisition #1",
                           grepl(pattern = "\\[3", x = name) ~ "Acquisition #2",
                           grepl(pattern = "\\[4", x = name) ~ "Extinction"),
         cue = case_when(grepl(pattern = "1]", x = name) ~ "CS+",
                         grepl(pattern = "2]", x = name) ~ "GS1",
                         grepl(pattern = "3]", x = name) ~ "GS2",
                         grepl(pattern = "4]", x = name) ~ "GS3")) %>% 
  mutate(block = factor(block, 
                        levels = c("Habituation",
                                   "Acquisition #1",
                                   "Acquisition #2",
                                   "Extinction")))


bcue_plot <- bcue_df %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0), linewidth = 2) +
  geom_density_ridges(aes(x = value, 
                          y = factor(block, levels = rev(levels(block))),
                          color = cue,
                          fill = cue),
                      scale = .95,
                      linewidth = 2,
                      alpha = .1) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP", breaks = seq(-.75, .75, by = .25)) +
  coord_cartesian(xlim = c(-.55,.55),
                  ylim = c(1,4.9),
                  expand = F) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Model 2: Effect of Cue by Block") + 
  theme_classic() +
  theme(text = element_text(size = fig_cue_scaling_text,family = "Arial"),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = -0.2, angle = 90),
        plot.title = element_text(hjust = .5),
        legend.position = "none")

bcue_plot
  
  
scaling_plot <-  model003_df %>%
  select(starts_with("scaling[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1.5) +
  geom_density(aes(x = value, 
                   color = name,
                   fill = name),
               linewidth = 2,
               alpha = .1) +
  annotate("text",
           x = 1.05, 
           y = 4.5,
           label = "CS+",
           family = "Arial",
           color = cue_color[1],
           size = 13) +
  annotate("text",
           x = 1.05, 
           y = 4.5 - .9,
           label = "GS1",
           family = "Arial",
           color = cue_color[2],
           size = 13) +
  annotate("text",
           x = 1.05, 
           y = 4.5 - .9 - .9,
           label = "GS2",
           family = "Arial",
           color = cue_color[3],
           size = 13) +
  annotate("text",
           x = 1.05, 
           y = 4.5 - .9 - .9 - .9,
           label = "GS3",
           family = "Arial",
           color = cue_color[4],
           size = 13) +
  coord_cartesian(ylim = c(0, 5.1),
                  xlim = c(-.575, 1.25),
                  expand = F) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP") +
  ggtitle("Model 3: Associate Value Effect on Cue ssVEP") + 
  theme_classic() +
  theme(text = element_text(size = fig_cue_scaling_text,family = "Arial"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

### CS+ scaling####
trials_to_plot <-1:176
number_of_samples_to_plot <- 100
participants_to_plot <- c(1, 5 , 14)

annotation_df <- tibble(
  participant = c(1, 5, 14),
  x = 17.5,
  y = 0.5,
  label = c("Participant 1", "Participant 5", "Participant 14")
)

set.seed(0)
samples_to_plot <- sample(CSP_plot_df$.draw, # further subsetting the 4,000 samples left
                          size = number_of_samples_to_plot,
                          replace = F)

# take a while to run
scaling_par_trial_df <- CSP_plot_df %>%
  filter(.draw %in% samples_to_plot) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(CSP_pred = value * `scaling[1]`) %>%
  mutate(GS1_pred = value * `scaling[2]`) %>%
  mutate(GS2_pred = value * `scaling[3]`) %>%
  mutate(GS4_pred = value * `scaling[4]`) %>%
  select(contains("pred"), ".draw", "trial","participant") %>%
  pivot_longer(contains("pred")) %>%
  mutate(name = factor(name, levels = unique(name)))

scaling_par_trial_avg_df <- CSP_df %>%
  pivot_longer(starts_with("CSP_ass")) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(participant = rep(gaborgen_stan_list$participant,
                           nrow(CSP_df)),
         trial = rep(gaborgen_stan_list$trial,
                     nrow(CSP_df)),
         cue = rep(gaborgen_stan_list$cue,
                   nrow(CSP_df))) %>%
  group_by(participant, trial) %>%
  mutate(mean_value = mean(value),
         median_value = median(value),.before = 1) %>%
  ungroup() %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(CSP_pred = value * `scaling[1]`) %>%
  mutate(GS1_pred = value * `scaling[2]`) %>%
  mutate(GS2_pred = value * `scaling[3]`) %>%
  mutate(GS4_pred = value * `scaling[4]`) %>%
  select(contains("pred"), ".draw", "trial","participant") %>%
  pivot_longer(contains("pred")) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  group_by(participant,trial, name) %>%
  reframe(mean_value = mean(value))

CSP_cue_trial_plot <- scaling_par_trial_df %>% 
  filter(participant %in% participants_to_plot) %>% 
  # Old data management, this didn't treat the matched samples right
  # CSP_df %>% 
  # # filter(participant %in% participants_to_plot) %>% doesn't work yet
  # pivot_longer(starts_with("CSP_ass")) %>% 
  # mutate(name = factor(name, levels = unique(name))) %>% 
  # mutate(participant = rep(gaborgen_stan_list$participant,
  #                          nrow(CSP_df)),
  #        trial = rep(gaborgen_stan_list$trial, 
  #                    nrow(CSP_df)),
  #        cue = rep(gaborgen_stan_list$cue, 
  #                  nrow(CSP_df))) %>% 
  # group_by(.draw, name) %>% 
  # mutate(CSP_pred = value * model003_df[model003_df$.draw == .draw,]$`scaling[1]`) %>% 
  # mutate(GS1_pred = value * model003_df[.draw,]$`scaling[2]`) %>% 
  # mutate(GS2_pred = value * model003_df[.draw,]$`scaling[3]`) %>% 
  # mutate(GS3_pred = value * model003_df[.draw,]$`scaling[4]`) %>% 
  # select(contains("pred"), ".draw", "trial","participant") %>% 
  # pivot_longer(contains("pred")) %>% 
  # mutate(name = factor(name, levels = unique(name))) %>% 
  # filter(.draw %in% samples_to_plot,
  #        participant %in% participants_to_plot) %>% 
  ggplot(aes(x = trial, 
             y = value,
             color = name)) +
  geom_line(aes(
    group = interaction(.draw, name)),
    linewidth = .75,
    alpha = .1) +
  geom_line(data = scaling_par_trial_avg_df %>% 
              filter(participant %in% participants_to_plot),
            aes(x = trial, 
                y = mean_value,
                group = name),
            color = "black",
            linewidth = 1.85) +
  geom_line(data = scaling_par_trial_avg_df %>% 
              filter(participant %in% participants_to_plot),
            aes(x = trial, 
                y = mean_value, 
                color = name),
            linewidth = 1.2) +
  # stat_summary(aes(group = name),
  #              fun = median,
  #              color = "black",
  #              geom = "line",
  #              linewidth = 1.85) +
  # stat_summary(aes(group = name),
  #              fun = median,
  #              geom = "line",
  #              linewidth = 1.2) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 7,
            family = "Arial",
            color = "black") +
  # annotate("text",
  #          x = 20, 
  #          y = .5,
  #          participant = 1,
  #          label = "Participant 1",
  #          family = "Arial",
  #          color = "black",
  #          size = 8) +
  # annotate("text",
  #          x = 20, 
  #          y = .5,
  #          participant = 6,
  #          label = "Participant 6",
  #          family = "Arial",
  #          color = "black",
  #          size = 8) +
  # annotate("text",
  #          x = 20, 
  #          y = .5,
  #          participant = 14,
  #          label = "Participant 14",
  #          family = "Arial",
  #          color = "black",
  #          size = 8) +
  scale_x_continuous(breaks = seq(0,170, by = 20),
                     name = "Trial") +
  scale_y_continuous(name = "Δ Z-Scored ssVEP",
                     breaks = seq(-.25, .75,by = .5)) +
  scale_color_manual(values = cue_color) +
  coord_cartesian(xlim = c(0, 176), 
                  expand = F,
                  ylim = c(-.51,1.1)) +
  theme_classic() +
  facet_grid(participant ~.,
             # scales = "free_y",
             space = "free") +
  ggtitle("Model 3: Associate Value Cue Effect Over Trials") +
  theme(text = element_text(family = "Arial",
                            size = fig_cue_scaling_text),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        legend.position = "none",
        # axis.line.y = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5, vjust = 0))







layout_grid <- c('
abb
acc
acc
')

bcue_plot +
scaling_plot+
CSP_cue_trial_plot +
  plot_annotation(title = "Model 2 vs Model 3 Cue Predictions",
                  theme = theme(
                    plot.title = element_text(family = "Arial",
                                              size = 30,
                                              color = "black",
                                              hjust = 0.5,
                                              vjust = -0.1,
                                              face = "bold")
                  ))   +
  plot_layout(
    design = layout_grid)


ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure4.png",
       dpi = 300, units = "in", 
       height = 3.5, width = 6,
       scale = 3)




# GM heat map, data recovery####
color_z_breaks <- seq(-.75,.75,by = .1)
color_z_breaks_minus_adapt <- seq(-.125,.125,by = .025)
heat_text_size <- 15
label_text_size <- 3.2
vline_linewidth <- .75
vline_linetype <- "dashed"

box_ymax <- 7.6
box_ymin <- 7


Oz_fft_df <- Oz_fft_df %>% 
  group_by(participant) %>% 
  mutate(zamp = scale(amplitude_15Hz_fft) %>% as.vector()) %>% 
  ungroup()

GM_cue_trial_df <- Oz_fft_df %>% 
  group_by(trial, cue) %>%
  reframe(mean_zamp = mean(zamp))

median(model008_df$intercept_average + (model008_df$fatigue_average * 1))

mu_pred_labels <- paste0("mu_pred[",1:gaborgen_stan_list$n_observations,"]")

intercept_labels <- paste0("intercept[",1:gaborgen_stan_list$n_participants,"]")

fatigue_labels <- paste0("fatigue[",1:gaborgen_stan_list$n_participants,"]")

fatigue_dfmod1 <- model008_fit$draws(variables = fatigue_labels,
                                 format = "df")
intercept_dfmod1 <- model008_fit$draws(variables = intercept_labels,
                                 format = "df")

fatigue_dfmod2 <- model001_fit$draws(variables = fatigue_labels,
                                 format = "df")
intercept_dfmod2 <- model001_fit$draws(variables = intercept_labels,
                                 format = "df")

fatigue_dfmod3 <- model003_fit$draws(variables = fatigue_labels,
                                 format = "df")
intercept_dfmod3 <- model003_fit$draws(variables = intercept_labels,
                                 format = "df")


mu_pred_dfmod3 <- model003_fit$draws(variables = mu_pred_labels,
                                 format = "df")


mu_pred_long_avg_infomod3 <- mu_pred_dfmod3 %>% 
  pivot_longer(starts_with("mu_p")) %>% 
  group_by(name) %>% 
  reframe(mean_zamp = mean(value)) %>%
  mutate(index = as.numeric(str_extract(name, "\\d+"))) %>%  # Extract numeric part
  arrange(index) %>%  # Arrange by extracted numeric part
  select(-index) %>%   # Remove temporary column if not needed
  mutate(trial_in_study = gaborgen_stan_list$cue_trial_count[
                            gaborgen_stan_list$indices_observed],
         cue = gaborgen_stan_list$cue[
                 gaborgen_stan_list$indices_observed],
         trial = gaborgen_stan_list$trial[
           gaborgen_stan_list$indices_observed],
         participant = gaborgen_stan_list$participant[
           gaborgen_stan_list$indices_observed]) %>% 
  mutate(mean_zamp_minus_adaptation = mean_zamp - 
           ((intercept_dfmod3[,paste0("intercept[",participant,"]")] %>% 
              unlist() %>% mean()) +
              ((fatigue_dfmod3[,paste0("fatigue[",participant,"]")] %>% 
                 unlist() %>% mean())*trial)
           ))

mu_pred_dfmod2 <- model001_fit$draws(variables = mu_pred_labels,
                                 format = "df")


mu_pred_long_avg_infomod2 <- mu_pred_dfmod2 %>% 
  pivot_longer(starts_with("mu_p")) %>% 
  group_by(name) %>% 
  reframe(mean_zamp = mean(value)) %>%
  mutate(index = as.numeric(str_extract(name, "\\d+"))) %>%  # Extract numeric part
  arrange(index) %>%  # Arrange by extracted numeric part
  select(-index) %>%   # Remove temporary column if not needed
  mutate(trial_in_study = gaborgen_stan_list$cue_trial_count[
                            gaborgen_stan_list$indices_observed],
         cue = gaborgen_stan_list$cue[
                 gaborgen_stan_list$indices_observed],
         trial = gaborgen_stan_list$trial[
           gaborgen_stan_list$indices_observed],
         participant = gaborgen_stan_list$participant[
           gaborgen_stan_list$indices_observed]) %>% 
  mutate(mean_zamp_minus_adaptation = mean_zamp - 
           ((intercept_dfmod2[,paste0("intercept[",participant,"]")] %>% 
               unlist() %>% mean()) +
              ((fatigue_dfmod2[,paste0("fatigue[",participant,"]")] %>% 
                  unlist() %>% mean())*trial)
           ))

mu_pred_dfmod1 <- model008_fit$draws(variables = mu_pred_labels,
                                 format = "df")


mu_pred_long_avg_infomod1 <- mu_pred_dfmod1 %>% 
  pivot_longer(starts_with("mu_p")) %>% 
  group_by(name) %>% 
  reframe(mean_zamp = mean(value)) %>%
  mutate(index = as.numeric(str_extract(name, "\\d+"))) %>%  # Extract numeric part
  arrange(index) %>%  # Arrange by extracted numeric part
  select(-index) %>%   # Remove temporary column if not needed
  mutate(trial_in_study = gaborgen_stan_list$cue_trial_count[
                            gaborgen_stan_list$indices_observed],
         cue = gaborgen_stan_list$cue[
                 gaborgen_stan_list$indices_observed],
         trial = gaborgen_stan_list$trial[
           gaborgen_stan_list$indices_observed],
         participant = gaborgen_stan_list$participant[
           gaborgen_stan_list$indices_observed]) %>% 
  mutate(mean_zamp_minus_adaptation = mean_zamp - 
           ((intercept_dfmod1[,paste0("intercept[",participant,"]")] %>% 
               unlist() %>% mean()) +
              ((fatigue_dfmod1[,paste0("fatigue[",participant,"]")] %>% 
                  unlist() %>% mean())*trial)
           ))

GM_mu_pred_holdmod2 <- mu_pred_long_avg_infomod2 %>% 
  group_by(cue,trial_in_study) %>%
  reframe(mean_zamp = mean(mean_zamp),
          mean_zamp_minus_adaptation = mean(mean_zamp_minus_adaptation)) %>% 
  mutate(cue = cue %>% as.character())

GM_mu_pred_duplicatesmod2 <- GM_mu_pred_holdmod2 %>% 
  filter(cue %in% c(2:4)) %>% 
  mutate(cue = paste0(cue, "_duplicate"))

GM_mu_pred_holdmod3 <- mu_pred_long_avg_infomod3 %>% 
  group_by(cue,trial_in_study) %>%
  reframe(mean_zamp = mean(mean_zamp),
          mean_zamp_minus_adaptation = mean(mean_zamp_minus_adaptation)) %>% 
  mutate(cue = cue %>% as.character())

GM_mu_pred_duplicatesmod3 <- GM_mu_pred_holdmod3 %>% 
  filter(cue %in% c(2:4)) %>% 
  mutate(cue = paste0(cue, "_duplicate"))


hold <- Oz_fft_df  %>% 
  mutate(mean_zamp_minus_adaptation = zamp - 
           ((intercept_dfmod1[,paste0("intercept[",stan_par_id,"]")] %>% 
               unlist() %>% mean()) +
              ((fatigue_dfmod1[,paste0("fatigue[",stan_par_id,"]")] %>% 
                  unlist() %>% mean())*trial)
           )) %>% 
  mutate(trial_in_study = gaborgen_stan_list$cue_trial_count[
    gaborgen_stan_list$indices_observed
  ]) %>% 
  group_by(cue,trial_in_study) %>%
  reframe(mean_zamp = mean(zamp),
          mean_zamp_minus_adaptation = mean(mean_zamp_minus_adaptation)) %>% 
  mutate(cue = cue %>% as.character()) 

hold_2 <- hold %>% 
  filter(cue %in% c("GS1", "GS2", "GS3")) %>%
  mutate(cue = paste0(cue, "_duplicate"))



### plot calls ####
raw_GM_heat <- bind_rows(hold,
                         hold_2) %>% 
  mutate(cue = factor(cue, levels = c("GS3",
                                      "GS2",
                                      "GS1",
                                      "CSP",
                                      "GS1_duplicate",
                                      "GS2_duplicate",
                                      "GS3_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  # geom_point(aes(color = mean_zamp), alpha = 0) + #just for color bar
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_contour_filled(aes(z = scales::squish(mean_zamp, 
  #                                            range = c(min(color_z_breaks) +.001, 
  #                                                      max(color_z_breaks) -.001))), 
  #              breaks = color_z_breaks) +
  # geom_contour(aes(z = scales::squish(mean_zamp, 
  #                                     range = c(min(color_z_breaks) +.001, 
  #                                               max(color_z_breaks) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks) +
  geom_contour_filled(aes(z = scales::squish(rolling_mean,
                                             range = c(min(color_z_breaks) +.001,
                                                       max(color_z_breaks) -.001))),
                      breaks = color_z_breaks) +
  geom_contour(aes(z = scales::squish(rolling_mean,
                                      range = c(min(color_z_breaks) +.001,
                                                max(color_z_breaks) -.001))),
               color = "black",
               breaks = color_z_breaks) +
  geom_vline(xintercept = c(8, 8+12,8+12+12), 
             linewidth = vline_linewidth,
             linetype = vline_linetype) +
  # annotate(geom = "text", 
  #          label = "Habituation", 
  #          x = 9/2, 
  #          y = mean(c(box_ymax,box_ymin)),
  #          family = "Arial",
  #          fontface = "bold",
  #          color = "black",
  #          size = label_text_size) +
  # annotate(geom = "text", 
  #          label = "Acquisition #1", 
  #          x = median(9:20), 
  #          y = mean(c(box_ymax,box_ymin)),
  #          family = "Arial",
  #          fontface = "bold",
  #          color = "black",
  #          size = label_text_size) +
  # annotate(geom = "text", 
  #          label = "Acquisition #2", 
  #          x = median(20:32), 
  #          y = mean(c(box_ymax,box_ymin)),
  #          family = "Arial",
  #          fontface = "bold",
  #          color = "black",
  #          size = label_text_size) +
  # annotate(geom = "text", 
  #          label = "Extinction", 
  #          x = median(32:44), 
  #          y = mean(c(box_ymax,box_ymin)),
  #          family = "Arial",
  #          fontface = "bold",
  #          color = "black",
  #          size = label_text_size) +
  annotate(geom = "rect",
           xmin = 1,
           xmax = 44,
           ymin = box_ymin,
           ymax = box_ymax,
           fill = "black") +
  annotate(geom = "text",
           label = "Habituation",
           x = 9/2,
           y = mean(c(box_ymax,box_ymin)),
           family = "Arial",
           fontface = "bold",
           color = "white",
           size = label_text_size) +
  # annotate(geom = "rect",
  #          xmin = median(9:20) - 4.3,
  #          xmax = median(9:20) + 4.3,
  #          ymin = box_ymin,
  #          ymax = box_ymax,
  #          fill = "black") +
  annotate(geom = "text",
           label = "Acquisition #1",
           x = median(9:20),
           y = mean(c(box_ymax,box_ymin)),
           family = "Arial",
           fontface = "bold",
           color = "white",
           size = label_text_size) +
  # annotate(geom = "rect",
  #          xmin = median(20:32) - 4.3,
  #          xmax = median(20:32) + 4.3,
  #          ymin = box_ymin,
  #          ymax = box_ymax,
  #          fill = "black") +
  annotate(geom = "text",
           label = "Acquisition #2",
           x = median(20:32),
           y = mean(c(box_ymax,box_ymin)),
           family = "Arial",
           fontface = "bold",
           color = "white",
           size = label_text_size) +
  # annotate(geom = "rect",
  #          xmin = median(32:44) - 4.3,
  #          xmax = median(32:44) + 4.3,
  #          ymin = box_ymin,
  #          ymax = box_ymax,
  #          fill = "black") +
  annotate(geom = "text",
           label = "Extinction",
           x = median(32:44),
           y = mean(c(box_ymax,box_ymin)),
           family = "Arial",
           fontface = "bold",
           color = "white",
           size = label_text_size) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_x_continuous(name = "Trial",
                     breaks = c(8,
                                8+12,
                                8+12+12,
                                8+12+12+12)) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks), # Gradient breaks
    limits = c(min(color_z_breaks), max(color_z_breaks)),
    breaks = seq(-.75, .75, by = .25),
    # breaks = color_z_breaks,
    # breaks = seq(-.5, 5, by = .5),
    name = "Z-scored ssVEP"
  ) +
  coord_cartesian(ylim = c(1, 7.6), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Grand Mean Per Cue and Trial") +
  theme_classic() +
  theme(text =element_text(size = heat_text_size, family = "Arial"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5,face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right",
        legend.direction = "horizontal")



# 
# raw_GM_heat & 
#   theme(
#     legend.title.position = "top",
#     legend.position = 'bottom',
#     legend.text = element_text(hjust = .5,),
#     legend.ticks = element_blank(),
#     legend.title = element_text(hjust = .5),
#     # legend.margin = margin(t = .5, r = .5, b = .5, l = .5),  # Reduce margin around the legend
#     legend.box.spacing = unit(0.01, "cm"),
#     legend.key.width = unit(2,"in")
#   )


mod2_GM_heat <- bind_rows(GM_mu_pred_holdmod2,
                          GM_mu_pred_duplicatesmod2) %>% 
  mutate(cue = factor(cue, levels = c("4",
                                      "3",
                                      "2",
                                      "1",
                                      "2_duplicate",
                                      "3_duplicate",
                                      "4_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_point(aes(color = mean_zamp), alpha = 0) + #just for color bar
  geom_contour_filled(aes(z = scales::squish(rolling_mean, 
                                             range = c(min(color_z_breaks) +.001, 
                                                       max(color_z_breaks) -.001))), 
                      breaks = color_z_breaks) +
  geom_contour(aes(z = scales::squish(rolling_mean, 
                                      range = c(min(color_z_breaks) +.001, 
                                                max(color_z_breaks) -.001))), 
               color = "black",
               breaks = color_z_breaks) +
  # geom_contour_filled(aes(z = scales::squish(mean_zamp, 
  #                                            range = c(min(color_z_breaks) +.001, 
  #                                                      max(color_z_breaks) -.001))), 
  #                     breaks = color_z_breaks) +
  # geom_contour(aes(z = scales::squish(mean_zamp, 
  #                                     range = c(min(color_z_breaks) +.001, 
  #                                               max(color_z_breaks) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks) +
  geom_vline(xintercept = c(8, 8+12,8+12+12), 
             linewidth = vline_linewidth,
             linetype = vline_linetype) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_x_continuous(name = "Trial",
                     breaks = c(8,
                                8+12,
                                8+12+12,
                                8+12+12+12)) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks), # Gradient breaks
    limits = c(min(color_z_breaks), max(color_z_breaks)),
    breaks = seq(-.75, .75, by = .25),
    # breaks = color_z_breaks,
    # breaks = seq(-.5, 5, by = .5),
    name = "Z-scored ssVEP"
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Model 2: Grand Mean Prediction") +
  theme_classic() +
  theme(text =element_text(size = heat_text_size, family = "Arial"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5,face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right",
        legend.direction = "horizontal")






mod3_GM_heat <- bind_rows(GM_mu_pred_holdmod3,
                          GM_mu_pred_duplicatesmod3) %>% 
  mutate(cue = factor(cue, levels = c("4",
                                      "3",
                                      "2",
                                      "1",
                                      "2_duplicate",
                                      "3_duplicate",
                                      "4_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_point(aes(color = mean_zamp), alpha = 0) + #just for color bar
  geom_contour_filled(aes(z = scales::squish(rolling_mean, 
                                             range = c(min(color_z_breaks) +.001, 
                                                       max(color_z_breaks) -.001))), 
                      breaks = color_z_breaks) +
  geom_contour(aes(z = scales::squish(rolling_mean, 
                                      range = c(min(color_z_breaks) +.001, 
                                                max(color_z_breaks) -.001))), 
               color = "black",
               breaks = color_z_breaks) +
  # geom_contour_filled(aes(z = scales::squish(mean_zamp, 
  #                                            range = c(min(color_z_breaks) +.001, 
  #                                                      max(color_z_breaks) -.001))), 
  #                     breaks = color_z_breaks) +
  # geom_contour(aes(z = scales::squish(mean_zamp, 
  #                                     range = c(min(color_z_breaks) +.001, 
  #                                               max(color_z_breaks) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks) +
  geom_vline(xintercept = c(8, 8+12,8+12+12), 
             linewidth = vline_linewidth,
             linetype = vline_linetype) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_x_continuous(name = "Trial",
                     breaks = c(8,
                                8+12,
                                8+12+12,
                                8+12+12+12)) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks), # Gradient breaks
    limits = c(min(color_z_breaks), max(color_z_breaks)),
    breaks = seq(-.75, .75, by = .25),
    # breaks = color_z_breaks,
    # breaks = seq(-.5, 5, by = .5),
    name = "Z-scored ssVEP"
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Model 3: Grand Mean Prediction") +
  theme_classic() +
  theme(text =element_text(size = heat_text_size, family = "Arial"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5,face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.direction = "horizontal")


## Heat combine save####
layout_grid <- '
aaaaaaaaaaaaaa
aaaaaaaaaaaaaa
aaaaaaaaaaaaaa
aaaaaaaaaaaaaa
aaaaaaaaaaaaaa
bbbbbbbbbbbbbb
bbbbbbbbbbbbbb
bbbbbbbbbbbbbb
bbbbbbbbbbbbbb
bbbbbbbbbbbbbb
cccccccccccccc
cccccccccccccc
cccccccccccccc
cccccccccccccc
cccccccccccccc
ddddddddddddd#
'

raw_GM_heat +
  mod2_GM_heat +
  mod3_GM_heat +
  guide_area() +
  plot_layout(
    design = layout_grid,
    guides = 'collect' # Combine legends from all plots
  ) &
  theme(axis.text = element_text(colour = "black"),
        plot.title = element_text(size = 13.5),
    legend.title.position = "top",
    legend.position = 'bottom',
    legend.justification = "center",
    legend.box.just = "center",
    legend.text = element_text(hjust = .5),
    legend.ticks = element_blank(),
    legend.title = element_text(hjust = .5, face = "bold"),
    legend.margin = margin(t = .5, r = .5, b = .5, l = .5),  # Reduce margin around the legend
    legend.box.spacing = unit(0.01, "cm"),
    legend.key.width = unit(.9,"in")
  )


ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure5.png",
       dpi = 300, units = "in", 
       height = 4, width = 2.5,
       scale = 2)


# missing trial interp figure####

missing_trial_posterior_samples <- 100
set.seed(0)
missing_trial_samples_to_plot <- sample(x = 1:(8*5e3), size = missing_trial_posterior_samples,replace = F)

participants_to_plot <- c(1,14)
# participants_to_plot <- c(1,5,14)

participant_trials <- (1:gaborgen_stan_list$n)[
  gaborgen_stan_list$participant %in% participants_to_plot]

participant_trial_observed <- participant_trials[
  participant_trials %in% gaborgen_stan_list$indices_observed]

participant_trial_missing <- participant_trials[
  participant_trials %in% gaborgen_stan_list$indices_missing]


cue_indices <- gaborgen_stan_list$cue[participant_trials]

participant_amp_array_labels <- paste0("amplitude_all[",
                                       participant_trials,
                                       "]")

amplitude_array_df <- model003_fit$draws(variables = participant_amp_array_labels,
                                         format = "df")



amplitude_array_df_long <- amplitude_array_df %>% 
  filter(.draw %in% missing_trial_samples_to_plot) %>% 
  pivot_longer(starts_with("ampl")) %>% 
  mutate(participant = rep(rep(x = participants_to_plot, 
                               each =  176),
                           missing_trial_posterior_samples)) %>% 
  mutate(trial = rep(rep(x = 1:gaborgen_stan_list$n_trials, 
                         length(participants_to_plot)),
                     missing_trial_posterior_samples)) %>% 
  mutate(cue = rep(x = cue_indices,
                   missing_trial_posterior_samples)) %>% 
  mutate(cue = case_when(cue == 1 ~ "CS+",
                         cue == 2 ~ "GS1",
                         cue == 3 ~ "GS2",
                         cue == 4 ~ "GS3")) %>% 
  mutate(missing = name %in% paste0("amplitude_all[",
                                    participant_trial_missing,
                                    "]")) %>% 
  mutate(single_dot_to_plot = F)

amplitude_array_df_long[
  1:(length(participants_to_plot)*gaborgen_stan_list$n_trials),]$single_dot_to_plot <- T


# participant_mu_array_labels <- paste0("amplitude_all[",
#                                       participant_trials,
#                                       "]")


observed_dot_size <- 3
quasi_spread <- .25
trial_limits <- c(110,125)
anno_spacing <-1.75
# trial_limits <- c(96,125)
# trial_limits <- c(96,145)
# trial_limits <- c(20,176)

amplitude_array_df_long %>% 
  # filter(trial >= trial_limits[1], trial <= trial_limits[2]) %>%
  ggplot() +
  geom_line(data = mu_pred_long_avg_infomod3 %>% 
              filter(participant %in% participants_to_plot) %>% 
              mutate(cue = case_when(cue == "1" ~ "CS+",
                                     cue == "2" ~ "GS1",
                                     cue == "3" ~ "GS2",
                                     cue == "4" ~ "GS3")),
            aes(x = trial, 
                y = mean_zamp,
                group = cue),
            color = "black",
            linewidth = 1.85) +
  geom_line(data = mu_pred_long_avg_infomod3 %>% 
              filter(participant %in% participants_to_plot) %>% 
              mutate(cue = case_when(cue == "1" ~ "CS+",
                                     cue == "2" ~ "GS1",
                                     cue == "3" ~ "GS2",
                                     cue == "4" ~ "GS3")),
            aes(x = trial, 
                y = mean_zamp, 
                color = cue),
            linewidth = 1.2) +
  geom_line(data = . %>% filter(missing == F, single_dot_to_plot),
             aes(x = trial, y = value)) +
  geom_point(data = . %>% filter(missing == F, single_dot_to_plot),
             aes(x = trial, y = value), 
             size = observed_dot_size + .75) +
  geom_point(data = . %>% filter(missing == F, single_dot_to_plot),
             aes(x = trial, y = value, color = cue),
             size = observed_dot_size) +
  geom_quasirandom(data = . %>% filter(missing == T),
                   aes(x = trial, y = value, group = cue),
                   size = .8,
                   color = "black",
                   width = quasi_spread,
                   alpha = 1) +
  geom_quasirandom(data = . %>% filter(missing == T),
                   aes(x = trial, y = value, color = cue),
                   size = .6,
                   width = quasi_spread,
                   alpha = 1) +
  geom_text(data = data.frame(trial = c(110.25,
                                        110.25+anno_spacing,
                                        110.25+anno_spacing+anno_spacing,
                                        110.25+anno_spacing+anno_spacing+anno_spacing),
                              value = c(1.75),
                              participant = c(1,1,1,1),
                              anno_label = c("CS+","GS1","GS2","GS3")),
            aes(x = trial,
                y = value,
                label = anno_label),
            inherit.aes = FALSE,
            color = cue_color,
            family = "Arial",
            size = 9) +
  scale_color_manual(values = cue_color) +
  scale_y_continuous(name = "Z-scored ssVEP",
                     breaks = seq(-2,2, by = 2)) +
  scale_x_continuous(name = "Trial",
                     breaks = seq(trial_limits[1], 175, by = 1)) +
  coord_cartesian(xlim = trial_limits, 
                  ylim = c(-2.5,2.5)
                  ) +
  theme_classic() +
  facet_wrap(~participant, ncol = 1,
             labeller = labeller(
               participant = c("1" = "Participant 1",
                               "14" = "Participant 14")
             )) +
  # facet_grid(.~participant,
  # facet_grid(participant ~.,
  #            # scales = "free_y",
  #            space = "free") +
  ggtitle("Model 3 Missing Trial Interpolation") +
  theme(text = element_text(family = "Arial",
                            size = fig_cue_scaling_text),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        # axis.line.y = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = .5, vjust = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)) # Add this line



ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/Figure6.png",
       dpi = 300, units = "in", 
       height = 2, width = 3,
       scale = 3)







mod2_GM_heat_minus_adapt <- bind_rows(GM_mu_pred_holdmod2,
          GM_mu_pred_duplicatesmod2) %>% 
  mutate(cue = factor(cue, levels = c("4",
                                      "3",
                                      "2",
                                      "1",
                                      "2_duplicate",
                                      "3_duplicate",
                                      "4_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp_minus_adaptation, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_point(aes(color = mean_zamp_minus_adaptation), alpha = 0) + #just for color bar
  geom_contour_filled(aes(z = scales::squish(rolling_mean, 
                                             range = c(min(color_z_breaks_minus_adapt) +.001, 
                                                       max(color_z_breaks_minus_adapt) -.001))), 
                      breaks = color_z_breaks_minus_adapt) +
  geom_contour(aes(z = scales::squish(rolling_mean, 
                                      range = c(min(color_z_breaks_minus_adapt) +.001, 
                                                max(color_z_breaks_minus_adapt) -.001))), 
               color = "black",
               breaks = color_z_breaks_minus_adapt) +
  # geom_contour_filled(aes(z = scales::squish(mean_zamp_minus_adaptation, 
  #                                            range = c(min(color_z_breaks_minus_adapt) +.001, 
  #                                                      max(color_z_breaks_minus_adapt) -.001))), 
  #                     breaks = color_z_breaks_minus_adapt) +
  # geom_contour(aes(z = scales::squish(mean_zamp_minus_adaptation, 
  #                                     range = c(min(color_z_breaks_minus_adapt) +.001, 
  #                                               max(color_z_breaks_minus_adapt) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks_minus_adapt) +
  geom_vline(xintercept = c(8, 8+12,8+12+12)) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks_minus_adapt), # Gradient breaks
    limits = c(min(color_z_breaks_minus_adapt), max(color_z_breaks_minus_adapt)),
    breaks = seq(-.1, .1, by = .1)
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Model 2: Cue by Block") +
  theme_classic()

mod3_GM_heat <- bind_rows(GM_mu_pred_holdmod3,
          GM_mu_pred_duplicatesmod3) %>% 
  mutate(cue = factor(cue, levels = c("4",
                                      "3",
                                      "2",
                                      "1",
                                      "2_duplicate",
                                      "3_duplicate",
                                      "4_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_point(aes(color = mean_zamp), alpha = 0) + #just for color bar
  geom_contour_filled(aes(z = scales::squish(rolling_mean, 
                                             range = c(min(color_z_breaks) +.001, 
                                                       max(color_z_breaks) -.001))), 
                      breaks = color_z_breaks) +
  geom_contour(aes(z = scales::squish(rolling_mean, 
                                      range = c(min(color_z_breaks) +.001, 
                                                max(color_z_breaks) -.001))), 
               color = "black",
               breaks = color_z_breaks) +
  # geom_contour_filled(aes(z = scales::squish(mean_zamp, 
  #                                            range = c(min(color_z_breaks) +.001, 
  #                                                      max(color_z_breaks) -.001))), 
  #                     breaks = color_z_breaks) +
  # geom_contour(aes(z = scales::squish(mean_zamp, 
  #                                     range = c(min(color_z_breaks) +.001, 
  #                                               max(color_z_breaks) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks) +
  geom_vline(xintercept = c(8, 8+12,8+12+12)) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks), # Gradient breaks
    limits = c(min(color_z_breaks), max(color_z_breaks)),
    breaks = seq(-.5, .5, by = .5)
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Model 3: Learning model") +
  theme_classic()

mod3_GM_heat_minus_adapt <- bind_rows(GM_mu_pred_holdmod3,
          GM_mu_pred_duplicatesmod3) %>% 
  mutate(cue = factor(cue, levels = c("4",
                                      "3",
                                      "2",
                                      "1",
                                      "2_duplicate",
                                      "3_duplicate",
                                      "4_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp_minus_adaptation, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_point(aes(color = mean_zamp_minus_adaptation), alpha = 0) + #just for color bar
  geom_contour_filled(aes(z = scales::squish(rolling_mean, 
                                             range = c(min(color_z_breaks_minus_adapt) +.001, 
                                                       max(color_z_breaks_minus_adapt) -.001))), 
                      breaks = color_z_breaks_minus_adapt) +
  geom_contour(aes(z = scales::squish(rolling_mean, 
                                      range = c(min(color_z_breaks_minus_adapt) +.001, 
                                                max(color_z_breaks_minus_adapt) -.001))), 
               color = "black",
               breaks = color_z_breaks_minus_adapt) +
  # geom_contour_filled(aes(z = scales::squish(mean_zamp_minus_adaptation, 
  #                                            range = c(min(color_z_breaks_minus_adapt) +.001, 
  #                                                      max(color_z_breaks_minus_adapt) -.001))), 
  #                     breaks = color_z_breaks_minus_adapt) +
  # geom_contour(aes(z = scales::squish(mean_zamp_minus_adaptation, 
  #                                     range = c(min(color_z_breaks_minus_adapt) +.001, 
  #                                               max(color_z_breaks_minus_adapt) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks_minus_adapt) +
  geom_vline(xintercept = c(8, 8+12,8+12+12)) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks_minus_adapt), # Gradient breaks
    limits = c(min(color_z_breaks_minus_adapt), max(color_z_breaks_minus_adapt)),
    breaks = seq(-.1, .1, by = .1)
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Model 3: Learning model") +
  theme_classic()



hold <- Oz_fft_df  %>% 
  mutate(mean_zamp_minus_adaptation = zamp - 
           ((intercept_dfmod1[,paste0("intercept[",stan_par_id,"]")] %>% 
               unlist() %>% mean()) +
              ((fatigue_dfmod1[,paste0("fatigue[",stan_par_id,"]")] %>% 
                  unlist() %>% mean())*trial)
           )) %>% 
  mutate(trial_in_study = gaborgen_stan_list$cue_trial_count[
    gaborgen_stan_list$indices_observed
  ]) %>% 
  group_by(cue,trial_in_study) %>%
  reframe(mean_zamp = mean(zamp),
          mean_zamp_minus_adaptation = mean(mean_zamp_minus_adaptation)) %>% 
  mutate(cue = cue %>% as.character()) 

hold_2 <- hold %>% 
  filter(cue %in% c("GS1", "GS2", "GS3")) %>%
  mutate(cue = paste0(cue, "_duplicate"))




raw_GM_heat <- bind_rows(hold,
          hold_2) %>% 
  mutate(cue = factor(cue, levels = c("GS3",
                                      "GS2",
                                      "GS1",
                                      "CSP",
                                      "GS1_duplicate",
                                      "GS2_duplicate",
                                      "GS3_duplicate"))) %>% 
  group_by(cue) %>%
  mutate(rolling_mean = zoo::rollapply(mean_zamp, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>%
  ungroup() %>%
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  # geom_point(aes(color = mean_zamp), alpha = 0) + #just for color bar
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  # geom_contour_filled(aes(z = scales::squish(mean_zamp, 
  #                                            range = c(min(color_z_breaks) +.001, 
  #                                                      max(color_z_breaks) -.001))), 
  #              breaks = color_z_breaks) +
  # geom_contour(aes(z = scales::squish(mean_zamp, 
  #                                     range = c(min(color_z_breaks) +.001, 
  #                                               max(color_z_breaks) -.001))), 
  #              color = "black",
  #              breaks = color_z_breaks) +
  geom_contour_filled(aes(z = scales::squish(rolling_mean,
                                             range = c(min(color_z_breaks) +.001,
                                                       max(color_z_breaks) -.001))),
               breaks = color_z_breaks) +
  geom_contour(aes(z = scales::squish(rolling_mean,
                                      range = c(min(color_z_breaks) +.001,
                                                max(color_z_breaks) -.001))),
               color = "black",
               breaks = color_z_breaks) +
  geom_vline(xintercept = c(8, 8+12,8+12+12)) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks), # Gradient breaks
    limits = c(min(color_z_breaks), max(color_z_breaks)),
    breaks = seq(-.5, 5, by = .5)
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Raw data averaged down to trial by cue") +
  theme_classic()

raw_GM_heat_minus_adapt <- bind_rows(hold,
          hold_2) %>% 
  mutate(cue = factor(cue, levels = c("GS3",
                                      "GS2",
                                      "GS1",
                                      "CSP",
                                      "GS1_duplicate",
                                      "GS2_duplicate",
                                      "GS3_duplicate"))) %>%
  group_by(cue) %>% 
  mutate(rolling_mean = zoo::rollapply(mean_zamp_minus_adaptation, width = 2, FUN = mean, fill = NA, align = "center", partial = TRUE)) %>% 
  ungroup() %>% 
  mutate(cue_int = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial_in_study, y = cue_int)) +
  geom_point(aes(color = rolling_mean), alpha = 0) + #just for color bar
  geom_contour_filled(aes(z = scales::squish(rolling_mean, 
                                             range = c(min(color_z_breaks_minus_adapt) +.001, 
                                                       max(color_z_breaks_minus_adapt) -.001))), 
                      breaks = color_z_breaks_minus_adapt) +
  geom_contour(aes(z = scales::squish(rolling_mean, 
                                      range = c(min(color_z_breaks_minus_adapt) +.001, 
                                                max(color_z_breaks_minus_adapt) -.001))), 
               color = "black",
               breaks = color_z_breaks_minus_adapt) +
  geom_vline(xintercept = c(8, 8+12,8+12+12)) +
  scale_fill_manual(
    values = colorRampPalette(c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"))(length(color_z_breaks)), # Generate colors automatically
    name = "Mean ZAmp",
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = 1:7,
    labels = c("GS3", "GS2", "GS1", "CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  scale_color_gradientn(
    colors = c("darkblue", "blue", "cyan", "green", "yellow", "red", "darkred"), 
    values = scales::rescale(color_z_breaks_minus_adapt), # Gradient breaks
    limits = c(min(color_z_breaks_minus_adapt), max(color_z_breaks_minus_adapt)),
    breaks = seq(-.1, .1, by = .1)
  ) +
  coord_cartesian(ylim = c(1, 7), 
                  xlim = c(1, 44), 
                  expand = F) +
  ggtitle("Raw data averaged down to trial by cue") +
  theme_classic()

raw_GM_heat /
  mod2_GM_heat /
  mod3_GM_heat

raw_GM_heat_minus_adapt /
  mod2_GM_heat_minus_adapt /
  mod3_GM_heat_minus_adapt


#data recovery 
mu_pred_data_rec <- Oz_fft_df %>% 
  mutate(mu_pred_mod001 = mu_pred_long_avg_infomod1$mean_zamp,
         mu_pred_mod002 = mu_pred_long_avg_infomod2$mean_zamp,
         mu_pred_mod003 = mu_pred_long_avg_infomod3$mean_zamp)

cor.test(mu_pred_data_rec$zamp, mu_pred_data_rec$mu_pred_mod001)
cor.test(mu_pred_data_rec$zamp, mu_pred_data_rec$mu_pred_mod002)
cor.test(mu_pred_data_rec$zamp, mu_pred_data_rec$mu_pred_mod003)

temp_dot_size <- .8
temp_dot_alpha <- 1

(mu_pred_data_rec %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(aes(x = zamp, 
                 y = mu_pred_mod001,
                 color = cue),
             size = temp_dot_size,
             alpha = temp_dot_alpha) +
  scale_color_manual(values = cue_color) +
  coord_cartesian(ylim = c(-1.1,1.1)) +
  theme_classic()) /

(mu_pred_data_rec %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(aes(x = zamp, 
                 y = mu_pred_mod002,
                 color = cue),
             size = temp_dot_size,
             alpha = temp_dot_alpha) +
  scale_color_manual(values = cue_color) +
   coord_cartesian(ylim = c(-1.1,1.1)) +
  theme_classic() )/

(mu_pred_data_rec %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(aes(x = zamp, 
                 y = mu_pred_mod003,
                 color = cue),
             size = temp_dot_size,
             alpha = temp_dot_alpha) +
  scale_color_manual(values = cue_color) +
   coord_cartesian(ylim = c(-1.1,1.1)) +
  theme_classic() )








  
Oz_fft_df %>% 
  mutate(trial_in_block = gaborgen_stan_list$cue_trial_count[
    gaborgen_stan_list$indices_observed
  ]) %>% 
  group_by(cue,trial_in_block) %>%
  reframe(mean_zamp = mean(zamp)) %>% 
  group_by(trial) %>% 
  ungroup() %>% 
  complete(trial, cue) %>%
  group_by(trial) %>% 
  mutate(mean_zamp_minus_adaptation = mean_zamp - median(model008_df$intercept_average + (model008_df$fatigue_average * trial))) %>% 
  group_by(cue) %>% 
  mutate(mean_zamp_m_adapt_moving_avg = zoo::rollapply(
      mean_zamp_minus_adaptation, 
      width = 20, # Use a window of 5 for smoothing
      FUN = mean,
      fill = NA,
      align = "center"
    )
  ) %>% 
  mutate(cue = cue %>% as.integer()) %>% 
  ggplot(aes(x = trial, y = cue, fill = mean_zamp_m_adapt_moving_avg)) +
  geom_raster(interpolate = TRUE) + # Smooth interpolation over both axes
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "green", "yellow", "red"), 
    values = scales::rescale(c(-.2, -0.05, 0, 0.05, .2)), # Gradient breaks
    limits = c(-.2, .2), # Color scale limits
    oob = scales::squish # Squish out-of-bounds values
  ) +
  scale_y_continuous(
    breaks = 1:4, 
    labels = c("CS+", "GS1", "GS2", "GS3") # Map y-axis labels to cues
  ) +
  coord_cartesian(expand = F) +
  theme_classic()



Oz_fft_df %>% 
  group_by(trial, cue) %>%
  reframe(mean_zamp = mean(zamp)) %>% 
  group_by(trial) %>% 
  ungroup() %>% 
  complete(trial, cue) %>%
  group_by(trial) %>% 
  mutate(mean_zamp_minus_adaptation = mean_zamp - median(model008_df$intercept_average + (model008_df$fatigue_average * trial))) %>% 
  mutate(cue = cue %>% as.integer()) %>% 
  # ggplot() +
  ggplot(aes(x = trial, y = cue)) +
  # stat_summary_2d(aes(z = mean_zamp, fill = after_stat(value)), bins = 20) +
  stat_summary_2d(aes(z = mean_zamp_minus_adaptation, fill = after_stat(value)), bins = 25) +
  # geom_raster(aes(x = trial, y = cue, fill = mean_zamp_minus_adaptation),interpolate = TRUE) +
  # geom_raster(aes(x = trial, y = cue, fill = mean_zamp),interpolate = TRUE) +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "green", "yellow", "red"), 
    values = scales::rescale(c(-.2, -0.1, 0, 0.1, .2)), # Rescales the values to 0-1
    limits = c(-.2, .2), # Ensures the color gradient covers the full range
    oob = scales::squish # Handles values outside the limits
  ) +
  # scale_fill_gradient2(low = "blue", mid = "green", high = "red") +
  scale_y_continuous(labels = c("CS+","GS1","GS2","GS3")) +
  coord_cartesian(expand = F) +
  theme_classic()
  # stat_density2d(aes(x = trial,
  #                    y = cue,
  #                    fill = mean_zamp), 
  #                geom = "tile", 
  #                contour = F)

Oz_fft_df %>% 
  group_by(trial, cue) %>%
  reframe(mean_zamp = mean(zamp)) %>% 
  group_by(trial) %>% 
  ungroup() %>% 
  complete(trial, cue) %>%
  group_by(trial) %>% 
  mutate(mean_zamp_minus_adaptation = mean_zamp - median(model008_df$intercept_average + (model008_df$fatigue_average * trial))) %>% 
  mutate(cue = cue %>% as.integer()) %>% 
  # ggplot() +
  ggplot(aes(x = trial, y = cue)) +
  stat_summary_2d(aes(z = mean_zamp, fill = after_stat(value)), bins = 20) +
  # stat_summary_2d(aes(z = mean_zamp_minus_adaptation, fill = after_stat(value)), bins = 25) +
  # geom_raster(aes(x = trial, y = cue, fill = mean_zamp_minus_adaptation),interpolate = TRUE) +
  # geom_raster(aes(x = trial, y = cue, fill = mean_zamp),interpolate = TRUE) +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "green", "yellow", "red"), 
    values = scales::rescale(c(-.5, -0.1, 0, 0.1, .5)), # Rescales the values to 0-1
    limits = c(-.5, .5), # Ensures the color gradient covers the full range
    oob = scales::squish # Handles values outside the limits
  ) +
  # scale_fill_gradient2(low = "blue", mid = "green", high = "red") +
  scale_y_continuous(labels = c("CS+","GS1","GS2","GS3")) +
  coord_cartesian(expand = F) +
  theme_classic()
  # stat_density2d(aes(x = trial,
  #                    y = cue,
  #                    fill = mean_zamp), 
  #                geom = "tile", 
  #                contour = F)
  






bcue_plot <- bcue_df %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, 
                          y = factor(block, levels = rev(levels(block))),
                          color = cue,
                          fill = cue),
                      scale = 1,
                      linewidth = 2,
                      alpha = .1) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP", breaks = seq(-.75, .75, by = .25)) +
  coord_cartesian(xlim = c(-.5,.5),
                  ylim = c(1.55,4.25)) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Effect of Cue Controlling for Adaptation") + 
  theme_classic() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())


scaling_plot <- model003_df %>%
  select(starts_with("scaling[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1.5) +
  geom_density(aes(x = value, 
                   color = name,
                   fill = name),
               linewidth = 2,
               alpha = .1) +
  annotate("text",
           x = 1, 
           y = 4.85,
           label = "CS+",
           family = "Arial",
           color = cue_color[1],
           size = 18) +
  annotate("text",
           x = 1, 
           y = 4.85 - .4,
           label = "GS1",
           family = "Arial",
           color = cue_color[2],
           size = 18) +
  annotate("text",
           x = 1, 
           y = 4.85 - .4 - .4,
           label = "GS2",
           family = "Arial",
           color = cue_color[3],
           size = 18) +
  annotate("text",
           x = 1, 
           y = 4.85 - .4 - .4 - .4,
           label = "GS3",
           family = "Arial",
           color = cue_color[4],
           size = 18) +
  coord_cartesian(ylim = c(0, 5.1),
                  xlim = c(-.575, 1.25),
                  expand = F) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP") +
  ggtitle("Model 2\nAssociate Value Effect on Cue ssVEP") + 
  theme_classic() +
  theme(text = element_text(size = 22,family = "Arial"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")











bayes_R2_residuals <- function(observed_data, 
                               predicted_means, 
                               correct_for_sample_by_obs = NULL) {
  
  if(!is.null(correct_for_sample_by_obs)){
    
    corrected_observed_data <- -1 * sweep(correct_for_sample_by_obs,
                                          2, observed_data)
    
    bayes_residuals <- corrected_observed_data - predicted_means
    
  } else {
    
    bayes_residuals <- -1 * sweep(predicted_means, 2, 
                                  observed_data)
    
  }
  variance_of_predicted_means <- apply(predicted_means, 1, var)
  
  variance_of_residuals <- apply(bayes_residuals, 1, var)
  
  
  bayesian_R_squared <- variance_of_predicted_means / 
    (variance_of_predicted_means + variance_of_residuals)
  
  return(bayesian_R_squared)
}

model001_fit_meta_data <- model001_fit$metadata()

model001_mu_pred <- model001_fit_meta_data$model_params[
  str_detect(model001_fit_meta_data$model_params, "mu_pred")]

model001_mu_pred_mat <- model001_fit$draws(variables = model001_mu_pred,
                                          format = "df") %>% 
  select(starts_with("mu_pred[")) %>% #necessary to remove hidden columns .draw and such
  as.matrix()

model003_fit_meta_data <- model003_fit$metadata()

model003__mu_pred <- model003_fit_meta_data$model_params[
  str_detect(model003_fit_meta_data$model_params, "mu_pred")]

model003_mu_pred_mat <- model003_fit$draws(
  variables = model003__mu_pred,
  format = "df") %>% 
  select(starts_with("mu_pred[")) %>% #necessary to remove hidden columns .draw and such
  as.matrix()

model008_fit_meta_data <- model008_fit$metadata()

model008_fit_mu_pred <- model008_fit_meta_data$model_params[
  str_detect(model008_fit_meta_data$model_params, "mu_pred")]

model008_mu_pred_mat <- model008_fit$draws(variables = model008_fit_mu_pred,
                                  format = "df") %>% 
  select(starts_with("mu_pred[")) %>% #necessary to remove hidden columns .draw and such
  as.matrix()

#call
dim(model001_mu_pred_mat)
dim(model003_mu_pred_mat)
dim(model008_mu_pred_mat)
length(gaborgen_stan_list$amplitude)

# participant filter
current_par <- "1"
current_par <- "11"
current_par <- "12"
current_par <- "19"
current_par <- "21"
model001_mu_pred_current_par <- model001_mu_pred_mat[,current_par ==
gaborgen_stan_list$participant[gaborgen_stan_list$indices_observed]
]
model003_mu_pred_current_par <- model003_mu_pred_mat[,current_par ==
gaborgen_stan_list$participant[gaborgen_stan_list$indices_observed]
]
model008_mu_pred_current_par <- model008_mu_pred_mat[,current_par ==
gaborgen_stan_list$participant[gaborgen_stan_list$indices_observed]
]

amp_current_par <- gaborgen_stan_list$amplitude[
  current_par ==
    gaborgen_stan_list$participant[gaborgen_stan_list$indices_observed]
]

model001_R2_current_par <- bayes_R2_residuals(
  observed_data = amp_current_par,
  predicted_means = model001_mu_pred_current_par)

model001_R2_current_par %>% density() %>% plot()

model003_R2_current_par <- bayes_R2_residuals(
  observed_data = amp_current_par,
  predicted_means = model003_mu_pred_current_par)

model003_R2_current_par %>% density() %>% plot()

(model003_R2_current_par - model001_R2_current_par) %>% 
  hist(breaks = 30)
  # density() %>% 
  # plot()

((model003_R2_current_par - model001_R2_current_par) > 0) %>% 
  sum() / 40000







model003_R2<- bayes_R2_residuals(observed_data = gaborgen_stan_list$amplitude,
                   predicted_means = model003_mu_pred_mat)

model003_R2 %>% density() %>% plot()

model008_R2<- bayes_R2_residuals(observed_data = gaborgen_stan_list$amplitude,
                   predicted_means = model008_mu_pred_mat)




model001_R2<- bayes_R2_residuals(observed_data = gaborgen_stan_list$amplitude,
                   predicted_means = model001_mu_pred_mat)

model001_R2 %>% density() %>% plot()

model003_R2<- bayes_R2_residuals(observed_data = gaborgen_stan_list$amplitude,
                   predicted_means = model003_mu_pred_mat)

model003_R2 %>% density() %>% plot()

model008_R2<- bayes_R2_residuals(observed_data = gaborgen_stan_list$amplitude,
                   predicted_means = model008_mu_pred_mat)

model008_R2 %>% density() %>% plot()

((model001_R2 - model008_R2) > 0) %>% 
  sum() / nrow(model001_mu_pred_mat)

density(model001_R2 - model008_R2) %>% plot()

((model003_R2 - model008_R2) > 0) %>% 
  sum() / nrow(model001_mu_pred_mat)

density(model003_R2 - model008_R2) %>% plot()

((model003_R2 - model001_R2) > 0) %>% 
  sum() / nrow(model001_mu_pred_mat)

density(model003_R2 - model001_R2) %>% plot()

# Loo R-squared####

loo_R2_cmdstanr <- function(
    observed_data,     # numeric vector of length N with observed outcomes
    predicted_means,   # numeric matrix (draws x N) of posterior E[y] for each obs
    log_lik_matrix     # numeric matrix (draws x N) of log-likelihood for each obs
) {
  # Make sure your inputs align:
  #   - predicted_means and log_lik_matrix each have the same number of rows = total MCMC draws
  #   - both have the same number of columns = N (length(observed_data))
  #   - observed_data is length N
  browser()
  if (!is.vector(observed_data)) {
    stop("`observed_data` must be a numeric vector.")
  }
  if (!(is.matrix(predicted_means) && is.matrix(log_lik_matrix))) {
    stop("`predicted_means` and `log_lik_matrix` must be numeric matrices.")
  }
  if (nrow(predicted_means) != nrow(log_lik_matrix)) {
    stop("Row counts of `predicted_means` and `log_lik_matrix` differ.")
  }
  if (ncol(predicted_means) != length(observed_data)) {
    stop("Number of columns in `predicted_means` != length(observed_data).")
  }
  if (any(dim(predicted_means) != dim(log_lik_matrix))) {
    stop("`predicted_means` and `log_lik_matrix` must have the same dimensions.")
  }
  
  # Number of draws (M) and observations (N)
  M <- nrow(predicted_means)
  N <- length(observed_data)
  
  # Convert log_lik_matrix to an exponentiated form for relative_eff()
  # You also need a vector of chain IDs. If your total draws M is from, say, 4 chains,
  # each with (M/4) post-warmup draws, you can do something like:
  chains <- 8  # <-- replace as needed
  draws_per_chain <- M / chains
  chain_id <- rep(seq_len(chains), each = draws_per_chain)
  
  # r_eff requires the 'loo' package
  r_eff <- loo::relative_eff(
    x        = exp(log_lik_matrix),
    chain_id = chain_id
  )
  
  # Build a PSIS object
  psis_object <- loo::psis(
    log_ratios = -log_lik_matrix,
    r_eff      = r_eff
  )
  
  # E_loo() gives the PSIS-LOO expectation of ypred
  #   i.e. the "leave-one-out prediction" for each observation.
  ypred_loo <- loo::E_loo(
    x           = predicted_means,
    psis_object = psis_object,
    log_ratios  = -log_lik_matrix
  )$value
  
  # Compute the LOO residuals
  eloo <- ypred_loo - observed_data
  
  # Now we do the Dirichlet resampling approach to get a distribution for R²
  # For details, see: https://avehtari.github.io/bayes_R2/bayes_R2.html
  # (You need `rdirichlet` or `rudirichlet` from the 'MCMCpack' or 'extraDistr' packages.)
  
  # Number of resamples:
  nsim <- 4000
  
  # Sample Dirichlet weights for N items
  # e.g. from library("extraDistr") => `rdiri(n, alpha)`
  # or library("MCMCpack") => `rdirichlet(n, alpha)`
  # For uniform alpha=1 => each draw is Dirichlet(1,...,1)
  
  if (!requireNamespace("extraDistr", quietly = TRUE)) {
    stop("Package 'extraDistr' must be installed for Dirichlet draws (or use another method).")
  }
  
  rd <- extraDistr::rdirichlet(nsim, rep(1, N))  # shape: (nsim x N)
  
  # We compute the variance of the observed data, y:
  #   var(y) ~ sum(weights * y^2) - [sum(weights * y)]^2
  # multiplied by (N/(N-1)) correction
  # rowSums() is for each row of rd
  vary <- (
    rowSums(sweep(rd, 2, observed_data^2, FUN = "*")) -
      rowSums(sweep(rd, 2, observed_data, FUN = "*"))^2
  ) * (N / (N - 1))
  
  # The variance of the LOO residuals, e_loo
  vareloo <- (
    rowSums(sweep(rd, 2, eloo^2, FUN = "*")) -
      rowSums(sweep(rd, 2, eloo, FUN = "*"))^2
  ) * (N / (N - 1))
  
  # The LOO-based R²
  looR2 <- 1 - (vareloo / vary)
  
  # Clip to [-1, 1]
  looR2[looR2 < -1] <- -1
  looR2[looR2 > 1]  <- 1
  
  return(looR2)
}

# Fit your model with cmdstanr
# mod <- cmdstan_model("my_model.stan")
# fit <- mod$sample(data = your_data_list, ...)

# 1) Extract your observed_data (y) from your original data or from the same data list
observed_data <- gaborgen_stan_list$amplitude  # or df$y, etc.

# 2) Extract the log_lik draws:
#    For this, your Stan model must have a 'log_lik[...]' in the parameters or generated quantities
model001_log_lik_array <- model001_fit$draws("log_lik", format = "draws_matrix") 
model003_log_lik_array <- model003_fit$draws("log_lik", format = "draws_matrix") 
model008_log_lik_array <- model008_fit$draws("log_lik", format = "draws_matrix") 
# shape: (#draws) x (#observations) 
# or possibly (#draws) x something else. Check dimension.

# 3) Extract the posterior prediction draws:
#    Suppose your Stan model has a parameter or generated quantity called 'y_hat' or 'mu_pred'
model001_predicted_means <- model001_fit$draws("mu_pred", format = "draws_matrix")
model003_predicted_means <- model003_fit$draws("mu_pred", format = "draws_matrix")
model008_predicted_means <- model008_fit$draws("mu_pred", format = "draws_matrix")
# shape: (#draws) x (#observations)

# 4) Call your LOO-based R² function
model001_loo_r2 <- loo_R2_cmdstanr(
  observed_data   = observed_data,
  predicted_means = model001_predicted_means,
  log_lik_matrix  = model001_log_lik_array
)
model003_loo_r2 <- loo_R2_cmdstanr(
  observed_data   = observed_data,
  predicted_means = model003_predicted_means,
  log_lik_matrix  = model003_log_lik_array
)
model008_loo_r2 <- loo_R2_cmdstanr(
  observed_data   = observed_data,
  predicted_means = model008_predicted_means,
  log_lik_matrix  = model008_log_lik_array
)

# 5) Summarize or visualize
summary(model001_R2)
hist(model001_R2, breaks = 30)
mean(model001_R2)
summary(model001_loo_r2)
hist(model001_loo_r2, breaks = 30)
mean(model001_loo_r2)

summary(model003_R2)
hist(model003_R2, breaks = 30)
mean(model003_R2)
summary(model003_loo_r2)
hist(model003_loo_r2, breaks = 30)
mean(model003_loo_r2)


summary(model003_R2 - model001_R2)
hist(model003_R2 - model001_R2, breaks = 30)
((model003_R2 - model001_R2) > 0) %>% 
  sum()/40000
  

summary(model003_loo_r2 - model001_loo_r2)
hist(model003_loo_r2 - model001_loo_r2, breaks = 30)
((model003_loo_r2 - model001_loo_r2)> 0) %>% 
  sum()/4000

summary(model003_R2 - model008_R2)
hist(model003_R2 - model008_R2, breaks = 30)
((model003_R2 - model008_R2) > 0) %>% 
  sum()/40000
  

summary(model003_loo_r2 - model008_loo_r2)
hist(model003_loo_r2 - model008_loo_r2, breaks = 30)
((model003_loo_r2 - model008_loo_r2)> 0) %>% 
  sum()/4000


# ELPD play####
loo_comapre_df <- rbind(
  data.frame(
  model008_fit_loo[["pointwise"]],
  observation_number = 1:gaborgen_stan_list$n_observations,
  model = 1),
  data.frame(
  model001_fit_loo[["pointwise"]],
  observation_number = 1:gaborgen_stan_list$n_observations,
  model = 2),
  data.frame(
  model003_fit_loo[["pointwise"]],
  observation_number = 1:gaborgen_stan_list$n_observations,
  model = 3)) %>% 
  pivot_wider(
    names_from = model,
    values_from = -observation_number
  ) %>% 
  mutate(.before = 1, 
         participant = gaborgen_stan_list$participant[gaborgen_stan_list$indices_observed],
         block = gaborgen_stan_list$block[gaborgen_stan_list$indices_observed],
         cue = gaborgen_stan_list$cue[gaborgen_stan_list$indices_observed],
         bad_k = model003_fit_loo[["pointwise"]][,5] > .7)

loo_comapre_df %>% 
  filter(bad_k == T)

loo_comapre_df %>% 
  # filter(participant %in% c(1,2,6,9,11,14,21)) %>% 
  group_by(participant) %>%
  reframe(mean_1 = mean(elpd_loo_1),
          mean_2 = mean(elpd_loo_2),
          mean_3 = mean(elpd_loo_3),
          sum_1 = sum(elpd_loo_1),
          sum_2 = sum(elpd_loo_2),
          sum_3 = sum(elpd_loo_3)) %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  # geom_text(aes(x = sum_2,
  #               y = sum_3,
  #               label = participant),check_overlap = T) +
  geom_text(aes(x = mean_2,
                y = mean_3,
                label = participant),check_overlap = T) +
  theme_classic()

loo_comapre_df %>% 
  # filter(participant %in% c(1,2,6,9,11,14,21)) %>% 
  group_by(block, cue) %>%
  reframe(mean_1 = mean(elpd_loo_1),
          mean_2 = mean(elpd_loo_2),
          mean_3 = mean(elpd_loo_3),
          sum_1 = sum(elpd_loo_1),
          sum_2 = sum(elpd_loo_2),
          sum_3 = sum(elpd_loo_3)) %>%
  mutate(cue = factor(cue, levels = unique(cue))) %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  # geom_text(aes(x = sum_2,
  #               y = sum_3,
  #               label = participant),check_overlap = T) +
  geom_text(aes(x = sum_2,
                y = sum_3,
                label = block,
                color = cue),
            check_overlap = T) +
  scale_color_manual(values = cue_color) +
  theme_classic()
  
  

# Generate data
x <- seq(-4, 4, length.out = 100)  # Values for the variable
mu <- 0                            # Mean of the standard normal
sigma <- 1                         # Standard deviation of the standard normal

# Compute log-likelihood
log_likelihood_std_df <- data.frame(x,
                                    log_likelihood = dnorm(x, mean = mu, sd = sigma, log = TRUE))

# Plot the log-likelihood
log_likelihood_std_df %>% 
  ggplot() +
  geom_line(aes(x = x, y = log_likelihood)) +
  scale_y_continuous(breaks = seq(-10,0,by = .5)) +
  theme_bw()

ggsave(filename = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/stdnormal_log_like.png",
       dpi = 300, units = "in", 
       height = 2, width = 3,
       scale = 2)

# Step 1: Generate random samples
set.seed(123)                     # For reproducibility
n_samples <- 100000                 # Number of samples
samples <- rnorm(n_samples, mean = 0, sd = 1)  # Standard normal distribution

# Step 2: Compute log-likelihoods
log_likelihood_samples_std_df <- data.frame(samples,
  log_likelihoods = dnorm(samples, mean = 0, sd = 1, log = TRUE))

log_likelihood_samples_std_df$log_likelihoods %>% mean()

(loo_comapre_df %>% 
  ggplot() +
  # geom_density(aes(x = elpd_loo_1)) +
  geom_histogram(aes(x = elpd_loo_1),bins = 1000) +
  coord_cartesian(xlim = c(-3, -0.75)) +
  theme_bw() )/
  (ggplot() +
  # geom_density(data = log_likelihood_samples_std_df,
  #                aes(x = log_likelihoods)) +
  geom_histogram(data = log_likelihood_samples_std_df,
                 aes(x = log_likelihoods),
                 bins = 1000) +
  coord_cartesian(xlim = c(-3, -0.75)) +
  theme_bw() )

loo_comapre_df %>% 
  ggplot() +
  geom_density(aes(x = elpd_loo_1)) +
  geom_density(data = log_likelihood_samples_std_df,
                 aes(x = log_likelihoods),
               color = "red") +
  coord_cartesian(xlim = c(-3, -0.75)) +
  theme_bw()

loo_comapre_df %>% 
  ggplot() +
  geom_density(aes(x = elpd_loo_1, y = ..count.. / sum(..count..) * 100)) + # Adjust density to percentages
  geom_density(data = log_likelihood_samples_std_df,
               aes(x = log_likelihoods, y = ..count.. / sum(..count..) * 100),
               color = "red") +
  coord_cartesian(xlim = c(-3, -0.75)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Use percentage labels
  theme_bw()

loo_comapre_df %>% 
  ggplot() +
  geom_density(aes(x = elpd_loo_1)) +
  geom_density(aes(x = elpd_loo_2),
               color = "red") +
  geom_density(aes(x = elpd_loo_3),
               color = "blue") +
  geom_density(data = log_likelihood_samples_std_df,
               aes(x = log_likelihoods),
               color = "red") +
  # coord_cartesian(xlim = c(-3, -0.75)) +
  theme_bw()


## Adaptation ####
# plot specifics
number_of_samples_to_plot <- 4000
number_of_chains_prior <- 4
prior_samples_per_chain <- number_of_samples_to_plot / number_of_chains_prior
line_width_adapt <- .1
line_alpha_adapt <- .1

## Median participant adaptation

median_model001_adaptation_df <- model001_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))

median_model003_adaptation_df <- model003_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))

median_model008_adaptation_df <- model008_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))



model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/priors_visualization_2.stan'

model_priors <- cmdstanr::cmdstan_model(model_path, 
                                        force_recompile = T)

#Model source code
# model_priors$print()

model_priors_fit <- model_priors$sample(refresh = 1000,
                                        seed = 4,
                                        iter_warmup = prior_samples_per_chain, 
                                        iter_sampling = prior_samples_per_chain, 
                                        save_warmup = F, 
                                        show_messages = T,
                                        output_dir = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains",
                                        chains = number_of_chains_prior,
                                        parallel_chains = 4)

model_priors_fit_df <- model_priors_fit$draws(format = "df")

number_of_samples_to_plot <- 4000

set.seed(0)
samples_to_plot <- sample(1:nrow(model001_df),
                          size = number_of_samples_to_plot,
                          replace = F)


adaptation_prior_plot <-
  model_priors_fit_df %>% 
  ggplot() +
  geom_abline(aes(intercept = intercept_average - (fatigue_average * (gaborgen_stan_list$n_trials/2)), 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  # geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  # geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  # geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
  # scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1),
  scale_y_continuous(limits =c(-2,2), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = c(1,seq(10, 176, by = 10)),
                     expand = c(0,0),
                     name = "Trial") +
  ggtitle("Prior for average adaptation over trials") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 15),
        axis.ticks.y = element_blank())

adaptation_model1_avg_plot <- model008_df[samples_to_plot,] %>%
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  # geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  # geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  # geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = c(1,seq(10, 176, by = 10)),
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma") +
  ggtitle("Model001: Only Adaptation") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 15),
        axis.ticks.y = element_blank())


cue_block_average_adapt_plot <- model001_df[samples_to_plot,] %>%
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  # geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  # geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  # geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = c(1,seq(10, 176, by = 10)),
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma") +
  ggtitle("Model002: Cue by Block") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 15),
        axis.ticks.y = element_blank())

learning_average_adapt_plot <- model003_df[samples_to_plot,] %>%
  ggplot() +
  geom_hline(yintercept = -1, color = "black", linetype = "dotted") +
  geom_hline(yintercept = -.5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_hline(yintercept = .5, color = "black", linetype = "dotted") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted") +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              linewidth = line_width_adapt, 
              alpha = line_alpha_adapt) +
  # geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  # geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  # geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
  scale_y_continuous(limits =c(-1,1), 
                     breaks = seq(-4, 4, by = 1),
                     name = "Z-scored ssVEP") +
  scale_x_continuous(limits =c(1,176), 
                     breaks = c(1,seq(10, 176, by = 10)),
                     expand = c(0,0),
                     name = "Trial") +
  scale_color_viridis_d(option = "plasma") +
  ggtitle("Model003: Learning") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 15),
        axis.ticks.y = element_blank())


layout_grid <- c('
A
B
C
D
')


adaptation_prior_plot +
  adaptation_model1_avg_plot +
  cue_block_average_adapt_plot +
  learning_average_adapt_plot +
  plot_annotation(title = "hold",
                  theme = theme(
                    plot.title = element_text(family = "Arial",
                                              size = 27,
                                              color = "black",
                                              hjust = 0.5,
                                              face = "bold")
                  )) +
  plot_layout(design = layout_grid)



geom_abline(data = median_model001_adaptation_df,
            aes(intercept = intercept, 
                slope = fatigue,
                color = index),
            linewidth = .4, 
            alpha = 1) 




(model_priors_fit_df %>% 
    # as.data.frame() %>% 
    ggplot() +
    geom_abline(aes(intercept = intercept_average - (fatigue_average * (gaborgen_stan_list$n_trials/2)), 
                    slope = fatigue_average),
                linewidth = line_width_adapt, 
                alpha = line_alpha_adapt) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
    geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
    # scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1),
    scale_y_continuous(limits =c(-2,2), 
                       breaks = seq(-4, 4, by = 1),
                       name = "Z-scored ssVEP") +
    scale_x_continuous(limits =c(1,176), 
                       breaks = c(1,seq(10, 176, by = 10)),
                       expand = c(0,0),
                       name = "Trial") +
    ggtitle("Prior for average adaptation over trials") +
    theme_bw() +
    theme(text = element_text(family = "arial", size = 15),
          axis.ticks.y = element_blank())) / 
  
  (model001_df[samples_to_plot,] %>%
     # (model001_df %>% 
     # select(intercept_average,fatigue_average) %>% 
     # pivot_longer(everything())
     ggplot() +
     geom_abline(aes(intercept = intercept_average, 
                     slope = fatigue_average),
                 linewidth = line_width_adapt, 
                 alpha = line_alpha_adapt) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     ggtitle("Model001 average adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank())) /
  (model003_df[samples_to_plot,] %>% 
     # (model003_df %>% 
     # select(intercept_average,fatigue_average) %>% 
     # pivot_longer(everything())
     ggplot() +
     geom_abline(aes(intercept = intercept_average, 
                     slope = fatigue_average),
                 # color = "blue",
                 linewidth = line_width_adapt, 
                 alpha = line_alpha_adapt) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     ggtitle("Model003 average adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank())) /
  (model005_df[samples_to_plot,] %>% 
     # select(intercept_average,fatigue_average) %>% 
     # pivot_longer(everything())
     ggplot() +
     geom_abline(aes(intercept = intercept_average, 
                     slope = fatigue_average),
                 # color = "blue",
                 linewidth = line_width_adapt, 
                 alpha = line_alpha_adapt) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     ggtitle("Model005 average adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank())) 


## Median participant adaptation ####

median_model001_adaptation_df <- model001_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))

median_model003_adaptation_df <- model003_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))


(median_model001_adaptation_df %>% 
    ggplot() +
    geom_abline(aes(intercept = intercept, 
                    slope = fatigue,
                    color = index),
                linewidth = .5, 
                alpha = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
    geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
    scale_y_continuous(limits =c(-2,2), 
                       breaks = seq(-4, 4, by = 1),
                       name = "Z-scored ssVEP") +
    scale_x_continuous(limits =c(1,176), 
                       breaks = c(1,seq(10, 176, by = 10)),
                       expand = c(0,0),
                       name = "Trial") +
    # scale_color_viridis_d(option = "cividis") +
    # scale_color_viridis_d(option = "viridis") +
    scale_color_viridis_d(option = "plasma") +
    ggtitle("Model001 Participant Median Adaptation") +
    theme_bw() +
    theme(text = element_text(family = "arial", size = 15),
          axis.ticks.y = element_blank())) /
  (median_model003_adaptation_df %>% 
     ggplot() +
     geom_abline(aes(intercept = intercept, 
                     slope = fatigue,
                     color = index),
                 linewidth = .5, 
                 alpha = 1) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     # scale_color_viridis_d(option = "cividis") +
     # scale_color_viridis_d(option = "viridis") +
     scale_color_viridis_d(option = "plasma") +
     ggtitle("Model003 Participant Median Adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank()))




### old loo below####


Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(participant, cue) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(x = participant,
                    ymin = sum_elpd_1_minus_3 - sum_se,
                    ymax = sum_elpd_1_minus_3 + sum_se,
                    y = sum_elpd_1_minus_3,
                    color = cue),
                width = 0,
                position = position_quasirandom(width = .1)) + 
  geom_quasirandom(aes(x = participant,
                       y = sum_elpd_1_minus_3,
                       color = cue),width = .1,
                   alpha = 1) +
  scale_color_manual(values = cue_color) +
  theme_classic() +
  theme(text = element_text(family = text_font,
                            size = text_size,
                            color = "black"),
        axis.line = element_line(linewidth = axis_line_thickness,
                                 lineend = "square"),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )






Oz_fft_df %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(cue) %>%  # Ensure rolling median is calculated per cue
  arrange(trial) %>%  # Sort by trial for correct rolling median calculation
  mutate(moving_sum = zoo::rollsum(elpd_1_minus_3,k = 75, fill = NA,align = "center")) %>%
  mutate(moving_median = zoo::rollmedian(elpd_1_minus_3, k = 25, fill = NA, align = "center")) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 36)) +
  geom_vline(aes(xintercept = 36 + 48)) +
  geom_vline(aes(xintercept = 36 + 48 + 48)) +
  # geom_quasirandom(aes(x = trial,
  #                      y = elpd_1_minus_3,
  #                      color = cue),
  #                  size = loo_dot_size,
  #                  alpha = loo_dot_alpha) +
  # geom_line(aes(x = trial, y = moving_median, color = cue)) +
  geom_line(aes(x = trial,
                y = moving_sum,
                color = cue),
            linewidth = 1.5) +
  scale_color_manual(values = cue_color) +
  theme_classic()





Oz_fft_df %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1))+
  geom_point(aes(x = model001_elpd,
                 y = model003_elpd,
                 color = cue)) +
  scale_color_manual(values = cue_color) +
  scale_y_continuous(breaks = seq(-10,0,by = .5)) +
  scale_x_continuous(breaks = seq(-10,0,by = .5)) +
  theme_bw()

Oz_fft_df %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1))+
  geom_point(aes(x = model001_elpd,
                 y = model005_elpd,
                 color = cue)) +
  scale_color_manual(values = cue_color) +
  theme_classic()
  
Oz_fft_df %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1))+
  geom_point(aes(x = model001_elpd,
                 y = model003_elpd,
                 color = cue),
             alpha = .5) +
  scale_color_manual(values = cue_color) +
  theme_classic() +
  facet_wrap(~block)


Oz_fft_df %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"],
         model005_elpd = model005_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  ggplot() +
  # geom_abline(aes(intercept = 0, slope = 1))+
  geom_point(aes(x = trial,
                 y = elpd_1_minus_3,
                 color = cue), alpha = .5) +
  scale_color_manual(values = cue_color) +
  theme_bw()

Oz_fft_df %>% 
    mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
           model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"],
           model005_elpd = model005_fit_loo$pointwise[, "elpd_loo"]) %>% 
    mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
    ggplot() +
    # geom_abline(aes(intercept = 0, slope = 1))+
    geom_point(aes(x = trial,
                   y = elpd_1_minus_3,
                   color = cue), alpha = .5) +
    scale_color_manual(values = cue_color) +
    theme_bw() +
  facet_wrap(~cue)

Oz_fft_df %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"],
         model005_elpd = model005_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(trial, cue) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3)) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  # geom_abline(aes(intercept = 0, slope = 1))+
  geom_line(aes(x = trial,
                y = mean_elpd_1_minus_3,
                color = cue), alpha = .5) +
  scale_color_manual(values = cue_color) +
    theme_bw()
  
Oz_fft_df %>% 
    mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
           model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"],
           model005_elpd = model005_fit_loo$pointwise[, "elpd_loo"]) %>% 
    mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
    mutate(elpd_5_minus_3 = model005_elpd - model003_elpd) %>% 
    ggplot() +
    # geom_abline(aes(intercept = 0, slope = 1))+
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 36)) +
    geom_vline(aes(xintercept = 36 + 48)) +
    geom_vline(aes(xintercept = 36 + 48 + 48)) +
    geom_smooth(aes(x = trial,
                   y = elpd_1_minus_3,
                   # y = elpd_5_minus_3,
                   color = cue), 
                , span = 0.25,
                alpha = .5,se = 0) +
    # geom_point(aes(x = trial,
    #                y = elpd_1_minus_3,
    #                color = cue), alpha = .5) +
    scale_color_manual(values = cue_color) +
    theme_bw()

Oz_fft_df %>% 
    mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
           model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"],
           model005_elpd = model005_fit_loo$pointwise[, "elpd_loo"]) %>% 
    mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
    mutate(elpd_5_minus_3 = model005_elpd - model003_elpd) %>% 
    mutate(elpd_1_minus_5 = model001_elpd - model005_elpd) %>% 
    ggplot() +
    # geom_abline(aes(intercept = 0, slope = 1))+
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 36)) +
    geom_vline(aes(xintercept = 36 + 48)) +
    geom_vline(aes(xintercept = 36 + 48 + 48)) +
    geom_smooth(aes(x = trial,
                   y = elpd_1_minus_5,
                   # y = elpd_5_minus_3,
                   color = cue), 
                alpha = .5,se = 0) +
    # geom_point(aes(x = trial,
    #                y = elpd_1_minus_3,
    #                color = cue), alpha = .5) +
    scale_color_manual(values = cue_color) +
    theme_bw()

## LOO by participant ####

Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(participant) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  mutate(sum_better_for_mod3 = sum_elpd_1_minus_3 < 0) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  # geom_pointrange(aes(x = participant,
  #                     y = mean_elpd_1_minus_3,
  #                     ymin = mean_elpd_1_minus_3 - mean_se,
  #                     ymax = mean_elpd_1_minus_3 + mean_se), 
  #                 alpha = 1) +
  geom_pointrange(aes(x = participant,
                      y = sum_elpd_1_minus_3,
                      ymin = sum_elpd_1_minus_3 - sum_se,
                      ymax = sum_elpd_1_minus_3 + sum_se,
                      color = sum_better_for_mod3),
                  alpha = 1) +
  theme_bw()

Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(participant, cue) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  # mutate(sum_better_for_mod3 = sum_elpd_1_minus_3 < 0) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  # geom_pointrange(aes(x = participant,
  #                     y = mean_elpd_1_minus_3,
  #                     ymin = mean_elpd_1_minus_3 - mean_se,
  #                     ymax = mean_elpd_1_minus_3 + mean_se), 
  #                 alpha = 1) +
  geom_pointrange(aes(x = participant,
                      y = sum_elpd_1_minus_3,
                      ymin = sum_elpd_1_minus_3 - sum_se,
                      ymax = sum_elpd_1_minus_3 + sum_se,
                      color = cue),
                  position = position_dodge(width = .1),
                  alpha = 1) +
  scale_color_manual(values = cue_color) +
  theme_bw()

Oz_fft_df %>% 
  mutate(participant = participant %>% as.factor() %>% as.integer() %>% as.factor) %>% 
  mutate(model001_elpd = model001_fit_loo$pointwise[, "elpd_loo"],
         model003_elpd = model003_fit_loo$pointwise[, "elpd_loo"]) %>% 
  mutate(elpd_1_minus_3 = model001_elpd - model003_elpd) %>% 
  group_by(participant, cue, block) %>% 
  summarise(mean_elpd_1_minus_3 = mean(elpd_1_minus_3),
            mean_se = plotrix::std.error(elpd_1_minus_3),
            sum_elpd_1_minus_3 = sum(elpd_1_minus_3),
            sum_se = var(elpd_1_minus_3) * sqrt(n())) %>% 
  # mutate(sum_better_for_mod3 = sum_elpd_1_minus_3 < 0) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0)) +
  # geom_pointrange(aes(x = participant,
  #                     y = mean_elpd_1_minus_3,
  #                     ymin = mean_elpd_1_minus_3 - mean_se,
  #                     ymax = mean_elpd_1_minus_3 + mean_se), 
  #                 alpha = 1) +
  geom_pointrange(aes(x = participant,
                      y = sum_elpd_1_minus_3,
                      ymin = sum_elpd_1_minus_3 - sum_se,
                      ymax = sum_elpd_1_minus_3 + sum_se,
                      color = cue),
                  position = position_dodge(width = .1),
                  alpha = 1) +
  scale_color_manual(values = cue_color) +
  facet_wrap(~block) +
  theme_bw()



# Compare models on their average prediction over trials ####
model001_fit_meta_data <- model001_fit$metadata()

model001_fit_relevant_parameters <- model001_fit_meta_data$model_params[
  !str_detect(model001_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model001_df <- posterior::as_draws_df(
  model001_fit$draws(variables = model001_fit_relevant_parameters))

model003_fit_meta_data <- model003_fit$metadata()

model003_fit_relevant_parameters <- model003_fit_meta_data$model_params[
  !str_detect(model003_fit_meta_data$model_params, "log_lik|mu_pred|amplitude|CSP_|S|raw")]

model003_df <- posterior::as_draws_df(
  model003_fit$draws(variables = model003_fit_relevant_parameters))

model005_fit_meta_data <- model005_fit$metadata()

model005_fit_relevant_parameters <- model005_fit_meta_data$model_params[
  !str_detect(model005_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]

model005_df <- model005_fit$draws(variables = model005_fit_relevant_parameters,
                                  format = "df")
# model006_fit_meta_data <- model006_fit$metadata()
# 
# model006_fit_relevant_parameters <- model006_fit_meta_data$model_params[
#   !str_detect(model006_fit_meta_data$model_params, "log_lik|mu|amplitude|raw|S")]
# 
# 
# model006_df <- model006_fit$draws(variables = model006_fit_relevant_parameters,
#                                   format = "df")


## Sigma ####
data.frame(model001_average = model001_df$sigma_average,
           model003_average = model003_df$sigma_average) %>% 
  # mutate(diff_mod1_mod3 = model001_average - model003_average) %>% 
pivot_longer(everything()) %>% 
  ggplot() +
  geom_vline(xintercept = 1, linewidth = 2) +
  geom_density(aes(x = value, 
                   color = name,
                   fill = name),
               linewidth = 2,
               alpha = .3) +
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(0, 33)) +
  theme_bw()

data.frame(model001_average = model001_df$sigma_average,
           model003_average = model003_df$sigma_average) %>% 
  summarise(diff_mod1_mod3 = sum(model001_average > model003_average)/(8*5000))


sigma_par_df <- rbind(model001_df %>% 
                        select(starts_with("sigma[")) %>% 
                        pivot_longer(starts_with("sigma[")) %>% 
                        mutate(name = factor(name, levels = unique(name)),
                               model = factor("model001", 
                                              levels = c("model001",
                                                         "model003"))),
                      model003_df %>% 
                        select(starts_with("sigma[")) %>% 
                        pivot_longer(starts_with("sigma[")) %>% 
                        mutate(name = factor(name, levels = unique(name)),
                               model = factor("model003", 
                                              levels = c("model001",
                                                         "model003"))))


sigma_par_df %>% 
  ggplot() +
  geom_density_ridges(aes(x = value, 
                          y = factor(name, levels = rev(levels(name))),
                          color = model,
                          fill = model),
                      linewidth = 1.25,
                      scale = 1,
                      alpha = 0.3) +
  coord_cartesian(xlim = c(0.8, 1.01),
                  ylim = c(1.4,24.75)) +
  theme_bw()

normalized_sigma_par_densities <- sigma_par_df %>% 
  ggplot(aes(x = value, 
             color = model, 
             fill = model)) +
  geom_density_ridges(aes(height = after_stat(scaled),
                          y = factor(name, levels = rev(levels(name)))),
                      linewidth = 1.25,
                      stat = "density",
                      scale = .9,
                      alpha = .3) +
  coord_cartesian(xlim = c(0.8, 1.01),
                  expand = 0,
                  ylim = c(1,25)) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank())

normalized_sigma_par_densities

# sigma_par_df %>% 
#   ggplot() +
#   geom_line(aes(x = value, 
#                 y = after_stat(scaled),
#                 color = model),
#             linewidth = 1.5,
#             stat = "density",
#             alpha = 1) +
#   coord_cartesian(xlim = c(0.8, 1.01)) +
#   theme_bw() 




## Adaptation ####
# plot specifics
number_of_samples_to_plot <- 4000
number_of_chains_prior <- 4
prior_samples_per_chain <- number_of_samples_to_plot / number_of_chains_prior
line_width_adapt <- .1
line_alpha_adapt <- .1


model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/priors_visualization_2.stan'

model_priors <- cmdstanr::cmdstan_model(model_path, 
                                        force_recompile = T)

#Model source code
# model_priors$print()

model_priors_fit <- model_priors$sample(refresh = 1000,
                                        seed = 4,
                                        iter_warmup = prior_samples_per_chain, 
                                        iter_sampling = prior_samples_per_chain, 
                                        save_warmup = F, 
                                        show_messages = T,
                                        output_dir = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains",
                                        chains = number_of_chains_prior,
                                        parallel_chains = 4)

model_priors_fit_df <- model_priors_fit$draws(format = "df")

number_of_samples_to_plot <- 4000

set.seed(0)
samples_to_plot <- sample(1:nrow(model001_df),
                          size = number_of_samples_to_plot,
                          replace = F)

(model_priors_fit_df %>% 
  # as.data.frame() %>% 
   ggplot() +
   geom_abline(aes(intercept = intercept_average - (fatigue_average * (gaborgen_stan_list$n_trials/2)), 
                   slope = fatigue_average),
               linewidth = line_width_adapt, 
               alpha = line_alpha_adapt) +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
   geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
 # scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1),
 scale_y_continuous(limits =c(-2,2), 
                    breaks = seq(-4, 4, by = 1),
                    name = "Z-scored ssVEP") +
   scale_x_continuous(limits =c(1,176), 
                      breaks = c(1,seq(10, 176, by = 10)),
                      expand = c(0,0),
                      name = "Trial") +
   ggtitle("Prior for average adaptation over trials") +
   theme_bw() +
   theme(text = element_text(family = "arial", size = 15),
         axis.ticks.y = element_blank())) / 
  
  (model001_df[samples_to_plot,] %>%
  # (model001_df %>% 
     # select(intercept_average,fatigue_average) %>% 
     # pivot_longer(everything())
     ggplot() +
     geom_abline(aes(intercept = intercept_average, 
                     slope = fatigue_average),
                 linewidth = line_width_adapt, 
                 alpha = line_alpha_adapt) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     ggtitle("Model001 average adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank())) /
  (model003_df[samples_to_plot,] %>% 
  # (model003_df %>% 
     # select(intercept_average,fatigue_average) %>% 
     # pivot_longer(everything())
     ggplot() +
     geom_abline(aes(intercept = intercept_average, 
                     slope = fatigue_average),
                 # color = "blue",
                 linewidth = line_width_adapt, 
                 alpha = line_alpha_adapt) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     ggtitle("Model003 average adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank())) /
  (model005_df[samples_to_plot,] %>% 
     # select(intercept_average,fatigue_average) %>% 
     # pivot_longer(everything())
     ggplot() +
     geom_abline(aes(intercept = intercept_average, 
                     slope = fatigue_average),
                 # color = "blue",
                 linewidth = line_width_adapt, 
                 alpha = line_alpha_adapt) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     ggtitle("Model005 average adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank())) 


## Median participant adaptation ####

median_model001_adaptation_df <- model001_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))

median_model003_adaptation_df <- model003_df %>% 
  select(starts_with("intercept["), starts_with("fatigue[")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("parameter_type", "index"),
    names_pattern = "(intercept|fatigue)\\[(.*)\\]",
    values_to = "value") %>%
  group_by(parameter_type, index) %>% 
  reframe(median_value = median(value)) %>% 
  pivot_wider(
    names_from = parameter_type,
    values_from = median_value) %>% 
  mutate(index = as.numeric(index)) %>% 
  arrange(index) %>% 
  mutate(index = factor(index, 
                        levels = 1:gaborgen_stan_list$n_participants))


(median_model001_adaptation_df %>% 
    ggplot() +
    geom_abline(aes(intercept = intercept, 
                    slope = fatigue,
                    color = index),
                linewidth = .5, 
                alpha = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
    geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
    scale_y_continuous(limits =c(-2,2), 
                       breaks = seq(-4, 4, by = 1),
                       name = "Z-scored ssVEP") +
    scale_x_continuous(limits =c(1,176), 
                       breaks = c(1,seq(10, 176, by = 10)),
                       expand = c(0,0),
                       name = "Trial") +
  # scale_color_viridis_d(option = "cividis") +
  # scale_color_viridis_d(option = "viridis") +
  scale_color_viridis_d(option = "plasma") +
    ggtitle("Model001 Participant Median Adaptation") +
    theme_bw() +
    theme(text = element_text(family = "arial", size = 15),
          axis.ticks.y = element_blank())) /
  (median_model003_adaptation_df %>% 
     ggplot() +
     geom_abline(aes(intercept = intercept, 
                     slope = fatigue,
                     color = index),
                 linewidth = .5, 
                 alpha = 1) +
     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
     geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
     geom_hline(yintercept = -1, color = "red", linetype = "dotted") +
     scale_y_continuous(limits =c(-2,2), 
                        breaks = seq(-4, 4, by = 1),
                        name = "Z-scored ssVEP") +
     scale_x_continuous(limits =c(1,176), 
                        breaks = c(1,seq(10, 176, by = 10)),
                        expand = c(0,0),
                        name = "Trial") +
     # scale_color_viridis_d(option = "cividis") +
     # scale_color_viridis_d(option = "viridis") +
     scale_color_viridis_d(option = "plasma") +
     ggtitle("Model003 Participant Median Adaptation") +
     theme_bw() +
     theme(text = element_text(family = "arial", size = 15),
           axis.ticks.y = element_blank()))



## Average Cue ####

# not sure this makes sense for all models

# Model 1 Bcue posteriors ####
bcue_df <- model001_df %>%
  select(starts_with("bcue")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(block = case_when(grepl(pattern = "\\[1", x = name) ~ "Habituation",
                         grepl(pattern = "\\[2", x = name) ~ "Acquisition #1",
                         grepl(pattern = "\\[3", x = name) ~ "Acquisition #2",
                         grepl(pattern = "\\[4", x = name) ~ "Extinction"),
         cue = case_when(grepl(pattern = "1]", x = name) ~ "CS+",
                           grepl(pattern = "2]", x = name) ~ "GS1",
                           grepl(pattern = "3]", x = name) ~ "GS2",
                           grepl(pattern = "4]", x = name) ~ "GS3")) %>% 
  mutate(block = factor(block, 
                        levels = c("Habituation",
                                   "Acquisition #1",
                                   "Acquisition #2",
                                   "Extinction")))



bcue_df %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = value, 
                          y = factor(block, levels = rev(levels(block))),
                          color = cue,
                          fill = cue),
                      scale = 1,
                      linewidth = 2,
                      alpha = .1) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP", breaks = seq(-.75, .75, by = .25)) +
  coord_cartesian(xlim = c(-.5,.5),
                  ylim = c(1.55,4.25)) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Effect of Cue Controlling for Adaptation") + 
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.y = element_blank())



# Model 2 Learning Rates ####

# pars_that_learned <- gaborgen_stan_list$participant[
#   gaborgen_stan_list$learned_at_end == 1] %>% 
#   unique()
# double checking flat lines aren't predictive of conditioning effect
# (model003_df %>%
#   select(starts_with("intercept[")) %>% 
#   pivot_longer(cols = everything()) %>% 
#   mutate(name = factor(name,
#                        levels = unique(name))) %>% 
#   ggplot(aes(x = value, 
#              y = factor(name, 
#                         levels = rev(levels(name)))
#              # ,fill = learned
#   )) +
#   geom_density_ridges() +
#   ggtitle("intercept") +
#   theme_bw() +
#   theme(text = element_text(size = 12),
#         axis.title = element_blank(),
#         axis.ticks = element_blank())) +
#     
#   model003_df %>%
#     select(starts_with("fatigue[")) %>% 
#     pivot_longer(cols = everything()) %>% 
#     mutate(name = factor(name,
#                          levels = unique(name))) %>% 
#     ggplot(aes(x = value, y = factor(name, levels = rev(levels(name))))) +
#     geom_density_ridges() +
#     # scale_y_discrete(labels = paste0("Participant ", 
#     #                                  gaborgen_stan_list$n_participants:1)) +
#     # coord_cartesian(xlim = c(-0.01,1.01), 
#     #                 expand = 0,
#     #                 ylim = c(1,25)) +
#     ggtitle("fatigue") +
#     theme_bw() +
#     theme(text = element_text(size = 12),
#           axis.title = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks = element_blank()) 
    


learning_rate_plot <-
  
  model003_df %>%
  select(starts_with("learning_paired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  # mutate(learned = if_else(
  #   grepl(
  #     pattern = (paste0("\\[", 
  #                     pars_that_learned,
  #                     collapse = "|")), 
  #     x = name), T, F)) %>% 
  ggplot(aes(x = value, 
             y = factor(name, 
                        levels = rev(levels(name)))
             # ,fill = learned
             )) +
  stat_density_ridges(fill = "black",
                      scale = 0.9,
                      rel_min_height = 0.01,
                      # panel_scaling = FALSE,
                      from = 0,
                      to = 1) +
                      coord_cartesian(xlim = c(0, 1)) +
  scale_y_discrete(labels = paste0("Participant ", 
                                   gaborgen_stan_list$n_participants:1)) +
  coord_cartesian(xlim = c(-0.01,1.01), expand = 0,
                  ylim = c(1,25))+
  ggtitle("Learning Rate Paired") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  
  model003_df %>%
  select(starts_with("learning_unpaired[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot(aes(x = value, y = factor(name, levels = rev(levels(name))))) +
  stat_density_ridges(fill = "black",
                      scale = 0.9,
                      rel_min_height = 0.01,
                      # panel_scaling = FALSE,
                      from = 0,
                      to = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_y_discrete(labels = paste0("Participant ", 
                                   gaborgen_stan_list$n_participants:1)) +
  coord_cartesian(xlim = c(-0.01,1.01), 
                  expand = 0,
                  ylim = c(1,25))+
  ggtitle("Learning Rate Unpaired") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) 

# sigma_par_df %>% 
#   filter(model == "model001") %>% 
#   ggplot(aes(x = value)) +
#   geom_density_ridges(aes(height = after_stat(scaled),
#                           y = factor(name, levels = rev(levels(name)))),
#                       linewidth = 1.25,
#                       stat = "density",
#                       scale = .9,
#                       alpha = .3) +
#   coord_cartesian(xlim = c(0.85, 1.01),
#                   expand = 0,
#                   ylim = c(1,25)) +
#   theme_bw() +
#   theme(axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text.y = element_blank())

learning_rate_plot +
  normalized_sigma_par_densities +
  geom_vline(xintercept = 1) +
  ggtitle("Error per model") +
  theme(legend.position = "none")
  
# CS+ associative value
line_width_csp <- .2
line_alpha_csp <- .2
line_width_csp_average <- 1.5

model003_fit_meta_data <- model003_fit$metadata()

CSP_parameters <- model003_fit_meta_data$model_params[
  str_detect(model003_fit_meta_data$model_params, "CSP_ass")] 

CSP_df <- posterior::as_draws_df(
  model003_fit$draws(variables = CSP_parameters))

number_of_samples_to_plot <- 4000

set.seed(0)
samples_to_plot <- sample(1:nrow(CSP_df),
                          size = number_of_samples_to_plot,
                          replace = F)



CSP_plot_df <- CSP_df %>% 
  pivot_longer(starts_with("CSP_ass")) %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  mutate(participant = rep(gaborgen_stan_list$participant, 
                           nrow(CSP_df)),
         trial = rep(gaborgen_stan_list$trial, 
                     nrow(CSP_df))) %>% 
  group_by(participant, trial) %>% 
  mutate(mean_value = mean(value),
         median_value = median(value)) %>% 
  ungroup() %>% 
  filter(.draw %in% samples_to_plot)
  
CSP_plot <- CSP_plot_df %>% 
  # group_by(participant, trial) %>% 
  # mutate(mean_value = mean(value)) %>% 
  # ungroup() %>% 
  filter(#participant %in% c(1,2,3),
         trial >= 30) %>% 
  ggplot() +
  geom_line(aes(x = trial, 
                y = value, 
                group = .draw),
            linewidth = line_width_csp,
            alpha = line_alpha_csp) +
  geom_line(aes(x = trial, 
                y = median_value),
            linewidth = line_width_csp_average,
            # alpha = line_alpha_csp,
            color = "red"
            ) +
  scale_x_continuous(breaks = seq(30,180, by = 10),
                     name = "Trial") +
  coord_cartesian(ylim = c(-0.05,1.05),expand = F) +
  theme_classic() +
  facet_grid(participant ~., 
             # scales = "free_y",
             space = "free") +
  ggtitle("Associate Value Over Trials (Black = Posterior Draws, Red = Median)") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = 12),)

layout_grid <- '
ABCCD
'
learning_rate_plot +
  CSP_plot +
  (normalized_sigma_par_densities +   
     geom_vline(xintercept = 1) +
     ggtitle("Error per Model (Green = RW)") +
     theme(legend.position = "none",
           text = element_text(size = 12))) +
    plot_layout(design = layout_grid)
     

## Scaling ####
bcue_df %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0), linewidth = 2) +
  geom_density_ridges(aes(x = value, 
                          y = factor(block, levels = rev(levels(block))),
                          color = cue,
                          fill = cue),
                      scale = .95,
                      linewidth = 2,
                      alpha = .1) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP", breaks = seq(-.75, .75, by = .25)) +
  coord_cartesian(xlim = c(-.55,.55),
                  ylim = c(1,4.9),
                  expand = F) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Model 1\nEffect of Cue by Block") + 
  theme_classic() +
  theme(text = element_text(size = 22,family = "Arial"),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +


model003_df %>%
  select(starts_with("scaling[")) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name,
                       levels = unique(name))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1.5) +
  geom_density(aes(x = value, 
                   color = name,
                   fill = name),
               linewidth = 2,
               alpha = .1) +
  annotate("text",
           x = 1, 
           y = 4.85,
           label = "CS+",
           family = "Arial",
           color = cue_color[1],
           size = 18) +
  annotate("text",
           x = 1, 
           y = 4.85 - .4,
           label = "GS1",
           family = "Arial",
           color = cue_color[2],
           size = 18) +
  annotate("text",
           x = 1, 
           y = 4.85 - .4 - .4,
           label = "GS2",
           family = "Arial",
           color = cue_color[3],
           size = 18) +
  annotate("text",
           x = 1, 
           y = 4.85 - .4 - .4 - .4,
           label = "GS3",
           family = "Arial",
           color = cue_color[4],
           size = 18) +
  coord_cartesian(ylim = c(0, 5.1),
                  xlim = c(-.575, 1.25),
                  expand = F) +
  scale_color_manual(values = cue_color) +
  scale_fill_manual(values = cue_color) +
  scale_x_continuous(name = "Δ Z-Scored ssVEP") +
  ggtitle("Model 2\nAssociate Value Effect on Cue ssVEP") + 
  theme_classic() +
  theme(text = element_text(size = 22,family = "Arial"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")



## Scaling times CS+strength ####
trials_to_plot <-1:176
number_of_samples_to_plot <- 200

set.seed(0)
samples_to_plot <- sample(1:nrow(CSP_df),
                          size = number_of_samples_to_plot,
                          replace = F)


CSP_df %>% 
  # select(all_of(paste0("CSP_associative_strength[",trials_to_plot,"]")), ".draw") %>% 
  pivot_longer(starts_with("CSP_ass")) %>% 
    mutate(name = factor(name, levels = unique(name))) %>% 
    mutate(participant = rep(gaborgen_stan_list$participant,
                             nrow(CSP_df)),
           trial = rep(gaborgen_stan_list$trial, 
                       nrow(CSP_df)),
           cue = rep(gaborgen_stan_list$cue, 
                       nrow(CSP_df))) %>% 
    # mutate(participant = rep(gaborgen_stan_list$participant[trials_to_plot],
    #                          nrow(CSP_df)),
    #        trial = rep(gaborgen_stan_list$trial[trials_to_plot], 
    #                    nrow(CSP_df)),
    #        cue = rep(gaborgen_stan_list$cue[trials_to_plot], 
    #                    nrow(CSP_df))) %>% 
  mutate(CSP_pred = value * model003_df$`scaling[1]`) %>% 
  mutate(GS1_pred = value * model003_df$`scaling[2]`) %>% 
  mutate(GS2_pred = value * model003_df$`scaling[3]`) %>% 
  mutate(GS3_pred = value * model003_df$`scaling[4]`) %>% 
  select(contains("pred"), ".draw", "trial","participant") %>% 
  pivot_longer(contains("pred")) %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  filter(.draw %in% samples_to_plot,
         participant %in% c(1, 6 , 11, 14, 21)) %>% 
  ggplot(aes(x = trial, 
             y = value,
             color = name)) +
  geom_line(aes(
    # x = trial, 
    #             y = value, 
    #             color = name,
                group = interaction(.draw, name)),
            linewidth = .5,
            alpha = .1) +
            # linewidth = line_width_csp,
            # alpha = line_alpha_csp) +
  stat_summary(aes(group = name),
               fun = median,
               color = "black",
               geom = "line",
               linewidth = 1.25) +
  stat_summary(aes(group = name),
               fun = median,
               geom = "line",
               linewidth = 1) +
            # linewidth = line_width_csp,
            # alpha = line_alpha_csp) +
  # geom_smooth(aes(x = trial, 
  #                 y = value, 
  #                 group = name),
  #           linewidth = 1.25,
  #           span = 0.1,
  #           color = "black",
  #           alpha = .2, se = 0) +
  #           # linewidth = line_width_csp,
  #           # alpha = line_alpha_csp) +
  # geom_smooth(aes(x = trial, 
  #                 y = value, 
  #                 color = name),
  #           linewidth = 1,
  #           span = 0.1,
  #           alpha = .2, se = 0) +
  #           # linewidth = line_width_csp,
  #           # alpha = line_alpha_csp) +
  scale_x_continuous(breaks = seq(0,180, by = 10),
                     name = "Trial") +
  scale_color_manual(values = cue_color) +
  # coord_cartesian(ylim = c(-0.05,1.05),expand = F) +
  theme_classic() +
  facet_grid(participant ~.,
             # scales = "free_y",
             space = "free") +
  ggtitle("Associate Value Over Trials (Black = Posterior Draws, Red = Median)") +
  theme(#strip.text.y = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = 12),)
  
  


## Show missing trial interpolation####
missing_trial_posterior_samples <- 100
set.seed(0)
missing_trial_samples_to_plot <- sample(x = 1:(8*5e3), size = missing_trial_posterior_samples,replace = F)

amplitude_array_to_plot_labels <- paste0("amplitude_all[",1:176,"]")
amplitude_array_to_plot_labels_less <- paste0("amplitude_all[",10:50,"]")



amplitude_array_df <- model001_fit$draws(variables = amplitude_array_to_plot_labels,
                                         format = "df")


amplitude_array_df %>% 
  # mutate()
  filter(.draw %in% missing_trial_samples_to_plot) %>% 
  select(all_of(paste0("amplitude_all[",75:100,"]"))) %>% 
  pivot_longer(starts_with("ampl")) %>% 
  mutate(name = factor(name, unique(name))) %>% 
  ggplot() +
  geom_quasirandom(aes(x = name, y = value))


## Prior visualization ####

### Fatigue line prior####

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/priors_visualization_2.stan'

#force_recompile = T is sometimes helpful
model_priors <- cmdstanr::cmdstan_model(model_path, 
                                        force_recompile = T)

#Model source code
# model_priors$print()

# Clear previous chains
# list.files(pattern = "priors",
#            path = where_to_save_chains, 
#            full.names = T) %>% 
#   file.remove()

# prior_list <- list()
# prior_list$intecept_prior_sd <- 0.75
# prior_list$fatigue_prior_sd <- 0.01
model_priors_fit <- model_priors$sample(refresh = 50,
                                        seed = 3,
                                        iter_warmup = 1000, 
                                        iter_sampling = 1000, 
                                        save_warmup = F, 
                                        show_messages = T,
                                        output_dir = where_to_save_chains,
                                        chains = number_of_chains,
                                        parallel_chains = number_of_parallel_chains)


model_priors_fit_summary <- 
  model_priors_fit$summary()


model_priors_fit_df <- model_priors_fit$draws(format = "df")


number_of_samples_to_plot <- 1000

set.seed(0)
samples_to_plot <- sample(1:nrow(model_priors_fit_df),
                          size = number_of_samples_to_plot,
                          replace = F)

model_priors_fit_df[samples_to_plot,] %>% 
  # as.data.frame() %>% 
  ggplot() +
  geom_abline(aes(intercept = intercept_average, 
                  slope = fatigue_average),
              size = .2, 
              alpha = .2) +
  scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1)) +
  scale_x_continuous(limits =c(1,176), breaks = seq(0, 176, by = 10)) +
  ggtitle("Fatigue Prior Predictive regression lines") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 15),
        axis.ticks.y = element_blank())

model_priors_fit_df[samples_to_plot,] %>% 
  # as.data.frame() %>% 
  ggplot() +
  geom_abline(aes(intercept = intercept, 
                  slope = fatigue),
              size = .2, 
              alpha = .2) +
  scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1)) +
  # scale_y_continuous(limits =c(-10,10), breaks = seq(-4, 4, by = 1)) +
  scale_x_continuous(limits =c(1,176), breaks = seq(0, 176, by = 10)) +
  ggtitle("Fatigue Prior Predictive regression lines") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 15),
        axis.ticks.y = element_blank())

### SD priors ####
set.seed(0)
number_of_samples_plot <- 5000
  
# X ~ exp(Z)
# Z ~ normal(m,s)
# median X is exp(m)
# mean X is exp(m+s^2/2)
# mode X is exp(m - s^2)

#intercept
rnorm(number_of_samples_plot,-1,.75) %>%
  exp() %>% # transform back to correct scale
  # median()
  data.frame(value = .) %>% 
  ggplot() +
  coord_cartesian(xlim = c(0,5)) +
  geom_histogram(aes(x = value), bins = 200)



#fatigue
rnorm(number_of_samples_plot,-3,.5) %>%
  exp() %>% # transform back to correct scale
  # median()
  # min()
  data.frame(value = .) %>% 
  ggplot() +
  coord_cartesian(xlim = c(0,1)) +
  geom_histogram(aes(x = value), bins = 100)


# intercept sd for each participant
rnorm(number_of_samples_plot,-1,1) %>% # raw unbounded prior to improve sampling
  exp() %>% # transform back to correct scale
  # median()
  hist(breaks = 10000,xlim = c(0,2))




text_size <- 10
intercept_average_mean <- 0
intercept_average_sd <- 0.75
fatigue_average_mean <- 0
fatigue_average_sd <- 0.01
intercept_fatigue_average_cov <- rbeta(number_of_samples_plot, 1.1,1.1) * -1

cov_mat <- matrix(intercept_average_sd, intercept_fatigue_average_cov,
                  intercept_fatigue_average_cov, fatigue_average_sd)

MASS::mvrnorm()


rbeta(number_of_samples_plot, 1.1,1.1) %>% 
  hist(breaks = 100)

### Deviation prior
set.seed(0)

# fatigue sd for each participant
rnorm(number_of_samples_plot,-3,1) %>% 
  exp() %>% 
  hist(breaks = 10000,xlim = c(0,1))
  # median()



# raw_deviation_samples <- rnorm(n = number_of_samples_plot, -3.5, .55) %>%
# # raw_deviation_samples <- rnorm(n = number_of_samples_plot, -1.25, .75) %>%
# # raw_deviation_samples <- rnorm(n = number_of_samples_plot, -0.5,1) %>%
#   data.frame(value = .) %>% 
#   mutate(value_exp = exp(value),
#          value_fatigue_sd = rexp(number_of_samples_plot, rate = 2))

# raw deviations that get prior
raw_deviation_samples %>% 
  ggplot() +
  geom_density(aes(x = value), fill = "gray") +
  # geom_histogram(aes(x = value), bins = 500) +
   # coord_cartesian(xlim = c(0,3)) +
  ggtitle("Prior for Raw Deviation (e.g., σRaw)") +
  theme_classic()+
  theme(text = element_text(family = "arial", size = text_size),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank())

# transformed to normal deviation with exp
raw_deviation_samples %>% 
  ggplot() +
  # geom_density(aes(x = value_exp), fill = "gray", bw = .8) +
  geom_histogram(aes(x = value_exp), fill = "gray", bins = 500) +
   coord_cartesian(xlim = c(0,1)) +
   # coord_cartesian(xlim = c(0,5)) +
  ggtitle("Prior for Deviation After Exp (e.g., σAverage)") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = text_size),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank())

raw_deviation_samples %>% 
  ggplot() +
  geom_histogram(aes(x = value_fatigue_sd), fill = "gray", bins = 1000) +
  # geom_histogram(aes(x = value), bins = 500) +
   coord_cartesian(xlim = c(0,5)) +
  ggtitle("Prior for Deviation After Exp (e.g., σAverage)") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = text_size),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank())
  
# intercept fatigue prior lines
set.seed(0)
number_of_samples_fatigue_plot <- 1000
fatigue_prior_df <- data.frame(intercept = rnorm(n = number_of_samples_fatigue_plot, 0, .75),
                               fatigue = rnorm(n = number_of_samples_fatigue_plot, 0, .01),
                               sample = 1:number_of_samples_fatigue_plot)

fatigue_prior_df %>% 
  ggplot() +
  geom_abline(aes(intercept = intercept, 
                  slope = fatigue),
              size = .2, 
              alpha = .2) +
  scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1)) +
  scale_x_continuous(limits =c(1,176), breaks = seq(0, 176, by = 10)) +
  ggtitle("Fatigue Prior Predictive regression lines") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = text_size),
        axis.ticks.y = element_blank())

fatigue_prior_per_par <- data.frame(intercept = rnorm(n = number_of_samples_fatigue_plot,
                                                      0,
                                                      sample(raw_deviation_samples$value_exp,
                                                             number_of_samples_fatigue_plot,
                                                             replace = F)),
                                    fatigue = rnorm(n = number_of_samples_fatigue_plot,
                                                    0,
                                                    sample(raw_deviation_samples$value_exp,
                                                           number_of_samples_fatigue_plot,
                                                           replace = F)))


fatigue_prior_per_par %>% 
  ggplot() +
  geom_density(aes(x = intercept))
  # geom_abline(aes(intercept = intercept, 
  #                 slope = fatigue),
  #             size = .2, 
  #             alpha = .2) +
  # scale_y_continuous(limits =c(-4,4), breaks = seq(-4, 4, by = 1)) +
  # scale_x_continuous(limits =c(1,176), breaks = seq(0, 176, by = 10)) +
  # ggtitle("Fatigue Prior Predictive regression lines") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = text_size),
        axis.ticks.y = element_blank())

## summaries ####

model001_fit_meta_data <- model001_fit$metadata()

model001_fit_relevant_parameters <- model001_fit_meta_data$model_params[
  !str_detect(model001_fit_meta_data$model_params, "log_lik|mu|amplitude")]
  # !str_detect(model001_fit_meta_data$model_params, "log_lik|mu|amplitude|raw")]

model001_fit_summary <- 
  model001_fit$summary(variables = model001_fit_relevant_parameters)

model005_fit_meta_data <- model005_fit$metadata()

model005_fit_relevant_parameters <- model005_fit_meta_data$model_params[
  !str_detect(model005_fit_meta_data$model_params, "log_lik|mu|amplitude|raw")]

model005_fit_summary <- 
  model005_fit$summary(variables = model005_fit_relevant_parameters)


model005_fit_meta_data <- model005_fit$metadata()

model005_fit_relevant_parameters <- model005_fit_meta_data$model_params[
  !str_detect(model005_fit_meta_data$model_params, "log_lik|mu|amplitude|raw")]

model005_df <- model005_fit$draws(variables = model005_fit_relevant_parameters,
                                  format = "df")
