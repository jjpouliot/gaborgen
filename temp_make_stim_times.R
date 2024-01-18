
dat_file_path <- "/Users/andrewfarkas/research_data/gaborgen/raw_data/GABORGEN24_120/DAT/gaborgen24_fMRI_Day1_120_logfile.dat"


log_file <- read.delim(dat_file_path,header = T, sep = ",")

stim_onset_and_duration <- log_file$timeSinceFirstTR
  
write_delim(stim_onset_and_duration, 
      file = "/Users/andrewfarkas/research_data/gaborgen/raw_data/GABORGEN24_120/stim_times.1D",
      col_names = F)
