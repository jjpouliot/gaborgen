# Make stim onset files
for (participant_index in 1:length(participants_to_preprocess)) {
  
  dat_file_path <- paste0(participant_directories[participant_index],
                          "/DAT/gaborgen24_fMRI_Day1_", 
                          participants_to_preprocess[participant_index],
                          "_logfile.dat")
  
  if(!file.exists(dat_file_path)){
    stop(paste0(dat_file_path, " does not exist. You must have the directory in the correct structure"))
  }
  
  log_file <- read.delim(dat_file_path, header = T, sep = ",")
  
  stim_onset <- log_file$timeSinceFirstTR
  #start here
  write(x = stim_onset, 
        file = paste0(participant_directories[participant_index], "/stim_times.1D"),
        ncolumns = 1)
  
}