
# used sapply so that dat files are in the same order as participant directories
found_dat_log_filepaths <- sapply(participant_directories, function(x){
  list.files(path = paste0(x, '/DAT'),
             pattern = 'logfile.dat$',
             full.names = T)
})


if(length(found_dat_log_filepaths) != length(participant_directories)){
  stop("There is not a dat log file for each participant directory. Directories are not in the proper structure or dat file is missing.")
}

for (directory_index in 1:length(participant_directories)) {
  
  log_file <- read.delim(found_dat_log_filepaths[directory_index], 
                         header = T, sep = ",")
  
  stim_onset <- log_file$timeSinceFirstTR
  
  write(x = stim_onset, 
        file = paste0(participant_directories[directory_index], "/stim_times.1D"),
        ncolumns = 1)
  
  phase_ids <- unique(log_file$phase)
  
  stim_ids <- unique(log_file$stim)
  
  for(phase_index in 1:length(phase_ids)) {
    for(stim_index in 1:length(stim_ids)){
      
      current_phase_stim_times <- subset(log_file, 
                                         phase == phase_ids[phase_index] & 
                                           stim == stim_ids[stim_index],
                                         select = timeSinceFirstTR,
                                         drop = T)
      
      write(x = current_phase_stim_times, 
            ncolumns = 1,
            file = paste0(participant_directories[directory_index], 
                          paste0("/phase_", 
                                 phase_ids[phase_index], 
                                 "_stim_",
                                 stim_ids[stim_index],
                                 "_stim_times.1D")))
      
    }
  }
}
