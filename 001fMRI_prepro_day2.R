# Start of user input ####

# The user may have to add paths to their local github repository and data here
# The code will force you to have it in the correct format

## Add where your local git repository is ####
known_gaborgen_local_gits <- c('~/Documents/GitHub/gaborgen',
                               '/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/code_repository/gaborgen')

## Add where your data is, it should be the same format and organization as it was in the dropbox ####
known_data_locations <- c('~/research_data/gaborgen/raw_data',
                          '/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/raw_data')

## Where the results should be saved
# Andrew Farkas Mac
where_results_should_be_saved <- "/Users/andrewfarkas/research_data/gaborgen/results"
# Andrew Farkas hipergator
# where_results_should_be_saved <- "/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/results"

if(!file.exists(where_results_should_be_saved)){
  stop("Choose a results folder that exists")
}

## Put the participant ID numbers that you would like to preprocess (eg c(118,119))
participants_to_preprocess <- c(123)

# End of user input ####


# Begin checking that everything is correctly organized ####
existing_git_locations <- file.exists(known_gaborgen_local_gits)

number_of_found_gits <- 0

for(path_index in 1:length(known_gaborgen_local_gits)){
  if (file.exists(known_gaborgen_local_gits[path_index])) {
    number_of_found_gits <- number_of_found_gits + 1
    local_git_directory <- known_gaborgen_local_gits[path_index]
  }
}

if (number_of_found_gits == 0) {
  stop('Add your data folder to known_data_locations vector above')
}

if (number_of_found_gits > 1) {
  stop('Multiple data locations found on the current computer')
}

setwd(local_git_directory)


existing_data_locations <- file.exists(known_data_locations)

number_of_found_locations <- 0

for(path_index in 1:length(known_data_locations)){
  if (file.exists(known_data_locations[path_index])) {
    number_of_found_locations <- number_of_found_locations + 1
    data_directory <- known_data_locations[path_index]
  }
}

if (number_of_found_locations == 0) {
  stop('Add your data folder to known_data_locations vector above')
}

if (number_of_found_locations > 1) {
  stop('Multiple data locations found on the current computer')
}


participant_directories <- paste0(data_directory,'/GABORGEN24_DAY2_', 
                                  participants_to_preprocess,"_2")

if (!all(file.exists(participant_directories))){
  stop('Not all the participant directories exist in the correct format')
}

## Check that the correct terminal is available ####
tryCatch({
  # Execute the command and check the exit status
  tcsh_status <- system('tcsh -c "echo tcsh is available"', 
                        intern = TRUE, 
                        ignore.stderr = F)
}, error = function(e) {
  tcsh_status <- "tcsh is unavailable"
  stop("The recommend terminal for afni (tcsh) is not available")
})

# Make stim onset file for every stimulus as prepro sanity check and to get
# nice quality control features
for (participant_index in 1:length(participants_to_preprocess)) {
  
  dat_file_path <- list.files(participant_directories[participant_index],
                              recursive = T,
                              pattern = "logfile.dat$", 
                              full.names = T)
  
  if(!file.exists(dat_file_path)){
    stop(paste0(dat_file_path, " does not exist. You must have the directory in the correct structure"))
  }
  
  log_file <- read.delim(dat_file_path, header = T, sep = ",")
  
  stim_onset <- log_file$timeSinceFirstTR
  
  write(x = stim_onset, 
        file = paste0(participant_directories[participant_index], "/stim_times.1D"),
        ncolumns = 1)
  
  # More stim files for later analyses
  # phase 1 = habituation, 2 = acquisition, 3 = extinction
  # stim 1 = CS+ to 4 = GS3
  # This would be much prettier in an apply function
  phase_4_stim_1_stim_times <- subset(log_file, 
                                      phase == 4 & stim == 1,
                                      select = timeSinceFirstTR,
                                      drop = T)
  phase_4_stim_2_stim_times <- subset(log_file, 
                                      phase == 4 & stim == 2,
                                      select = timeSinceFirstTR,
                                      drop = T)
  phase_4_stim_3_stim_times <- subset(log_file, 
                                      phase == 4 & stim == 3,
                                      select = timeSinceFirstTR,
                                      drop = T)
  phase_4_stim_4_stim_times <- subset(log_file, 
                                      phase == 4 & stim == 4,
                                      select = timeSinceFirstTR,
                                      drop = T)
  
  write(x = phase_4_stim_1_stim_times, 
        file = paste0(participant_directories[participant_index], "/phase_4_stim_1_stim_times.1D"),
        ncolumns = 1)
  write(x = phase_4_stim_2_stim_times, 
        file = paste0(participant_directories[participant_index], "/phase_4_stim_2_stim_times.1D"),
        ncolumns = 1)
  write(x = phase_4_stim_3_stim_times, 
        file = paste0(participant_directories[participant_index], "/phase_4_stim_3_stim_times.1D"),
        ncolumns = 1)
  write(x = phase_4_stim_4_stim_times, 
        file = paste0(participant_directories[participant_index], "/phase_4_stim_4_stim_times.1D"),
        ncolumns = 1)
  
}

# Begin data preprocessing ####

## Convert structural and functional images to nifti fomat, deoblique structural ####
## functional should not need to be deobliqued
source("format_mri.R", local = T, echo = T)

## Preprocess MRI per participant ####
source("afni_proc_py_with_regression_multi_stim.R", local = T, echo = T)

