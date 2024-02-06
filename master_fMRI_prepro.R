# Start of user input ####

# The user may have to add paths to their local github repository and data here
# The code will force you to have it in the correct format

## Add where your local git repository is ####
known_gaborgen_local_gits <- c('~/Documents/GitHub/gaborgen',
                               '/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/code_repository/gaborgen')

## Add where your data is, it should be the same format and organization as it was in the dropbox ####
known_data_locations <- c('~/research_data/gaborgen/raw_data',
                          '/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/raw_data')

## Put the participant ID numbers that you would like to preprocess (eg c(118,119))
participants_to_preprocess <- c(122)



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


participant_directories <- paste0(data_directory,'/GABORGEN24_', 
                                  participants_to_preprocess)

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




if(
  system2('tcsh', 
        args = c('-c', '"echo tcsh is available"'), 
        stdout = TRUE, 
        stderr = TRUE) != 
  "tcsh is available") {
  stop('The tcsh terminal is not available. This is the prefered terminal for AFNI')
}

# Begin data preprocessing ####

## Convert structural and functional images to nifti fomat, deoblique structural ####
## functional should not need to be deobliqued
source("format_mri.R", local = T, echo = T)

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

## Preprocess MRI per participant ####
where_results_should_be_saved <- "/Users/andrewfarkas/research_data/gaborgen/results"

source("afni_proc_py_prepro.R", local = T, echo = T)
