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
# participants_to_preprocess <- c(101:103,106:130, 135:136)
participants_to_preprocess <- c(145)

## Which days should be processed
# days_to_preprocess <- c(1)
# days_to_preprocess <- c(2)
days_to_preprocess <- c(1,2)

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

participant_directories <- c()
for (participant_index in 1:length(participants_to_preprocess)) {
  
  if (participants_to_preprocess[participant_index] < 123 & 1 %in% days_to_preprocess) {
    
    participant_directories <- c(participant_directories,
                                 paste0(data_directory,
                                        '/GABORGEN24_', 
                                        participants_to_preprocess[
                                          participant_index]))
    
  } else if (participants_to_preprocess[participant_index] >= 123) {
    
    
    
    participant_directories <- c(participant_directories,
                                 paste0(data_directory,
                                        '/GABORGEN24_DAY',
                                        days_to_preprocess,
                                        '_',
                                        participants_to_preprocess[
                                          participant_index]))
    
  }
}



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

# Make stim onset files
source("make_stim_times_files.R", local = T, echo = T)

# Begin data preprocessing ####

## Convert structural and functional images to nifti fomat, deoblique structural ####
## functional should not need to be deobliqued
source("format_mri.R", local = T, echo = T)

## Preprocess MRI per participant ####
source("afni_proc_py_with_regression_multi_stim.R", local = T, echo = T)

