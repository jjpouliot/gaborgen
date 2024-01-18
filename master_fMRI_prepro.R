# Add where your local git repository is
known_gaborgen_local_gits <- c('~/Documents/GitHub/gaborgen')

# Add where you data is, it should be the same format and organization as it was in the dropbox
known_data_locations <- c('~/research_data/gaborgen/raw_data/')

# Put the participant ID numbers that you would like to preprocess (eg c(118,119))
participants_to_preprocess <- c(119)


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




participant_directories <- paste0(data_directory,'GABORGEN24_', 
                                  participants_to_preprocess)

if (!all(file.exists(participant_directories))){
  stop('Not all the participant directories exist in the correct format')
}

# Check that the correct terminal is available
if(
  system2('tcsh', 
        args = c('-c', '"echo tcsh is available"'), 
        stdout = TRUE, 
        stderr = TRUE) != 
  "tcsh is available") {
  stop('The tcsh terminal is not available. This is the prefered terminal for AFNI')
}

# Convert structural and functional images to nifti fomat, deoblique structural
# functional should not need to be deobliqued
source("format_mri.R", local = T, echo = T)

# Make stim onset files
for (participant_index in 1:length(participants_to_preprocess)) {
  
  dat_file_path <- paste0(participant_directories[participant_index],
                          "/DAT/gaborgen24_fMRI_Day1_", 
                          participants_to_preprocess[participant_index],
                          "_logfile.dat")
  
  file.exists(dat_file_path)
  
  log_file <- read.delim(dat_file_path, header = T, sep = ",")
  
  stim_onset_and_duration <- log_file$timeSinceFirstTR
  #start here
  write_delim(stim_onset_and_duration, 
              file = "/Users/andrewfarkas/research_data/gaborgen/raw_data/GABORGEN24_120/stim_times.1D",
              col_names = F)
  
}

# Preprocess MRI per participant



