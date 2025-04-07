data_location <- '/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results'
where_to_copy <- '/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder'

# copy design matrices
full_design_matrix_paths <- list.files(path = data_location, pattern = 'X.nocensor.xmat.1D$', recursive = T)

full_design_matrix_new_names <- sub(replacement = "_", pattern = "/", x = full_design_matrix_paths)

file.copy(from = full_design_matrix_paths,to = paste0(where_to_copy, "/", full_design_matrix_new_names))

# copy censor info
censor_info_paths <- list.files(path = data_location, pattern = 'combined_2.1D$', recursive = T)

file.copy(from = censor_info_paths, to = paste0(where_to_copy, "/", censor_info_paths))


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


### tests
system2('tcsh', 
        args = c('-c', shQuote(paste('dcm2niix', current_functional_directory))))




