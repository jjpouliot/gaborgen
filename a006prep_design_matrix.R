# The plan for this script is to create a design matrix for the hemodynamic response per each stimulus and each shock.
# So lots of columns that are going to will be grouped by some multilevelness.
# I am going to organize the columns by block, cue, then trial (eg. hab cs+ trial_1, hab cs+ trial_2, ...)

# data_location <- '/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results'
# where_to_copy <- '/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder'

# data_location <- '/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results'
# where_to_copy <- '/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder'

# borrowed from a001 so that the stim onsets are made and environmental variables are set
known_gaborgen_local_gits <- c(
  '~/Documents/GitHub/gaborgen',
  '/home/andrewfarkas/Repositories/gaborgen',
  '/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/code_repository/gaborgen'
)

## Add where your data is, it should be the same format and organization as it was in the dropbox ####
known_data_locations <- c(
  '~/research_data/gaborgen/raw_data',
  '/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data',
  '/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/raw_data',
  '/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data'
)

## Where the results should be saved
# Andrew Farkas Mac
where_results_should_be_saved <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/fmri_preprocessed"
# where_results_should_be_saved <- "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/fmri_preprocessed"
# where_results_should_be_saved <- "/Users/andrewfarkas/research_data/gaborgen/results"
# Andrew Farkas hipergator
# where_results_should_be_saved <- "/blue/akeil/andrew.farkas/gaborgen24_eeg_fmri/results"

if (!file.exists(where_results_should_be_saved)) {
  stop("Choose a results folder that exists")
}

## Put the participant ID numbers that you would like to preprocess (eg c(118,119))
# participants_to_preprocess <- c(101:103,106:130, 135:136)
participants_to_preprocess <- c(101:155)
participants_to_preprocess <- c(101:103, 105:155, 157:158)

## Which days should be processed
days_to_preprocess <- c(1)
# days_to_preprocess <- c(2)
#days_to_preprocess <- c(1,2)

# End of user input ####

# Begin checking that everything is correctly organized ####
existing_git_locations <- file.exists(known_gaborgen_local_gits)

number_of_found_gits <- 0

for (path_index in 1:length(known_gaborgen_local_gits)) {
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

for (path_index in 1:length(known_data_locations)) {
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
  if (
    participants_to_preprocess[participant_index] < 123 &
      1 %in% days_to_preprocess
  ) {
    participant_directories <- c(
      participant_directories,
      paste0(
        data_directory,
        '/GABORGEN24_',
        participants_to_preprocess[
          participant_index
        ]
      )
    )
  } else if (participants_to_preprocess[participant_index] >= 123) {
    participant_directories <- c(
      participant_directories,
      paste0(
        data_directory,
        '/GABORGEN24_DAY',
        days_to_preprocess,
        '_',
        participants_to_preprocess[
          participant_index
        ]
      )
    )
  }
}


if (!all(file.exists(participant_directories))) {
  paste0(
    "Directories that don't exist: ",
    participant_directories[!file.exists(participant_directories)]
  )
  stop('Not all the participant directories exist in the correct format')
}

## Check that the correct terminal is available ####
tryCatch(
  {
    # Execute the command and check the exit status
    tcsh_status <- system(
      'tcsh -c "echo tcsh is available"',
      intern = TRUE,
      ignore.stderr = F
    )
  },
  error = function(e) {
    tcsh_status <- "tcsh is unavailable"
    stop("The recommend terminal for afni (tcsh) is not available")
  }
)

# Make stim onset files
by_block <- T
shock_times <- T
source("make_stim_times_files.R", local = T, echo = T)

# Make design matrix minus movement
for (participant_index in 1:length(participant_directories)) {
  setwd(participant_directories[participant_index])

  design_matrix_script <- paste0(
    '3dDeconvolve
    -nodata 1070 2
    -num_stimts 17
    -stim_times_IM 1 ',
    paste0(
      participant_directories[participant_index],
      "/habituation_CS+_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 1 habituation_CS+
    -stim_times_IM 2 ',
    paste0(
      participant_directories[participant_index],
      "/habituation_GS1_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 2 habituation_GS1
    -stim_times_IM 3 ',
    paste0(
      participant_directories[participant_index],
      "/habituation_GS2_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 3 habituation_GS2
    -stim_times_IM 4 ',
    paste0(
      participant_directories[participant_index],
      "/habituation_GS3_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 4 habituation_GS3
    -stim_times_IM 5 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_1_CS+_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 5 acquisition_block_1_CS+
    -stim_times_IM 6 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_1_GS1_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 6 acquisition_block_1_GS1
    -stim_times_IM 7 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_1_GS2_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 7 acquisition_block_1_GS2
    -stim_times_IM 8 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_1_GS3_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 8 acquisition_block_1_GS3
    -stim_times_IM 9 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_2_CS+_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 9 acquisition_block_2_CS+
    -stim_times_IM 10 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_2_GS1_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 10 acquisition_block_2_GS1
    -stim_times_IM 11 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_2_GS2_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 11 acquisition_block_2_GS2
    -stim_times_IM 12 ',
    paste0(
      participant_directories[participant_index],
      "/acquisition_block_2_GS3_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 12 acquisition_block_2_GS3
    -stim_times_IM 13 ',
    paste0(
      participant_directories[participant_index],
      "/extinction_CS+_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 13 extinction_CS+
    -stim_times_IM 14 ',
    paste0(
      participant_directories[participant_index],
      "/extinction_GS1_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 14 extinction_GS1
    -stim_times_IM 15 ',
    paste0(
      participant_directories[participant_index],
      "/extinction_GS2_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 15 extinction_GS2
    -stim_times_IM 16 ',
    paste0(
      participant_directories[participant_index],
      "/extinction_GS3_stim_times_for_stan.1D "
    ),
    '"BLOCK(2,1)"
    -stim_label 16 extinction_GS3
    -stim_times_IM 17 ',
    paste0(
      participant_directories[participant_index],
      "/shock_times_for_stan.1D "
    ),
    '"BLOCK(.1,1)"
    -stim_label 17 shock
    -x1D X_IM.xmat.1D -xjpeg X_IM.jpg -x1D_uncensored X_IM.nocensor.xmat.1D'
  )
  design_matrix_script <- gsub(
    pattern = "\n",
    replacement = "",
    x = design_matrix_script
  )

  system2('tcsh', args = c('-c', shQuote('source $HOME/.tcshrc')))

  tryCatch(
    {
      system2('tcsh', args = c('-c', shQuote(design_matrix_script)))
    },
    error = function(e) {
      # Display the actual error message
      message("Error encountered: ", e$message)

      # Custom stop message
      stop(
        "If python error, then you likely need to change your terminal configuration files by running this code in the R console: "
      )
    }
  )
}
