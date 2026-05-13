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
# participants_to_preprocess <- c(101:155)
# participants_to_preprocess <- c(101:103, 105:155, 157:158)
participants_to_preprocess <- c(
  "101", # made
  "102", # 22% censored# made
  "103", # made
  "106", # made
  "107", # made
  "108", # made
  "109", # made
  "113", # made
  "114", # 23% censored # made
  "115", # made
  "116", # made
  "117", # made
  "119", # made
  "120", # 16% censored # moved over 2 mm # made
  "121", # made
  "122", # made
  "123", # made
  #"124", # 34% censored moved head to every cue
  "125", # made
  "126", # made
  "127", # made
  "128", # made
  "129", # made
  "131", # made
  "132", # made
  "133", # made
  "134", # made
  "135", # made
  # "136", # 50% censored
  "137", # made
  "138", # made
  # "139", # 31% censored, movement over 3mm
  "140", # 23.7% censored # made
  "141", # made
  # double-check these
  #"142", #probably asleep # made
  "143", #probably asleep, not convinced enough to keep them out # made
  #"144", # 29% censored severe TSNR warnings, large pitch shifts above 5 degrees on way out of acquisition
  #
  "145", # made
  #"146", # messed up alignment
  #"147", #asleep , not aligned right?
  #"148", #messed up alignment
  "149", # made
  "150", # made
  "151", # made
  "152", # made
  "153", # made
  "154", # made
  "155", # made
  "158", # made
  "159", # made
  "160", # made
  "161" # made
)

## Which days should be processed
days_to_preprocess <- c(1)
# days_to_preprocess <- c(2)
#days_to_preprocess <- c(1,2)

## Toggle: TRUE = one parameter per trial (IM solution)
##         FALSE = one parameter per cue per block (standard solution)
by_trial <- FALSE

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

  stim_flag <- if (by_trial) "-stim_times_IM" else "-stim_times"

  make_stim_entry <- function(index, label, filename, hrf) {
    paste0(
      "    ",
      stim_flag,
      " ",
      index,
      " ",
      participant_directories[participant_index],
      "/",
      filename,
      " ",
      '"',
      hrf,
      '"',
      "\n",
      "    -stim_label ",
      index,
      " ",
      label
    )
  }

  stim_entries <- paste(
    c(
      make_stim_entry(
        1,
        "habituation_CS+",
        "habituation_CS+_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        2,
        "habituation_GS1",
        "habituation_GS1_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        3,
        "habituation_GS2",
        "habituation_GS2_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        4,
        "habituation_GS3",
        "habituation_GS3_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        5,
        "acquisition_block_1_CS+",
        "acquisition_block_1_CS+_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        6,
        "acquisition_block_1_GS1",
        "acquisition_block_1_GS1_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        7,
        "acquisition_block_1_GS2",
        "acquisition_block_1_GS2_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        8,
        "acquisition_block_1_GS3",
        "acquisition_block_1_GS3_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        9,
        "acquisition_block_2_CS+",
        "acquisition_block_2_CS+_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        10,
        "acquisition_block_2_GS1",
        "acquisition_block_2_GS1_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        11,
        "acquisition_block_2_GS2",
        "acquisition_block_2_GS2_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        12,
        "acquisition_block_2_GS3",
        "acquisition_block_2_GS3_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        13,
        "extinction_CS+",
        "extinction_CS+_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        14,
        "extinction_GS1",
        "extinction_GS1_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        15,
        "extinction_GS2",
        "extinction_GS2_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        16,
        "extinction_GS3",
        "extinction_GS3_stim_times_for_stan.1D",
        "BLOCK(2,1)"
      ),
      make_stim_entry(17, "shock", "shock_times_for_stan.1D", "BLOCK(.1,1)")
    ),
    collapse = "\n"
  )

  xmat_prefix <- if (by_trial) "X_IM" else "X"

  design_matrix_script <- paste0(
    "3dDeconvolve",
    " -nodata 1070 2",
    " -num_stimts 17\n",
    stim_entries,
    "\n",
    "    -x1D ",
    xmat_prefix,
    ".xmat.1D",
    " -xjpeg ",
    xmat_prefix,
    ".jpg",
    " -x1D_uncensored ",
    xmat_prefix,
    ".nocensor.xmat.1D"
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
