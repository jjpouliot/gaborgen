## rerun_afni_stats.R
## Parses a stats.REML_cmd file and reruns 3dDeconvolve + 3dREMLfit.
##
## Usage: Rscript rerun_afni_stats.R
## ---------------------------------------------------------------

# Rerun make stim files from a001fMRI_prepro.R
# This is a special version that does a slightly different regression with the cues by block and with shock as a predictor
# To be used as a whole-brain comparison to the Bayesian GP models

# Start of user input ####

# The user may have to add paths to their local github repository and data here
# The code will force you to have it in the correct format

## Add where your local git repository is ####
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
# participants_to_preprocess <- c(101:103, 105:158)

# test day 2
participants_to_preprocess <- c(
  # "123" #,
  "125",
  "126",
  "127",
  "128",
  "129",
  "130",
  "132",
  "133",
  "135",
  "136",
  "137",
  "138",
  "139",
  "141",
  "142",
  "143",
  "144",
  "145",
  "147",
  "150",
  "151",
  "153",
  "154",
  "159",
  "160"
)

#day 1 used
# participants_to_preprocess <- c(
#   "101", # made
#   "102", # 22% censored# made
#   "103", # made
#   "106", # made
#   "107", # made
#   "108", # made
#   "109", # made
#   "113", # made
#   "114", # 23% censored # made
#   "115", # made
#   "116", # made
#   "117", # made
#   "119", # made
#   "120", # 16% censored # moved over 2 mm # made
#   "121", # made
#   "122", # made
#   "123", # made
#   #"124", # 34% censored moved head to every cue
#   "125", # made
#   "126", # made
#   "127", # made
#   "128", # made
#   "129", # made
#   "131", # made
#   "132", # made
#   "133", # made
#   "134", # made
#   "135", # made
#   # "136", # 50% censored
#   "137", # made
#   "138", # made
#   # "139", # 31% censored, movement over 3mm
#   "140", # 23.7% censored # made
#   "141", # made
#   # double-check these
#   #"142", #probably asleep # made
#   "143", #probably asleep, not convinced enough to keep them out # made
#   #"144", # 29% censored severe TSNR warnings, large pitch shifts above 5 degrees on way out of acquisition
#   #
#   "145", # made
#   #"146", # messed up alignment
#   #"147", #asleep , not aligned right?
#   #"148", #messed up alignment
#   "149", # made
#   "150", # made
#   "151", # made
#   "152", # made
#   "153", # made
#   "154", # made
#   "155", # made
#   "158", # made
#   "159", # made
#   "160", # made
#   "161" # made
# )

## Which days should be processed
# days_to_preprocess <- c(1)
days_to_preprocess <- c(2)
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
source("make_stim_times_files.R", local = T, echo = T)


# Rerun regression for each participant

for (participant_index in 1:length(participant_directories)) {
  # get results directory
  # not necessary because I am not going into each results folder, just one big folder
  # current_results_folder <- dir(
  #   where_results_should_be_saved,
  #   pattern = paste0(
  #     basename(participant_directories[participant_index]),
  #     ".results"
  #   ),
  #   full.names = T
  # )

  # # replace stimuli files in results folder from those in raw_data folder
  # file.copy(
  #   from = list.files(
  #     path = participant_directories[participant_index],
  #     pattern = "phase",
  #     full.names = T
  #   ),
  #   to = paste0(current_results_folder, "/stimuli"),
  #   overwrite = T
  # )

  # move to results directory and rerun regression
  # current_results_folder <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info"
  setwd(where_results_should_be_saved)
  #day 2 not done
  if (grepl("DAY2", participant_directories[participant_index])) {
    regression_terminal_script = paste0(
      '3dDeconvolve -input ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        ".results/",
        'pb05.',
        basename(participant_directories[participant_index]),
        '.r01.scale+tlrc.HEAD',
        " "
      ),
      '-censor ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        ".results/",
        'censor_',
        basename(participant_directories[participant_index]),
        '_combined_2.1D',
        " "
      ),
      '-ortvec ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        ".results/",
        'mot_demean.r01.1D mot_demean_r01'
      ),
      ' -polort 15 -num_stimts 4 
  -stim_times 1 ',
      paste0(
        participant_directories[participant_index],
        '/phase_4_stim_1_stim_times.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 1 csp
  -stim_times 2 ',
      paste0(
        participant_directories[participant_index],
        '/phase_4_stim_2_stim_times.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 2 gs1
  -stim_times 3 ',
      paste0(
        participant_directories[participant_index],
        '/phase_4_stim_3_stim_times.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 3 gs2
  -stim_times 4 ',
      paste0(
        participant_directories[participant_index],
        '/phase_4_stim_4_stim_times.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 4 gs3 
  -fout -tout -x1D ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        '_X.xmat_fixed_T1_off.1D'
      ),
      ' -xjpeg ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        '_X_fixed_T1_off.jpg'
      ),
      ' -x1D_uncensored ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        '_X.nocensor.xmat.1D'
      ),
      ' -errts ',
      paste0(
        where_results_should_be_saved,
        "/",
        'errts.',
        basename(participant_directories[participant_index]),
        '_fixed_T1_off',
        " "
      ),
      '-bucket ',
      paste0(
        where_results_should_be_saved,
        "/",
        'stats.',
        basename(participant_directories[participant_index]),
        '_fixed_T1_off'
      )
    )

    regression_terminal_script <- gsub(
      pattern = "\n",
      replacement = "",
      x = regression_terminal_script
    )

    system2('tcsh', args = c('-c', shQuote('source $HOME/.tcshrc')))

    tryCatch(
      {
        system2('tcsh', args = c('-c', shQuote(regression_terminal_script)))
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
  } else {
    regression_terminal_script = paste0(
      '3dDeconvolve -input ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        "/",
        'pb05.',
        basename(participant_directories[participant_index]),
        '.r01.scale+tlrc.HEAD',
        " "
      ),
      '-censor ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        "/",
        'censor_',
        basename(participant_directories[participant_index]),
        '_combined_2.1D',
        " "
      ),
      '-ortvec ',
      paste0(
        where_results_should_be_saved,
        "/",
        basename(participant_directories[participant_index]),
        "/",
        'mot_demean.r01.1D mot_demean_r01'
      ),
      ' -polort 15 -num_stimts 17 
  -stim_times 1 ',
      paste0(
        participant_directories[participant_index],
        '/habituation_CS+_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 1 hab_csp
  -stim_times 2 ',
      paste0(
        participant_directories[participant_index],
        '/habituation_GS1_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 2 hab_gs1
  -stim_times 3 ',
      paste0(
        participant_directories[participant_index],
        '/habituation_GS2_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 3 hab_gs2
  -stim_times 4 ',
      paste0(
        participant_directories[participant_index],
        '/habituation_GS3_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 4 hab_gs3 
  -stim_times 5 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_1_CS+_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 5 acq_1_csp 
  -stim_times 6 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_1_GS1_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 6 acq_1_gs1 
  -stim_times 7 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_1_GS2_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 7 acq_1_gs2 
  -stim_times 8 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_1_GS3_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 8 acq_1_gs3 
  -stim_times 9 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_2_CS+_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 9 acq_2_csp 
  -stim_times 10 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_2_GS1_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 10 acq_2_gs1 
  -stim_times 11 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_2_GS2_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 11 acq_2_gs2
  -stim_times 12 ',
      paste0(
        participant_directories[participant_index],
        '/acquisition_block_2_GS3_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 12 acq_2_gs3 
  -stim_times 13 ',
      paste0(
        participant_directories[participant_index],
        '/extinction_CS+_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 13 ext_csp
  -stim_times 14 ',
      paste0(
        participant_directories[participant_index],
        '/extinction_GS1_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 14 ext_gs1
  -stim_times 15 ',
      paste0(
        participant_directories[participant_index],
        '/extinction_GS2_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 15 ext_gs2
  -stim_times 16 ',
      paste0(
        participant_directories[participant_index],
        '/extinction_GS3_stim_times_for_stan.1D'
      ),
      ' "BLOCK(2,1)" -stim_label 16 ext_gs3
  -stim_times 17 ',
      paste0(
        participant_directories[participant_index],
        '/shock_times_for_stan.1D'
      ),
      ' "BLOCK(0.1,1)" -stim_label 17 shock  
  -fout -tout -x1D ',
      paste0(
        basename(participant_directories[participant_index]),
        '_X.xmat_fixed_T1_off.1D'
      ),
      ' -xjpeg ',
      paste0(
        basename(participant_directories[participant_index]),
        '_X_fixed_T1_off.jpg'
      ),
      ' -x1D_uncensored ',
      paste0(
        basename(participant_directories[participant_index]),
        '_X.nocensor.xmat.1D'
      ),
      ' -errts ',
      paste0(
        'errts.',
        basename(participant_directories[participant_index]),
        '_fixed_T1_off',
        " "
      ),
      '-bucket ',
      paste0(
        'stats.',
        basename(participant_directories[participant_index]),
        '_fixed_T1_off'
      )
    )

    regression_terminal_script <- gsub(
      pattern = "\n",
      replacement = "",
      x = regression_terminal_script
    )

    system2('tcsh', args = c('-c', shQuote('source $HOME/.tcshrc')))

    tryCatch(
      {
        system2('tcsh', args = c('-c', shQuote(regression_terminal_script)))
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

  REML_regression_terminal_script = paste0(
    '3dREMLfit -matrix ',
    paste0(
      where_results_should_be_saved,
      "/",
      basename(participant_directories[participant_index]),
      '_X.xmat_fixed_T1_off.1D'
    ),
    ' -input ',
    paste0(
      where_results_should_be_saved,
      "/",
      basename(participant_directories[participant_index]),
      ".results/",
      'pb05.',
      basename(participant_directories[participant_index]),
      '.r01.scale+tlrc.HEAD',
      " "
    ),
    '-fout -tout -Rbuck ',
    paste0(
      where_results_should_be_saved,
      "/",
      'stats.',
      basename(participant_directories[participant_index]),
      '_fixed_T1_off_REML',
      " "
    ),
    '-Rvar ',
    paste0(
      where_results_should_be_saved,
      "/",
      'stats.',
      basename(participant_directories[participant_index]),
      '_fixed_T1_off_REMLvar',
      " "
    ),
    '-Rerrts ',
    paste0(
      where_results_should_be_saved,
      "/",
      'errts.',
      basename(participant_directories[participant_index]),
      '_fixed_T1_off_REML',
      " "
    ),
    '-verb $*'
  )

  REML_regression_terminal_script <- gsub(
    pattern = "\n",
    replacement = "",
    x = REML_regression_terminal_script
  )

  system2('tcsh', args = c('-c', shQuote('source $HOME/.tcshrc')))

  tryCatch(
    {
      system2('tcsh', args = c('-c', shQuote(REML_regression_terminal_script)))
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
