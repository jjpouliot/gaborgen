## rerun_afni_stats.R
## Parses a stats.REML_cmd file and reruns 3dDeconvolve + 3dREMLfit.
##
## Usage: Rscript rerun_afni_stats.R
## ---------------------------------------------------------------

# Rerun make stim files from a001fMRI_prepro.R
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
participants_to_preprocess <- c(133)

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
  current_results_folder <- dir(
    where_results_should_be_saved,
    pattern = paste0(
      basename(participant_directories[participant_index]),
      ".results"
    ),
    full.names = T
  )

  # replace stimuli files in results folder from those in raw_data folder
  file.copy(
    from = list.files(
      path = participant_directories[participant_index],
      pattern = "phase",
      full.names = T
    ),
    to = paste0(current_results_folder, "/stimuli"),
    overwrite = T
  )

  # move to results directory and rerun regression
  setwd(current_results_folder)

  if (grepl("DAY2", participant_directories[participant_index])) {
    regression_terminal_script = paste0(
      '3dDeconvolve -input ',
      paste0(
        'pb05.',
        basename(participant_directories[participant_index]),
        '.r01.scale+tlrc.HEAD',
        " "
      ),
      '-censor ',
      paste0(
        'censor_',
        basename(participant_directories[participant_index]),
        '_combined_2.1D',
        " "
      ),
      '-ortvec mot_demean.r01.1D mot_demean_r01 
  -polort 15 -num_stimts 4 
  -stim_times 1 stimuli/phase_4_stim_1_stim_times.1D "BLOCK(2,1)" -stim_label 1 phase_4_stim_1
  -stim_times 2 stimuli/phase_4_stim_2_stim_times.1D "BLOCK(2,1)" -stim_label 2 phase_4_stim_2
  -stim_times 3 stimuli/phase_4_stim_3_stim_times.1D "BLOCK(2,1)" -stim_label 3 phase_4_stim_3
  -stim_times 4 stimuli/phase_4_stim_4_stim_times.1D "BLOCK(2,1)" -stim_label 4 phase_4_stim_4
  -fout -tout -x1D X.xmat_fixed_T1_off.1D -xjpeg X_fixed_T1_off.jpg -x1D_uncensored X.nocensor.xmat.1D
  -errts ',
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
  } else {
    regression_terminal_script = paste0(
      '3dDeconvolve -input ',
      paste0(
        'pb05.',
        basename(participant_directories[participant_index]),
        '.r01.scale+tlrc.HEAD',
        " "
      ),
      '-censor ',
      paste0(
        'censor_',
        basename(participant_directories[participant_index]),
        '_combined_2.1D',
        " "
      ),
      '-ortvec mot_demean.r01.1D mot_demean_r01 
  -polort 15 -num_stimts 12 
  -stim_times 1 stimuli/phase_1_stim_1_stim_times.1D "BLOCK(2,1)" -stim_label 1 phase_1_stim_1 
  -stim_times 2 stimuli/phase_1_stim_2_stim_times.1D "BLOCK(2,1)" -stim_label 2 phase_1_stim_2 
  -stim_times 3 stimuli/phase_1_stim_3_stim_times.1D "BLOCK(2,1)" -stim_label 3 phase_1_stim_3 
  -stim_times 4 stimuli/phase_1_stim_4_stim_times.1D "BLOCK(2,1)" -stim_label 4 phase_1_stim_4 
  -stim_times 5 stimuli/phase_2_stim_1_stim_times.1D "BLOCK(2,1)" -stim_label 5 phase_2_stim_1 
  -stim_times 6 stimuli/phase_2_stim_2_stim_times.1D "BLOCK(2,1)" -stim_label 6 phase_2_stim_2 
  -stim_times 7 stimuli/phase_2_stim_3_stim_times.1D "BLOCK(2,1)" -stim_label 7 phase_2_stim_3 
  -stim_times 8 stimuli/phase_2_stim_4_stim_times.1D "BLOCK(2,1)" -stim_label 8 phase_2_stim_4 
  -stim_times 9 stimuli/phase_3_stim_1_stim_times.1D "BLOCK(2,1)" -stim_label 9 phase_3_stim_1 
  -stim_times 10 stimuli/phase_3_stim_2_stim_times.1D "BLOCK(2,1)" -stim_label 10 phase_3_stim_2 
  -stim_times 11 stimuli/phase_3_stim_3_stim_times.1D "BLOCK(2,1)" -stim_label 11 phase_3_stim_3 
  -stim_times 12 stimuli/phase_3_stim_4_stim_times.1D "BLOCK(2,1)" -stim_label 12 phase_3_stim_4 
  -fout -tout -x1D X.xmat_fixed_T1_off.1D -xjpeg X_fixed_T1_off.jpg -x1D_uncensored X.nocensor.xmat.1D 
  -errts ',
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
    '3dREMLfit -matrix X.xmat_fixed_T1_off.1D -input ',
    paste0(
      'pb05.',
      basename(participant_directories[participant_index]),
      '.r01.scale+tlrc.HEAD',
      " "
    ),
    '-fout -tout -Rbuck ',
    paste0(
      'stats.',
      basename(participant_directories[participant_index]),
      '_fixed_T1_off_REML',
      " "
    ),
    '-Rvar ',
    paste0(
      'stats.',
      basename(participant_directories[participant_index]),
      '_fixed_T1_off_REMLvar',
      " "
    ),
    '-Rerrts ',
    paste0(
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
