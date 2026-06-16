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

## Toggle (only applies when by_trial = FALSE):
##   TRUE  = find the earliest CS+ onset across all of acquisition and treat
##           all cue trials (for every cue) that occurred at or before that
##           onset as habituation. The rationale is that no conditioning can
##           have occurred until the first CS+/US pairing is possible, so any
##           GS or CS+ trials preceding that moment should behave like
##           naive habituation responses. Modified stim-times files are written
##           with a "_first_to_hab" suffix and the output xmat prefix also gets
##           "_first_to_hab" appended, so originals are preserved.
##   FALSE = standard treatment (no onsets are reassigned across phases)
first_acq_to_hab <- TRUE

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

# When first_acq_to_hab is requested, write modified stim-times files.
# The cutoff is the earliest CS+ onset anywhere in acquisition. Every cue
# onset (across both acquisition blocks) that falls at or before that time
# is moved into the habituation regressor for that cue. Separate
# "_first_to_hab" files are written for habituation, acquisition block 1,
# and acquisition block 2; all other files (extinction, shock) are unchanged.
if (!by_trial && first_acq_to_hab) {
  cues <- c("CS+", "GS1", "GS2", "GS3")

  # Small helper: write a stim-times file, using an asterisk placeholder if
  # there are no remaining onsets (AFNI requires a non-empty file).
  write_stim_file <- function(times, path) {
    if (length(times) == 0) {
      writeLines("*", path)
    } else {
      write(times, ncolumns = 1, file = path)
    }
  }

  for (participant_index in 1:length(participant_directories)) {
    pdir <- participant_directories[participant_index]

    # Temporal cutoff: the first CS+ presentation anywhere in acquisition
    acq1_csplus_times <- scan(
      file.path(pdir, "acquisition_block_1_CS+_stim_times_for_stan.1D"),
      quiet = TRUE
    )
    acq2_csplus_times <- scan(
      file.path(pdir, "acquisition_block_2_CS+_stim_times_for_stan.1D"),
      quiet = TRUE
    )
    first_csplus_time <- min(c(acq1_csplus_times, acq2_csplus_times))

    message(
      basename(pdir),
      ": first CS+ onset in acquisition = ",
      first_csplus_time,
      "s"
    )

    for (cue in cues) {
      hab_times <- scan(
        file.path(pdir, paste0("habituation_", cue, "_stim_times_for_stan.1D")),
        quiet = TRUE
      )
      acq1_times <- scan(
        file.path(
          pdir,
          paste0("acquisition_block_1_", cue, "_stim_times_for_stan.1D")
        ),
        quiet = TRUE
      )
      acq2_times <- scan(
        file.path(
          pdir,
          paste0("acquisition_block_2_", cue, "_stim_times_for_stan.1D")
        ),
        quiet = TRUE
      )

      # Any acquisition onset at or before the first CS+ goes to habituation
      hab_modified <- sort(c(
        hab_times,
        acq1_times[acq1_times <= first_csplus_time],
        acq2_times[acq2_times <= first_csplus_time]
      ))
      acq1_modified <- acq1_times[acq1_times > first_csplus_time]
      acq2_modified <- acq2_times[acq2_times > first_csplus_time]

      write_stim_file(
        hab_modified,
        file.path(
          pdir,
          paste0("habituation_", cue, "_stim_times_for_stan_first_to_hab.1D")
        )
      )
      write_stim_file(
        acq1_modified,
        file.path(
          pdir,
          paste0(
            "acquisition_block_1_",
            cue,
            "_stim_times_for_stan_first_to_hab.1D"
          )
        )
      )
      write_stim_file(
        acq2_modified,
        file.path(
          pdir,
          paste0(
            "acquisition_block_2_",
            cue,
            "_stim_times_for_stan_first_to_hab.1D"
          )
        )
      )
    }
  }
}

# Make design matrix minus movement
for (participant_index in 1:length(participant_directories)) {
  setwd(participant_directories[participant_index])

  stim_flag <- if (by_trial) "-stim_times_IM" else "-stim_times"

  # Suffix for stim-times files that are modified by first_acq_to_hab.
  # Extinction and shock files are never modified so they always use "".
  fth_suffix <- if (!by_trial && first_acq_to_hab) "_first_to_hab" else ""

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
        paste0("habituation_CS+_stim_times_for_stan", fth_suffix, ".1D"),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        2,
        "habituation_GS1",
        paste0("habituation_GS1_stim_times_for_stan", fth_suffix, ".1D"),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        3,
        "habituation_GS2",
        paste0("habituation_GS2_stim_times_for_stan", fth_suffix, ".1D"),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        4,
        "habituation_GS3",
        paste0("habituation_GS3_stim_times_for_stan", fth_suffix, ".1D"),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        5,
        "acquisition_block_1_CS+",
        paste0(
          "acquisition_block_1_CS+_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        6,
        "acquisition_block_1_GS1",
        paste0(
          "acquisition_block_1_GS1_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        7,
        "acquisition_block_1_GS2",
        paste0(
          "acquisition_block_1_GS2_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        8,
        "acquisition_block_1_GS3",
        paste0(
          "acquisition_block_1_GS3_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        9,
        "acquisition_block_2_CS+",
        paste0(
          "acquisition_block_2_CS+_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        10,
        "acquisition_block_2_GS1",
        paste0(
          "acquisition_block_2_GS1_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        11,
        "acquisition_block_2_GS2",
        paste0(
          "acquisition_block_2_GS2_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
        "BLOCK(2,1)"
      ),
      make_stim_entry(
        12,
        "acquisition_block_2_GS3",
        paste0(
          "acquisition_block_2_GS3_stim_times_for_stan",
          fth_suffix,
          ".1D"
        ),
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

  xmat_prefix <- paste0(
    basename(participant_directories[participant_index]),
    "_",
    if (by_trial) "X_IM" else "X",
    fth_suffix
  )

  design_matrix_script <- paste0(
    "3dDeconvolve",
    " -nodata 1070 2",
    " -polort 15",
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

  design_matrices_dir <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/design_matrices"
  nocensor_file <- file.path(
    participant_directories[participant_index],
    paste0(xmat_prefix, ".nocensor.xmat.1D")
  )
  file.copy(nocensor_file, design_matrices_dir, overwrite = TRUE)
}
