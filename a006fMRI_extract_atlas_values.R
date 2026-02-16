# data_location <- '/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results'
# where_to_copy <- '/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder'
data_location <- '/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/fmri_preprocessed'
where_to_copy <- '/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/roi_timeseries'

# copy design matrices
full_design_matrix_paths <- list.files(
  path = data_location,
  pattern = 'X.nocensor.xmat.1D$',
  recursive = T
)

full_design_matrix_new_names <- sub(
  replacement = "_",
  pattern = "/",
  x = full_design_matrix_paths
)

file.copy(
  from = paste0(data_location, "/", full_design_matrix_paths),
  to = paste0(where_to_copy, "/", full_design_matrix_new_names)
)

# copy censor info
censor_info_paths <- list.files(
  path = data_location,
  pattern = 'combined_2.1D$',
  recursive = T
)

file.copy(
  from = paste0(data_location, "/", censor_info_paths),
  to = paste0(where_to_copy, "/", basename(censor_info_paths))
)


# copy functional data

processed_functional_paths <- list.files(
  path = data_location,
  pattern = "pb05",
  recursive = T
)

file.copy(
  from = paste0(data_location, "/", processed_functional_paths),
  to = paste0(where_to_copy, "/", basename(processed_functional_paths))
)


# calculate ROI times series per participant, just day 1 right now

participants_to_process <- c(120, 137, 138, 139, 144, 145, 151, 153, 155, 158)

for (i in 1:length(participants_to_process)) {
  setwd(where_to_copy) # go to ROI folder

  # run this terminal command, expect per participant, gets timeseries per ROI into text file
  #   3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_149.r01.scale+tlrc > 149_roi_stats.txt

  if (participants_to_process[i] < 123) {
    system2(
      'tcsh',
      args = c(
        '-c',
        shQuote(paste(
          '3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_',
          participants_to_process[i],
          ".r01.scale+tlrc > ",
          participants_to_process[i],
          "_roi_stats.txt",
          sep = ""
        ))
      )
    )
  } else {
    system2(
      'tcsh',
      args = c(
        '-c',
        shQuote(paste(
          '3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_',
          participants_to_process[i],
          ".r01.scale+tlrc > ",
          participants_to_process[i],
          "_roi_stats.txt",
          sep = ""
        ))
      )
    )
  }
}

# Make design matrices with shock, BLOCK

participants_to_find <- c(101:154)

# data_directory <- "~/research_data/gaborgen/raw_data"
data_directory <- "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data"

log_file_paths <- list.files(
  path = data_directory,
  pattern = "logfile.dat$",
  recursive = T,
  full.names = T
)

# where_to_save <- "/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder"
where_to_save <- "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/shock_design_matrix_column_to_add/"

# write shock times

for (i in 1:length(log_file_paths)) {
  log_file <- read.delim(log_file_paths[i], header = T, sep = ",")

  name_to_save <- sub(
    basename(log_file_paths[i]),
    pattern = "logfile.dat$",
    replacement = "shock_onsets.1D"
  )

  # approximate shock onset,
  shock_onsets <- subset(
    log_file,
    paired == 1,
    select = timeSinceFirstTR,
    drop = T
  ) +
    2 # shock comes two seconds after the paired stimulus (approximately, double check needed)

  write(
    x = shock_onsets,
    file = paste0(where_to_save, name_to_save),
    ncolumns = 1
  )
}


#BLOCK(d, p) ; d = duration seconds, p = amplitude, 1 is easiest to understand and should be default unless you know what you're doing

shock_onset_stim_files <- list.files(
  path = where_to_save,
  pattern = "shock_onsets.1D$",
  recursive = T,
  full.names = T
)

setwd(where_to_save)


#enters(\n) will cause command to enter early

for (i in 1:length(shock_onset_stim_files)) {
  sub_name <- sub(
    basename(shock_onset_stim_files[i]),
    pattern = "shock_onsets.1D$",
    replacement = ""
  )

  system2(
    'tcsh',
    args = c(
      '-c',
      shQuote(paste(
        '3dDeconvolve -nodata 1070 2.0 -num_stimts 1 -stim_times 1',
        paste0(shock_onset_stim_files[i], " "),
        '"BLOCK(.1,1)" -x1D',
        paste0(sub_name, "design_matrix.xmat.1D "),
        '-xjpeg',
        paste0(sub_name, "design_matrix.jpg "),
        '-x1D_stop'
      ))
    )
  )
}

# 3dDeconvolve -nodata 1070 2.0 \
# -num_stimts 1 \
# -stim_times 1 /Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder/gaborgen24_fMRI_Day1_121_shock_onsets.1D "BLOCK(.1,1)" \
# -x1D design_matrix.xmat.1D \
# -xjpeg design_matrix.jpg \
# -x1D_stop
#
#
#
#
#
# 3dDeconvolve -nodata 200 2.0 \
# -num_stimts 1 \
# -stim_times 1 stimfile.1D 'BLOCK(2,1)' \
# -x1D design_matrix.xmat.1D \
# -xjpeg design_matrix.jpg \
# -x1D_stop
