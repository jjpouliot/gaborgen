motion_parent_directory <- "/home/andrewfarkas/Downloads/motion/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results"

participant_directories <- list.dirs(motion_parent_directory)[grepl(
  pattern = "GABOR",
  x = list.dirs(motion_parent_directory)
)]

for (i in 1:length(participant_directories)) {
  current_directory <- participant_directories[i]
  current_motion_file <- paste0(current_directory, "/mot_demean.r01.1D")
  file.copy(
    from = current_motion_file,
    to = paste0(
      "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/default_afni_movement_censor/",
      basename(current_directory),
      ".mot_demean.r01.1D"
    )
  )
}
