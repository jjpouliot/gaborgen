# used sapply so that dat files are in the same order as participant directories
found_dat_log_filepaths <- unlist(lapply(participant_directories, function(x) {
  list.files(path = paste0(x, '/DAT'), pattern = 'logfile.dat$', full.names = T)
}))


if (length(found_dat_log_filepaths) != length(participant_directories)) {
  stop(
    "There is not a dat log file for each participant directory. Directories are not in the proper structure or dat file is missing."
  )
}

T1_off_participant_regex <- "101|103|106|110|111|114|117|119|122|DAY1_123|DAY1_127|DAY1_130|DAY1_131|DAY1_132|DAY1_134|DAY1_135|DAY1_137|DAY1_144|DAY1_146|DAY1_148|DAY1_149|DAY1_150|DAY1_154|DAY1_158|DAY1_159|DAY1_160|DAY2_133"

for (directory_index in 1:length(participant_directories)) {
  log_file <- read.delim(
    found_dat_log_filepaths[directory_index],
    header = T,
    sep = ","
  )

  stim_onset <- log_file$timeSinceFirstTR

  # participants that had a T1 off trigger start and need to have stim time adjusted
  if (
    grepl(T1_off_participant_regex, participant_directories[directory_index])
  ) {
    stim_onset <- log_file$timeSinceFirstTR + 2
  }

  write(
    x = stim_onset,
    file = paste0(participant_directories[directory_index], "/stim_times.1D"),
    ncolumns = 1
  )

  phase_ids <- unique(log_file$phase)

  stim_ids <- unique(log_file$stim)

  for (phase_index in 1:length(phase_ids)) {
    for (stim_index in 1:length(stim_ids)) {
      current_phase_stim_times <- subset(
        log_file,
        phase == phase_ids[phase_index] &
          stim == stim_ids[stim_index],
        select = timeSinceFirstTR,
        drop = T
      )

      # participants that had a T1 off trigger start and need to have stim time adjusted
      if (
        grepl(
          T1_off_participant_regex,
          participant_directories[directory_index]
        )
      ) {
        current_phase_stim_times <- current_phase_stim_times + 2
      }

      write(
        x = current_phase_stim_times,
        ncolumns = 1,
        file = paste0(
          participant_directories[directory_index],
          paste0(
            "/phase_",
            phase_ids[phase_index],
            "_stim_",
            stim_ids[stim_index],
            "_stim_times.1D"
          )
        )
      )
    }
  }
}


# by_block can be created to get sitmulus onsets by block instead of phase
if (exists("by_block")) {
  for (directory_index in 1:length(participant_directories)) {
    log_file <- read.delim(
      found_dat_log_filepaths[directory_index],
      header = T,
      sep = ","
    )

    stim_onset <- log_file$timeSinceFirstTR

    # participants that had a T1 off trigger start and need to have stim time adjusted
    if (
      grepl(T1_off_participant_regex, participant_directories[directory_index])
    ) {
      stim_onset <- log_file$timeSinceFirstTR + 2
    }

    write(
      x = stim_onset,
      file = paste0(participant_directories[directory_index], "/stim_times.1D"),
      ncolumns = 1
    )

    # no for loops because there isn't a proper block (1 to 4) in log file
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 1 & log_file$stim == 1],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/habituation_CS+_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 1 & log_file$stim == 2],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/habituation_GS1_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 1 & log_file$stim == 3],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/habituation_GS2_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 1 & log_file$stim == 4],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/habituation_GS3_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 1 & log_file$stim == 1
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_1_CS+_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 1 & log_file$stim == 2
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_1_GS1_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 1 & log_file$stim == 3
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_1_GS2_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 1 & log_file$stim == 4
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_1_GS3_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 2 & log_file$stim == 1
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_2_CS+_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 2 & log_file$stim == 2
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_2_GS1_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 2 & log_file$stim == 3
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_2_GS2_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[
        log_file$phase == 2 & log_file$block == 2 & log_file$stim == 4
      ],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/acquisition_block_2_GS3_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 3 & log_file$stim == 1],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/extinction_CS+_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 3 & log_file$stim == 2],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/extinction_GS1_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 3 & log_file$stim == 3],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/extinction_GS2_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
    write(
      x = log_file$timeSinceFirstTR[log_file$phase == 3 & log_file$stim == 4],
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        "/extinction_GS3_stim_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
      )
    )
  }
}


# write out shock times if shock_times exists
if (exists("shock_times")) {
  for (directory_index in 1:length(participant_directories)) {
    log_file <- read.delim(
      found_dat_log_filepaths[directory_index],
      header = T,
      sep = ","
    )

    shock_times_since_first_TR <- (log_file$timeSinceFirstTR[
      log_file$paired == 1
    ]) +
      2

    # participants that had a T1 off trigger start and need to have stim time adjusted
    if (
      grepl(T1_off_participant_regex, participant_directories[directory_index])
    ) {
      shock_times_since_first_TR <- shock_times_since_first_TR + 2
    }

    write(
      x = shock_times_since_first_TR,
      ncolumns = 1,
      file = paste0(
        participant_directories[directory_index],
        paste0(
          "/shock_times_for_stan.1D" # added "for_stan" so it doesn't get used by the afni_proc_py regression
        )
      )
    )
  }
}
