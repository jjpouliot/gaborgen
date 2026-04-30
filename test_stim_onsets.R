log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_145/DAT/gaborgen24_fMRI_Day1_145_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_145",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_145",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_158/DAT/gaborgen24_fMRI_Day1_158_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_158",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_158",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_137/DAT/gaborgen24_fMRI_Day1_137_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_137",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_137",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)


log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_138/DAT/gaborgen24_fMRI_Day1_138_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_138",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_138",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_144/DAT/gaborgen24_fMRI_Day1_144_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_144",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_144",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)


log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_101/DAT/gaborgen24_fMRI_Day1_101_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_101",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_101",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_102/DAT/gaborgen24_fMRI_Day1_102_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_102",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_102",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_103/DAT/gaborgen24_fMRI_Day1_103_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_103",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_103",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_106/DAT/gaborgen24_fMRI_Day1_106_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_106",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_106",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)
log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_107/DAT/gaborgen24_fMRI_Day1_107_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_107",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_107",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_108/DAT/gaborgen24_fMRI_Day1_108_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_108",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_108",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_109/DAT/gaborgen24_fMRI_Day1_109_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_109",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_109",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_111/DAT/gaborgen24_fMRI_Day1_111_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_111",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_111",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_151/DAT/gaborgen24_fMRI_Day1_151_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_151",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY1_151",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

# day 2

log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_137/DAT/gaborgen24_fMRI_Day2_137_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_137",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_137",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)


#####
log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_138/DAT/gaborgen24_fMRI_Day2_138_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_138",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_138",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

####
log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_144/DAT/gaborgen24_fMRI_Day2_144_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_144",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_144",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)

####
log_file <- read.delim(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_145/DAT/gaborgen24_fMRI_Day2_145_logfile.dat",
  header = T,
  sep = ","
)


stim_onset <- log_file$timeSinceFirstTR


write(
  x = stim_onset,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_145",
    "/stim_times.1D"
  ),
  ncolumns = 1
)


stim_onset_t1off <- log_file$timeSinceFirstTR + 2


write(
  x = stim_onset_t1off,
  file = paste0(
    "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data/GABORGEN24_DAY2_145",
    "/stim_times_t1off.1D"
  ),
  ncolumns = 1
)
