library(tidyverse)
library(cmdstanr)
# library(signal)

# we recommend running this in a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# Filter used, same notation as matlab
# fifth order (n = 5)
# cut off frequency = 0.026 * .25 (nyquist) = .0065 Hz
highpass <- signal::butter(n = 5, W = 0.026, type = "high")


# load data from text files ####

data_dir <- c(
  # "/home/andrew/Downloads/roi_data_and_info"
  # "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info"
  # "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info"
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI"
)

roi_key <- read_delim(paste0(
  data_dir,
  '/roi_data_and_info/HCPex_SUIT_labels.txt'
))

colnames(roi_key) <- c("roi_id", "roi")

# useable_participants <- c("145")
# used for a lot of the early testing
# useable_participants <- c(
#   "116",
#   "117",
#   "119",
#   #"120", no ROI?
#   "121",
#   "122",
#   "125",
#   "126",
#   "127",
#   "128",
#   "129",
#   "131",
#   "132",
#   "133",
#   "145",
#   "149"
# )

# useable_participants <- c(
#   "101",
#   #"102", # 22% censored
#   "103",
#   "106",
#   "107",
#   "108",
#   "109",
#   "113",
#   #"114", alot censored
#   "115",
#   "116",
#   "117",
#   "119",
#   #"120", no ROI? # lots censored
#   "121",
#   "122",
#   "123",
#   #"124", # lots censored
#   "125",
#   "126",
#   "127",
#   "128",
#   "129",
#   #"131", #missing qc
#   "132",
#   #"133", #missing qc
#   #"134",#missing qc
#   #"135",#missing qc
#   #"136",#lots censored
#   #"137", make roi
#   #"138", make roi
#   #"139", #lots censored
#   #"140", #lots censored
#   "141",
#   #"142", #probably asleep
#   #"143", #probably asleep
#   #"144", #lots censored
#   "145",
#   "149"
# )

#124 has too much censored, 157 doesn't have blip scans, could add them back in
useable_participants <- c(
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
  "158" # made
)

# participants used for HBM tech report
# useable_participants <- c(
#   "101",
#   "102", # 22% censored
#   "103",
#   "106",
#   "107",
#   "108",
#   "109",
#   "113",
#   #"114", alot censored
#   "115",
#   "116",
#   "117",
#   "119",
#   #"120", no ROI? # lots censored
#   "121",
#   "122",
#   "123",
#   #"124", # lots censored
#   "125",
#   "126",
#   "127",
#   "128",
#   "129",
#   #"131", excess movement
#   "132",
#   "133",
#   #"134",#asleep?
#   "135",
#   #"136",#lots censored
#   #"137", make roi
#   #"138", make roi
#   #"139", #lots censored
#   #"140", #lots censored
#   "141"
#   #"142", #probably asleep
#   #"143", #probably asleep
#   #"144", #lots censored
#   # "145",
#   # "149"
# )

add_shock <- T

bold_per_roi_df <- data.frame(
  "par" = numeric(),
  "roi" = character(),
  "roi_id" = numeric(),
  "seconds" = numeric(),
  "uncensored" = numeric(),
  "bold" = numeric()
)


for (i in 1:length(useable_participants)) {
  bold_text_file <- list.files(
    path = paste0(data_dir, "/roi_data_and_info/roi_timeseries"),
    pattern = paste0(useable_participants[i], "_roi_stats.txt"),
    recursive = T,
    full.names = T,
    include.dirs = F
  )

  all_bold <- read_delim(
    bold_text_file,
    trim_ws = T
  )

  all_bold <- all_bold %>%
    select(starts_with("NZMEAN")) %>%
    mutate(across(everything(), ~ . - 100)) %>%
    mutate(across(everything(), ~ signal::filtfilt(highpass, .)))

  # why are there extra areas that are not in the key?
  long_mean_bold <- all_bold %>%
    mutate(time = seq(0, 1069 * 2, by = 2)) %>%
    pivot_longer(cols = contains("NZMEAN_")) %>%
    select(name, time, value) %>%
    mutate(
      roi_id = sub(pattern = "NZMean_", replacement = "", x = name),
      par = useable_participants[i]
    ) %>%
    rename("bold" = value) %>%
    merge(x = ., y = roi_key, by.x = "roi_id", by.y = "roi_id")

  if (useable_participants[i] < 123) {
    censor_info <- read_delim(
      paste0(
        paste0(data_dir, "/roi_data_and_info/default_afni_movement_censor"),
        '/censor_GABORGEN24_',
        useable_participants[i],
        '_combined_2.1D'
      ),
      col_names = F
    )[[2]]
  } else {
    censor_info <- read_delim(
      paste0(
        paste0(data_dir, "/roi_data_and_info/default_afni_movement_censor"),
        '/censor_GABORGEN24_DAY1_',
        useable_participants[i],
        '_combined_2.1D'
      ),
      col_names = F
    )[[2]]
  }

  current_bold_per_roi_df <- data.frame(
    "par" = long_mean_bold$par,
    "roi" = long_mean_bold$roi,
    "roi_id" = as.numeric(long_mean_bold$roi_id),
    "time_sec" = long_mean_bold$time,
    "bold" = long_mean_bold$bold
  ) %>%
    arrange(roi_id, time_sec) %>%
    mutate("censor" = rep(censor_info, n() / 1070), .before = "bold") # so everything lines up right

  bold_per_roi_df <- rbind.data.frame(
    bold_per_roi_df,
    current_bold_per_roi_df
  )
}

# by participant design matrix
all_par_design_matrices <- tibble()

for (i in 1:length(useable_participants)) {
  if (useable_participants[i] < 123) {
    old_design_matrix <- read_delim(
      paste0(
        paste0(data_dir, "/roi_data_and_info/design_matrices"),
        '/GABORGEN24_',
        useable_participants[i],
        '.results_X.nocensor.xmat.1D'
      ),
      skip = 64,
      col_names = F
    )
    # just keep the last few columns that are for motion
    motion_design_matrix <- old_design_matrix[
      1:(nrow(old_design_matrix) - 1),
      (ncol(old_design_matrix) - 5):(ncol(old_design_matrix))
    ]

    # Add in per stimulus variables
    stimulus_design_matrix <- read_delim(
      paste0(
        paste0(data_dir, "/raw_data/"),
        '/GABORGEN24_',
        useable_participants[i],
        '/X_IM.nocensor.xmat.1D'
      ),
      skip = 84,
      col_names = F
    )
    stimulus_design_matrix <- stimulus_design_matrix[
      1:(nrow(stimulus_design_matrix) - 1),
      4:ncol(stimulus_design_matrix)
    ]

    current_design_matrix <- cbind(
      stimulus_design_matrix,
      motion_design_matrix
    )
  } else {
    old_design_matrix <- read_delim(
      paste0(
        paste0(data_dir, "/roi_data_and_info/design_matrices"),
        '/GABORGEN24_DAY1_',
        useable_participants[i],
        '.results_X.nocensor.xmat.1D'
      ),
      skip = 64,
      col_names = F
    )
    # just keep the last few columns that are for motion
    motion_design_matrix <- old_design_matrix[
      1:(nrow(old_design_matrix) - 1),
      (ncol(old_design_matrix) - 5):(ncol(old_design_matrix))
    ]

    # Add in per stimulus variables
    stimulus_design_matrix <- read_delim(
      paste0(
        paste0(data_dir, "/raw_data/"),
        '/GABORGEN24_DAY1_',
        useable_participants[i],
        '/X_IM.nocensor.xmat.1D'
      ),
      skip = 84,
      col_names = F
    )
    stimulus_design_matrix <- stimulus_design_matrix[
      1:(nrow(stimulus_design_matrix) - 1),
      4:ncol(stimulus_design_matrix)
    ]

    current_design_matrix <- cbind(
      stimulus_design_matrix,
      motion_design_matrix
    )
  }

  colnames(current_design_matrix) <- c(
    paste0("habituation_CS+_", c(1:8)),
    paste0("habituation_GS1_", c(1:8)),
    paste0("habituation_GS2_", c(1:8)),
    paste0("habituation_GS3_", c(1:8)),
    paste0("acquisition_block1_CS+_", c(1:12)),
    paste0("acquisition_block1_GS1_", c(1:12)),
    paste0("acquisition_block1_GS2_", c(1:12)),
    paste0("acquisition_block1_GS3_", c(1:12)),
    paste0("acquisition_block2_CS+_", c(1:12)),
    paste0("acquisition_block2_GS1_", c(1:12)),
    paste0("acquisition_block2_GS2_", c(1:12)),
    paste0("acquisition_block2_GS3_", c(1:12)),
    paste0("extinction_CS+_", c(1:12)),
    paste0("extinction_GS1_", c(1:12)),
    paste0("extinction_GS2_", c(1:12)),
    paste0("extinction_GS3_", c(1:12)),
    paste0("shock_", c(1:15)),
    "mot_demean_r01_0",
    "mot_demean_r01_1",
    "mot_demean_r01_2",
    "mot_demean_r01_3",
    "mot_demean_r01_4",
    "mot_demean_r01_5"
  )

  current_design_matrix <- current_design_matrix %>%
    mutate(time = seq(0, 1069 * 2, by = 2), .before = 1) %>%
    mutate(participant = useable_participants[i], .before = 1)

  all_par_design_matrices <- rbind(
    all_par_design_matrices,
    current_design_matrix
  )
}

#get which csp trials were paired with the US per participant in order.
#this is used in model 35 to not select those CS+ column from the design matrix
#Assuming US is the only driver of activity on paired trials by only using the HDF of the shock column
paired_csp_indices <- tibble(
  "participant" = character(length = 41 * 44),
  "trial" = integer(length = 41 * 44),
  "paired" = logical(length = 41 * 44)
)
paired_csp_matrix <- matrix(nrow = 41, ncol = 44)

all_log_files <- list.files(
  paste0(data_dir, "/raw_data/"),
  pattern = "logfile.dat$",
  recursive = T,
  full.names = T
)

useable_log_files <- all_log_files[!grepl("DAY2", all_log_files)]

useable_log_files <- useable_log_files[grepl(
  paste(useable_participants, collapse = "|"),
  useable_log_files
)]


row_index <- 1
for (i in 1:length(useable_participants)) {
  current_par <- useable_participants[i]

  current_log_file <- read.delim(useable_log_files[i], sep = ",")

  current_csp_trial <- 1
  for (j in 1:176) {
    if (current_log_file$stim[j] == 1 & current_log_file$paired[j] == 0) {
      paired_csp_indices[row_index, "participant"] <- useable_participants[i]
      paired_csp_indices[row_index, "trial"] <- current_csp_trial
      paired_csp_indices[row_index, "paired"] <- F
      paired_csp_matrix[i, current_csp_trial] <- 0
      row_index <- row_index + 1
      current_csp_trial <- current_csp_trial + 1
    } else if (
      current_log_file$stim[j] == 1 & current_log_file$paired[j] == 1
    ) {
      paired_csp_indices[csp_index, "participant"] <- useable_participants[i]
      paired_csp_indices[csp_index, "trial"] <- current_csp_trial
      paired_csp_indices[csp_index, "paired"] <- T
      paired_csp_matrix[i, current_csp_trial] <- 1
      row_index <- row_index + 1
      current_csp_trial <- current_csp_trial + 1
    }
  }
}

paired_csp_indices
paired_csp_matrix

rowMeans(paired_csp_matrix)

# create stan list ####
used_df <- bold_per_roi_df %>%
  # filter(roi_id %in% c(1, 2)) #
  filter(
    roi_id %in%
      c(
        1, # V1 L
        181, # V1 R
        4, # V4 L
        184, # V4 R
        22, # V5/MT L
        23, # V5/MT L
        202, # V5/MT R
        203, # V5/MT R
        8, # V6 L
        188, # V6 R
        89, # TE L
        90, # TE L
        91, # TE L
        269, # TE R
        270, # TE R
        271, # TE R
        98, # TPJ L
        99, # TPJ L
        278, # TPJ R
        279, # TPJ R
        144, # ACC L
        324, # ACC R
        387, # Amygdala L
        420, # Amygdala R
        157, # Orbital L
        337, # Orbital R
        68, # Anterior insula L
        69, # Anterior insula L
        73, # Anterior insula L
        248, # Anterior insula R
        249, # Anterior insula R
        253, # Anterior insula R
        384, # Nucleus Accubems L
        417, # Nucleus Accubems R
        80, # Hippocampus L
        260 # Hippocampus R
      )
  ) #
# filter(roi_id %in% c(1:3, 181:183)) # visual

used_df <- used_df %>%
  mutate(
    roi_merge = case_when(
      roi_id %in% c(22, 23) ~ "V5_MT_L",
      roi_id %in% c(202, 203) ~ "V5_MT_R",
      roi_id %in% c(89, 90, 91) ~ "TE_L",
      roi_id %in% c(269, 270, 271) ~ "TE_R",
      roi_id %in% c(98, 99) ~ "TPJ_L",
      roi_id %in% c(278, 279) ~ "TPJ_R",
      roi_id %in% c(68, 69, 73) ~ "Ant_Ins_L",
      roi_id %in% c(248, 249, 253) ~ "Ant_Ins_R",
      roi_id %in% c(144) ~ "ACC_L",
      roi_id %in% c(324) ~ "ACC_R",
      .default = roi
    )
  )

used_df$roi_merge %>% unique()


# Currently the merge id doesn't really get used in stan model, it just loops 1 to end.
# It would be simipler to sort the merge id so that the numbers within participant come in order
used_df <- used_df %>%
  mutate(
    roi_merge_id = as.integer(factor(
      roi_merge,
      levels = c(
        "Primary_Visual_Cortex_L",
        "Primary_Visual_Cortex_R",
        "Fourth_Visual_Area_L",
        "Fourth_Visual_Area_R",
        "V5_MT_L",
        "V5_MT_R",
        "Sixth_Visual_Area_L",
        "Sixth_Visual_Area_R",
        "TE_L",
        "TE_R",
        "TPJ_L",
        "TPJ_R",
        "Hippocampus_L",
        "Hippocampus_R",
        "Amygdala_L",
        "Amygdala_R",
        "Nucleus_Accumbens_L",
        "Nucleus_Accumbens_R",
        "Ant_Ins_L",
        "Ant_Ins_R",
        "ACC_L",
        "ACC_R",
        "Orbital_Frontal_Complex_L",
        "Orbital_Frontal_Complex_R"
      )
    ))
  ) %>%
  group_by(par, time_sec, censor, roi_merge, roi_merge_id) %>%
  reframe(bold = mean(bold))

# used_df <- used_df %>%
#   arrange(par, roi_merge_id)

# used_df %>%
#   mutate(par_int = as.integer(as.factor(par))) %>%
#   filter(
#     roi_merge %in% "Primary_Visual_Cortex_R",
#     par == "108"
#   ) %>%
#   pull(bold) %>%
#   plot(type = "l")

# used_df %>%
#   mutate(par_int = as.integer(as.factor(par))) %>%
#   filter(roi_merge %in% "Primary_Visual_Cortex_R") %>%
#   print(n = 6000)

# used_df$roi_merge %>% unique()
# used_df$par %>% unique()

# extra filtering for quality checking
# used_df <- used_df %>%
#   filter(
#     roi_merge %in%
#       c(
#         "Ant_Ins_L",
#         "Ant_Ins_R"
#       )
#     )
#
#
# used_df <- used_df %>%
#   filter(
#     roi_merge %in%
#       c(
#         "Primary_Visual_Cortex_L",
#         "Primary_Visual_Cortex_R"
#       )
#   )
#
used_df <- used_df %>%
  filter(
    roi_merge %in%
      c(
        # "Primary_Visual_Cortex_L",
        # "Primary_Visual_Cortex_R" ,
        # "Fourth_Visual_Area_L",
        # "Fourth_Visual_Area_R",
        "V5_MT_L",
        "V5_MT_R" #,
        # "Sixth_Visual_Area_L",
        # "Sixth_Visual_Area_R",
        # "TE_L",
        # "TE_R",
        # "TPJ_L",
        # "TPJ_R",
        # "Hippocampus_L",
        # "Hippocampus_R",
        # "Amygdala_L",
        # "Amygdala_R",
        # "Nucleus_Accumbens_L",
        # "Nucleus_Accumbens_R",
        # "Ant_Ins_L",
        # "Ant_Ins_R",
        # "ACC_L",
        # "ACC_R",
        # "Orbital_Frontal_Complex_L",
        # "Orbital_Frontal_Complex_R"
      )
  )


fmri_stan_list <- list()

fmri_stan_list$n_par <- used_df$par %>% unique() %>% length()

# fmri_stan_list$n_roi <- used_df$roi_id %>% unique() %>% length()
fmri_stan_list$n_roi <- used_df$roi_merge_id %>% unique() %>% length()

fmri_stan_list$n_bold <- 1070

fmri_stan_list$n <- fmri_stan_list$n_par *
  fmri_stan_list$n_roi *
  fmri_stan_list$n_bold

fmri_stan_list$par <- used_df$par %>% as.factor() %>% as.integer()

# fmri_stan_list$roi <- used_df$roi_id %>% as.factor() %>% as.integer()

fmri_stan_list$roi <- used_df$roi_merge_id

fmri_stan_list$paired_mat_one_true <- paired_csp_matrix

# just vector
# fmri_stan_list$bold <- used_df$bold
# [par, time]
# fmri_stan_list$bold <- used_df$bold %>%
#   array(
#     dim = c(
#       fmri_stan_list$n_bold * fmri_stan_list$n_roi,
#       # fmri_stan_list$n_roi,
#       fmri_stan_list$n_par
#     )
#   ) %>%
#   aperm(c(2, 1))
# [par, roi, time]
fmri_stan_list$bold <- used_df$bold %>%
  array(
    dim = c(
      1070,
      fmri_stan_list$n_roi,
      fmri_stan_list$n_par
    )
  ) %>%
  aperm(c(3, 2, 1))

# just vector
# fmri_stan_list$usable_bold_indices <- used_df$censor %>% as.integer()
# [par, roi, time] to [par, time]; same for each roi
fmri_stan_list$usable_bold_indices_one_is_true <- used_df$censor %>%
  array(
    dim = c(
      1070,
      fmri_stan_list$n_roi,
      fmri_stan_list$n_par
    )
  ) %>%
  aperm(c(3, 2, 1))

# We first subset with drop = FALSE to preserve all dimensions:
temp <- fmri_stan_list$usable_bold_indices_one_is_true[, 1, , drop = FALSE]
# temp now has dimensions: [n_par, 1, 1070]

# Then we manually reshape the result by dropping the singleton second dimension:
dim(temp) <- c(dim(temp)[1], dim(temp)[3])
# Now temp is a matrix of dimensions: [n_par, 1070]

# Finally, assign it back:
fmri_stan_list$usable_bold_indices_one_is_true <- temp

# fmri_stan_list$usable_bold_indices <- fmri_stan_list$usable_bold_indices[, 1, ] %>% array(dim =c(fmri_stan_list$n_par, 1070))

# array [n_par] int
censor_vec <- vector()
for (p in 1:fmri_stan_list$n_par) {
  censor_vec <- c(
    censor_vec,
    (1070 - sum(fmri_stan_list$usable_bold_indices_one_is_true[p, ] == 1))
  )
}
fmri_stan_list$n_censor <- censor_vec %>% as.integer()

# try without motion betas
# all_par_design_matrices <- all_par_design_matrices %>%
#   select(!starts_with("mot"))

fmri_stan_list$n_beta <- ncol(all_par_design_matrices) - 2

# array[n_par, 1070, n_beta] real design_array
design_list <- split(
  all_par_design_matrices,
  all_par_design_matrices$participant
)
design_matrices <- lapply(
  design_list,
  function(df) as.matrix(df[, setdiff(names(df), c("participant", "time"))])
)
tmp_array <- array(
  unlist(design_matrices),
  dim = c(
    1070,
    fmri_stan_list$n_beta,
    fmri_stan_list$n_par
  )
)
design_array <- aperm(tmp_array, c(3, 1, 2))

fmri_stan_list$design_array <- design_array

rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

csp_linespace_per_par <- tibble()

for (p in 1:fmri_stan_list$n_par) {
  current_csp_linespace_row <- c(1:44)[
    fmri_stan_list$paired_mat_one_true[p, ] == 0
  ]

  csp_linespace_per_par <- rbind(
    csp_linespace_per_par,
    current_csp_linespace_row
  )
}

fmri_stan_list$csp_linespace_per_par <- as.matrix(csp_linespace_per_par)


fmri_stan_list$useable_cue_betas_per_par


# fmri_stan_list$n_beta_stim <- 13

# cmdstanr::write_stan_json(
#   fmri_stan_list,
#   file = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_no_shock.json"
# )
# cmdstanr::write_stan_json(
#   fmri_stan_list,
#   file = "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_visual.json"
# )

# model 18 doesn't need the roi dimension... only works for 1 roi
# fmri_stan_list$bold <- fmri_stan_list$bold[, 1, ]

# save(
#   fmri_stan_list,
#   file = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_antinsula_24_single_trialDM.RData"
# )

# cmdstanr::write_stan_json(
#   fmri_stan_list,
#   file = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_antinsula_24_single_trialDM.json"
# )

# save(
#   fmri_stan_list,
#   file = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_single_trialDM_ant_insulas.RData"
# )

# cmdstanr::write_stan_json(
#   fmri_stan_list,
#   file = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_single_trialDM_ant_insulas.json"
# )

save(
  fmri_stan_list,
  file = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_single_trialDM_V5_2.RData"
)

cmdstanr::write_stan_json(
  fmri_stan_list,
  file = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_single_trialDM_V5_2.json"
)

# cmdstanr::write_stan_json(
#   fmri_stan_list,
#   file = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_no_mot_betas.json"
# )

# load(
#   "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_antinsula_24_single_trialDM.RData"
# )

load(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_single_trialDM.RData"
)

# Stan settings ####
number_of_chains <- 1
warmup_samples_per_chain <- 2000
posterior_samples_per_chain <- 2000
# where_to_save_chains <- '/home/andrew/Documents/stan_chains_ssd/'
# where_to_save_chains <- '/run/media/andrew/Barracuda_8tb/stan_chains/'
# where_to_save_chains <- '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains'
where_to_save_chains <- '/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains'

# First map_rect model ####

# model_path <- '/home/andrew/Documents/GitHub/gaborgen/stan_models/fMRI/Model012.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model019.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model021.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model018.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model018.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model018_no_ML.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model025.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model026.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model027.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model028.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model030.stan' # used this model, next few aren't necessary
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model033.stan'
# model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model034.stan'
model_path <- '/home/andrewfarkas/Repositories/gaborgen/stan_models/fMRI/Model035.stan' #editted model 30

# Fit models
model <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE)
  # cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)


model_opt <- model$optimize(
  data = fmri_stan_list,
  iter = 1e6,
  jacobian = T,
  init = .1,
  threads = 10
)

model_opt$metadata()

model_opt_summary <- model_opt$summary()

model$variables()$parameters

model_opt_2 <- model$optimize(
  data = fmri_stan_list,
  iter = 1e6,
  jacobian = T,
  init = .1,
  threads = 10
)

model_opt_2$metadata()

model_opt_summary_2 <- model_opt_2$summary()

model$variables()$parameters


model_fit <- model$sample(
  data = fmri_stan_list,
  refresh = 1,
  seed = 3,
  threads_per_chain = 16,
  # adapt_delta = .7, # not default
  init = .1,
  # metric = "dense_e", # not default
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  # opencl_ids = c(0, 0),
  output_dir = where_to_save_chains,
  chains = number_of_chains,
  parallel_chains = number_of_chains
)


model012_fit <- cmdstanr::as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/stan_chains/Model012-202504092356-1-23e85e.csv",
    "/home/andrewf/Downloads/stan_chains/Model012-202504092356-2-23e85e.csv",
    "/home/andrewf/Downloads/stan_chains/Model012-202504092356-3-23e85e.csv",
    "/home/andrewf/Downloads/stan_chains/Model012-202504092356-4-23e85e.csv"
  )
)

model012_draws <- model012_fit$draws(format = "df")


df_long <- model012_draws %>%
  select(starts_with("beta")) %>%
  pivot_longer(
    cols = starts_with("betas["),
    names_to = c("par", "roi", "coef"),
    names_pattern = "betas\\[(\\d+),(\\d+),(\\d+)\\]",
    values_to = "beta"
  ) %>%
  mutate(coef = factor(coef, levels = c(1:18)))

df_long %>%
  filter(
    par == 1,
    # coef %in% c(9:12)
    # coef %in% c(1:4)
    coef %in% c(5:8)
  ) %>%
  ggplot() +
  geom_density(aes(x = beta, color = coef))

df_long %>%
  filter(
    par == 2,
    coef %in% c(9:12)
    # coef %in% c(1:4)
    # coef %in% c(5:8)
  ) %>%
  ggplot() +
  geom_density(aes(x = beta, color = coef))

#

#

#

split()

split(
  all_par_design_matrices,
  all_par_design_matrices$participant,
  drop = c("participant")
)

fmri_stan_list$design_matrices[1, , ] <- all_par_design_matrices %>%
  filter(participant == 149) %>%
  select(3:ncol(all_par_design_matrices)) %>%
  as.matrix()

#no polort high-pass
fmri_stan_list$design_matrix <- design_matrix[, 17:ncol(design_matrix)]
# fmri_stan_list$design_matrix <- design_matrix

fmri_stan_list$n_DM_cols <- ncol(fmri_stan_list$design_matrix)


roi_df <- read_delim(
  '/home/andrew/Downloads/stan_data/145_first_attempt/roi_stats.txt',
  trim_ws = T
)

design_matrix <- design_matrix[
  1:(nrow(design_matrix) - 1),
  2:ncol(design_matrix)
]

design_matrix$X2 <- as.numeric(design_matrix$X2)

censor_info <- read_delim(
  '/home/andrew/Downloads/stan_data/145_first_attempt/censor_GABORGEN24_DAY1_145_combined_2.1D',
  col_names = F
)[[2]]

fmri_stan_list <- list()

fmri_stan_list$n_participants <- 1
fmri_stan_list$n_ROIs <- 1

# 69 = Anterior_Ventral_Insular_Area_L
fmri_stan_list$amplitude_no_censor <- roi_df %>%
  select(paste0("Mean_", 69)) %>%
  pull()

fmri_stan_list$censor <- c(1:1070)[censor_info == 1]

fmri_stan_list$n_amplitude <- fmri_stan_list$amplitude_no_censor %>%
  length()

fmri_stan_list$n_censor <- 1070 - length(fmri_stan_list$censor)

#no polort high-pass
fmri_stan_list$design_matrix <- design_matrix[, 17:ncol(design_matrix)]
# fmri_stan_list$design_matrix <- design_matrix

fmri_stan_list$n_DM_cols <- ncol(fmri_stan_list$design_matrix)


highpass <- butter(n = 5, W = 0.026, type = "high")


# 145
roi_df <- read_delim(
  '/home/andrew/Downloads/stan_data/145_first_attempt/roi_stats.txt',
  trim_ws = T
)
roi_key <- read_delim(
  '/home/andrew/Downloads/stan_data/145_first_attempt/HCPex_SUIT_labels.txt'
)

colnames(roi_key) <- c("roi_id", "roi")

design_matrix <- read_delim(
  '/home/andrew/Downloads/stan_data/145_first_attempt/X.nocensor.xmat.1D',
  skip = 64,
  col_names = F
)
design_matrix <- design_matrix[
  1:(nrow(design_matrix) - 1),
  2:ncol(design_matrix)
]

design_matrix$X2 <- as.numeric(design_matrix$X2)

censor_info <- read_delim(
  '/home/andrew/Downloads/stan_data/145_first_attempt/censor_GABORGEN24_DAY1_145_combined_2.1D',
  col_names = F
)[[2]]

fmri_stan_list <- list()

fmri_stan_list$n_participants <- 1
fmri_stan_list$n_ROIs <- 1

# 69 = Anterior_Ventral_Insular_Area_L
fmri_stan_list$amplitude_no_censor <- roi_df %>%
  select(paste0("Mean_", 69)) %>%
  pull()

fmri_stan_list$censor <- c(1:1070)[censor_info == 1]

fmri_stan_list$n_amplitude <- fmri_stan_list$amplitude_no_censor %>%
  length()

fmri_stan_list$n_censor <- 1070 - length(fmri_stan_list$censor)

# fmri_stan_list$design_matrix <- design_matrix

#no polort high-pass
fmri_stan_list$design_matrix <- design_matrix[, 17:ncol(design_matrix)]

fmri_stan_list$n_DM_cols <- ncol(fmri_stan_list$design_matrix)


# add more to list

participant_dir <- c("/home/andrew/Downloads/stan_data/")

participant_ids <- c("145", "149")

bold_per_roi_df <- data.frame(
  "par" = numeric(),
  "roi" = character(),
  "roi_id" = numeric(),
  "seconds" = numeric(),
  "censor" = numeric(),
  "bold" = numeric()
)


for (i in 1:length(participant_ids)) {
  participant_files <- list.files(
    path = participant_dir,
    pattern = participant_ids[i],
    recursive = T,
    full.names = T,
    include.dirs = F
  )

  bold_text_file <- participant_files[grepl(
    x = participant_files,
    pattern = "roi"
  )]

  all_bold <- read_delim(
    bold_text_file,
    trim_ws = T
  )

  censor_info <-
    # why are there extra areas that are not in the key?
    long_mean_bold <- all_bold %>%
      mutate(time = seq(0, 1069 * 2, by = 2)) %>%
      pivot_longer(cols = contains("NZMEAN_")) %>%
      select(name, time, value) %>%
      mutate(
        roi_id = sub(pattern = "NZMean_", replacement = "", x = name),
        par = participant_ids[i]
      ) %>%
      rename("raw_bold" = value) %>%
      merge(x = ., y = roi_key, by.x = "roi_id", by.y = "roi_id")

  current_bold_per_roi_df <- data.frame(
    "par" = long_mean_bold$par,
    "roi" = long_mean_bold$roi,
    "roi_id" = long_mean_bold$roi_id,
    "time_sec" = long_mean_bold$time,
    "censor" = long,
    "bold" = numeric()
  )
  current_bold_per_roi_df$par <- long_mean_bold$par
  current_bold_per_roi_df$roi <- long_mean_bold$roi
  current_bold_per_roi_df$roi_id <- long_mean_bold$roi_id

  bold_per_roi_df$

  roi_df <- read_delim(
    '/home/andrew/Downloads/stan_data/145_first_attempt/roi_stats.txt',
    trim_ws = T
  )

  design_matrix <- design_matrix[
    1:(nrow(design_matrix) - 1),
    2:ncol(design_matrix)
  ]

  design_matrix$X2 <- as.numeric(design_matrix$X2)

  censor_info <- read_delim(
    '/home/andrew/Downloads/stan_data/145_first_attempt/censor_GABORGEN24_DAY1_145_combined_2.1D',
    col_names = F
  )[[2]]

  fmri_stan_list <- list()

  fmri_stan_list$n_participants <- 1
  fmri_stan_list$n_ROIs <- 1

  # 69 = Anterior_Ventral_Insular_Area_L
  fmri_stan_list$amplitude_no_censor <- roi_df %>%
    select(paste0("Mean_", 69)) %>%
    pull()

  fmri_stan_list$censor <- c(1:1070)[censor_info == 1]

  fmri_stan_list$n_amplitude <- fmri_stan_list$amplitude_no_censor %>%
    length()

  fmri_stan_list$n_censor <- 1070 - length(fmri_stan_list$censor)

  #no polort high-pass
  fmri_stan_list$design_matrix <- design_matrix[, 17:ncol(design_matrix)]
  # fmri_stan_list$design_matrix <- design_matrix

  fmri_stan_list$n_DM_cols <- ncol(fmri_stan_list$design_matrix)
}


# Stan settings ####
number_of_chains <- 1
warmup_samples_per_chain <- 200
posterior_samples_per_chain <- 200
# where_to_save_chains <- '/home/andrew/Documents/stan_chains_ssd/'
# where_to_save_chains <- '/run/media/andrew/Barracuda_8tb/stan_chains/'
where_to_save_chains <- '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains'

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model014.stan'


# Fit models
model <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(
    stan_threads = TRUE #,
    # stan_opencl = TRUE
  )
)

model_fit <- model$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 5,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  output_dir = where_to_save_chains,
  chains = number_of_chains,
  parallel_chains = number_of_chains
)

model001_fit <- model001$optimize(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  thread = 8,
  jacobian = T,
  iter = 10000
  # iter_warmup = warmup_samples_per_chain,
  # iter_sampling = posterior_samples_per_chain,
  # save_warmup = F,
  # show_messages = T,
  # output_dir = where_to_save_chains,
  # chains = number_of_chains,
  # parallel_chains = number_of_chains
)


model001_fit_meta_data <- model001_fit$metadata()

model001_fit_relevant_parameters <- model001_fit_meta_data$model_params[
  !str_detect(model001_fit_meta_data$model_params, "log_lik|amp|Cov|Cor|D")
]

model001_fit_summary <-
  model001_fit$summary(variables = model001_fit_relevant_parameters)


sample_rate <- .5 #hz from a TR of 2seconds
nyquist <- sample_rate / 2
cutoff_frequency <- 0.0065 #hz what afni finds in 3dDeconvolve for (p - 2) / duration
w_parameter <- cutoff_frequency / nyquist

highpass <- butter(n = 5, W = 0.026, type = "high")


data.frame(
  time = seq(0, 2138, by = 2),
  no_filt = roi_df$Mean_69 - 100,
  filtered = filtfilt(highpass, roi_df$Mean_69 - 100)
) %>%
  pivot_longer(cols = contains("filt")) %>%
  ggplot() +
  geom_line(
    aes(
      x = time,
      color = name,
      y = value
    ),
    linewidth = .2
  ) +
  theme_classic()

## most current list ####
fmri_stan_list <- list()

fmri_stan_list$n_participants <- 1
fmri_stan_list$n_ROIs <- 1

# 69 = Anterior_Ventral_Insular_Area_L
fmri_stan_list$amplitude_no_censor <- roi_df %>%
  select(paste0("Mean_", 69)) %>%
  mutate(
    centered_amp = Mean_69 - 100,
    filtered_amp = filtfilt(highpass, centered_amp)
  ) %>%
  pull(filtered_amp)

fmri_stan_list$censor <- c(1:1070)[censor_info == 1]

fmri_stan_list$n_amplitude <- fmri_stan_list$amplitude_no_censor %>%
  length()

fmri_stan_list$n_censor <- 1070 - length(fmri_stan_list$censor)

#no polort high-pass
fmri_stan_list$design_matrix <- design_matrix[, 17:ncol(design_matrix)]

fmri_stan_list$n_DM_cols <- ncol(fmri_stan_list$design_matrix)

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model002.stan'


# Fit models
model002 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE)
)

model002_fit <- model002$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 2,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  output_dir = where_to_save_chains,
  chains = number_of_chains,
  parallel_chains = number_of_chains
)

model002_fit <- model002$optimize(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  thread = 8,
  jacobian = T,
  iter = 10000
  # iter_warmup = warmup_samples_per_chain,
  # iter_sampling = posterior_samples_per_chain,
  # save_warmup = F,
  # show_messages = T,
  # output_dir = where_to_save_chains,
  # chains = number_of_chains,
  # parallel_chains = number_of_chains
)


model002_fit_meta_data <- model002_fit$metadata()

model002_fit_relevant_parameters <- model002_fit_meta_data$model_params[
  !str_detect(model002_fit_meta_data$model_params, "log_lik|amp|Cov|Cor|D")
]

model002_fit_summary <-
  model002_fit$summary(variables = model002_fit_relevant_parameters)

# Model 3 ####
model_path <- '~/Documents/GitHub/gaborgen/stan_models/fMRI/Model003.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model003.stan'

# Fit models
model003 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)

model003_fit <- model003$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 1,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  opencl_ids = c(0, 0),
  output_dir = where_to_save_chains,
  chains = number_of_chains,
  parallel_chains = number_of_chains
)

# model003_fit <- model003$optimize(
#   data = fmri_stan_list,
#   refresh = 50,
#   seed = 3,
#   thread = 8,
#   jacobian = T,
#   iter = 10000
#   # iter_warmup = warmup_samples_per_chain,
#   # iter_sampling = posterior_samples_per_chain,
#   # save_warmup = F,
#   # show_messages = T,
#   # output_dir = where_to_save_chains,
#   # chains = number_of_chains,
#   # parallel_chains = number_of_chains
# )

# save(model003_fit, file = '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/model003_fmri_fit.RData')

model003_fit_meta_data <- model003_fit$metadata()

model003_fit_relevant_parameters <- model003_fit_meta_data$model_params[
  !str_detect(model003_fit_meta_data$model_params, "log_lik|amp|Cov|Cor|D")
]

model003_fit_summary <-
  model003_fit$summary(variables = model003_fit_relevant_parameters)


model003_draws_array <- model003_fit$draws()

# even though it converges anyway, now can get convergence with 15 samples
posterior::rhat_nested(
  posterior::extract_variable_matrix(
    model003_fit,
    variable = c("lp__")
  )[1:15, ],
  superchain_ids = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
)


model003_draws <- model003_fit$draws(format = "df")

model003_plot_betas_df <- model003_draws %>%
  pivot_longer(starts_with("Betas")) %>%
  dplyr::filter(
    name %in%
      c(
        "Betas[1]",
        "Betas[2]",
        "Betas[3]",
        "Betas[4]",
        "Betas[5]",
        "Betas[6]",
        "Betas[7]",
        "Betas[8]",
        "Betas[9]",
        "Betas[10]",
        "Betas[11]",
        "Betas[12]"
      )
  ) %>%
  mutate(
    phase = case_when(
      name %in%
        c(
          "Betas[1]",
          "Betas[2]",
          "Betas[3]",
          "Betas[4]"
        ) ~
        "Habituation",
      name %in%
        c(
          "Betas[5]",
          "Betas[6]",
          "Betas[7]",
          "Betas[8]"
        ) ~
        "Acquisition",
      name %in%
        c(
          "Betas[9]",
          "Betas[10]",
          "Betas[11]",
          "Betas[12]"
        ) ~
        "Extinction"
    ),
    cue = case_when(
      name %in%
        c(
          "Betas[1]",
          "Betas[5]",
          "Betas[9]"
        ) ~
        "CSP",
      name %in%
        c(
          "Betas[2]",
          "Betas[6]",
          "Betas[10]"
        ) ~
        "GS1",
      name %in%
        c(
          "Betas[3]",
          "Betas[7]",
          "Betas[11]"
        ) ~
        "GS2",
      name %in%
        c(
          "Betas[4]",
          "Betas[8]",
          "Betas[12]"
        ) ~
        "GS3"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "Habituation",
        "Acquisition",
        "Extinction"
      )
    )
  )

model003_plot_betas_df %>%
  group_by(phase, cue) %>%
  reframe(n())


model003_plot_betas_df %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(color = cue, x = value)) +
  facet_wrap(~phase, ncol = 1) +
  theme_classic()


# Model 6:3 but faster ####
model_path <- '~/Documents/GitHub/gaborgen/stan_models/fMRI/Model006.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model003.stan'

# Fit models
model006 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)

model006_fit <- model006$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 1,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  opencl_ids = c(0, 0),
  output_dir = where_to_save_chains,
  chains = 4,
  parallel_chains = 4
)

# Model 7:3 but faster ####
model_path <- '~/Documents/GitHub/gaborgen/stan_models/fMRI/Model007.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model003.stan'

# Fit models
model007 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)

model007_fit <- model007$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 1,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  opencl_ids = c(0, 0),
  output_dir = where_to_save_chains,
  chains = 4,
  parallel_chains = 4
)


model007_fit_meta_data <- model007_fit$metadata()

model007_fit_relevant_parameters <- model007_fit_meta_data$model_params[
  !str_detect(model007_fit_meta_data$model_params, "log_lik|amp|Cov|Cor|D|Mu")
]

model007_fit_summary <-
  model007_fit$summary(variables = model007_fit_relevant_parameters)

# Model 8:3 but f-aster ####
model_path <- '~/Documents/GitHub/gaborgen/stan_models/fMRI/Model008.stan'
# model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model003.stan'

# Fit models
model008 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)

model008_fit <- model008$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 1,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  opencl_ids = c(0, 0),
  output_dir = where_to_save_chains,
  chains = 4,
  parallel_chains = 4
)

# model003_fit <- model003$optimize(
#   data = fmri_stan_list,
#   refresh = 50,
#   seed = 3,
#   thread = 8,
#   jacobian = T,
#   iter = 10000
#   # iter_warmup = warmup_samples_per_chain,
#   # iter_sampling = posterior_samples_per_chain,
#   # save_warmup = F,
#   # show_messages = T,
#   # output_dir = where_to_save_chains,
#   # chains = number_of_chains,
#   # parallel_chains = number_of_chains
# )

# save(model003_fit, file = '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/model003_fmri_fit.RData')

model003_fit_meta_data <- model003_fit$metadata()

model003_fit_relevant_parameters <- model003_fit_meta_data$model_params[
  !str_detect(model003_fit_meta_data$model_params, "log_lik|amp|Cov|Cor|D")
]

model003_fit_summary <-
  model003_fit$summary(variables = model003_fit_relevant_parameters)


# Try faster latent correlation

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model004.stan'


# Fit models
model004 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE)
)

model004_fit <- model004$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 2,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  output_dir = where_to_save_chains,
  chains = number_of_chains,
  parallel_chains = number_of_chains
)

model004_fit <- model004$optimize(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  thread = 8,
  jacobian = T,
  iter = 10000,
  # iter_warmup = warmup_samples_per_chain,
  # iter_sampling = posterior_samples_per_chain,
  # save_warmup = F,
  # show_messages = T,
  output_dir = where_to_save_chains,
  # chains = number_of_chains,
  # parallel_chains = number_of_chains
)

model004_fit_meta_data <- model004_fit$metadata()

model004_fit_relevant_parameters <- model004_fit_meta_data$model_params[
  !str_detect(model004_fit_meta_data$model_params, "log_lik|amp|Cov|Cor|D")
]

model004_fit_summary <-
  model004_fit$summary(variables = model004_fit_relevant_parameters)
