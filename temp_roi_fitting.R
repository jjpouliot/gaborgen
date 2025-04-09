library(tidyverse)
library(cmdstanr)
# library(signal)

# we recommend running this in a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# Filter used, same notation as matlab
highpass <- signal::butter(n = 5, W = 0.026, type = "high")


# load data from text files ####

data_dir <- c(
  "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info"
)

roi_key <- read_delim(paste0(data_dir, '/HCPex_SUIT_labels.txt'))

colnames(roi_key) <- c("roi_id", "roi")

useable_participants <- c("131", "145", "149")

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
    path = data_dir,
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
        data_dir,
        '/censor_GABORGEN24_',
        useable_participants[i],
        '_combined_2.1D'
      ),
      col_names = F
    )[[2]]
  } else {
    censor_info <- read_delim(
      paste0(
        data_dir,
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
    current_design_matrix <- read_delim(
      paste0(
        data_dir,
        '/GABORGEN24_',
        useable_participants[i],
        '.results_X.nocensor.xmat.1D'
      ),
      skip = 64,
      col_names = F
    )
  } else {
    current_design_matrix <- read_delim(
      paste0(
        data_dir,
        '/GABORGEN24_DAY1_',
        useable_participants[i],
        '.results_X.nocensor.xmat.1D'
      ),
      skip = 64,
      col_names = F
    )
  }
  current_design_matrix <- current_design_matrix[
    1:(nrow(current_design_matrix) - 1),
    2:ncol(current_design_matrix)
  ]

  current_design_matrix$X2 <- as.numeric(current_design_matrix$X2)

  # don't need high pass variables
  current_design_matrix <- current_design_matrix[,
    17:ncol(current_design_matrix)
  ]

  colnames(current_design_matrix) <- c(
    "phase_1_stim_1",
    "phase_1_stim_2",
    "phase_1_stim_3",
    "phase_1_stim_4",
    "phase_2_stim_1",
    "phase_2_stim_2",
    "phase_2_stim_3",
    "phase_2_stim_4",
    "phase_3_stim_1",
    "phase_3_stim_2",
    "phase_3_stim_3",
    "phase_3_stim_4",
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

# create stan list ####
used_df <- bold_per_roi_df %>%
  filter(roi_id %in% c(1, 69))


fmri_stan_list <- list()

fmri_stan_list$n_par <- used_df$par %>% unique() %>% length()

fmri_stan_list$n_roi <- used_df$roi_id %>% unique() %>% length()

fmri_stan_list$n_bold <- 1070

fmri_stan_list$n <- fmri_stan_list$n_par *
  fmri_stan_list$n_roi *
  fmri_stan_list$n_bold

fmri_stan_list$par <- used_df$par %>% as.factor() %>% as.integer()

fmri_stan_list$roi <- used_df$roi_id %>% as.factor() %>% as.integer()

# just vector
# fmri_stan_list$bold <- used_df$bold
# [par, time]
fmri_stan_list$bold <- used_df$bold %>%
  array(
    dim = c(
      fmri_stan_list$n_bold * fmri_stan_list$n_roi,
      # fmri_stan_list$n_roi,
      fmri_stan_list$n_par
    )
  ) %>%
  aperm(c(2, 1))
# [par, roi, time]
# fmri_stan_list$bold <- used_df$bold %>%
#   array(
#     dim = c(
#       1070,
#       fmri_stan_list$n_roi,
#       fmri_stan_list$n_par
#     )
#   ) %>%
#   aperm(c(3, 2, 1))

# just vector
# fmri_stan_list$usable_bold_indices <- used_df$censor %>% as.integer()
# [par, roi, time] to [par, time]; same for each roi
fmri_stan_list$usable_bold_indices <- used_df$censor %>%
  array(
    dim = c(
      1070,
      fmri_stan_list$n_roi,
      fmri_stan_list$n_par
    )
  ) %>%
  aperm(c(3, 2, 1))
fmri_stan_list$usable_bold_indices <- fmri_stan_list$usable_bold_indices[, 1, ]


# array [n_par] int
censor_vec <- vector()
for (p in 1:fmri_stan_list$n_par) {
  censor_vec <- c(
    censor_vec,
    (1070 - sum(fmri_stan_list$usable_bold_indices[p, ] == 1))
  )
}
fmri_stan_list$n_censor <- censor_vec %>% as.integer()

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

# Stan settings ####
number_of_chains <- 2
warmup_samples_per_chain <- 200
posterior_samples_per_chain <- 200
# where_to_save_chains <- '/home/andrew/Documents/stan_chains_ssd/'
# where_to_save_chains <- '/run/media/andrew/Barracuda_8tb/stan_chains/'
where_to_save_chains <- '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains'

# First map_rect model ####

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model012.stan'


# Fit models
model012 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE)
  # cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)

model012_fit <- model012$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 3,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = T,
  show_messages = T,
  output_dir = where_to_save_chains,
  chains = number_of_chains,
  parallel_chains = number_of_chains
)

model012_fit$output()

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
number_of_chains <- 16
warmup_samples_per_chain <- 200
posterior_samples_per_chain <- 200
where_to_save_chains <- '/home/andrew/Documents/stan_chains_ssd/'
# where_to_save_chains <- '/run/media/andrew/Barracuda_8tb/stan_chains/'
# where_to_save_chains <- '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains'

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model001.stan'


# Fit models
model001 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE, stan_opencl = TRUE)
)

model001_fit <- model001$sample(
  data = fmri_stan_list,
  refresh = 50,
  seed = 3,
  threads_per_chain = 2,
  iter_warmup = warmup_samples_per_chain,
  iter_sampling = posterior_samples_per_chain,
  save_warmup = F,
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
