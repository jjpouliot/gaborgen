library(tidyverse)
library(cmdstanr)
library(signal)


roi_df <- read_delim(
  '/home/andrewf/Downloads/stan_data/roi_stats.txt',
  trim_ws = T
)
roi_key <- read_delim('/home/andrewf/Downloads/stan_data/HCPex_SUIT_labels.txt')
design_matrix <- read_delim(
  '/home/andrewf/Downloads/stan_data/X.nocensor.xmat.1D',
  skip = 64,
  col_names = F
)
design_matrix <- design_matrix[
  1:(nrow(design_matrix) - 1),
  2:ncol(design_matrix)
]

design_matrix$X2 <- as.numeric(design_matrix$X2)

censor_info <- read_delim(
  '/home/andrewf/Downloads/stan_data/censor_GABORGEN24_DAY1_145_combined_2.1D',
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

fmri_stan_list$design_matrix <- design_matrix

fmri_stan_list$n_DM_cols <- ncol(fmri_stan_list$design_matrix)


# Stan settings
number_of_chains <- 4
warmup_samples_per_chain <- 1000
posterior_samples_per_chain <- 1000
where_to_save_chains <- '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains'

model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model001.stan'


# Fit models
model001 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE)
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


model_path <- '/home/andrewf/Repositories/gaborgen/stan_models/fMRI/Model003.stan'


# Fit models
model003 <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  force_recompile = T,
  cpp_options = list(stan_threads = TRUE)
)

# model003_fit <- model003$sample(
#   data = fmri_stan_list,
#   refresh = 50,
#   seed = 3,
#   threads_per_chain = 2,
#   iter_warmup = warmup_samples_per_chain,
#   iter_sampling = posterior_samples_per_chain,
#   save_warmup = T,
#   show_messages = T,
#   output_dir = where_to_save_chains,
#   chains = number_of_chains,
#   parallel_chains = number_of_chains
# )

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
