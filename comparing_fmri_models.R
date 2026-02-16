library(tidyverse)
library(cmdstanr)
library(patchwork)
library(ggridges)

# we recommend running this in a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# load(
#   "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_HBM.RData"
# )
# load(
#   "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_HBM.RData"
# )

cue_color <- c("red1", "green1", "purple1", "blue1", "black")

# data_dir <- "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains"
data_dir <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains"

# 41 participants

#left right anterior insula rois: 68, 69, 73, 248, 249, 253

# Model 35: For overlapping CS+/US just uses shock coefficent no CS+ cue there
# The model actually still estimates that CS+ betas that are paired, but those
# betas don't get used for prediction. They the get created from the ML-structure
# based on the correlation of trial unpaired around them.

model35_V5_fit <- as_cmdstan_fit(
  files = c(
    # paste0(data_dir, "/model035_V5_chain_14588605_1.csv")
    paste0(data_dir, "/model035_V5_chain_14599422_1.csv")
  )
)

model35_V1_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model035_V1_chain_14588657_1.csv")
  )
)

model35_AntIns_fit <- as_cmdstan_fit(
  files = c(
    # paste0(data_dir, "/model035_V5_chain_14588605_1.csv")
    paste0(data_dir, "/model035_AntIns_chain_14598819_1.csv")
  )
)


# Model 30: tries to deconvolve overlapping CS+ and US
model30_V1_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030V1_chain_14415120_1.csv")
  )
)

model30_V4_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030V4_chain_14448457_1.csv")
  )
)

model30_V5_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030V5_chain_14415685_1.csv")
  )
)

model30_V6_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030V6_chain_14448485_1.csv")
  )
)

model30_TE_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030TE_chain_14448524_1.csv")
  )
)

model30_TPJ_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030TPJ_chain_14448526_1.csv")
  )
)

model30_HIP_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030HIP_chain_14448581_1.csv")
  )
)

model30_amyg_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030amyg_chain_14428904_1.csv")
  )
)

model30_NA_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030NA_chain_14448709_1.csv")
  )
)

model30_ant_ins_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030_chain_14414680_1.csv")
  )
)

model30_ACC_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030ACC_chain_14448724_1.csv")
  )
)

model30_ORB_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030ORB_chain_14448730_1.csv")
  )
)

prepare_data_and_plot_cue_per_trial <- function(model_fit, ROI_name_string) {
  model_fit_meta_data <- model_fit$metadata()
  model_fit_beta_parameters <- model_fit_meta_data$model_params[
    str_detect(
      model_fit_meta_data$model_params,
      "betas\\["
      # "log_lik|L_|K_"
    )
  ]

  model_fit_beta_draws <- model_fit$draws(
    variables = model_fit_beta_parameters,
    format = "df"
  )

  start_index <- 1
  stop_index <- 8
  csp_indices <- c(start_index:stop_index)
  start_index <- start_index + 8
  stop_index <- stop_index + 8
  gs1_indices <- c(start_index:stop_index)
  start_index <- start_index + 8
  stop_index <- stop_index + 8
  gs2_indices <- c(start_index:stop_index)
  start_index <- start_index + 8
  stop_index <- stop_index + 8
  gs3_indices <- c(start_index:stop_index)
  start_index <- start_index + 8
  stop_index <- stop_index + 12
  csp_indices <- c(csp_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs1_indices <- c(gs1_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs2_indices <- c(gs2_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs3_indices <- c(gs3_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  csp_indices <- c(csp_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs1_indices <- c(gs1_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs2_indices <- c(gs2_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs3_indices <- c(gs3_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  csp_indices <- c(csp_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs1_indices <- c(gs1_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs2_indices <- c(gs2_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 12
  gs3_indices <- c(gs3_indices, start_index:stop_index)
  start_index <- start_index + 12
  stop_index <- stop_index + 15
  shock_indices <- c(start_index:stop_index)

  beta_key <- data.frame(
    beta_index = 1:176,
    trial_per_cue = c(
      rep(1:8, 4),
      rep(9:20, 4),
      rep(21:32, 4),
      rep(33:44, 4)
    )
  )

  model_fit_beta_draws_long <- model_fit_beta_draws %>%
    pivot_longer(
      cols = starts_with("betas["),
      names_to = c(".unused", "roi", "trial_index"),
      names_sep = "\\[|,|\\]",
      values_to = "beta_value"
    ) %>%
    select(-.unused) %>%
    mutate(
      roi = as.integer(roi),
      trial_index = as.integer(trial_index)
    ) %>%
    mutate(
      cue = case_when(
        trial_index %in% csp_indices ~ "csp",
        trial_index %in% gs1_indices ~ "gs1",
        trial_index %in% gs2_indices ~ "gs2",
        trial_index %in% gs3_indices ~ "gs3",
        trial_index %in% shock_indices ~ "shock"
      )
    ) %>%
    mutate(
      block = case_when(
        trial_index %in% c(1:32) ~ "habituation",
        trial_index %in% c(33:80) ~ "acquisition #1",
        trial_index %in% c(81:128) ~ "acquisition #2",
        trial_index %in% c(129:176) ~ "extinction"
      )
    ) %>%
    mutate(
      side = case_when(
        roi == 1 ~ "left",
        roi == 2 ~ "right"
      )
    )

  model_fit_beta_draws_long <- merge(
    x = model_fit_beta_draws_long,
    y = beta_key,
    by.x = "trial_index",
    by.y = "beta_index",
    all.x = T
  )

  lower_bound <- .17
  upper_bound <- .83

  left_roi_plot <- model_fit_beta_draws_long %>%
    filter(cue != "shock", side == "left") %>%
    group_by(
      trial_per_cue,
      cue,
      trial_index
    ) %>%
    reframe(
      median_posterior = median(beta_value),
      lower_2_5 = quantile(beta_value, lower_bound),
      lower_97_5 = quantile(beta_value, upper_bound)
      # lower_2_5 = quantile(beta_value, .1),
      # lower_97_5 = quantile(beta_value, .9)
    ) %>%
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 8) +
    geom_vline(xintercept = 8 + 12) +
    geom_vline(xintercept = 8 + 12 + 12) +
    geom_ribbon(
      aes(
        x = trial_per_cue,
        y = median_posterior,
        ymin = lower_2_5,
        ymax = lower_97_5,
        fill = cue
      ),
      linewidth = .5,
      color = NA,
      alpha = .5
    ) +
    scale_fill_manual(values = cue_color) +
    ggtitle(paste0("Change in left ", ROI_name_string, " Over Trials")) +
    theme_bw() +
    theme(text = element_text(family = "Arial", size = 20))

  right_roi_plot <- model_fit_beta_draws_long %>%
    filter(cue != "shock", side == "right") %>%
    group_by(
      trial_per_cue,
      cue,
      trial_index
    ) %>%
    reframe(
      median_posterior = median(beta_value),
      lower_2_5 = quantile(beta_value, lower_bound),
      lower_97_5 = quantile(beta_value, upper_bound)
      # lower_2_5 = quantile(beta_value, .1),
      # lower_97_5 = quantile(beta_value, .9)
    ) %>%
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 8) +
    geom_vline(xintercept = 8 + 12) +
    geom_vline(xintercept = 8 + 12 + 12) +
    geom_ribbon(
      aes(
        x = trial_per_cue,
        y = median_posterior,
        ymin = lower_2_5,
        ymax = lower_97_5,
        fill = cue
      ),
      linewidth = .5,
      color = NA,
      alpha = .5
    ) +
    scale_fill_manual(values = cue_color) +
    ggtitle(paste0("Change in right ", ROI_name_string, " Over Trials")) +
    theme_bw() +
    theme(text = element_text(family = "Arial", size = 20))

  avg_roi_plot <- model_fit_beta_draws_long %>%
    filter(
      cue != "shock",
    ) %>%
    group_by(.draw, cue, trial_per_cue) %>%
    reframe(
      avg_roi_draw = mean(beta_value)
    ) %>%
    group_by(cue, trial_per_cue) %>%
    reframe(
      median_posterior = median(avg_roi_draw),
      lower_2_5 = quantile(avg_roi_draw, lower_bound),
      lower_97_5 = quantile(avg_roi_draw, upper_bound)
      # lower_2_5 = quantile(beta_value, .1),
      # lower_97_5 = quantile(beta_value, .9)
    ) %>%
    ggplot() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 8) +
    geom_vline(xintercept = 8 + 12) +
    geom_vline(xintercept = 8 + 12 + 12) +
    geom_ribbon(
      aes(
        x = trial_per_cue,
        y = median_posterior,
        ymin = lower_2_5,
        ymax = lower_97_5,
        fill = cue
      ),
      linewidth = .5,
      color = NA,
      alpha = .5
    ) +
    scale_fill_manual(values = cue_color) +
    ggtitle(paste0("Change in average ", ROI_name_string, " Over Trials")) +
    theme_bw() +
    theme(text = element_text(family = "Arial", size = 20))

  shock_posteriors_plot <- model_fit_beta_draws_long %>%
    filter(
      cue == "shock",
      trial_index %in% shock_indices
    ) %>%
    group_by(.draw, trial_index) %>%
    mutate(
      avg_roi_draw = mean(beta_value),
      trial_index = factor(trial_index, levels = unique(trial_index))
    ) %>%
    ggplot() +
    geom_vline(aes(xintercept = 0)) +
    geom_density_ridges(aes(x = avg_roi_draw, y = trial_index)) +
    # coord_cartesian(xlim = c(-.2, .4)) +
    theme_bw() +
    ggtitle("Bold Change to Shocks") +
    theme_bw() +
    theme(text = element_text(family = "Arial", size = 20))

  ((left_roi_plot / right_roi_plot / avg_roi_plot) | shock_posteriors_plot) +
    patchwork::plot_annotation(
      title = paste0(ROI_name_string, " (N = 41)"),
      theme = theme(plot.title = element_text(face = "bold", size = 30))
    )
}


prepare_data_and_plot_cue_per_trial(model30_V1_fit, "Primary Visual 1")

prepare_data_and_plot_cue_per_trial(model35_V1_fit, "Primary Visual 1")

prepare_data_and_plot_cue_per_trial(model30_V4_fit, "Visual 4")

prepare_data_and_plot_cue_per_trial(model30_V5_fit, "V5/MT")

prepare_data_and_plot_cue_per_trial(model35_V5_fit, "V5/MT")

prepare_data_and_plot_cue_per_trial(model30_V6_fit, "Visual 6")

prepare_data_and_plot_cue_per_trial(model30_TE_fit, "Temporal (TE)")

prepare_data_and_plot_cue_per_trial(
  model30_TPJ_fit,
  "Temporal Parietal Junction"
)

prepare_data_and_plot_cue_per_trial(model30_HIP_fit, "Hippocampus")

prepare_data_and_plot_cue_per_trial(model30_amyg_fit, "Amygdala")

prepare_data_and_plot_cue_per_trial(model30_NA_fit, "Nucleus Accumbens")

prepare_data_and_plot_cue_per_trial(model30_ant_ins_fit, "Anterior Insula")

prepare_data_and_plot_cue_per_trial(model35_AntIns_fit, "Anterior Insula")

prepare_data_and_plot_cue_per_trial(model30_NA_fit, "Anterior Cingulate")

prepare_data_and_plot_cue_per_trial(model30_ORB_fit, "Orbital")


model30_ant_ins_fit_meta_data <- model30_ant_ins_fit$metadata()

model30_ant_ins_fit_relevant_parameters <- model30_ant_ins_fit_meta_data$model_params[
  !str_detect(
    model30_ant_ins_fit_meta_data$model_params,
    "log_lik|L|K|delta_z|betas_motion|betas_z|rho_time_z|theta|sigma_z|sigma\\[|delta\\[|rho_time|raw"
  )
]

model30_ant_ins_summary <- model30_ant_ins_fit$summary(
  variables = model30_ant_ins_fit_relevant_parameters
)

model30_ant_ins_draws <- model30_ant_ins_fit$draws(
  variables = model30_ant_ins_fit_relevant_parameters,
  format = "df"
)

model30_ant_ins_beta_parameters <- model30_ant_ins_fit_meta_data$model_params[
  str_detect(
    model30_ant_ins_fit_meta_data$model_params,
    "betas\\["
    # "log_lik|L_|K_"
  )
]

model30_ant_ins_beta_draws <- model30_ant_ins_fit$draws(
  variables = model30_ant_ins_beta_parameters,
  format = "df"
)

model30_V5_fit_meta_data <- model30_V5_fit$metadata()

model30_V5_fit_relevant_parameters <- model30_V5_fit_meta_data$model_params[
  !str_detect(
    model30_V5_fit_meta_data$model_params,
    "log_lik|L|K|delta_z|betas_motion|betas_z|rho_time_z|theta|sigma_z|sigma\\[|delta\\[|rho_time|raw"
  )
]

model30_V5_summary <- model30_V5_fit$summary(
  variables = model30_V5_fit_relevant_parameters
)

model30_V5_draws <- model30_V5_fit$draws(
  variables = model30_V5_fit_relevant_parameters,
  format = "df"
)

model30_V5_beta_parameters <- model30_V5_fit_meta_data$model_params[
  str_detect(
    model30_V5_fit_meta_data$model_params,
    "betas\\["
    # "log_lik|L_|K_"
  )
]

model30_V5_beta_draws <- model30_V5_fit$draws(
  variables = model30_V5_beta_parameters,
  format = "df"
)


model30_amyg_fit_meta_data <- model30_amyg_fit$metadata()

model30_amyg_fit_relevant_parameters <- model30_amyg_fit_meta_data$model_params[
  !str_detect(
    model30_amyg_fit_meta_data$model_params,
    "log_lik|L|K|delta_z|betas_motion|betas_z|rho_time_z|theta|sigma_z|sigma\\[|delta\\[|rho_time|raw"
  )
]

model30_amyg_summary <- model30_amyg_fit$summary(
  variables = model30_amyg_fit_relevant_parameters
)

model30_amyg_draws <- model30_amyg_fit$draws(
  variables = model30_amyg_fit_relevant_parameters,
  format = "df"
)

model30_amyg_beta_parameters <- model30_amyg_fit_meta_data$model_params[
  str_detect(
    model30_amyg_fit_meta_data$model_params,
    "betas\\["
    # "log_lik|L_|K_"
  )
]

model30_amyg_beta_draws <- model30_amyg_fit$draws(
  variables = model30_amyg_beta_parameters,
  format = "df"
)


start_index <- 1
stop_index <- 8
csp_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 8
gs1_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 8
gs2_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 8
gs3_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 12
csp_indices <- c(csp_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs1_indices <- c(gs1_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs2_indices <- c(gs2_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs3_indices <- c(gs3_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
csp_indices <- c(csp_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs1_indices <- c(gs1_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs2_indices <- c(gs2_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs3_indices <- c(gs3_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
csp_indices <- c(csp_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs1_indices <- c(gs1_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs2_indices <- c(gs2_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs3_indices <- c(gs3_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 15
shock_indices <- c(start_index:stop_index)


beta_key <- data.frame(
  beta_index = 1:176,
  trial_per_cue = c(
    rep(1:8, 4),
    rep(9:20, 4),
    rep(21:32, 4),
    rep(33:44, 4)
  )
)


model30_ant_ins_beta_draws_long <- model30_ant_ins_beta_draws %>% # your tibble
  pivot_longer(
    cols = starts_with("betas["),
    names_to = c(".unused", "roi", "trial_index"),
    names_sep = "\\[|,|\\]",
    values_to = "beta_value"
  ) %>%
  select(-.unused) %>%
  mutate(
    roi = as.integer(roi),
    trial_index = as.integer(trial_index)
  ) %>%
  mutate(
    cue = case_when(
      trial_index %in% csp_indices ~ "csp",
      trial_index %in% gs1_indices ~ "gs1",
      trial_index %in% gs2_indices ~ "gs2",
      trial_index %in% gs3_indices ~ "gs3",
      trial_index %in% shock_indices ~ "shock"
    )
  ) %>%
  mutate(
    block = case_when(
      trial_index %in% c(1:32) ~ "habituation",
      trial_index %in% c(33:80) ~ "acquisition #1",
      trial_index %in% c(81:128) ~ "acquisition #2",
      trial_index %in% c(129:176) ~ "extinction"
    )
  ) %>%
  mutate(
    side = case_when(
      roi == 1 ~ "left",
      roi == 2 ~ "right"
    )
  )

model30_ant_ins_beta_draws_long <- merge(
  x = model30_ant_ins_beta_draws_long,
  y = beta_key,
  by.x = "trial_index",
  by.y = "beta_index",
  all.x = T
)

lower_bound <- .2
upper_bound <- .8

model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  group_by(
    trial_per_cue,
    cue,
    trial_index
  ) %>%
  reframe(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, lower_bound),
    lower_97_5 = quantile(beta_value, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Left Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock", side == "right") %>%
  group_by(
    trial_per_cue,
    cue,
    trial_index
  ) %>%
  reframe(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, lower_bound),
    lower_97_5 = quantile(beta_value, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Right Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_ant_ins_beta_draws_long %>%
  filter(
    cue != "shock",
  ) %>%
  group_by(.draw, cue, trial_per_cue) %>%
  reframe(
    avg_roi_draw = mean(beta_value)
  ) %>%
  group_by(cue, trial_per_cue) %>%
  reframe(
    median_posterior = median(avg_roi_draw),
    lower_2_5 = quantile(avg_roi_draw, lower_bound),
    lower_97_5 = quantile(avg_roi_draw, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  filter(trial_index %in% csp_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


# V5
model30_V5_beta_draws_long <- model30_V5_beta_draws %>% # your tibble
  pivot_longer(
    cols = starts_with("betas["),
    names_to = c(".unused", "roi", "trial_index"),
    names_sep = "\\[|,|\\]",
    values_to = "beta_value"
  ) %>%
  select(-.unused) %>%
  mutate(
    roi = as.integer(roi),
    trial_index = as.integer(trial_index)
  ) %>%
  mutate(
    cue = case_when(
      trial_index %in% csp_indices ~ "csp",
      trial_index %in% gs1_indices ~ "gs1",
      trial_index %in% gs2_indices ~ "gs2",
      trial_index %in% gs3_indices ~ "gs3",
      trial_index %in% shock_indices ~ "shock"
    )
  ) %>%
  mutate(
    block = case_when(
      trial_index %in% c(1:32) ~ "habituation",
      trial_index %in% c(33:80) ~ "acquisition #1",
      trial_index %in% c(81:128) ~ "acquisition #2",
      trial_index %in% c(129:176) ~ "extinction"
    )
  ) %>%
  mutate(
    side = case_when(
      roi == 1 ~ "left",
      roi == 2 ~ "right"
    )
  )

model30_V5_beta_draws_long <- merge(
  x = model30_V5_beta_draws_long,
  y = beta_key,
  by.x = "trial_index",
  by.y = "beta_index",
  all.x = T
)

lower_bound <- .2
upper_bound <- .8

model30_V5_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  group_by(
    trial_per_cue,
    cue,
    trial_index
  ) %>%
  reframe(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, lower_bound),
    lower_97_5 = quantile(beta_value, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Left Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_V5_beta_draws_long %>%
  filter(cue != "shock", side == "right") %>%
  group_by(
    trial_per_cue,
    cue,
    trial_index
  ) %>%
  reframe(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, lower_bound),
    lower_97_5 = quantile(beta_value, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Right Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_V5_beta_draws_long %>%
  filter(
    cue != "shock",
  ) %>%
  group_by(.draw, cue, trial_per_cue) %>%
  reframe(
    avg_roi_draw = mean(beta_value)
  ) %>%
  group_by(cue, trial_per_cue) %>%
  reframe(
    median_posterior = median(avg_roi_draw),
    lower_2_5 = quantile(avg_roi_draw, lower_bound),
    lower_97_5 = quantile(avg_roi_draw, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Bilateral V5 Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  filter(trial_index %in% csp_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


# amyg
model30_amyg_beta_draws_long <- model30_amyg_beta_draws %>% # your tibble
  pivot_longer(
    cols = starts_with("betas["),
    names_to = c(".unused", "roi", "trial_index"),
    names_sep = "\\[|,|\\]",
    values_to = "beta_value"
  ) %>%
  select(-.unused) %>%
  mutate(
    roi = as.integer(roi),
    trial_index = as.integer(trial_index)
  ) %>%
  mutate(
    cue = case_when(
      trial_index %in% csp_indices ~ "csp",
      trial_index %in% gs1_indices ~ "gs1",
      trial_index %in% gs2_indices ~ "gs2",
      trial_index %in% gs3_indices ~ "gs3",
      trial_index %in% shock_indices ~ "shock"
    )
  ) %>%
  mutate(
    block = case_when(
      trial_index %in% c(1:32) ~ "habituation",
      trial_index %in% c(33:80) ~ "acquisition #1",
      trial_index %in% c(81:128) ~ "acquisition #2",
      trial_index %in% c(129:176) ~ "extinction"
    )
  ) %>%
  mutate(
    side = case_when(
      roi == 1 ~ "left",
      roi == 2 ~ "right"
    )
  )

model30_amyg_beta_draws_long <- merge(
  x = model30_amyg_beta_draws_long,
  y = beta_key,
  by.x = "trial_index",
  by.y = "beta_index",
  all.x = T
)

lower_bound <- .2
upper_bound <- .8

model30_amyg_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  group_by(
    trial_per_cue,
    cue,
    trial_index
  ) %>%
  reframe(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, lower_bound),
    lower_97_5 = quantile(beta_value, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Left Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_amyg_beta_draws_long %>%
  filter(cue != "shock", side == "right") %>%
  group_by(
    trial_per_cue,
    cue,
    trial_index
  ) %>%
  reframe(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, lower_bound),
    lower_97_5 = quantile(beta_value, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Right Anterior Insula Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_amyg_beta_draws_long %>%
  filter(
    cue != "shock",
  ) %>%
  group_by(.draw, cue, trial_per_cue) %>%
  reframe(
    avg_roi_draw = mean(beta_value)
  ) %>%
  group_by(cue, trial_per_cue) %>%
  reframe(
    median_posterior = median(avg_roi_draw),
    lower_2_5 = quantile(avg_roi_draw, lower_bound),
    lower_97_5 = quantile(avg_roi_draw, upper_bound)
    # lower_2_5 = quantile(beta_value, .1),
    # lower_97_5 = quantile(beta_value, .9)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Bilateral amyg Over Trials (N =41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


model30_amyg_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  filter(trial_index %in% csp_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_amyg_beta_draws_long %>%
  filter(cue != "shock", side == "right") %>%
  filter(trial_index %in% csp_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock", side == "left") %>%
  filter(trial_index %in% gs1_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock", side == "right") %>%
  filter(trial_index %in% gs1_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_ant_ins_beta_draws_long %>%
  filter(cue != "shock") %>%
  filter(trial_index %in% gs1_indices) %>%
  mutate(
    trial_per_cue = factor(trial_per_cue, levels = unique(trial_per_cue))
  ) %>%
  group_by(.draw, trial_per_cue) %>%
  reframe(
    avg_roi_draw = mean(beta_value)
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = avg_roi_draw, y = trial_per_cue)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


model30_ant_ins_beta_draws_long %>%
  filter(
    cue == "shock",
    trial_index %in% shock_indices,
    side == "left"
  ) %>%
  group_by(.draw, trial_index) %>%
  mutate(
    avg_roi_draw = mean(beta_value),
    trial_index = factor(trial_index, levels = unique(trial_index))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = avg_roi_draw, y = trial_index)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model30_ant_ins_beta_draws_long %>%
  filter(
    cue == "shock",
    trial_index %in% shock_indices,
    side == "right"
  ) %>%
  group_by(.draw, trial_index) %>%
  mutate(
    avg_roi_draw = mean(beta_value),
    trial_index = factor(trial_index, levels = unique(trial_index))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = avg_roi_draw, y = trial_index)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


model30_ant_ins_beta_draws_long %>%
  filter(
    cue == "shock",
    trial_index %in% shock_indices
  ) %>%
  group_by(.draw, trial_index) %>%
  mutate(
    avg_roi_draw = mean(beta_value),
    trial_index = factor(trial_index, levels = unique(trial_index))
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = avg_roi_draw, y = trial_index)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=41)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


# old

model030_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model030_chain_13938457_1.csv")
  )
)

model030_fit_meta_data <- model030_fit$metadata()

model030_fit_relevant_parameters <- model030_fit_meta_data$model_params[
  !str_detect(
    model030_fit_meta_data$model_params,
    "log_lik|L|K|delta_z|betas_motion|betas_z|rho_time_z|theta|sigma_z|sigma\\[|delta\\[|rho_time|raw"
  )
]

model030_summary <- model030_fit$summary(
  variables = model030_fit_relevant_parameters
)

model030_draws <- model030_fit$draws(
  variables = model030_fit_relevant_parameters,
  format = "df"
)

# used roi key
used_df %>%
  group_by(m)


# first GP model
model026_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model026_chain_11922182_1.csv"),
    paste0(data_dir, "/model026_chain_11922182_2.csv"),
    paste0(data_dir, "/model026_chain_11922182_3.csv"),
    paste0(data_dir, "/model026_chain_11922182_4.csv"),
    paste0(data_dir, "/model026_chain_11922182_5.csv")
  )
)
# faster model 26?
model027_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model027_chain_12183289_1.csv"),
    paste0(data_dir, "/model027_chain_12183289_2.csv"),
    paste0(data_dir, "/model027_chain_12183289_3.csv"),
    paste0(data_dir, "/model027_chain_12183289_5.csv")
  )
)

model028_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model028_chain_12509072_1.csv"),
    paste0(data_dir, "/model028_chain_12509072_2.csv"),
    paste0(data_dir, "/model028_chain_12509072_3.csv"),
    paste0(data_dir, "/model028_chain_12509072_4.csv"),
    paste0(data_dir, "/model028_chain_12509072_5.csv")
  )
)


model026_fit_meta_data <- model026_fit$metadata()

model027_fit_meta_data <- model027_fit$metadata()

model028_fit_meta_data <- model028_fit$metadata()

model026_fit_relevant_parameters <- model026_fit_meta_data$model_params[
  str_detect(
    model026_fit_meta_data$model_params,
    "betas\\["
    # "betas\\[|lp_"
  )
]

model026_fit_relevant_parameters <- model026_fit_meta_data$model_params[
  !str_detect(
    model026_fit_meta_data$model_params,
    "betas\\["
    # "log_lik|L_|K_"
  )
]


model027_fit_relevant_parameters <- model027_fit_meta_data$model_params[
  str_detect(
    model027_fit_meta_data$model_params,
    "betas\\["
    # "betas\\[|lp_"
  )
]

model027_fit_relevant_parameters <- model027_fit_meta_data$model_params[
  !str_detect(
    model027_fit_meta_data$model_params,
    "log_lik|L|K|delta_z|betas_motion|betas_z|rho_time_z|theta"
  )
]

model028_fit_relevant_parameters <- model028_fit_meta_data$model_params[
  !str_detect(
    model028_fit_meta_data$model_params,
    "log_lik|L|K|delta_z|betas_motion|betas_z|rho_time_z|theta"
  )
]

model026_loo <- model026_fit$loo()

model027_loo <- model027_fit$loo()

model028_loo <- model028_fit$loo()

loo::loo_compare(
  model018_loo,
  model026_loo
)
loo::loo_compare(
  model026_loo,
  model027_loo,
  model028_loo
)

model026_summary <- model026_fit$summary(
  variables = model026_fit_relevant_parameters
)

model026_draws <- model026_fit$draws(
  variables = model026_fit_relevant_parameters,
  format = "df"
)

model027_summary <- model027_fit$summary(
  variables = model027_fit_relevant_parameters
)

model027_draws <- model027_fit$draws(
  variables = model027_fit_relevant_parameters,
  format = "df"
)

model028_summary <- model028_fit$summary(
  variables = model028_fit_relevant_parameters
)

model028_draws <- model028_fit$draws(
  variables = model027_fit_relevant_parameters,
  format = "df"
)

exp(-.5 * (1 / model028_draws$csp_rho)^2) %>% hist()
exp(-.5 * (10 / model028_draws$csp_rho)^2) %>% hist()

# csp median
exp(-.5 * (1 / 9.02)^2)
exp(-.5 * (2 / 9.02)^2)
exp(-.5 * (3 / 9.02)^2)
exp(-.5 * (4 / 9.02)^2)
exp(-.5 * (5 / 9.02)^2)
exp(-.5 * (6 / 9.02)^2)
exp(-.5 * (7 / 9.02)^2)
exp(-.5 * (8 / 9.02)^2)
exp(-.5 * (9 / 9.02)^2)
exp(-.5 * (10 / 9.02)^2)
exp(-.5 * (11 / 9.02)^2)
exp(-.5 * (12 / 9.02)^2)
exp(-.5 * (13 / 9.02)^2)
exp(-.5 * (14 / 9.02)^2)
exp(-.5 * (15 / 9.02)^2)
exp(-.5 * (16 / 9.02)^2)
exp(-.5 * (17 / 9.02)^2)
exp(-.5 * (18 / 9.02)^2)
# gs1 median
exp(-.5 * (1 / 8.05)^2)
exp(-.5 * (2 / 8.05)^2)
exp(-.5 * (3 / 8.05)^2)
exp(-.5 * (4 / 8.05)^2)
exp(-.5 * (5 / 8.05)^2)
exp(-.5 * (6 / 8.05)^2)
exp(-.5 * (7 / 8.05)^2)
exp(-.5 * (8 / 8.05)^2)
exp(-.5 * (9 / 8.05)^2)
exp(-.5 * (10 / 8.05)^2)
exp(-.5 * (11 / 8.05)^2)
exp(-.5 * (12 / 8.05)^2)
# gs2 median
exp(-.5 * (1 / 2.95)^2)
exp(-.5 * (2 / 2.95)^2)
exp(-.5 * (3 / 2.95)^2)
exp(-.5 * (4 / 2.95)^2)
exp(-.5 * (5 / 2.95)^2)
exp(-.5 * (6 / 2.95)^2)
exp(-.5 * (7 / 2.95)^2)
exp(-.5 * (8 / 2.95)^2)
exp(-.5 * (9 / 2.95)^2)
exp(-.5 * (10 / 2.95)^2)
exp(-.5 * (11 / 2.95)^2)
exp(-.5 * (12 / 2.95)^2)
# gs3 median
exp(-.5 * (1 / 3.09)^2)
exp(-.5 * (2 / 3.09)^2)
exp(-.5 * (3 / 3.09)^2)
exp(-.5 * (4 / 3.09)^2)
exp(-.5 * (5 / 3.09)^2)
exp(-.5 * (6 / 3.09)^2)
exp(-.5 * (7 / 3.09)^2)
exp(-.5 * (8 / 3.09)^2)
exp(-.5 * (9 / 3.09)^2)
exp(-.5 * (10 / 3.09)^2)
exp(-.5 * (11 / 3.09)^2)
exp(-.5 * (12 / 3.09)^2)


start_index <- 1
stop_index <- 8
csp_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 8
gs1_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 8
gs2_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 8
gs3_indices <- c(start_index:stop_index)
start_index <- start_index + 8
stop_index <- stop_index + 12
csp_indices <- c(csp_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs1_indices <- c(gs1_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs2_indices <- c(gs2_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs3_indices <- c(gs3_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
csp_indices <- c(csp_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs1_indices <- c(gs1_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs2_indices <- c(gs2_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs3_indices <- c(gs3_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
csp_indices <- c(csp_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs1_indices <- c(gs1_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs2_indices <- c(gs2_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 12
gs3_indices <- c(gs3_indices, start_index:stop_index)
start_index <- start_index + 12
stop_index <- stop_index + 15
shock_indices <- c(start_index:stop_index)


model026_draws_betas_long <- model026_draws %>% # your tibble
  pivot_longer(
    cols = matches("^betas\\[\\d+\\]$"), # only the betas[*] columns
    names_to = "beta_index",
    values_to = "beta_value",
    names_pattern = "betas\\[(\\d+)\\]" # capture the digits inside []
  ) %>%
  mutate(beta_index = as.integer(beta_index))

model027_draws_betas_long <- model027_draws %>% # your tibble
  pivot_longer(
    cols = matches("^betas\\[\\d+\\]$"), # only the betas[*] columns
    names_to = "beta_index",
    values_to = "beta_value",
    names_pattern = "betas\\[(\\d+)\\]" # capture the digits inside []
  ) %>%
  mutate(beta_index = as.integer(beta_index))

model028_draws_betas_long <- model028_draws %>% # your tibble
  pivot_longer(
    cols = matches("^betas\\[\\d+\\]$"), # only the betas[*] columns
    names_to = "beta_index",
    values_to = "beta_value",
    names_pattern = "betas\\[(\\d+)\\]" # capture the digits inside []
  ) %>%
  mutate(beta_index = as.integer(beta_index))

beta_key <- data.frame(
  beta_index = 1:176,
  trial_per_cue = c(
    rep(1:8, 4),
    rep(9:20, 4),
    rep(21:32, 4),
    rep(33:44, 4)
  )
)

model026_draws_betas_long <- model026_draws_betas_long %>%
  mutate(
    cue = case_when(
      beta_index %in% csp_indices ~ "csp",
      beta_index %in% gs1_indices ~ "gs1",
      beta_index %in% gs2_indices ~ "gs2",
      beta_index %in% gs3_indices ~ "gs3",
      beta_index %in% shock_indices ~ "shock"
    )
  ) %>%
  mutate(
    block = case_when(
      beta_index %in% c(1:32) ~ "habituation",
      beta_index %in% c(33:80) ~ "acquisition #1",
      beta_index %in% c(81:128) ~ "acquisition #2",
      beta_index %in% c(129:176) ~ "extinction"
    )
  )

model026_draws_betas_long <- merge(
  x = model026_draws_betas_long,
  y = beta_key,
  by.x = "beta_index",
  by.y = "beta_index",
  all.x = T
)

model028_draws_betas_long <- model028_draws_betas_long %>%
  mutate(
    cue = case_when(
      beta_index %in% csp_indices ~ "csp",
      beta_index %in% gs1_indices ~ "gs1",
      beta_index %in% gs2_indices ~ "gs2",
      beta_index %in% gs3_indices ~ "gs3",
      beta_index %in% shock_indices ~ "shock"
    )
  ) %>%
  mutate(
    block = case_when(
      beta_index %in% c(1:32) ~ "habituation",
      beta_index %in% c(33:80) ~ "acquisition #1",
      beta_index %in% c(81:128) ~ "acquisition #2",
      beta_index %in% c(129:176) ~ "extinction"
    )
  )

model028_draws_betas_long <- merge(
  x = model028_draws_betas_long,
  y = beta_key,
  by.x = "beta_index",
  by.y = "beta_index",
  all.x = T
)


model026_draws_betas_long %>%
  filter(cue != "shock") %>%
  group_by(
    trial_per_cue,
    cue,
    beta_index
  ) %>%
  summarise(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, .3),
    lower_97_5 = quantile(beta_value, .7)
    # lower_2_5 = quantile(beta_value, .5+.341),
    # lower_97_5 = quantile(beta_value, .5-.341)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Left Anterior Insula Over Trials (N =24)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model028_draws_betas_long %>%
  filter(cue != "shock") %>%
  group_by(
    trial_per_cue,
    cue,
    beta_index
  ) %>%
  summarise(
    median_posterior = median(beta_value),
    lower_2_5 = quantile(beta_value, .3),
    lower_97_5 = quantile(beta_value, .7)
    # lower_2_5 = quantile(beta_value, .5+.341),
    # lower_97_5 = quantile(beta_value, .5-.341)
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 8) +
  geom_vline(xintercept = 8 + 12) +
  geom_vline(xintercept = 8 + 12 + 12) +
  geom_ribbon(
    aes(
      x = trial_per_cue,
      y = median_posterior,
      ymin = lower_2_5,
      ymax = lower_97_5,
      fill = cue
    ),
    linewidth = .5,
    color = NA,
    alpha = .5
  ) +
  scale_fill_manual(values = cue_color) +
  ggtitle("Change in Left Anterior Insula Over Trials (N =24)") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))


model026_draws_betas_long %>%
  filter(beta_index %in% csp_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_density(aes(x = beta_value))

model026_draws_betas_long %>%
  filter(beta_index %in% csp_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = beta_index)) +
  coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full CS+ Posteriors (N=24") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model026_draws_betas_long %>%
  filter(beta_index %in% gs1_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = beta_index)) +
  # coord_cartesian(xlim = c(-.2, .4)) +
  theme_bw() +
  ggtitle("Full GS1 Posteriors (N=24") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 20))

model026_draws_betas_long %>%
  filter(beta_index %in% gs2_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = beta_index)) +
  theme_bw()

model026_draws_betas_long %>%
  filter(beta_index %in% gs3_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = beta_index)) +
  theme_bw()


model026_draws_betas_long %>%
  filter(beta_index %in% shock_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = beta_index)) +
  theme_bw() +
  ggtitle("Shock Posteriors (N=24)") +
  theme(text = element_text(family = "Arial", size = 20))

model028_draws_betas_long %>%
  filter(beta_index %in% shock_indices) %>%
  mutate(beta_index = factor(beta_index, levels = unique(beta_index))) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density_ridges(aes(x = beta_value, y = beta_index)) +
  theme_bw() +
  ggtitle("Shock Posteriors (N=24)") +
  theme(text = element_text(family = "Arial", size = 20))


## work that goes into tech report manuscript

# left anterior insula full multilevel priors
model018_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model018_chain_8872315_1.csv"),
    paste0(data_dir, "/model018_chain_8872315_2.csv"),
    paste0(data_dir, "/model018_chain_8872315_3.csv"),
    paste0(data_dir, "/model018_chain_8872315_4.csv"),
    paste0(data_dir, "/model018_chain_8872315_5.csv")
  )
)

model018_fit_meta_data <- model018_fit$metadata()

model018_fit_relevant_parameters <- model018_fit_meta_data$model_params[
  !str_detect(
    model018_fit_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta|shard_ind"
  )
]

model018_summary <- model018_fit$summary(
  variables = model018_fit_relevant_parameters
)

model018_draws <- model018_fit$draws(
  variables = model018_fit_relevant_parameters,
  format = "df"
)

model018_loo <- model018_fit$loo()

model018_no_ml_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model018_no_ML_chain_8955138_1.csv"),
    paste0(data_dir, "/model018_no_ML_chain_8955138_2.csv"),
    paste0(data_dir, "/model018_no_ML_chain_8955138_3.csv"),
    paste0(data_dir, "/model018_no_ML_chain_8955138_4.csv"),
    paste0(data_dir, "/model018_no_ML_chain_8955138_5.csv")
  )
)

model018_no_ml_fit_meta_data <- model018_no_ml_fit$metadata()

model018_no_ml_fit_relevant_parameters <- model018_no_ml_fit_meta_data$model_params[
  !str_detect(
    model018_no_ml_fit_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta|shard_ind"
  )
]

model018_no_ml_summary <- model018_no_ml_fit$summary(
  variables = model018_no_ml_fit_relevant_parameters
)

model018_no_ml_draws <- model018_no_ml_fit$draws(
  variables = model018_no_ml_fit_relevant_parameters,
  format = "df"
)

model018_no_ml_loo <- model018_no_ml_fit$loo()

loo::loo_compare(
  model018_loo,
  model018_no_ml_loo
)

loo::loo_model_weights(list(
  model018_loo,
  model018_no_ml_loo
))

bad_ks <- model018_no_ml_loo$pointwise[, 5] > .7

model018_no_ml_loo$pointwise[bad_ks, ] %>% sum()

model018_loo$pointwise[bad_ks, ] %>% sum()

fmri_stan_list$par[bad_ks]


model018_mu_beta_df <- model018_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^mu_betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "ROI", "condition"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
    values_to = "mu_beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(condition %in% c(1:13)) %>%
  mutate(condition = factor(condition, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      condition %in% c(1:4) ~ "habituation",
      condition %in% c(5:8, 13) ~ "acquisition",
      condition %in% c(9:12) ~ "extinction"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "habituation",
        "acquisition",
        "extinction"
      )
    )
  ) %>%
  mutate(
    cue = factor(
      case_when(
        condition %in% c(1, 5, 9) ~ "CS+",
        condition %in% c(2, 6, 10) ~ "GS1",
        condition %in% c(3, 7, 11) ~ "GS2",
        condition %in% c(4, 8, 12) ~ "GS3",
        condition %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  )

# there are no mu beatas for this model
# model018_no_ml_mu_beta_df <- model018_no_ml_draws %>%
#   # 1. select the .draw column plus all the mu_betas[...] cols
#   dplyr::select(
#     .draw,
#     tidyselect::matches("^mu_betas\\[")
#   ) %>%
#   # 2. pivot them longer, extracting the numbers inside the brackets
#   tidyr::pivot_longer(
#     cols = -.draw,
#     names_to = c("parameter", "ROI", "condition"),
#     # regex:
#     #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
#     #   \\[       → literal “[”
#     #   (\\d+)    → one or more digits  (we call that condition)
#     #   ,         → literal comma
#     #   (\\d+)    → one or more digits  (we call that ROI)
#     #   \\]       → literal “]”
#     names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
#     values_to = "mu_beta"
#   ) %>%
#   # 3. drop the now-redundant “parameter” column if you like
#   select(-parameter) %>%
#   filter(condition %in% c(1:13)) %>%
#   mutate(condition = factor(condition, levels = 1:13)) %>%
#   mutate(
#     phase = case_when(
#       condition %in% c(1:4) ~ "habituation",
#       condition %in% c(5:8, 13) ~ "acquisition",
#       condition %in% c(9:12) ~ "extinction"
#     )
#   ) %>%
#   mutate(
#     phase = factor(
#       phase,
#       levels = c(
#         "habituation",
#         "acquisition",
#         "extinction"
#       )
#     )
#   ) %>%
#   mutate(
#     cue = factor(
#       case_when(
#         condition %in% c(1, 5, 9) ~ "CS+",
#         condition %in% c(2, 6, 10) ~ "GS1",
#         condition %in% c(3, 7, 11) ~ "GS2",
#         condition %in% c(4, 8, 12) ~ "GS3",
#         condition %in% c(13) ~ "shock",
#       ),
#       levels = c("CS+", "GS1", "GS2", "GS3", "shock")
#     )
#   )

(model018_mu_beta_df %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = mu_beta, color = cue)) +
  scale_color_manual(values = cue_color) +
  facet_wrap(~phase, ncol = 1) +
  coord_cartesian(xlim = c(-.175, .35)) +
  ggtitle("Motion Multilevel Prior") +
  theme_bw())

model018_beta_df <- model018_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "par", "ROI", "coef"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+),(\\d+)\\]",
    values_to = "beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(coef %in% c(1:13)) %>%
  mutate(coef = factor(coef, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      coef %in% c(1:4) ~ "Habituation",
      coef %in% c(5:8, 13) ~ "Acquisition",
      coef %in% c(9:12) ~ "Extinction"
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
  ) %>%
  mutate(par = factor(par)) %>%
  mutate(
    cue = factor(
      case_when(
        coef %in% c(1, 5, 9) ~ "CS+",
        coef %in% c(2, 6, 10) ~ "GS1",
        coef %in% c(3, 7, 11) ~ "GS2",
        coef %in% c(4, 8, 12) ~ "GS3",
        coef %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  )

model018_no_ml_beta_df <- model018_no_ml_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "par", "ROI", "coef"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+),(\\d+)\\]",
    values_to = "beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(coef %in% c(1:13)) %>%
  mutate(coef = factor(coef, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      coef %in% c(1:4) ~ "Habituation",
      coef %in% c(5:8, 13) ~ "Acquisition",
      coef %in% c(9:12) ~ "Extinction"
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
  ) %>%
  mutate(par = factor(par)) %>%
  mutate(
    cue = factor(
      case_when(
        coef %in% c(1, 5, 9) ~ "CS+",
        coef %in% c(2, 6, 10) ~ "GS1",
        coef %in% c(3, 7, 11) ~ "GS2",
        coef %in% c(4, 8, 12) ~ "GS3",
        coef %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  )

# figure settings
text_font <- "Arial"
text_size <- 15
axis_line_thickness <- 1
participant_line_thickness <- .5
participant_line_alpha <- .5
xaxis_limits <- c(-.8, 1.5)
# yaxis_limits <- c(.0575, 1)
yaxis_limits <- c(0, 15)
vertical_line_thickness <- 1.25


fmri_multilevel_participant_posteriors_fig <-
  ggplot() +
  stat_density(
    data = model018_beta_df,
    aes(
      x = beta,
      # y = after_stat(scaled),
      color = cue,
      group = interaction(par, cue)
    ),
    geom = "line",
    # linetype = "dashed",
    position = "identity",
    alpha = participant_line_alpha,
    linewidth = participant_line_thickness,
  ) +
  # stat_density(
  #   data = model018_mu_beta_df,
  #   aes(x = mu_beta, y = after_stat(scaled), group = cue),
  #   geom = "line",
  #   position = "identity",
  #   alpha = 1,
  #   linewidth = 2,
  # ) +
  # stat_density(
  #   data = model018_mu_beta_df,
  #   aes(x = mu_beta, y = after_stat(scaled), color = cue),
  #   geom = "line",
  #   position = "identity",
  #   alpha = 1,
  #   linewidth = 1.5,
  # ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = vertical_line_thickness
  ) +
  scale_color_manual(values = cue_color) +
  scale_y_continuous(name = "Density", breaks = c(5, 10, 15)) +
  scale_x_continuous(
    name = "fMRI BOLD Percent Change",
    breaks = c(-.5, 0, 0.5, 1)
  ) +
  facet_wrap(~phase, ncol = 1) +
  coord_cartesian(
    xlim = xaxis_limits,
    expand = F,
    ylim = yaxis_limits
  ) +
  ggtitle("Multilevel Priors:\nParticipant Cue Effects by Phase") +
  theme_bw() +
  theme(
    text = element_text(family = text_font, size = text_size, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    panel.spacing = unit(0.15, "lines"),
    legend.position = "none",
    axis.line = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold"),
    # panel.border = element_rect(
    #   color = "black",
    #   fill = NA,
    #   linewidth = axis_line_thickness
    # ), # Add a border
    # axis.line = element_line(linewidth = axis_line_thickness,
    #                          lineend = "square"),
    axis.ticks = element_blank(),
    axis.text = element_text(face = "bold")
  )

fmri_multilevel_participant_posteriors_fig


ggsave(
  filename = "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/fMRI_multilevel_post.png",
  dpi = 300,
  units = "in",
  height = 1.5,
  width = 1.5,
  scale = 3.5
)

fmri_nonmultilevel_participant_posteriors_fig <-
  ggplot() +
  stat_density(
    data = model018_no_ml_beta_df,
    aes(
      x = beta,
      # y = after_stat(scaled),
      color = cue,
      group = interaction(par, cue)
    ),
    geom = "line",
    # linetype = "dashed",
    position = "identity",
    alpha = participant_line_alpha,
    linewidth = participant_line_thickness,
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = vertical_line_thickness
  ) +
  scale_color_manual(values = cue_color) +
  scale_y_continuous(name = "Density", breaks = c(5, 10, 15)) +
  scale_x_continuous(
    name = "fMRI BOLD Percent Change",
    breaks = c(-.5, 0, .5, 1)
  ) +
  facet_wrap(~phase, ncol = 1) +
  coord_cartesian(
    xlim = xaxis_limits,
    expand = F,
    ylim = yaxis_limits
  ) +
  ggtitle("Uninformative Priors:\nParticipant Cue Effects by Phase") +
  theme_bw() +
  theme(
    text = element_text(family = text_font, size = text_size, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    panel.spacing = unit(0.15, "lines"),
    legend.position = "none",
    axis.line = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold"),
    # panel.border = element_rect(
    #   color = "black",
    #   fill = NA,
    #   linewidth = axis_line_thickness
    # ), # Add a border
    # axis.line = element_line(linewidth = axis_line_thickness,
    #                          lineend = "square"),
    axis.ticks = element_blank(),
    axis.text = element_text(face = "bold")
  )

fmri_nonmultilevel_participant_posteriors_fig

ggsave(
  filename = "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/fMRI_nonmultilevel_post.png",
  dpi = 300,
  units = "in",
  height = 1.5,
  width = 1.5,
  scale = 3.5
)


# timeseries plot

par_index <- 1
posterior_indices <- sample(1:5000, 100)
beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')

model018_beta_draws_1 <- model018_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

model018_no_ml_beta_draws_1 <- model018_no_ml_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

uncensored_1 <- fmri_stan_list$usable_bold_indices_one_is_true[par_index, ] == 1

current_par_design_matrix_1 <- fmri_stan_list$design_array[
  par_index,
  ,
  # beta
]

mod18_posterior_mu_1 <- (current_par_design_matrix_1 %*%
  t(model018_beta_draws_1[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)) %>%
  pivot_longer(starts_with("V")) %>%
  mutate(
    value = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_1],
      value,
      NA
    )
  )

mod18_posterior_mu_CI_1 <- mod18_posterior_mu_1 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025, na.rm = T),
    ub_value = quantile(value, probs = .975, na.rm = T)
  )

mod18_no_ml_posterior_mu_1 <- (current_par_design_matrix_1 %*%
  t(model018_no_ml_beta_draws_1[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)) %>%
  pivot_longer(starts_with("V")) %>%
  mutate(
    value = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_1],
      value,
      NA
    )
  )

mod18_no_ml_posterior_mu_CI_1 <- mod18_no_ml_posterior_mu_1 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025, na.rm = T),
    ub_value = quantile(value, probs = .975, na.rm = T)
  )


bold_plot_df_1 <- data.frame(
  time_seconds = seq(0, 1069 * 2, by = 2),
  bold = fmri_stan_list$bold[par_index, ]
) %>%
  mutate(
    bold = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_1],
      bold,
      NA
    )
  )

fmri_plot_1 <- ggplot() +
  geom_line(
    data = bold_plot_df_1,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .5,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_no_ml_posterior_mu_CI_1,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "blue",
    color = NA,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_posterior_mu_CI_1,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "red",
    color = NA,
    alpha = .5
  ) +
  scale_y_continuous(name = "fMRI BOLD Percent Change") +
  scale_x_continuous(
    name = "Time (seconds)",
    breaks = seq(0, 2000, by = 250)
  ) +
  coord_cartesian(
    xlim = c(-50, 2190),
    # ylim = c(-1.25, 1.25),
    expand = F
  ) +
  ggtitle(paste0("Participant ", par_index)) +
  theme_bw() +
  theme(
    text = element_text(family = text_font, size = text_size, color = "black"),
    # strip.background = element_blank(),
    # strip.text = element_text(size = text_size),
    panel.spacing = unit(0.15, "lines"),
    legend.position = "none",
    axis.line = element_blank(),
    # panel.border = element_rect(
    #   color = "black",
    #   fill = NA,
    #   linewidth = axis_line_thickness
    # ), # Add a border
    # axis.line = element_line(linewidth = axis_line_thickness,
    #                          lineend = "square"),
    axis.ticks = element_blank(),
    axis.text = element_text(face = "bold"),
    axis.title = element_blank()
  )

fmri_plot_1

par_index <- 5
posterior_indices <- sample(1:5000, 100)
beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')

model018_beta_draws_2 <- model018_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

model018_no_ml_beta_draws_2 <- model018_no_ml_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

uncensored_2 <- fmri_stan_list$usable_bold_indices_one_is_true[par_index, ] == 1

current_par_design_matrix_2 <- fmri_stan_list$design_array[
  par_index,
  ,
  # beta
]

mod18_posterior_mu_2 <- (current_par_design_matrix_2 %*%
  t(model018_beta_draws_2[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)) %>%
  pivot_longer(starts_with("V")) %>%
  mutate(
    value = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_2],
      value,
      NA
    )
  )

mod18_posterior_mu_CI_2 <- mod18_posterior_mu_2 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025, na.rm = T),
    ub_value = quantile(value, probs = .975, na.rm = T)
  )

mod18_no_ml_posterior_mu_2 <- (current_par_design_matrix_2 %*%
  t(model018_no_ml_beta_draws_2[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)) %>%
  pivot_longer(starts_with("V")) %>%
  mutate(
    value = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_2],
      value,
      NA
    )
  )

mod18_no_ml_posterior_mu_CI_2 <- mod18_no_ml_posterior_mu_2 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025, na.rm = T),
    ub_value = quantile(value, probs = .975, na.rm = T)
  )


bold_plot_df_2 <- data.frame(
  time_seconds = seq(0, 1069 * 2, by = 2),
  bold = fmri_stan_list$bold[par_index, ]
) %>%
  mutate(
    bold = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_2],
      bold,
      NA
    )
  )

fmri_plot_2 <- ggplot() +
  geom_line(
    data = bold_plot_df_2,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .5,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_no_ml_posterior_mu_CI_2,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "blue",
    color = NA,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_posterior_mu_CI_2,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "red",
    color = NA,
    alpha = .5
  ) +
  scale_y_continuous(name = "fMRI BOLD Percent Change") +
  scale_x_continuous(
    name = "Time (seconds)",
    breaks = seq(0, 2000, by = 250)
  ) +
  coord_cartesian(
    xlim = c(-50, 2190),
    # ylim = c(-1.25, 1.25),
    expand = F
  ) +
  ggtitle(paste0("Participant ", par_index)) +
  theme_bw() +
  theme(
    text = element_text(family = text_font, size = text_size, color = "black"),
    # strip.background = element_blank(),
    # strip.text = element_text(size = text_size),
    panel.spacing = unit(0.15, "lines"),
    legend.position = "none",
    axis.line = element_blank(),
    # panel.border = element_rect(
    #   color = "black",
    #   fill = NA,
    #   linewidth = axis_line_thickness
    # ), # Add a border
    # axis.line = element_line(linewidth = axis_line_thickness,
    #                          lineend = "square"),
    axis.ticks = element_blank(),
    axis.text = element_text(face = "bold"),
    axis.title.x = element_blank()
  )

fmri_plot_2


par_index <- 10
posterior_indices <- sample(1:5000, 100)
beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')

model018_beta_draws_3 <- model018_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

model018_no_ml_beta_draws_3 <- model018_no_ml_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

uncensored_3 <- fmri_stan_list$usable_bold_indices_one_is_true[par_index, ] == 1

current_par_design_matrix_3 <- fmri_stan_list$design_array[
  par_index,
  ,
  # beta
]

mod18_posterior_mu_3 <- (current_par_design_matrix_3 %*%
  t(model018_beta_draws_3[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)) %>%
  pivot_longer(starts_with("V")) %>%
  mutate(
    value = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_3],
      value,
      NA
    )
  )

mod18_posterior_mu_CI_3 <- mod18_posterior_mu_3 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025, na.rm = T),
    ub_value = quantile(value, probs = .975, na.rm = T)
  )

mod18_no_ml_posterior_mu_3 <- (current_par_design_matrix_3 %*%
  t(model018_no_ml_beta_draws_3[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)) %>%
  pivot_longer(starts_with("V")) %>%
  mutate(
    value = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_3],
      value,
      NA
    )
  )

mod18_no_ml_posterior_mu_CI_3 <- mod18_no_ml_posterior_mu_3 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025, na.rm = T),
    ub_value = quantile(value, probs = .975, na.rm = T)
  )


bold_plot_df_3 <- data.frame(
  time_seconds = seq(0, 1069 * 2, by = 2),
  bold = fmri_stan_list$bold[par_index, ]
) %>%
  mutate(
    bold = if_else(
      time_seconds %in% seq(0, 1069 * 2, by = 2)[uncensored_3],
      bold,
      NA
    )
  )

fmri_plot_3 <- ggplot() +
  geom_line(
    data = bold_plot_df_3,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .5,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_no_ml_posterior_mu_CI_3,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "blue",
    color = NA,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_posterior_mu_CI_3,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "red",
    color = NA,
    alpha = .5
  ) +
  scale_y_continuous(name = "fMRI BOLD Percent Change") +
  scale_x_continuous(
    name = "Time (seconds)",
    breaks = seq(0, 2000, by = 250)
  ) +
  coord_cartesian(
    xlim = c(-50, 2190),
    # ylim = c(-1.25, 1.25),
    expand = F
  ) +
  ggtitle(paste0("Participant ", par_index)) +
  theme_bw() +
  theme(
    text = element_text(family = text_font, size = text_size, color = "black"),
    # strip.background = element_blank(),
    # strip.text = element_text(size = text_size),
    panel.spacing = unit(0.15, "lines"),
    legend.position = "none",
    axis.line = element_blank(),
    # panel.border = element_rect(
    #   color = "black",
    #   fill = NA,
    #   linewidth = axis_line_thickness
    # ), # Add a border
    # axis.line = element_line(linewidth = axis_line_thickness,
    #                          lineend = "square"),
    axis.ticks = element_blank(),
    axis.text = element_text(face = "bold"),
    axis.title.y = element_blank()
  )

fmri_plot_3


((fmri_plot_1) /
  (fmri_plot_2) /
  (fmri_plot_3)) +
  patchwork::plot_annotation(
    title = "Multilevel (Red) vs Non-Multilevel (Blue) Regression Mean",
    theme = theme(
      plot.title = element_text(
        hjust = .5,
        face = "bold",
        family = "Arial",
        size = 17.5,
      )
    )
  )

ggsave(
  filename = "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/fMRI_multilevel_2.png",
  dpi = 300,
  units = "in",
  height = 1.5 * 3,
  width = 3.75,
  scale = 2
)

layout_grid <- c(
  '
AB
AB
AB
AC
AC
AC
'
)

layout_grid <- c(
  '
AD
AD
BD
BE
CE
CE
'
)

fmri_timeseries_plot +
  fmri_multilevel_participant_posteriors_fig +
  fmri_nonmultilevel_participant_posteriors_fig +
  plot_annotation(
    title = "Effects of Multilevel Structure on Per Participant fMRI Deconvolution",
    theme = theme(
      plot.title = element_text(
        family = "Arial",
        size = 27,
        color = "black",
        hjust = 0.5,
        face = "bold"
      )
    )
  ) +
  plot_layout(design = layout_grid)


# old below
# left anterior insula full multilevel priors
model018_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model018_chain_66251546_1.csv"),
    paste0(data_dir, "/model018_chain_66251546_2.csv"),
    paste0(data_dir, "/model018_chain_66251546_3.csv"),
    paste0(data_dir, "/model018_chain_66251546_4.csv"),
    paste0(data_dir, "/model018_chain_66251546_5.csv")
  )
)

model018_fit_meta_data <- model018_fit$metadata()

model018_fit_relevant_parameters <- model018_fit_meta_data$model_params[
  !str_detect(
    model018_fit_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta|shard_ind"
  )
]

model018_summary <- model018_fit$summary(
  variables = model018_fit_relevant_parameters
)

model018_draws <- model018_fit$draws(
  variables = model018_fit_relevant_parameters,
  format = "df"
)

model018_loo <- model018_fit$loo()


model018_mu_beta_df <- model018_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^mu_betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "ROI", "condition"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
    values_to = "mu_beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(condition %in% c(1:13)) %>%
  mutate(condition = factor(condition, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      condition %in% c(1:4) ~ "habituation",
      condition %in% c(5:8, 13) ~ "acquisition",
      condition %in% c(9:12) ~ "extinction"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "habituation",
        "acquisition",
        "extinction"
      )
    )
  ) %>%
  mutate(
    cue = factor(
      case_when(
        condition %in% c(1, 5, 9) ~ "CS+",
        condition %in% c(2, 6, 10) ~ "GS1",
        condition %in% c(3, 7, 11) ~ "GS2",
        condition %in% c(4, 8, 12) ~ "GS3",
        condition %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  )


# anterior insula no multilevel on motion
model022_fit <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model022_chain_66382183_1.csv"),
    paste0(data_dir, "/model022_chain_66382183_2.csv"),
    paste0(data_dir, "/model022_chain_66382183_3.csv"),
    paste0(data_dir, "/model022_chain_66382183_4.csv"),
    paste0(data_dir, "/model022_chain_66382183_5.csv")
  )
)

model022_fit_meta_data <- model022_fit$metadata()

model022_fit_relevant_parameters <- model022_fit_meta_data$model_params[
  !str_detect(
    model022_fit_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta|shard_ind"
  )
]

model022_summary <- model022_fit$summary(
  variables = model022_fit_relevant_parameters
)

model022_loo <- model022_fit$loo()


model022_draws <- model022_fit$draws(
  variables = model022_fit_relevant_parameters,
  format = "df"
)

loo::loo_compare(
  model018_loo,
  model022_loo
)

loo::loo_model_weights(list(
  model018_loo,
  model022_loo
))

#
# save(
#   fmri_stan_list,
#   file = '/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/ant_insula_n15_stan_list.RData'
# )

model018_mu_beta_df <- model018_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^mu_betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "ROI", "condition"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
    values_to = "mu_beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(condition %in% c(1:13)) %>%
  mutate(condition = factor(condition, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      condition %in% c(1:4) ~ "habituation",
      condition %in% c(5:8, 13) ~ "acquisition",
      condition %in% c(9:12) ~ "extinction"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "habituation",
        "acquisition",
        "extinction"
      )
    )
  ) %>%
  mutate(
    cue = factor(
      case_when(
        condition %in% c(1, 5, 9) ~ "CS+",
        condition %in% c(2, 6, 10) ~ "GS1",
        condition %in% c(3, 7, 11) ~ "GS2",
        condition %in% c(4, 8, 12) ~ "GS3",
        condition %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  )

model022_mu_beta_df <- model022_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^mu_betas_stim\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "ROI", "condition"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
    values_to = "mu_beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(condition %in% c(1:13)) %>%
  mutate(condition = factor(condition, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      condition %in% c(1:4) ~ "habituation",
      condition %in% c(5:8, 13) ~ "acquisition",
      condition %in% c(9:12) ~ "extinction"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "habituation",
        "acquisition",
        "extinction"
      )
    )
  ) %>%
  mutate(
    cue = factor(
      case_when(
        condition %in% c(1, 5, 9) ~ "CS+",
        condition %in% c(2, 6, 10) ~ "GS1",
        condition %in% c(3, 7, 11) ~ "GS2",
        condition %in% c(4, 8, 12) ~ "GS3",
        condition %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  )

(model018_mu_beta_df %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = mu_beta, color = cue)) +
  scale_color_manual(values = cue_color) +
  facet_wrap(~phase, ncol = 1) +
  coord_cartesian(xlim = c(-.175, .35)) +
  ggtitle("Motion Multilevel Prior") +
  theme_bw()) /

  (model022_mu_beta_df %>%
    ggplot() +
    geom_vline(xintercept = 0) +
    geom_density(aes(x = mu_beta, color = cue)) +
    scale_color_manual(values = cue_color) +
    facet_wrap(~phase, ncol = 1) +
    coord_cartesian(xlim = c(-.175, .35)) +
    ggtitle("No Motion Multilevel Prior") +
    theme_bw())


par_index <- 1
posterior_indices <- sample(1:2000, 100)
beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')

model022_beta_draws_1 <- model022_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

model018_beta_draws_1 <- model018_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()


uncensored_1 <- fmri_stan_list$usable_bold_indices_one_is_true[par_index, ] == 1

current_par_design_matrix_1 <- fmri_stan_list$design_array[
  par_index,
  uncensored_1,
  # beta
]

mod18_posterior_mu_1 <- (current_par_design_matrix_1 %*%
  t(model018_beta_draws_1[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored_1]) %>%
  pivot_longer(starts_with("V"))

mod18_posterior_mu_CI_1 <- mod18_posterior_mu_1 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025),
    ub_value = quantile(value, probs = .975)
  )

mod22_posterior_mu_1 <- (current_par_design_matrix_1 %*%
  t(model022_beta_draws_1[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored_1]) %>%
  pivot_longer(starts_with("V"))

mod22_posterior_mu_CI_1 <- mod22_posterior_mu_1 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025),
    ub_value = quantile(value, probs = .975)
  )


bold_plot_df_1 <- data.frame(
  time_seconds = seq(0, 1069 * 2, by = 2)[uncensored_1],
  bold = fmri_stan_list$bold[par_index, 1, uncensored_1]
)

fmri_plot_1 <- ggplot() +
  geom_line(
    data = bold_plot_df_1,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .5,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod22_posterior_mu_CI_1,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "red",
    color = NA,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_posterior_mu_CI_1,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "blue",
    color = NA,
    alpha = .5
  ) +
  scale_y_continuous(name = "Percent Change fMRI BOLD") +
  scale_x_continuous(
    name = "Time (seconds)",
    breaks = seq(0, 2000, by = 250)
  ) +
  coord_cartesian(
    xlim = c(-50, 2190),
    # ylim = c(-1.25, 1.25),
    expand = F
  ) +
  ggtitle("Participant 1") +
  theme_classic() +
  theme()

par_index <- 4
posterior_indices <- sample(1:2000, 100)
beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')

model022_beta_draws_2 <- model022_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

model018_beta_draws_2 <- model018_fit$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()


uncensored_2 <- fmri_stan_list$usable_bold_indices_one_is_true[par_index, ] == 1

current_par_design_matrix_2 <- fmri_stan_list$design_array[
  par_index,
  uncensored_2,
  # beta
]

mod18_posterior_mu_2 <- (current_par_design_matrix_2 %*%
  t(model018_beta_draws_2[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored_2]) %>%
  pivot_longer(starts_with("V"))

mod18_posterior_mu_CI_2 <- mod18_posterior_mu_2 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025),
    ub_value = quantile(value, probs = .975)
  )

mod22_posterior_mu_2 <- (current_par_design_matrix_2 %*%
  t(model022_beta_draws_2[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored_2]) %>%
  pivot_longer(starts_with("V"))

mod22_posterior_mu_CI_2 <- mod22_posterior_mu_2 %>%
  group_by(time_seconds) %>%
  reframe(
    avg_value = mean(value),
    lb_value = quantile(value, probs = .025),
    ub_value = quantile(value, probs = .975)
  )


bold_plot_df_2 <- data.frame(
  time_seconds = seq(0, 1069 * 2, by = 2)[uncensored_2],
  bold = fmri_stan_list$bold[par_index, 1, uncensored_2]
)

fmri_plot_2 <- ggplot() +
  geom_line(
    data = bold_plot_df_2,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .5,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod22_posterior_mu_CI_2,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "red",
    color = NA,
    alpha = .5
  ) +
  geom_ribbon(
    data = mod18_posterior_mu_CI_2,
    aes(
      x = time_seconds,
      y = avg_value,
      ymin = lb_value,
      ymax = ub_value
    ),
    linewidth = .5,
    fill = "blue",
    color = NA,
    alpha = .5
  ) +
  scale_y_continuous(name = "Percent Change fMRI BOLD") +
  scale_x_continuous(
    name = "Time (seconds)",
    breaks = seq(0, 2000, by = 250)
  ) +
  coord_cartesian(
    xlim = c(-50, 2190),
    # ylim = c(-1.25, 1.25),
    expand = F
  ) +
  ggtitle("Participant 2") +
  theme_classic() +
  theme()

((fmri_plot_1) /
  (fmri_plot_2)) +
  patchwork::plot_annotation(
    title = "Multilevel Priors (Blue) Lead to Better fMRI Deconvolution"
  )

ggsave(
  filename = "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/fMRI_multilevel.png",
  dpi = 300,
  units = "in",
  height = 1.5 * 2,
  width = 3.75,
  scale = 2
)

# geom_line(
#   data = mod22_posterior_mu_CI,
#   aes(x = time_seconds, y = avg_value),
#   color = "black",
#   linewidth = 1.1
# ) +
# geom_line(
#   data = mod22_posterior_mu_CI,
#   aes(x = time_seconds, y = avg_value),
#   color = "red",
#   linewidth = 1
# )

geom_errorbar(
  data = mod18_posterior_mu_CI,
  aes(x = time_seconds, y = avg_value, ymin = lb_value, ymax = ub_value),
  color = "blue",
  alpha = .4
) +
  # theme_bw()

  ggplot() +
  geom_line(
    data = bold_plot_df,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .2
  ) +
  # geom_line(
  #   data = mod22_posterior_mu[mod22_posterior_mu$used, ],
  geom_line(
    data = mod22_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .2,
    alpha = .2,
    color = "blue"
  ) +
  geom_line(
    data = mod18_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .2,
    alpha = .2,
    color = "red"
  ) +
  theme_classic()


## old below

model025_fit_mpi <- as_cmdstan_fit(
  files = c(
    paste0(data_dir, "/model025_chain_3729566_1.csv"),
    paste0(data_dir, "/model025_chain_3729566_2.csv"),
    paste0(data_dir, "/model025_chain_3729566_3.csv"),
    paste0(data_dir, "/model025_chain_3729566_4.csv")
  )
)

model025_fit_mpi_meta_data <- model025_fit_mpi$metadata()

model025_fit_mpi_relevant_parameters <- model025_fit_mpi_meta_data$model_params[
  !str_detect(
    model025_fit_mpi_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta"
  )
]


model025_fit_mpi_summary <- model025_fit_mpi$summary(
  variables = model025_fit_mpi_relevant_parameters
)

model025_fit_draws <- model025_fit_mpi$draws(
  variables = model025_fit_mpi_relevant_parameters,
  format = "df"
)


model025_fit_loo <- model025_fit_mpi$loo()

model025_fit_draws$`rho_z_betas_z[5,1]` %>% tanh() %>% density() %>% plot()
model025_fit_draws$`rho_z_betas_z[5,1]` %>% tanh() %>% hist()


model025_fit_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^mu_betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "condition", "ROI"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
    values_to = "mu_beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(condition %in% c(1:13)) %>%
  mutate(condition = factor(condition, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      condition %in% c(1:4) ~ "habituation",
      condition %in% c(5:8, 13) ~ "acquisition",
      condition %in% c(9:12) ~ "extinction"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "habituation",
        "acquisition",
        "extinction"
      )
    )
  ) %>%
  mutate(
    cue = factor(
      case_when(
        condition %in% c(1, 5, 9) ~ "CS+",
        condition %in% c(2, 6, 10) ~ "GS1",
        condition %in% c(3, 7, 11) ~ "GS2",
        condition %in% c(4, 8, 12) ~ "GS3",
        condition %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  ) %>%
  mutate(
    roi_name = factor(
      case_when(
        ROI == 1 ~ "V1",
        ROI == 2 ~ "V2",
      ),
      levels = c("V1", "V2")
    )
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), linewidth = 1.5) +
  geom_density(aes(x = mu_beta, color = cue), linewidth = 1.5) +
  scale_color_manual(values = cue_color) +
  scale_x_continuous(name = "% BOLD Change") +
  facet_wrap(~ phase * roi_name, ncol = 2) +
  theme_bw() +
  theme(text = element_text(size = 20))


model024_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewfarkas/Downloads/model024_chain_3294503_1.csv",
    "/home/andrewfarkas/Downloads/model024_chain_3294503_2.csv",
    "/home/andrewfarkas/Downloads/model024_chain_3294503_4.csv",
    "/home/andrewfarkas/Downloads/model024_chain_3294503_5.csv"
  )
)

model024_fit_mpi_meta_data <- model024_fit_mpi$metadata()

model024_fit_mpi_relevant_parameters <- model024_fit_mpi_meta_data$model_params[
  !str_detect(
    model024_fit_mpi_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta"
  )
]


model024_fit_mpi_summary <- model024_fit_mpi$summary(
  variables = model024_fit_mpi_relevant_parameters
)


model024_fit_loo <- model024_fit_mpi$loo()

loo::loo_compare(
  model024_fit_loo,
  model025_fit_loo
)


model024_fit_draws <- model024_fit_mpi$draws(
  variables = model024_fit_mpi_relevant_parameters
)

# tau_betas should have had a lower bound, because it didn't it could reach a negative or postive value and reach the same results
# so some beta_zs are reversed, but the betas are the same, so most parameters are interpretable
bayesplot::mcmc_trace(model024_fit_draws[,, "betas_z[1,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas_z[2,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas_z[4,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas[4,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas_z[5,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas[5,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas_z[7,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas[7,1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "betas_z[4,3,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[1,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[1,2]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[5,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[6,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[7,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[8,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[9,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[10,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[12,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[13,1]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "mu_betas[5,2]"])
bayesplot::mcmc_trace(model024_fit_draws[,, "tau_betas[1,1]"])


model024_fit_draws <- model024_fit_mpi$draws(
  variables = model024_fit_mpi_relevant_parameters,
  format = "df"
)


# assume your tibble is called `draws_df`
model024_fit_draws %>%
  # 1. select the .draw column plus all the mu_betas[...] cols
  dplyr::select(
    .draw,
    tidyselect::matches("^mu_betas\\[")
  ) %>%
  # 2. pivot them longer, extracting the numbers inside the brackets
  tidyr::pivot_longer(
    cols = -.draw,
    names_to = c("parameter", "condition", "ROI"),
    # regex:
    #   ([^[]+)   → everything before the “[”  (we call that “parameter”)
    #   \\[       → literal “[”
    #   (\\d+)    → one or more digits  (we call that condition)
    #   ,         → literal comma
    #   (\\d+)    → one or more digits  (we call that ROI)
    #   \\]       → literal “]”
    names_pattern = "([^\\[]+)\\[(\\d+),(\\d+)\\]",
    values_to = "mu_beta"
  ) %>%
  # 3. drop the now-redundant “parameter” column if you like
  select(-parameter) %>%
  filter(condition %in% c(1:13)) %>%
  mutate(condition = factor(condition, levels = 1:13)) %>%
  mutate(
    phase = case_when(
      condition %in% c(1:4) ~ "habituation",
      condition %in% c(5:8, 13) ~ "acquisition",
      condition %in% c(9:12) ~ "extinction"
    )
  ) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "habituation",
        "acquisition",
        "extinction"
      )
    )
  ) %>%
  mutate(
    cue = factor(
      case_when(
        condition %in% c(1, 5, 9) ~ "CS+",
        condition %in% c(2, 6, 10) ~ "GS1",
        condition %in% c(3, 7, 11) ~ "GS2",
        condition %in% c(4, 8, 12) ~ "GS3",
        condition %in% c(13) ~ "shock",
      ),
      levels = c("CS+", "GS1", "GS2", "GS3", "shock")
    )
  ) %>%
  mutate(
    roi_name = factor(
      case_when(
        ROI == 1 ~ "V1",
        ROI == 2 ~ "V2",
      ),
      levels = c("V1", "V2")
    )
  ) %>%
  ggplot() +
  geom_vline(aes(xintercept = 0), linewidth = 1.5) +
  geom_density(aes(x = mu_beta, color = cue), linewidth = 1.5) +
  scale_color_manual(values = cue_color) +
  scale_x_continuous(name = "% BOLD Change") +
  facet_wrap(~ phase * roi_name, ncol = 2) +
  theme_bw() +
  theme(text = element_text(size = 20))


model018_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewfarkas/Downloads/model018_chain_66251546_1.csv",
    "/home/andrewfarkas/Downloads/model018_chain_66251546_2.csv",
    "/home/andrewfarkas/Downloads/model018_chain_66251546_3.csv",
    "/home/andrewfarkas/Downloads/model018_chain_66251546_4.csv"
    # "/home/andrewfarkas/Downloads/model018_chain_66251546_5.csv"
  )
)

model018_fit_mpi_meta_data <- model018_fit_mpi$metadata()

model018_fit_mpi_relevant_parameters <- model018_fit_mpi_meta_data$model_params[
  !str_detect(
    model018_fit_mpi_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta"
  )
]


model018_fit_mpi_summary <- model018_fit_mpi$summary(
  variables = model018_fit_mpi_relevant_parameters
)


model018_fit_loo <- model018_fit_mpi$loo()


model018_fit_betas_names <- (model018_fit_mpi_meta_data$model_params[
  str_detect(
    model018_fit_mpi_meta_data$model_params,
    '^betas\\['
  )
])

model018_fit_mu_betas_namecspacq <- (model018_fit_mpi_meta_data$model_params[
  str_detect(
    model018_fit_mpi_meta_data$model_params,
    'mu_betas\\[1,5]'
  )
])

model018_fit_mu_betas_shock <- (model018_fit_mpi_meta_data$model_params[
  str_detect(
    model018_fit_mpi_meta_data$model_params,
    'mu_betas\\[1,13]'
  )
])

model018_fit_mu_betas_csp_hab <- (model018_fit_mpi_meta_data$model_params[
  str_detect(
    model018_fit_mpi_meta_data$model_params,
    'mu_betas\\[1,1]'
  )
])

# model018_fit_beta_csp_acq <- model018_fit_betas_names[
#   str_detect(
#     model018_fit_betas_names,
#     '1,5]$'
#   )
# ]

model018_beta_draws <- model018_fit_mpi$draws(
  format = "df",
  variables = model018_fit_betas_names
)

model018_mu_beta_cspacq_draws <- model018_fit_mpi$draws(
  format = "df",
  variables = model018_fit_mu_betas_namecspacq
) %>%
  select(starts_with("mu_bet")) %>%
  pivot_longer(everything())

model018_mu_beta_csphab_draws <- model018_fit_mpi$draws(
  format = "df",
  variables = model018_fit_mu_betas_csp_hab
) %>%
  select(starts_with("mu_bet")) %>%
  pivot_longer(everything())

model018_mu_beta_shock_draws <- model018_fit_mpi$draws(
  format = "df",
  variables = model018_fit_mu_betas_shock
) %>%
  select(starts_with("mu_bet")) %>%
  pivot_longer(everything())

model018_beta_csp_acq_draws <- model018_beta_draws %>%
  select(ends_with('1,5]'))

model018_beta_shock_draws <- model018_beta_draws %>%
  select(ends_with('1,13]'))

model018_beta_csp_hab_draws <- model018_beta_draws %>%
  select(ends_with('1,1]'))

model018_beta_csp_acq_draws %>%
  pivot_longer(everything()) %>%
  mutate(name_factor = factor(name, levels = unique(name))) %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = value, color = name_factor)) +
  geom_density(
    data = model018_mu_beta_cspacq_draws,
    aes(x = value),
    linewidth = 2
  ) +
  scale_color_discrete(name = "Participant", labels = c(1:15)) +
  scale_x_continuous(
    name = "Percent Change BOLD",
    breaks = seq(-.4, .8, by = .2)
  ) +
  ggtitle("CS+ Acquisition Anterior Insula") +
  theme_classic() +
  theme(text = element_text(size = 25))


model018_beta_shock_draws %>%
  pivot_longer(everything()) %>%
  mutate(name_factor = factor(name, levels = unique(name))) %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = value, color = name_factor)) +
  geom_density(
    data = model018_mu_beta_shock_draws,
    aes(x = value),
    linewidth = 2
  ) +
  scale_color_discrete(name = "Participant", labels = c(1:15)) +
  scale_x_continuous(
    name = "Percent Change BOLD",
    breaks = seq(-.4, .8, by = .2)
  ) +
  ggtitle("Shock Acquisition Anterior Insula") +
  theme_classic() +
  theme(text = element_text(size = 25))


model018_beta_csp_hab_draws %>%
  pivot_longer(everything()) %>%
  mutate(name_factor = factor(name, levels = unique(name))) %>%
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = value, color = name_factor)) +
  geom_density(
    data = model018_mu_beta_csphab_draws,
    aes(x = value),
    linewidth = 2
  ) +
  scale_color_discrete(name = "Participant", labels = c(1:15)) +
  scale_x_continuous(
    name = "Percent Change BOLD",
    breaks = seq(-.4, .8, by = .2)
  ) +
  ggtitle("CS+ Habituation Anterior Insula") +
  theme_classic() +
  theme(text = element_text(size = 25))


model018_fit_no_mot <- as_cmdstan_fit(
  files = c(
    "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model018_chain_66662361_1.csv",
    "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model018_chain_66662361_2.csv" #,
    # "/home/andrewf/Downloads/model022_chain_66382183_3.csv"
  )
)

model018_fit_no_mot_meta_data <- model018_fit_no_mot$metadata()

model018_fit_no_mot_relevant_parameters <- model018_fit_no_mot_meta_data$model_params[
  !str_detect(
    model018_fit_no_mot_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta"
  )
]

model018_fit_no_mot_loo <- model018_fit_no_mot$loo()


model018_fit_no_mot_summary <- model018_fit_no_mot$summary(
  variables = model018_fit_no_mot_relevant_parameters
)

model022_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model022_chain_66382183_1.csv",
    "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model022_chain_66382183_2.csv" #,
    # "/home/andrewf/Downloads/model022_chain_66382183_3.csv"
  )
)

model022_fit_mpi_meta_data <- model022_fit_mpi$metadata()

model022_fit_mpi_relevant_parameters <- model022_fit_mpi_meta_data$model_params[
  !str_detect(
    model022_fit_mpi_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta"
  )
]


model022_fit_mpi_summary <- model022_fit_mpi$summary(
  variables = model022_fit_mpi_relevant_parameters
)


model021_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/model021_chain_66466225_1.csv",
    "/home/andrewf/Downloads/model021_chain_66466225_2.csv" #,
    # "/home/andrewf/Downloads/model021_chain_66466225_3.csv"
  )
)

model021_fit_mpi_meta_data <- model021_fit_mpi$metadata()

model021_fit_mpi_relevant_parameters <- model021_fit_mpi_meta_data$model_params[
  !str_detect(
    model021_fit_mpi_meta_data$model_params,
    "log_lik|mu_pred|amplitude|S|theta"
  )
]


model021_fit_mpi_summary <- model021_fit_mpi$summary(
  variables = model021_fit_mpi_relevant_parameters
)


# model018_fit_mpi <- as_cmdstan_fit(
#   files = c(
#     "/home/andrewf/Downloads/model018_chain_66251546_1.csv",
#     "/home/andrewf/Downloads/model018_chain_66251546_2.csv",
#     "/home/andrewf/Downloads/model018_chain_66251546_3.csv"
#   )
# )

model022_loo <- model022_fit_mpi$loo()
# model021_loo <- model021_fit_mpi$loo()
model018_loo <- model018_fit_mpi$loo()
model018_fit_no_mot_loo <- model018_fit_no_mot$loo()


loo::loo_compare(
  model022_loo,
  # model021_loo,
  model018_loo
  # model018_fit_no_mot_loo
)

loo::loo_model_weights(
  list(
    model022_loo,
    # model021_loo,
    model018_loo
  )
)

# look at the mu recreation of data

par_index <- 5
posterior_indices <- sample(1:2000, 100)

beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')
model022_beta_draws <- model022_fit_mpi$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

beta_names <- paste0('betas[', par_index, ',1,', c(1:19), ']')
model018_beta_draws <- model018_fit_mpi$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

beta_names <- paste0('betas[', par_index, ',1,', c(1:13), ']')
model018_no_mot_beta_draws <- model018_fit_no_mot$draws(
  format = "df",
  variables = beta_names
) %>%
  select(starts_with("betas")) %>%
  as.matrix()

uncensored <- fmri_stan_list$usable_bold_indices_one_is_true[par_index, ] == 1

current_par_design_matrix <- fmri_stan_list$design_array[
  par_index,
  uncensored,
  # beta
]

mod18_no_mot_posterior_mu <- (current_par_design_matrix[, 1:13] %*%
  t(model018_no_mot_beta_draws[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored]) %>%
  pivot_longer(starts_with("V"))

mod22_posterior_mu <- (current_par_design_matrix %*%
  t(model022_beta_draws[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored]) %>%
  pivot_longer(starts_with("V"))

mod18_posterior_mu <- (current_par_design_matrix %*%
  t(model018_beta_draws[posterior_indices, ])) %>%
  as_tibble() %>%
  mutate(time_seconds = seq(0, 1069 * 2, by = 2)[uncensored]) %>%
  pivot_longer(starts_with("V"))

bold_plot_df <- data.frame(
  time_seconds = seq(0, 1069 * 2, by = 2)[uncensored],
  bold = fmri_stan_list$bold[par_index, uncensored]
)

ggplot() +
  geom_line(
    data = bold_plot_df,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .2
  ) +
  geom_line(
    data = mod18_no_mot_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .5,
    alpha = .2,
    color = "red"
  ) +
  theme_classic()

ggplot() +
  geom_line(
    data = bold_plot_df,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .2
  ) +
  # geom_line(
  #   data = mod22_posterior_mu[mod22_posterior_mu$used, ],
  geom_line(
    data = mod22_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .2,
    alpha = .2,
    color = "red"
  ) +
  theme_classic()

ggplot() +
  geom_line(
    data = bold_plot_df,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .2
  ) +
  # geom_line(
  #   data = mod22_posterior_mu[mod22_posterior_mu$used, ],
  geom_line(
    data = mod18_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .2,
    alpha = .2,
    color = "red"
  ) +
  theme_classic()

ggplot() +
  geom_line(
    data = bold_plot_df,
    aes(
      x = time_seconds,
      y = bold
    ),
    linewidth = .2
  ) +
  # geom_line(
  #   data = mod22_posterior_mu[mod22_posterior_mu$used, ],
  geom_line(
    data = mod22_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .2,
    alpha = .2,
    color = "blue"
  ) +
  geom_line(
    data = mod18_posterior_mu,
    aes(
      x = time_seconds,
      y = value,
      group = name
    ),
    linewidth = .2,
    alpha = .2,
    color = "red"
  ) +
  theme_classic()


# model018_fit_mpi_summary <- model018_fit_mpi$summary()

model018_fit_mpi_relevant_parameters <- model018_fit_mpi_meta_data$model_params[
  str_detect(
    model018_fit_mpi_meta_data$model_params,
    "log_lik"
  )
]

model018_mpi_draws <- model018_fit_mpi$draws(
  format = "df",
  variables = model018_fit_mpi_relevant_parameters
)

model018_mpi_draws %>%
  select(contains("log_lik[")) %>%
  # unlist() %>%
  # sum()
  rowSums() %>%
  density() %>%
  plot()

model018_loo <- model018_fit_mpi$loo()


model019_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/model019_chain_66250437_1.csv"
  )
)

model019_fit_mpi_meta_data <- model019_fit_mpi$metadata()

model019_fit_mpi_relevant_parameters <- model019_fit_mpi_meta_data$model_params[
  !str_detect(
    model019_fit_mpi_meta_data$model_params,
    "mu|amplitude|raw|S|theta|bold_z"
  )
]


model019_fit_mpi_summary <- model019_fit_mpi$summary(
  variables = model019_fit_mpi_relevant_parameters
)

model019_fit_mpi_summary <- model019_fit_mpi$summary()


model019_mpi_draws <- model019_fit_mpi$draws(format = "df")

model019_loo <- model019_fit_mpi$loo()


model013_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_1.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_2.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_3.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_4.csv",
    "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains/model013_chain_65201387_5.csv"
  )
)

model013_fit_mpi_summary <- model013_fit_mpi$summary()


model013_mpi_draws <- model013_fit_mpi$draws(format = "df")

model013_loo <- model013_fit_mpi$loo()


model016_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/model016_chain_65887620_1.csv"
  )
)

model016_draws <- model016_fit_mpi$draws(format = "df")

model018_fit_mpi <- as_cmdstan_fit(
  files = c(
    "/home/andrewf/Downloads/model018_chain_65957320_1.csv",
    "/home/andrewf/Downloads/model018_chain_65957320_3.csv",
    "/home/andrewf/Downloads/model018_chain_65986961_1.csv"
  )
)

model018_fit_mpi_summary <- model018_fit_mpi$summary()

model018_draws <- model018_fit_mpi$draws(format = "df")

names(model016_draws)
names(model018_draws)

model016_draws %>%
  select(starts_with('betas[')) %>%
  select(ends_with(c(",1,5]", ",1,13]"))) %>%
  # pivot_longer(cols = ends_with(",1,5]")) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model018_draws %>%
  select(starts_with('betas[')) %>%
  select(ends_with(c(",1,5]", ",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,5]")) %>%
  # pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model016_draws %>%
  select(starts_with(c(
    'betas[8,1,5]',
    'betas[13,1,5]',
    'betas[14,1,5]',
    'betas[15,1,5]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,5]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model018_draws %>%
  select(starts_with(c(
    'betas[8,1,5]',
    'betas[13,1,5]',
    'betas[14,1,5]',
    'betas[15,1,5]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,5]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model016_draws %>%
  select(starts_with(c(
    'betas[8,1,13]',
    'betas[13,1,13]',
    'betas[14,1,13]',
    'betas[15,1,13]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))

model018_draws %>%
  select(starts_with(c(
    'betas[8,1,13]',
    'betas[13,1,13]',
    'betas[14,1,13]',
    'betas[15,1,13]'
  ))) %>%
  # select(ends_with(c(",1,5]",",1,13]"))) %>%
  pivot_longer(cols = ends_with(",1,13]")) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))


# can I get the per element log_lik

model018_draws$`log_lik[1]`[1]

posterior_indices <- sample(1:3000, size = 10)

fmri_stan_list$usable_bold_indices_one_is_true[1, ]

cov_mat <- matrix(nrow = 1070, ncol = 1070)

for (i in 1:1070) {
  for (j in 1:1070) {
    if (i == j) {
      cov_mat[i, j] = (model018_draws$`sigma[1,1]`[1])^2
    } else {
      cov_mat[i, j] = model018_draws$`delta[1,1]`[1] *
        ((model018_draws$`rho_time[1,1]`[1])^abs(i - j)) *
        ((model018_draws$`sigma[1,1]`[1])^2)
    }
  }
}

l_cov_mat <- t(chol(cov_mat[
  fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1,
  fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1
]))

dm <- fmri_stan_list$design_array[1, , ] %>% matrix(nrow = 1070)

beta_vector <- model018_draws %>%
  select(starts_with("betas[1,")) %>%
  slice(1) %>%
  unlist() %>%
  matrix(ncol = 1)

# visualize
(dm %*% beta_vector) %>% plot(type = "l")
fmri_stan_list$bold[1, ] %>% plot(type = "l")

mu <- dm %*% beta_vector
bold <- fmri_stan_list$bold[1, ]

mu_censor <- mu[fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1]
bold_censor <- bold[fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1]

mu_censor %>% plot(type = "l")
bold_censor %>% plot(type = "l")

data.frame(
  time = seq(0, 1069 * 2, by = 2),
  usable = fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1,
  mu,
  bold
) %>%
  # filter(usable ==T) %>%
  ggplot() +
  geom_line(aes(x = time, y = bold)) +
  geom_line(aes(x = time, y = mu), color = "red") +
  theme_bw()

bold_censor_z <- solve(l_cov_mat, (bold_censor - mu_censor))

bold_censor_z %>% plot(type = "l")

log_liks <- vector(length = length(bold_censor))

for (i in 1:length(bold_censor)) {
  log_liks[i] <- dnorm(bold_censor_z[i], log = T) - log(l_cov_mat[i, i])
}

# matches perfectly
sum(log_liks)
model018_draws$`log_lik[1]`[1]


# visualize posterior draws
posterior_indices <- sample(1:3000, size = 10)
cov_list <- list()

for (i in 1:length(posterior_indices)) {
  # vectorized cov_mat for fun
  # pull out your parameters
  sigma_val <- model018_draws$`sigma[1,1]`[i]
  delta_val <- model018_draws$`delta[1,1]`[i]
  rho_val <- model018_draws$`rho_time[1,1]`[i]
  s2 <- sigma_val^2

  # index vector
  n <- 1070
  idx <- seq_len(n)

  # 1) build the |i-j| distance matrix
  dist_mat <- abs(outer(idx, idx, FUN = "-"))

  # 2) fill in covariances: delta*rho^|i-j|*sigma^2 everywhere
  cov_mat <- s2 * delta_val * rho_val^dist_mat

  # 3) override the diagonal to be sigma^2
  diag(cov_mat) <- s2

  cov_list[[i]] <- cov_mat

  dm <- fmri_stan_list$design_array[1, , ] %>% matrix(nrow = 1070)

  #must censor dm before multiplying with beta
  dm_censor <- dm[fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1, ]

  beta_vector <- model018_draws %>%
    select(starts_with("betas[1,")) %>%
    slice(posterior_indices) %>%
    unlist() %>%
    matrix(ncol = length(posterior_indices))

  mu_censor <- dm_censor %*% beta_vector
}


bold <- fmri_stan_list$bold[1, ]

bold_censor <- bold[fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1]

# visualize
str(mu)

bold_plot_df <- data.frame(
  time = seq(0, 1069 * 2, by = 2),
  usable = fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1,
  bold
) %>%
  filter(usable == T)

mu_plot_df <- data.frame(
  time = seq(0, 1069 * 2, by = 2)[
    fmri_stan_list$usable_bold_indices_one_is_true[1, ] == 1
  ],
  mu_censor
) %>%
  pivot_longer(starts_with("X"))


ggplot() +
  geom_line(
    data = bold_plot_df,
    aes(x = time, y = bold)
  ) +
  geom_line(
    data = mu_plot_df,
    aes(x = time, y = value, group = name),
    color = "red",
    line_width = 0.2,
    alpha = .5
  ) +
  # coord_cartesian(xlim = c(0,500)) +
  theme_bw()


(dm %*% beta_vector) %>% plot(type = "l")
fmri_stan_list$bold[1, ] %>% plot(type = "l")


mu_censor %>% plot(type = "l")
bold_censor %>% plot(type = "l")
