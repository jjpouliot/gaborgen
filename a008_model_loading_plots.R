# ============================================================================
#  Gaborgen24 fMRI — Stan model loading + plotting
# ----------------------------------------------------------------------------
#  Refactored from the original misc/model_loading_plots.R, which read each
#  ROI's chains and built each figure inline. Here everything is driven from
#  a registry tibble + a handful of small functions.
#
#  Layout
#  ------
#    1. Libraries
#    2. Configuration         (paths, palettes, intervals, y-limits)
#    3. Beta metadata         (mu_betas[j] -> phase/cue
#                              betas[trial_idx] -> cue/block/trial_per_cue)
#    4. Chain registry        (one row per model x ROI)
#    5. Path / loader funs    (read_mu_betas_chains, read_trial_betas_chains)
#    6. Tidiers               (tidy_mu_betas, tidy_trial_betas,
#                              compute_cue_metrics)
#    7. Plotters              (plot_phase_cue_posterior,
#                              plot_trial_timeseries,
#                              plot_metrics_timeseries,
#                              plot_roi_summary, save_roi_summary)
#    8. Optional: LOO         (compute_loo)
#    9. Example driver        (a single ROI end-to-end at the bottom)
#
#  To use as a library, source() the file and ignore section 9.
# ============================================================================

# ---- 1. Libraries ----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(patchwork)
library(posterior)
library(loo)
library(tidybayes)


# ---- 2. Configuration ------------------------------------------------------

DATA_DIR <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains"
PLOT_DIR <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/misc"
parent_folder <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI"

# Cue palettes — pretty labels for mu_betas plots, lowercase for the
# trial-level data coming out of model038.
cue_palette_pretty <- c(
  "CS+" = "red1",
  "GS1" = "green1",
  "GS2" = "purple1",
  "GS3" = "blue1",
  "Shock" = "black"
)

cue_palette_raw <- c(
  csp = "red1",
  gs1 = "green1",
  gs2 = "purple1",
  gs3 = "blue1",
  shock = "black"
)

# 68% interval (~ 1 SD) — used for the time-series ribbon
band_lower <- 0.5 - 0.341
band_upper <- 0.5 + 0.341

# Vertical lines that mark block boundaries on the trial_per_cue axis
trial_block_breaks <- c(9, 21, 32)

# Per-ROI y-axis limits (NULL = let ggplot pick). Used by plot_roi_summary().
roi_ylim <- list(
  "V4" = c(-0.25, 0.7),
  "V5" = c(-0.025, 0.26),
  "V5" = c(-0.1, 0.3),
  "TE" = c(-0.11, 0.2),
  "TPJ" = c(-0.05, 0.2),
  "ACC" = c(-0.1, 0.25),
  "Amygdala" = c(-0.075, 0.33),
  "Hippocampus" = c(-0.05, 0.15),
  "Anterior Insula" = c(-0.05, 0.33)
  # "OFC" = c(-0.05, 0.33)
  # add more as needed
)


# ---- 3. Beta metadata ------------------------------------------------------

#  mu_betas[roi, j]  for j = 1..17  (model039 / model040)
#    j  1- 4 : Habituation    x (CS+, GS1, GS2, GS3)
#    j  5- 8 : Acquisition #1 x (CS+, GS1, GS2, GS3)
#    j  9-12 : Acquisition #2 x (CS+, GS1, GS2, GS3)
#    j 13-16 : Extinction     x (CS+, GS1, GS2, GS3)
#    j 17    : Shock
mu_beta_meta <- tibble(
  j = 1:17,
  phase = factor(
    c(
      rep("Habituation", 4),
      rep("Acquisition #1", 4),
      rep("Acquisition #2", 4),
      rep("Extinction", 4),
      "Shock"
    ),
    levels = c(
      "Habituation",
      "Acquisition #1",
      "Acquisition #2",
      "Shock",
      "Extinction"
    ) # Shock between Acq2 & Ext
  ),
  cue = factor(
    c(rep(c("CS+", "GS1", "GS2", "GS3"), 4), "Shock"),
    levels = c("CS+", "GS1", "GS2", "GS3", "Shock")
  )
)

# betas[roi, n] for n = 1..191  (model038, single-trial)
#    1- 32  habituation     :  4 cues x  8 trials
#   33- 80  acquisition #1  :  4 cues x 12 trials
#   81-128  acquisition #2  :  4 cues x 12 trials
#  129-176  extinction      :  4 cues x 12 trials
#  177-191  shock           : 15 trials, x-positions chosen so they line up
#                              with acquisition trial spacing
trial_beta_meta <- bind_rows(
  tibble(
    beta_index = 1:32,
    block = "habituation",
    cue = rep(c("csp", "gs1", "gs2", "gs3"), each = 8),
    trial_per_cue = rep(1:8, times = 4)
  ),
  tibble(
    beta_index = 33:80,
    block = "acquisition #1",
    cue = rep(c("csp", "gs1", "gs2", "gs3"), each = 12),
    trial_per_cue = rep(9:20, times = 4)
  ),
  tibble(
    beta_index = 81:128,
    block = "acquisition #2",
    cue = rep(c("csp", "gs1", "gs2", "gs3"), each = 12),
    trial_per_cue = rep(21:32, times = 4)
  ),
  tibble(
    beta_index = 129:176,
    block = "extinction",
    cue = rep(c("csp", "gs1", "gs2", "gs3"), each = 12),
    trial_per_cue = rep(33:44, times = 4)
  ),
  tibble(
    beta_index = 177:191,
    block = "shock",
    cue = "shock",
    trial_per_cue = c(9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32)
  )
)

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
  "158", # made
  "159", # made
  "160", # made
  "161" # made
)

# ---- 4. Chain registry -----------------------------------------------------
# One row per (model_id, roi_id). `chain_nums` is a comma-separated string
# of chain numbers to load — usually "1,2", but some fits only have one
# chain (or only chain 2, in V6's case) and that's encoded here.

chain_registry <- tibble::tribble(
  ~model_id  , ~roi_id                   , ~roi_label          , ~suffix    , ~job_id    , ~chain_nums ,
  # --- model042 (hierarchical: mu_betas + per-subject betas with t-distribution regularization) ---------------
  "model042" , "V1F"                     , "V1F"               , "pb03_res" , "33726894" , "1,2"       ,
  "model042" , "V4"                      , "V4"                , "pb03_res" , "33218289" , "1,2"       ,
  "model042" , "V5"                      , "V5"                , "pb03_res" , "33233099" , "1,2"       ,
  "model042" , "V6"                      , "V6"                , "pb03_res" , "33233105" , "1,2"       ,
  "model042" , "TE"                      , "TE"                , "pb03_res" , "33233115" , "1,2"       ,
  "model042" , "TPJ"                     , "TPJ"               , "pb03_res" , "33233116" , "1,2"       ,
  "model042" , "ACC"                     , "ACC"               , "pb03_res" , "33233118" , "1,2"       ,
  "model042" , "NA"                      , "Nucleus Accumbens" , "pb03_res" , "33233130" , "1,2"       ,
  "model042" , "AMY"                     , "Amygdala"          , "pb03_res" , "33099488" , "1,2"       ,
  "model042" , "HIP"                     , "Hippocampus"       , "pb03_res" , "33233131" , "1,2"       ,
  "model042" , "ant_insula"              , "Anterior Insula"   , "pb03_res" , "33066991" , "1,2"       ,
  "model042" , "OFC"                     , "OFC"               , "pb03_res" , "33349421" , "1"         ,

  # --- model039 (hierarchical: mu_betas + per-subject betas) ---------------
  "model039" , "V4"                      , "V4"                , "pb03_res" , "32674652" , "1,2"       ,
  "model039" , "V5"                      , "V5"                , "pb03_res" , "32674687" , "1,2"       ,
  "model039" , "V6"                      , "V6"                , "pb03_res" , "32674694" , "1,2"       ,
  "model039" , "TE"                      , "TE"                , "pb03_res" , "32674818" , "1,2"       ,
  "model039" , "TPJ"                     , "TPJ"               , "pb03_res" , "32674863" , "1,2"       ,
  "model039" , "ACC"                     , "ACC"               , "pb03_res" , "32675094" , "1,2"       ,
  "model039" , "NA"                      , "Nucleus Accumbens" , "pb03_res" , "32675096" , "1,2"       ,
  "model039" , "AMY"                     , "Amygdala"          , "pb03_res" , "32675097" , "1,2"       ,
  "model039" , "HIP"                     , "Hippocampus"       , "pb03_res" , "32675098" , "1,2"       ,
  "model039" , "ant_insula"              , "Anterior Insula"   , "pb03_res" , "32658847" , "1,2"       ,
  "model039" , "OFC"                     , "OFC"               , "pb03_res" , "32675123" , "1,2"       ,
  # "model040",  "V4",                        "V4 (m040)",         "pb03_res", "32816165", "1",

  # --- model038 (single-trial betas) --------------------------------------
  "model038" , "V1F"                     , "V1F"               , "pb03_res" , "31602869" , "1"         ,
  "model038" , "V4"                      , "V4"                , "pb03_res" , "31836148" , "1,2"       ,
  "model038" , "V5"                      , "V5"                , "pb03_res" , "32007843" , "1,2"       ,
  "model038" , "V6"                      , "V6"                , "pb03"     , "32236921" , "1,2"       ,
  "model038" , "TE"                      , "TE"                , "pb03_res" , "32236914" , "1,2"       ,
  "model038" , "TPJ"                     , "TPJ"               , "pb03_res" , "32236922" , "1,2"       ,
  "model038" , "ACC"                     , "ACC"               , "pb03_res" , "32236924" , "1,2"       ,
  "model038" , "Amygdala"                , "Amygdala"          , "pb03_res" , "31811588" , "1,2"       ,
  "model038" , "Hippocampus"             , "Hippocampus"       , "pb03_res" , "31835580" , "1,2"       ,
  "model038" , "ant_insula"              , "Anterior Insula"   , "pb03_res" , "31603369" , "1,2"       ,
  "model038" , "Orbital_Frontal_Complex" , "OFC"               , "pb03_res" , "32315249" , "1,2"       ,
)


# ---- 5. Path / loader functions --------------------------------------------

# Build the chain CSV paths from a registry row.
chain_paths <- function(
  model_id,
  roi_id,
  suffix,
  job_id,
  chain_nums = "1,2",
  data_dir = DATA_DIR
) {
  nums <- as.integer(strsplit(chain_nums, ",", fixed = TRUE)[[1]])
  fname <- function(i) {
    sprintf("%s_%s_%s_chain_%s_%d.csv", model_id, roi_id, suffix, job_id, i)
  }
  file.path(data_dir, vapply(nums, fname, character(1)))
}

# Look up a (model_id, roi_id) row in the registry and return CSV paths.
chain_paths_lookup <- function(
  model_id,
  roi_id,
  registry = chain_registry,
  data_dir = DATA_DIR
) {
  row <- registry %>%
    filter(.data$model_id == .env$model_id, .data$roi_id == .env$roi_id)
  if (nrow(row) != 1) {
    stop(
      "chain_paths_lookup(): expected exactly 1 row for ",
      model_id,
      " / ",
      roi_id,
      ", got ",
      nrow(row)
    )
  }
  chain_paths(
    row$model_id,
    row$roi_id,
    row$suffix,
    row$job_id,
    row$chain_nums,
    data_dir = data_dir
  )
}

# Loader for hierarchical fits (model039 / model040). Defaults to reading
# only `mu_betas` — much faster than pulling all per-subject betas. Pass
# `variables = c("mu_betas", "betas", "log_lik")` to get everything.
read_chains <- function(
  model_id = NULL,
  roi_id = NULL,
  variables = "mu_betas",
  files = NULL,
  registry = chain_registry,
  data_dir = DATA_DIR
) {
  if (is.null(files)) {
    files <- chain_paths_lookup(model_id, roi_id, registry, data_dir)
  }
  fit <- read_cmdstan_csv(files = files, variables = variables)
  fit$post_warmup_draws
}

# Loader for single-trial fits (model038). Trial-level betas live in
# columns betas[roi, n] for roi in 1..n_roi and n in 1..n_trial_betas
# (i.e. total betas minus motion regressors). n_trial_betas = 191 covers
# all 176 cue trials plus 15 shock trials.
read_trial_betas_chains <- function(
  model_id = NULL,
  roi_id = NULL,
  n_roi = 2,
  n_trial_betas = 191,
  files = NULL,
  registry = chain_registry,
  data_dir = DATA_DIR
) {
  vars <- paste0(
    "betas[",
    rep(seq_len(n_roi), each = n_trial_betas),
    ",",
    rep(seq_len(n_trial_betas), times = n_roi),
    "]"
  )
  if (is.null(files)) {
    files <- chain_paths_lookup(model_id, roi_id, registry, data_dir)
  }
  fit <- read_cmdstan_csv(files = files, variables = vars)
  fit$post_warmup_draws
}


# ---- 6. Tidiers ------------------------------------------------------------

# Long form: average over ROI dimension (i) within each draw, then attach
# phase/cue labels. One row per (draw, j).
tidy_mu_betas <- function(draws, meta = mu_beta_meta) {
  draws %>%
    spread_draws(mu_betas[i, j]) %>%
    group_by(.draw, j) %>%
    summarise(mu_beta = mean(mu_betas), .groups = "drop") %>%
    left_join(meta, by = "j") %>%
    filter(!is.na(phase))
}

# Long form for single-trial betas: one row per (draw, roi, trial), with
# cue / block / trial_per_cue attached.
tidy_trial_betas <- function(draws, meta = trial_beta_meta) {
  posterior::as_draws_df(draws) %>%
    pivot_longer(
      cols = starts_with("betas["),
      names_to = c(".unused", "roi", "beta_index"),
      names_sep = "\\[|,|\\]",
      values_to = "beta_value"
    ) %>%
    select(-.unused) %>%
    mutate(
      roi = as.integer(roi),
      beta_index = as.integer(beta_index)
    ) %>%
    left_join(meta, by = "beta_index")
}

# Per-draw generalization / selectivity / sharpening over the cue trials.
compute_cue_metrics <- function(trial_long) {
  trial_long %>%
    filter(cue != "shock") %>%
    group_by(.draw, trial_per_cue) %>%
    summarise(
      csp = mean(beta_value[cue == "csp"]),
      gs1 = mean(beta_value[cue == "gs1"]),
      gs2 = mean(beta_value[cue == "gs2"]),
      gs3 = mean(beta_value[cue == "gs3"]),
      .groups = "drop"
    ) %>%
    mutate(
      generalization = (csp + gs1 + gs2 + gs3) / 4,
      selectivity = csp - (gs1 + gs2 + gs3) / 3,
      sharpening = csp - gs1 - gs2 + gs3
    ) %>%
    pivot_longer(
      c(generalization, selectivity, sharpening),
      names_to = "metric",
      values_to = "value"
    )
}

# Per-participant tidied betas from model039: averages over ROI within
# each draw and joins phase/cue labels. One row per (.draw, participant,
# j). Used by compute_participant_metrics() and plot_participant_cue_betas().
tidy_participant_betas <- function(draws, meta = mu_beta_meta) {
  draws %>%
    spread_draws(betas[participant, roi, j]) %>%
    group_by(.draw, participant, j) %>%
    summarise(beta = mean(betas), .groups = "drop") %>% # avg over ROI
    left_join(meta, by = "j")
}


# Per-participant, per-phase generalization / selectivity / sharpening from
# model039's per-subject `betas[participant, roi, j]`. ROIs are averaged
# first; the same three formulas as compute_cue_metrics() are then applied
# within each (draw, participant, phase). Shock is dropped since the
# metrics need all four cues.
#
# Returns long form: one row per (.draw, participant, phase, metric).
compute_participant_metrics <- function(draws, meta = mu_beta_meta) {
  draws %>%
    spread_draws(betas[participant, roi, j]) %>%
    group_by(.draw, participant, j) %>%
    summarise(beta = mean(betas), .groups = "drop") %>% # avg over ROI
    left_join(meta, by = "j") %>%
    filter(cue != "Shock") %>%
    select(.draw, participant, phase, cue, beta) %>%
    pivot_wider(names_from = cue, values_from = beta) %>%
    mutate(
      avg_gen = (`CS+` + GS1 + GS2 + GS3) / 4, # was "generalization": overall level
      linear_gen = -1 * ((-3 * `CS+` - 1 * GS1 + 1 * GS2 + 3 * GS3) / 10), # linear slope per similarity step
      selectivity = `CS+` - (GS1 + GS2 + GS3) / 3,
      sharpening = `CS+` - GS1 - GS2 + GS3 # = quadratic contrast (CSP + GS3) - (GS1 + GS2)
    ) %>%
    select(
      .draw,
      participant,
      phase,
      avg_gen,
      linear_gen,
      selectivity,
      sharpening
    ) %>%
    pivot_longer(
      c(avg_gen, linear_gen, selectivity, sharpening),
      names_to = "metric",
      values_to = "value"
    )
}


# ---- 7. Plotters -----------------------------------------------------------

# Phase x cue posterior (right-hand panel of the summary figure).
plot_phase_cue_posterior <- function(
  mu_long,
  roi_label,
  ylim = NULL,
  palette = cue_palette_pretty
) {
  p <- mu_long %>%
    ggplot(aes(x = cue, y = mu_beta, fill = cue, color = cue)) +
    stat_halfeye(
      .width = .341,
      slab_alpha = 0.55,
      point_colour = "black",
      point_size = 3,
      interval_colour = "black",
      interval_size_range = c(1, 2)
    ) +
    geom_hline(yintercept = 0) +
    facet_grid(. ~ phase, scales = "free_x", space = "free_x", switch = "x") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(x = NULL, y = "Posterior of mean β (avg. over ROI)") +
    ggtitle(paste0(roi_label, " — cue x block posteriors")) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0.2, "lines"),
      strip.placement = "outside",
      legend.position = "right"
    )

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  p
}

# Time-series across trial_per_cue (left-hand panel).
plot_trial_timeseries <- function(
  trial_long,
  roi_label,
  ylim = NULL,
  palette = cue_palette_raw,
  ci_lower = band_lower,
  ci_upper = band_upper,
  breaks = trial_block_breaks
) {
  summarised <- trial_long %>%
    group_by(.draw, cue, trial_per_cue, beta_index) %>%
    summarise(avg_roi_draw = mean(beta_value), .groups = "drop") %>%
    group_by(cue, trial_per_cue, beta_index) %>%
    summarise(
      median_posterior = median(avg_roi_draw),
      lower = quantile(avg_roi_draw, ci_lower),
      upper = quantile(avg_roi_draw, ci_upper),
      .groups = "drop"
    )

  p <- summarised %>%
    ggplot(aes(x = trial_per_cue, color = cue, fill = cue)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = breaks) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(aes(y = median_posterior), linewidth = 0.4, alpha = 0.35) +
    geom_smooth(
      aes(y = median_posterior),
      se = FALSE,
      linewidth = 1.3,
      span = 0.4
    ) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    ggtitle(paste0("Change in average ", roi_label, " over trials")) +
    theme_bw(base_size = 20)

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  p
}

# Metrics (generalization / selectivity / sharpening) over trials.
plot_metrics_timeseries <- function(metrics_long, breaks = trial_block_breaks) {
  metrics_long %>%
    ggplot(aes(x = trial_per_cue, y = value, fill = metric, color = metric)) +
    geom_vline(xintercept = breaks) +
    tidybayes::stat_lineribbon(.width = .341, alpha = 0.25) +
    geom_hline(yintercept = 0) +
    theme_bw(base_size = 20)
}

# Per-participant posteriors of generalization / selectivity / sharpening,
# laid out as point-ranges with beeswarm-style horizontal jitter inside
# each phase. One pointrange per participant; the central credible
# interval (`.width`) shows the inner density of each participant's
# posterior, and ggbeeswarm spreads overlapping participants apart so the
# group-level distribution is also legible.
#
# `.width` is the symmetric credible interval (default 68% ~ 1 SD).
# `jitter_width` controls how wide each phase's beeswarm can spread.
plot_participant_metrics <- function(
  participant_metrics,
  .width = 0.68,
  point_size = 0.3,
  fatten = 1.8,
  alpha = 0.75,
  free_y = F,
  show_labels = FALSE
) {
  half <- .width / 2

  per_participant <- participant_metrics %>%
    mutate(participant = factor(participant)) %>%
    group_by(participant, phase, metric) %>%
    summarise(
      median = median(value),
      lower = quantile(value, 0.5 - half),
      upper = quantile(value, 0.5 + half),
      .groups = "drop"
    )

  ggplot(
    per_participant,
    aes(x = participant, y = median, ymin = lower, ymax = upper)
  ) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_pointrange(size = point_size, fatten = fatten, alpha = alpha) +
    facet_grid(
      metric ~ phase,
      scales = if (free_y) {
        "free_y"
      } else {
        NULL
      },
      switch = "y"
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.placement = "outside",
      axis.text.x = if (show_labels) {
        element_text(angle = 90, hjust = 1, size = 8)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_labels) {
        element_line()
      } else {
        element_blank()
      },
      panel.spacing.x = unit(0.4, "lines")
    ) +
    labs(
      x = "Participant",
      y = paste0(
        "Per-participant posterior (median +/- ",
        round(.width * 100),
        "% CI)"
      )
    )
}

# Per-participant beta per cue x block. Same horizontal layout as
# plot_participant_metrics (participants ordered by ID, same slot across
# all blocks) but with 4 pointranges per participant in each non-shock
# phase -- one per cue, dodged side-by-side and colored using
# `cue_palette_pretty`. Shock appears as its own phase column with one
# pointrange per participant, matching the mu_betas headline figure.
#
# With 4 cues and 44 participants the dodge can get tight; bump `dodge_w`
# up to spread cues further apart, or shrink `point_size` / `fatten` if
# the points look heavy.
plot_participant_cue_betas <- function(
  participant_beta_long,
  .width = 0.68,
  point_size = 0.2,
  fatten = 1.5,
  alpha = 0.8,
  dodge_w = 0.75,
  show_labels = FALSE,
  palette = cue_palette_pretty
) {
  half <- .width / 2

  per_participant <- participant_beta_long %>%
    mutate(participant = factor(participant)) %>%
    group_by(participant, phase, cue) %>%
    summarise(
      median = median(beta),
      lower = quantile(beta, 0.5 - half),
      upper = quantile(beta, 0.5 + half),
      .groups = "drop"
    )

  ggplot(
    per_participant,
    aes(x = participant, y = median, ymin = lower, ymax = upper, color = cue)
  ) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_pointrange(
      position = position_dodge(width = dodge_w),
      size = point_size,
      fatten = fatten,
      alpha = alpha
    ) +
    facet_grid(. ~ phase, scales = "free_x", space = "free_x", switch = "x") +
    scale_color_manual(values = palette) +
    theme_bw(base_size = 14) +
    theme(
      strip.placement = "outside",
      axis.text.x = if (show_labels) {
        element_text(angle = 90, hjust = 1, size = 6)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_labels) {
        element_line()
      } else {
        element_blank()
      },
      panel.spacing.x = unit(0.4, "lines"),
      legend.position = "right"
    ) +
    labs(
      x = "Participant",
      y = paste0("Per-participant β (median +/- ", round(.width * 100), "% CI)")
    )
}


# Headline figure: time-series on the left, phase x cue posterior on the
# right. `ylim` defaults to roi_ylim[[roi_label]] when available, NULL
# otherwise.
plot_roi_summary <- function(
  trial_draws,
  mu_draws,
  roi_label,
  ylim = roi_ylim[[roi_label]]
) {
  mu_long <- tidy_mu_betas(mu_draws)
  trial_long <- tidy_trial_betas(trial_draws)

  ts_plot <- plot_trial_timeseries(trial_long, roi_label, ylim = ylim)
  posterior_plot <- plot_phase_cue_posterior(mu_long, roi_label, ylim = ylim)

  posterior_plot +
    ts_plot +
    theme(legend.position = "right")
}

# Save a combined ggplot/patchwork object to PLOT_DIR.
save_roi_summary <- function(
  combined_plot,
  roi_label,
  dir = PLOT_DIR,
  prefix = "mod38_39_",
  width = 24,
  height = 12,
  dpi = 100
) {
  fname <- file.path(
    dir,
    paste0(prefix, gsub(" ", "_", roi_label), ".png")
  )
  ggsave(
    filename = fname,
    plot = combined_plot,
    width = width,
    height = height,
    dpi = dpi
  )
  invisible(fname)
}


# ---- 8. Optional: LOO ------------------------------------------------------

# Approximate leave-one-out CV. Requires log_lik to have been read in with
# the draws (pass variables = c("mu_betas", "log_lik") to the loader).
compute_loo <- function(draws) {
  ll <- draws %>%
    posterior::subset_draws(variable = "log_lik") %>%
    posterior::as_draws_array()
  r_eff <- loo::relative_eff(exp(ll), cores = 10)
  loo::loo(ll, r_eff = r_eff, cores = 10)
}

# -- temp ratings

Day1_ratings_paths <- list.files(
  path = parent_folder,
  pattern = "Day1.*ratings.dat$",
  recursive = T,
  full.names = T
)


for (i in 1:length(Day1_ratings_paths)) {
  current_ratings <- read.csv(Day1_ratings_paths[i])
  if (i == 1) {
    ratings_df <- current_ratings
  } else {
    ratings_df <- rbind.data.frame(ratings_df, current_ratings)
  }
}

ratings_df$partInd <- as.character(ratings_df$partInd)

# experiment mislabels a rating
names(ratings_df)[9] <- "ar_gs2"

ratings_df <- ratings_df %>%
  select(!contains("Dur")) %>%
  pivot_longer(
    cols = c(
      val_csp,
      val_gs1,
      val_gs2,
      val_gs3,
      ar_csp,
      ar_gs1,
      ar_gs2,
      ar_gs3,
      exp_csp,
      exp_gs1,
      exp_gs2,
      exp_gs3
    ),
    names_to = c(".value", "condition"),
    names_sep = "_"
  )

ratings2 <- ratings_df %>%
  filter(partInd %in% useable_participants) %>%
  mutate(
    block = recode(
      ratInd,
      `1` = "hab",
      `2` = "acq_1",
      `3` = "acq_2",
      `4` = "ext"
    )
  ) %>% # ADJUST to your mapping
  group_by(partInd) %>%
  mutate(val_btw = mean(val), val_wth = val - mean(val)) %>% # keep the within/between split
  mutate(ar_btw = mean(ar), ar_wth = ar - mean(ar)) %>% # keep the within/between split
  mutate(exp_btw = mean(exp), exp_wth = exp - mean(exp)) %>% # keep the within/between split
  ungroup() %>%
  mutate(participant = as.integer(as.factor(partInd)))


ai_mu_full <- read_chains(
  "model042",
  "ant_insula",
  variables = c("mu_betas", "betas")
)

ai_pp_cue_betas <- tidy_participant_betas(ai_mu_full) %>%
  # filter(j < 18)
  filter(j < 17)

amy_mu_full <- read_chains(
  "model042",
  "AMY",
  variables = c("mu_betas", "betas")
)

amy_pp_cue_betas <- tidy_participant_betas(amy_mu_full) %>%
  # filter(j < 18)
  filter(j < 17)


# 1. Brain: collapse the 4000 draws to mean + SD per cell (insula only)
brain <- ai_pp_cue_betas %>%
  # brain <- amy_pp_cue_betas %>%
  group_by(participant, phase, cue) %>%
  summarise(brain_mean = mean(beta), brain_sd = sd(beta), .groups = "drop")

# 2. Reconcile labels so the keys match ratings2 exactly
cue_map <- c("CS+" = "csp", "GS1" = "gs1", "GS2" = "gs2", "GS3" = "gs3")
phase_map <- c(
  "Habituation" = "hab",
  "Acquisition #1" = "acq_1",
  "Acquisition #2" = "acq_2",
  "Extinction" = "ext"
)

brain <- brain %>%
  mutate(
    condition = cue_map[as.character(cue)],
    block = phase_map[as.character(phase)]
  ) #%>%
# left_join(part_key, by = "participant") # you supply: participant(int) <-> partInd(chr)

# 3. Join ratings (already carries val, val_wth, val_btw)
df <- brain %>%
  # df_amy <- brain %>%
  inner_join(
    ratings2 %>%
      select(
        participant,
        block,
        condition,
        val,
        val_wth,
        val_btw,
        ar,
        ar_wth,
        ar_btw,
        exp,
        exp_wth,
        exp_btw,
      ),
    by = c("participant", "block", "condition")
  )

nrow(df) # sanity check; chase down any cells dropped to label mismatch or NA val

gm <- mean(df$exp) # 0–100 grand mean (cell-weighted; fine for centering)
df <- df %>%
  mutate(
    exp_abs = (exp - gm) / 10, # absolute, grand-mean-centered  (= exp_wth10 + exp_btw10)
    exp_wth10 = exp_wth / 10, # within-person deviation
    exp_btw10 = (exp_btw - gm) / 10 # person mean, centered
  )

#extra models for exp testing
library(brms)

ctrl <- list(adapt_delta = 0.95)
common <- function(form) {
  brm(
    bf(form),
    data = df,
    chains = 4,
    cores = 4,
    save_pars = save_pars(all = TRUE),
    backend = "cmdstanr",
    control = ctrl
  )
}

# M1  current: within/between split (reference)
m_wb <- common(
  brain_mean | se(brain_sd, sigma = TRUE) ~
    exp_wth10 * block + exp_btw10 + (1 + exp_wth10 | participant)
)

# M2  absolute scale as ONE quantity — no within/between distinction
m_abs <- common(
  brain_mean | se(brain_sd, sigma = TRUE) ~
    exp_abs * block + (1 + exp_abs | participant)
)

# M3  absolute main slope + person mean as a confound guard (best of both)
m_abs_btw <- common(
  brain_mean | se(brain_sd, sigma = TRUE) ~
    exp_abs * block + exp_btw10 + (1 + exp_abs | participant)
)

# M4  nonlinear absolute shape (threshold / saturation), per block
m_abs_s <- common(
  brain_mean | se(brain_sd, sigma = TRUE) ~
    s(exp_abs, by = block) + block + exp_btw10 + (1 + exp_abs | participant)
)

m_ctx <- common(
  brain_mean | se(brain_sd, sigma = TRUE) ~
    exp_abs + exp_btw10 + block + (1 + exp_abs | participant)
)
hypothesis(m_ctx, "exp_btw10 = 0") # coef(exp_btw10) = between − within (the contextual effect)

m_wb <- add_criterion(m_wb, "loo", moment_match = TRUE)
m_abs <- add_criterion(m_abs, "loo", moment_match = TRUE)
m_abs_btw <- add_criterion(m_abs_btw, "loo", moment_match = TRUE)
m_abs_s <- add_criterion(m_abs_s, "loo", moment_match = TRUE)
loo_compare(m_wb, m_abs, m_abs_btw, m_abs_s)

m_wb$criteria$loo # prints ELPD + the Pareto-k table
loo::pareto_k_table(m_wb$criteria$loo) # counts binned good / ok / bad / very bad
ids <- pareto_k_ids(m_wb$criteria$loo, 0.7) # WHICH rows are unreliable
df[ids, c("participant", "block", "cue", "brain_mean", "brain_sd", "exp")]
plot(m_wb$criteria$loo) # k value per observation

m_wb <- add_criterion(m_wb, "loo", reloo = TRUE, overwrite = TRUE)
m_abs <- add_criterion(m_abs, "loo", reloo = TRUE, overwrite = TRUE)
m_abs_btw <- add_criterion(m_abs_btw, "loo", reloo = TRUE, overwrite = TRUE)
m_abs_s <- add_criterion(m_abs_s, "loo", reloo = TRUE, overwrite = TRUE)

m_wb$criteria$loo
m_abs$criteria$loo
m_abs_btw$criteria$loo
m_abs_s$criteria$loo
loo_compare(m_wb, m_abs, m_abs_btw, m_abs_s)

posterior_summary(m_ctx, variable = c("b_exp_abs", "b_exp_btw10"))

conditional_effects(m_abs_btw, effects = "exp_abs:block")

fixef(m_abs_btw)

summary(m_abs_btw)

# Same grand mean used to build exp_abs — must match
block_labels_named <- c(
  hab = "Habituation",
  acq_1 = "Acquisition 1",
  acq_2 = "Acquisition 2",
  ext = "Extinction"
)

cue_colors_val <- c(
  "CS+" = "#E63946",
  "GS1" = "#2DC653",
  "GS2" = "#9B5DE5",
  "GS3" = "#457B9D"
)

gm <- mean(df$exp)

df_cues <- df |> filter(cue %in% c("CS+", "GS1", "GS2", "GS3"))

# ---- population conditional effect over the ABSOLUTE expectancy axis ----
cond_grid_exp <- df |>
  group_by(block) |>
  summarise(
    exp_abs = list(seq(min(exp_abs), max(exp_abs), length.out = 50)),
    exp_btw10 = mean(df$exp_btw10),
    brain_sd = mean(df$brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(exp_abs)

cond_line_exp <- cond_grid_exp |>
  add_epred_draws(m_abs_btw, re_formula = NA) |>
  group_by(block, exp_abs) |>
  summarise(
    epred = median(.epred),
    lo = quantile(.epred, 0.025),
    hi = quantile(.epred, 0.975),
    .groups = "drop"
  )

# replicate the population line into each cue row for Plot B
cond_line_exp_by_cue <- cond_line_exp |>
  tidyr::crossing(
    cue = factor(
      c("CS+", "GS1", "GS2", "GS3"),
      levels = c("CS+", "GS1", "GS2", "GS3")
    )
  )

# ---- fading reference lines on the TRUE 0–100 grid (thickest at grand mean) ----
exp_label_vals <- seq(0, 100, by = 10)
vlines_exp <- tibble(
  label = exp_label_vals,
  x = (exp_label_vals - gm) / 10, # position in exp_abs units
  dist = abs(exp_label_vals - gm) / 10
) |>
  mutate(
    line_alpha = pmax(0.05, 0.55 * (0.72^dist)),
    line_lwd = pmax(0.10, 0.70 * (0.72^dist))
  )

# axis re-referenced to real shock-probability %
exp_x_scale <- scale_x_continuous(
  breaks = (seq(0, 100, by = 20) - gm) / 10,
  labels = seq(0, 100, by = 20),
  name = "Shock expectancy (%)"
)

# ---- Plot A: faceted by block ----
plot_a_exp <- ggplot() +
  geom_vline(
    data = vlines_exp,
    aes(xintercept = x, linewidth = line_lwd, alpha = line_alpha),
    colour = "grey30"
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(
    data = cond_line_exp,
    aes(x = exp_abs, ymin = lo, ymax = hi),
    alpha = 0.25,
    fill = "grey40"
  ) +
  geom_line(
    data = cond_line_exp,
    aes(x = exp_abs, y = epred),
    linewidth = 1.2,
    colour = "grey20"
  ) +
  # geom_point(data = df_cues,
  #            aes(x = exp_abs, y = brain_mean, colour = cue),
  #            size = 2.5, alpha = 0.8) +
  ggbeeswarm::geom_quasirandom(
    data = df_cues,
    aes(x = exp_abs, y = brain_mean, colour = cue),
    method = "quasirandom", # offset ∝ local density → fans out the dense stacks
    orientation = "x", # group by expectancy value, spread horizontally
    width = 0.30, # spread in exp_abs units (0.30 ≈ 3 percentage points)
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ factor(block, levels = names(block_labels_named)),
    labeller = as_labeller(block_labels_named),
    nrow = 2
  ) +
  scale_colour_manual(values = cue_colors_val) +
  exp_x_scale +
  guides(alpha = "none", linewidth = "none") +
  labs(
    y = "Anterior insula \u03b2",
    colour = "Cue",
    title = "Expectancy \u2192 anterior insula: conditional effect by block",
    subtitle = "Population-level line \u00b1 95% CI (m_abs_btw)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 14, colour = "grey50")
  )

ggsave(
  file.path(PLOT_DIR, "ai_exp_m_abs_btw_by_block.png"),
  plot = plot_a_exp,
  width = 20,
  height = 14,
  dpi = 150
)

# ---- Plot B: faceted cue × block ----
plot_b_exp <- ggplot() +
  # geom_vline(
  #   data = vlines_exp,
  #   aes(xintercept = x, linewidth = line_lwd, alpha = line_alpha),
  #   colour = "grey30"
  # ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(
    data = cond_line_exp_by_cue,
    aes(x = exp_abs, ymin = lo, ymax = hi),
    alpha = 0.25,
    fill = "grey40"
  ) +
  geom_line(
    data = cond_line_exp_by_cue,
    aes(x = exp_abs, y = epred),
    linewidth = 0.9,
    colour = "grey20"
  ) +
  # geom_point(data = df_cues,
  #            aes(x = exp_abs, y = brain_mean, colour = cue),
  #            size = 2.5, alpha = 0.8) +
  ggbeeswarm::geom_quasirandom(
    data = df_cues,
    aes(x = exp_abs, y = brain_mean, colour = cue),
    method = "quasirandom", # offset ∝ local density → fans out the dense stacks
    orientation = "x", # group by expectancy value, spread horizontally
    width = 0.30, # spread in exp_abs units (0.30 ≈ 3 percentage points)
    size = 2,
    alpha = 0.8
  ) +
  facet_grid(
    cue ~ factor(block, levels = names(block_labels_named)),
    labeller = labeller(block = as_labeller(block_labels_named))
  ) +
  scale_colour_manual(values = cue_colors_val, guide = "none") +
  exp_x_scale +
  guides(alpha = "none", linewidth = "none") +
  labs(
    y = "Anterior insula \u03b2",
    title = "Expectancy \u2192 anterior insula: conditional effect by cue \u00d7 block",
    subtitle = "Population-level line \u00b1 95% CI (m_abs_btw)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.subtitle = element_text(size = 14, colour = "grey50")
  )

ggsave(
  file.path(PLOT_DIR, "ai_exp_m_abs_btw_by_cue_block.png"),
  plot = plot_b_exp,
  width = 20,
  height = 14,
  dpi = 150
)

#
fit_ai_val <- brms::brm(
  # fit_amy_val <- brms::brm(
  brms::bf(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      val_wth + val_btw + (1 + val_wth | participant)
  ),
  data = df,
  chains = 4,
  cores = 4,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr"
)

fit_ai_ar <- brms::brm(
  # fit_amy_ar <- brms::brm(
  brms::bf(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      ar_wth + ar_btw + (1 + ar_wth | participant)
  ),
  data = df,
  chains = 4,
  cores = 4,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr"
)

fit_ai_exp <- brms::brm(
  # fit_amy_exp <- brms::brm(
  brms::bf(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      exp_wth + exp_btw + (1 + exp_wth | participant)
  ),
  data = df,
  chains = 4,
  cores = 4,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr"
)

# fit_amy_val_by_block <- brms::brm(
#   brms::bf(
#     brain_mean | se(brain_sd, sigma = TRUE) ~
#       val_wth + val_btw + (1 + val_wth | participant)
#   ),
#   data = df,
#   chains = 4,
#   cores = 4,
#   backend = "cmdstanr"
# )

# fit_amy_ar_by_block <- brms::brm(
#   brms::bf(
#     brain_mean | se(brain_sd, sigma = TRUE) ~
#       ar_wth + ar_btw + (1 + ar_wth | participant)
#   ),
#   data = df,
#   chains = 4,
#   cores = 4,
#   backend = "cmdstanr"
# )

fit_ai_val_by_block <- brms::brm(
  # fit_amy_val <- brms::brm(
  brms::bf(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      val_wth * block + val_btw + (1 + val_wth | participant)
  ),
  data = df,
  chains = 4,
  cores = 4,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr"
)

fit_ai_ar_by_block <- brms::brm(
  # fit_amy_ar <- brms::brm(
  brms::bf(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      ar_wth * block + ar_btw + (1 + ar_wth | participant)
  ),
  data = df,
  chains = 4,
  cores = 4,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr"
)

fit_ai_exp_by_block <- brms::brm(
  # fit_amy_exp_by_block <- brms::brm(
  brms::bf(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      exp_wth * block + exp_btw + (1 + exp_wth | participant)
  ),
  data = df,
  chains = 4,
  cores = 4,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr"
)

library(posterior)
summary(fit_amy_val)
fit_amy_val_summary <- summarise_draws(as_draws_df(fit_amy_val))
summary(fit_amy_ar)
fit_amy_ar_summary <- summarise_draws(as_draws_df(fit_amy_ar))
summary(fit_amy_exp)
fit_amy_exp_summary <- summarise_draws(as_draws_df(fit_amy_exp))
summary(fit_amy_exp_by_block)
fit_amy_exp_by_block_summary <- summarise_draws(as_draws_df(
  fit_amy_exp_by_block
))

summary(fit_ai_val)
fit_ai_val_summary <- summarise_draws(as_draws_df(fit_ai_val))
summary(fit_ai_ar)
fit_ai_ar_summary <- summarise_draws(as_draws_df(fit_ai_ar))
summary(fit_ai_exp)
fit_ai_exp_summary <- summarise_draws(as_draws_df(fit_ai_exp))
summary(fit_ai_val_by_block)
fit_ai_val_by_block_summary <- summarise_draws(as_draws_df(
  fit_ai_val_by_block
))
summary(fit_ai_ar_by_block)
fit_ai_ar_by_block_summary <- summarise_draws(as_draws_df(
  fit_ai_ar_by_block
))
summary(fit_ai_exp_by_block)
fit_ai_exp_by_block_summary <- summarise_draws(as_draws_df(
  fit_ai_exp_by_block
))

emmeans::emtrends(fit_amy_exp_by_block, ~block, var = "exp_wth")
emm <- emmeans::emtrends(fit_amy_exp_by_block, ~block, var = "exp_wth")
pairs(emm) # acq_1 vs acq_2, acq_1 vs ext, etc., with HPDs

emmeans::emtrends(fit_ai_val_by_block, ~block, var = "val_wth")
emm <- emmeans::emtrends(fit_ai_val_by_block, ~block, var = "val_wth")
pairs(emm) # acq_1 vs acq_2, acq_1 vs ext, etc., with HPDs

emmeans::emtrends(fit_ai_exp_by_block, ~block, var = "exp_wth")
emm <- emmeans::emtrends(fit_ai_exp_by_block, ~block, var = "exp_wth")
pairs(emm) # acq_1 vs acq_2, acq_1 vs ext, etc., with HPDs

emmeans::emtrends(fit_ai_exp_by_block, ~block, var = "exp_wth")
emm <- emmeans::emtrends(fit_ai_exp_by_block, ~block, var = "exp_wth")
pairs(emm) # acq_1 vs acq_2, acq_1 vs ext, etc., with HPDs


# brms::posterior_summary(fit_amy_val)
brms::hypothesis(fit_amy_val, "val_wth > 0") # gives posterior prob + an evidence ratio
brms::hypothesis(fit_amy_ar, "ar_wth > 0") # gives posterior prob + an evidence ratio
brms::hypothesis(fit_amy_exp, "exp_wth > 0") # gives posterior prob + an evidence ratio
brms::hypothesis(fit_amy_exp_by_block, "exp_wth > 0") # gives posterior prob + an evidence ratio
brms::hypothesis(fit_ai_val_by_block, "val_wth > 0") # gives posterior prob + an evidence ratio
brms::hypothesis(fit_ai_ar_by_block, "ar_wth > 0") # gives posterior prob + an evidence ratio
brms::hypothesis(fit_ai_exp_by_block, "exp_wth > 0") # gives posterior prob + an evidence ratio


plot(fit_amy_val) # trace + density per parameter
plot(fit_amy_ar) # trace + density per parameter
plot(fit_amy_exp) # trace + density per parameter
plot(fit_amy_exp_by_block) # trace + density per parameter
plot(fit_ai_val_by_block) # trace + density per parameter
plot(fit_ai_ar_by_block) # trace + density per parameter
plot(fit_ai_exp_by_block) # trace + density per parameter


brms::mcmc_plot(fit_amy_val, variable = c("b_val_wth", "b_val_btw")) # interval plot of the two
brms::mcmc_plot(fit_amy_ar, variable = c("b_ar_wth", "b_ar_btw")) # interval plot of the two
brms::mcmc_plot(fit_amy_exp, variable = c("b_exp_wth", "b_exp_btw")) # interval plot of the two
brms::mcmc_plot(fit_amy_exp_by_block, variable = c("b_exp_wth", "b_exp_btw")) # interval plot of the two
brms::mcmc_plot(fit_ai_val_by_block, variable = c("b_val_wth", "b_val_btw")) # interval plot of the two
brms::mcmc_plot(fit_ai_ar_by_block, variable = c("b_ar_wth", "b_ar_btw")) # interval plot of the two
brms::mcmc_plot(fit_ai_exp_by_block, variable = c("b_exp_wth", "b_exp_btw")) # interval plot of the two


brms::conditional_effects(fit_amy_val) # predicted across valence
brms::conditional_effects(fit_amy_ar) # predicted across arousal
brms::conditional_effects(fit_amy_exp) # predicted across expectancy
brms::conditional_effects(fit_ai_val) # predicted across valence
brms::conditional_effects(fit_ai_ar) # predicted across arousal
brms::conditional_effects(fit_ai_exp) # predicted across expectancy
brms::conditional_effects(fit_amy_exp_by_block) # predicted across expectancy
brms::conditional_effects(fit_ai_val_by_block) # predicted across expectancy
brms::conditional_effects(fit_ai_ar_by_block) # predicted across expectancy
brms::conditional_effects(fit_ai_exp_by_block) # predicted across expectancy

brms::ranef(fit_amy_val)$participant[,, "val_wth"] # each person's deviation from the average slope

brms::ranef(fit_amy_exp)$participant[,, "exp_wth"] # each person's deviation from the average slope

library(tidybayes)

slope_draws_ai <- fit_val %>%
  spread_draws(b_val_wth, r_participant[participant, term]) %>%
  filter(term == "val_wth") %>%
  mutate(slope = b_val_wth + r_participant) # fixed + random = total slope

slope_draws_ai %>%
  ggplot(aes(x = slope, y = reorder(factor(participant), slope))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  stat_halfeye(.width = c(.66, .95), normalize = "xy", fill = "#4477AA") +
  labs(
    x = "Valence → anterior insula slope",
    y = "Participant",
    title = "Per-participant valence coupling (anterior insula)"
  ) +
  theme_minimal()

slope_draws_ai <- fit_ar %>%
  spread_draws(b_ar_wth, r_participant[participant, term]) %>%
  filter(term == "ar_wth") %>%
  mutate(slope = b_ar_wth + r_participant) # fixed + random = total slope

slope_draws_ai %>%
  ggplot(aes(x = slope, y = reorder(factor(participant), slope))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  stat_halfeye(.width = c(.66, .95), normalize = "xy", fill = "#4477AA") +
  labs(
    x = "Arousal → anterior insula slope",
    y = "Participant",
    title = "Per-participant arousal coupling (anterior insula)"
  ) +
  theme_minimal()

slope_draws_ai <- fit_exp %>%
  spread_draws(b_exp_wth, r_participant[participant, term]) %>%
  filter(term == "exp_wth") %>%
  mutate(slope = b_exp_wth + r_participant) # fixed + random = total slope

slope_draws_ai %>%
  ggplot(aes(x = slope, y = reorder(factor(participant), slope))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  stat_halfeye(.width = c(.66, .95), normalize = "xy", fill = "#4477AA") +
  labs(
    x = "Expectancy → anterior insula slope",
    y = "Participant",
    title = "Per-participant expectancy coupling (anterior insula)"
  ) +
  theme_minimal()

slope_draws <- fit_amy_exp %>%
  spread_draws(b_exp_wth, r_participant[participant, term]) %>%
  filter(term == "exp_wth") %>%
  mutate(slope = b_exp_wth + r_participant) # fixed + random = total slope

slope_draws %>%
  ggplot(aes(x = slope, y = reorder(factor(participant), slope))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  stat_halfeye(.width = c(.66, .95), normalize = "xy", fill = "#4477AA") +
  labs(
    x = "Expectancy → amygdala slope",
    y = "Participant",
    title = "Per-participant expectancy coupling (amygdala)"
  ) +
  theme_minimal()

df_amy <- df

library(modelr)

fan <- df_amy %>% # your amygdala-joined frame
  group_by(participant) %>%
  data_grid(exp_wth = seq_range(exp_wth, n = 30)) %>% # each person's own range
  mutate(exp_btw = mean(df_amy$exp_btw), brain_sd = 1) %>% # brain_sd is a dummy here
  add_epred_draws(fit_amy_exp, ndraws = 200, re_formula = NULL) %>%
  group_by(participant, exp_wth) %>%
  summarise(epred = median(.epred), .groups = "drop")

grp <- tibble(
  exp_wth = seq_range(df_amy$exp_wth, 30),
  exp_btw = mean(df_amy$exp_btw),
  brain_sd = 1
) %>%
  add_epred_draws(fit_amy_exp, re_formula = NA) %>% # population line only
  group_by(exp_wth) %>%
  summarise(epred = median(.epred), .groups = "drop")

ggplot(fan, aes(exp_wth, epred, group = participant)) +
  geom_line(alpha = .35, colour = "#4477AA") +
  geom_line(data = grp, aes(group = NULL), colour = "black", linewidth = 1.2) +
  labs(
    x = "Expectancy (within-person)",
    y = "Predicted amygdala",
    title = "Individual expectancy→amygdala lines, group in black"
  ) +
  theme_minimal()

fan <- df %>% # your amygdala-joined frame
  group_by(participant) %>%
  data_grid(exp_wth = seq_range(exp_wth, n = 30)) %>% # each person's own range
  mutate(exp_btw = mean(df$exp_btw), brain_sd = 1) %>% # brain_sd is a dummy here
  add_epred_draws(fit_exp, ndraws = 200, re_formula = NULL) %>%
  group_by(participant, exp_wth) %>%
  summarise(epred = median(.epred), .groups = "drop")

grp <- tibble(
  exp_wth = seq_range(df$exp_wth, 30),
  exp_btw = mean(df$exp_btw),
  brain_sd = 1
) %>%
  add_epred_draws(fit_exp, re_formula = NA) %>% # population line only
  group_by(exp_wth) %>%
  summarise(epred = median(.epred), .groups = "drop")

ggplot(fan, aes(exp_wth, epred, group = participant)) +
  geom_line(alpha = .35, colour = "#4477AA") +
  geom_line(data = grp, aes(group = NULL), colour = "black", linewidth = 1.2) +
  labs(
    x = "Expectancy (within-person)",
    y = "Predicted anterior insula",
    title = "Individual expectancy→anterior insula lines, group in black"
  ) +
  theme_minimal()

re <- coef(fit_amy_exp)$participant
tibble(
  participant = rownames(re),
  baseline = re[, "Estimate", "Intercept"],
  slope = re[, "Estimate", "exp_wth"]
) %>%
  ggplot(aes(baseline, slope)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "#4477AA") +
  labs(x = "Per-participant amygdala baseline", y = "Expectancy slope") +
  theme_minimal()

re <- coef(fit_exp)$participant
tibble(
  participant = rownames(re),
  baseline = re[, "Estimate", "Intercept"],
  slope = re[, "Estimate", "exp_wth"]
) %>%
  ggplot(aes(baseline, slope)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "#4477AA") +
  labs(x = "Per-participant anterior insula baseline", y = "Expectancy slope") +
  theme_minimal()

cor_row <- function(fit, term, label) {
  par <- paste0("cor_participant__Intercept__", term)
  brms::posterior_summary(fit, variable = par) |>
    as_tibble() |>
    mutate(model = label, .before = 1)
}

bind_rows(
  cor_row(fit_exp, "exp_wth", "Insula  × expectancy"),
  cor_row(fit_amy_exp, "exp_wth", "Amygdala × expectancy")
)

# ---------------------------------------------------------------------------
# Conditional effect: does val_wth predict anterior insula equally across cues?
# Per-participant regression lines from fit_ai_val_by_block, dots coloured by
# cue (CS+, GS1, GS2, GS3), faceted by block.
# ---------------------------------------------------------------------------

# Build per-participant prediction grid (one line per participant × block)
pred_grid_val <- df |>
  filter(cue %in% c("CS+", "GS1", "GS2", "GS3")) |>
  group_by(participant, block) |>
  summarise(
    val_wth = list(seq(min(val_wth), max(val_wth), length.out = 20)),
    val_btw = mean(val_btw),
    brain_sd = mean(brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(val_wth)

# Median posterior prediction line per participant × block
fan_val <- pred_grid_val |>
  add_epred_draws(fit_ai_val_by_block, ndraws = 100, re_formula = NULL) |>
  group_by(participant, block, val_wth) |>
  summarise(epred = median(.epred), .groups = "drop")

# Build per-participant prediction grid (one line per participant × block)
pred_grid_ar <- df |>
  filter(cue %in% c("CS+", "GS1", "GS2", "GS3")) |>
  group_by(participant, block) |>
  summarise(
    ar_wth = list(seq(min(ar_wth), max(ar_wth), length.out = 20)),
    ar_btw = mean(ar_btw),
    brain_sd = mean(brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(ar_wth)

# Median posterior prediction line per participant × block
fan_ar <- pred_grid_ar |>
  add_epred_draws(fit_ai_ar_by_block, ndraws = 100, re_formula = NULL) |>
  group_by(participant, block, ar_wth) |>
  summarise(epred = median(.epred), .groups = "drop")


# Build per-participant prediction grid (one line per participant × block)
pred_grid_exp <- df |>
  filter(cue %in% c("CS+", "GS1", "GS2", "GS3")) |>
  group_by(participant, block) |>
  summarise(
    exp_wth = list(seq(min(exp_wth), max(exp_wth), length.out = 20)),
    exp_btw = mean(exp_btw),
    brain_sd = mean(brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(exp_wth)

# Median posterior prediction line per participant × block
fan_exp <- pred_grid_exp |>
  add_epred_draws(fit_ai_exp_by_block, ndraws = 100, re_formula = NULL) |>
  group_by(participant, block, exp_wth) |>
  summarise(epred = median(.epred), .groups = "drop")

block_labels_named <- c(
  hab = "Habituation",
  acq_1 = "Acquisition 1",
  acq_2 = "Acquisition 2",
  ext = "Extinction"
)

cue_colors_val <- c(
  "CS+" = "#E63946",
  "GS1" = "#2DC653",
  "GS2" = "#9B5DE5",
  "GS3" = "#457B9D"
)
cue_colors_ar <- cue_colors_val
cue_colors_exp <- cue_colors_val

ggplot() +
  geom_line(
    data = fan_val,
    aes(x = val_wth, y = epred, group = participant),
    alpha = 0.3,
    linewidth = 0.5,
    colour = "grey40"
  ) +
  geom_point(
    data = df |> filter(cue %in% c("CS+", "GS1", "GS2", "GS3")),
    aes(x = val_wth, y = brain_mean, colour = cue),
    alpha = 0.7,
    size = 1.8
  ) +
  facet_wrap(
    ~ factor(block, levels = names(block_labels_named)),
    labeller = as_labeller(block_labels_named),
    nrow = 2
  ) +
  scale_colour_manual(values = cue_colors_val) +
  labs(
    x = "Valence rating (within-person centered)",
    y = "Anterior insula \u03b2",
    colour = "Cue",
    title = "Does valence predict anterior insula equally across cues?",
    subtitle = "Grey lines: per-participant model slopes (fit_ai_val_by_block) \u00b7 Points: observed data by cue"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 20),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 9, colour = "grey50")
  )

# ---------------------------------------------------------------------------
# Population-level conditional effect (line + 95% CI ribbon), dots by cue
# Plot A: faceted by block
# Plot B: faceted by cue × block
# ---------------------------------------------------------------------------

# Population-level conditional effect: val_wth across its full range, per block
cond_grid_val <- df |>
  group_by(block) |>
  summarise(
    val_wth = list(seq(min(val_wth), max(val_wth), length.out = 50)),
    val_btw = mean(df$val_btw),
    brain_sd = mean(df$brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(val_wth)

cond_line_val <- cond_grid_val |>
  add_epred_draws(fit_ai_val_by_block, re_formula = NA) |>
  group_by(block, val_wth) |>
  summarise(
    epred = median(.epred),
    lo = quantile(.epred, 0.025),
    hi = quantile(.epred, 0.975),
    .groups = "drop"
  )

cond_grid_ar <- df |>
  group_by(block) |>
  summarise(
    ar_wth = list(seq(min(ar_wth), max(ar_wth), length.out = 50)),
    ar_btw = mean(df$ar_btw),
    brain_sd = mean(df$brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(ar_wth)

cond_line_ar <- cond_grid_ar |>
  add_epred_draws(fit_ai_ar_by_block, re_formula = NA) |>
  group_by(block, ar_wth) |>
  summarise(
    epred = median(.epred),
    lo = quantile(.epred, 0.025),
    hi = quantile(.epred, 0.975),
    .groups = "drop"
  )

cond_grid_exp <- df |>
  group_by(block) |>
  summarise(
    exp_wth = list(seq(min(exp_wth), max(exp_wth), length.out = 50)),
    exp_btw = mean(df$exp_btw),
    brain_sd = mean(df$brain_sd),
    .groups = "drop"
  ) |>
  tidyr::unnest(exp_wth)

cond_line_exp <- cond_grid_exp |>
  add_epred_draws(fit_ai_exp_by_block, re_formula = NA) |>
  group_by(block, exp_wth) |>
  summarise(
    epred = median(.epred),
    lo = quantile(.epred, 0.025),
    hi = quantile(.epred, 0.975),
    .groups = "drop"
  )

df_cues <- df |> filter(cue %in% c("CS+", "GS1", "GS2", "GS3"))

# Plot A: participant IDs as labels, population line + ribbon, faceted by block
# -- geom_text version (simpler, more cluttered) --
# ggplot() +
#   geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
#   geom_ribbon(
#     data = cond_line_val,
#     aes(x = val_wth, ymin = lo, ymax = hi),
#     alpha = 0.25,
#     fill = "grey40"
#   ) +
#   geom_line(
#     data = cond_line_val,
#     aes(x = val_wth, y = epred),
#     linewidth = 1,
#     colour = "grey20"
#   ) +
#   geom_text(
#     data = df_cues,
#     aes(x = val_wth, y = brain_mean, label = participant, colour = cue),
#     size = 4,
#     alpha = 0.85
#   ) +
#   facet_wrap(
#     ~ factor(block, levels = names(block_labels_named)),
#     labeller = as_labeller(block_labels_named),
#     nrow = 2
#   ) +
#   scale_colour_manual(values = cue_colors_val) +
#   labs(
#     x = "Valence rating (within-person centered)",
#     y = "Anterior insula \u03b2",
#     colour = "Cue",
#     title = "Valence \u2192 anterior insula: conditional effect by block",
#     subtitle = "Population-level line \u00b1 95% CI (fit_ai_val_by_block) \u00b7 labels = participant ID"
#   ) +
#   theme_minimal(base_size = 20) +
#   theme(
#     strip.text = element_text(face = "bold"),
#     legend.position = "bottom",
#     plot.subtitle = element_text(size = 14, colour = "grey50")
#   )

# -- geom_label_repel version with fading vlines + saved at large size --
library(ggrepel)

# Fading vertical reference lines: darkest/thickest at x=0, decay outward
vlines_df <- tibble(
  x = seq(-5L, 7L, by = 1L),
  dist = abs(x)
) |>
  mutate(
    line_alpha = pmax(0.05, 0.55 * (0.72^dist)),
    line_lwd = pmax(0.10, 0.70 * (0.72^dist))
  )

plot_a <- ggplot() +
  geom_vline(
    data = vlines_df,
    aes(xintercept = x, linewidth = line_lwd, alpha = line_alpha),
    colour = "grey30"
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(
    data = cond_line_val,
    aes(x = val_wth, ymin = lo, ymax = hi),
    alpha = 0.25,
    fill = "grey40"
  ) +
  geom_line(
    data = cond_line_val,
    aes(x = val_wth, y = epred),
    linewidth = 1.2,
    colour = "grey20"
  ) +
  geom_point(
    data = df_cues,
    aes(x = val_wth, y = brain_mean, colour = cue),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_label_repel(
    data = df_cues,
    aes(x = val_wth, y = brain_mean, label = participant, colour = cue),
    fill = alpha("white", 0.8),
    size = 5,
    fontface = "bold",
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    segment.size = 0.7,
    segment.alpha = 0.6,
    max.overlaps = 30,
    force = 3,
    seed = 812,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ factor(block, levels = names(block_labels_named)),
    labeller = as_labeller(block_labels_named),
    nrow = 2
  ) +
  scale_colour_manual(values = cue_colors_val) +
  scale_x_continuous(breaks = seq(-5, 7, by = 1)) +
  guides(alpha = "none", linewidth = "none") +
  labs(
    x = "Valence rating (within-person centered)",
    y = "Anterior insula \u03b2",
    colour = "Cue",
    title = "Valence \u2192 anterior insula: conditional effect by block",
    subtitle = "Population-level line \u00b1 95% CI (fit_ai_val_by_block) \u00b7 labels = participant ID"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 14, colour = "grey50")
  )

plot_a
ggsave(
  file.path(PLOT_DIR, "ai_val_by_block_conditional_effect.png"),
  plot = plot_a,
  width = 20,
  height = 24,
  dpi = 150
)

plot_a_ar <- ggplot() +
  geom_vline(
    data = vlines_df,
    aes(xintercept = x, linewidth = line_lwd, alpha = line_alpha),
    colour = "grey30"
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(
    data = cond_line_ar,
    aes(x = ar_wth, ymin = lo, ymax = hi),
    alpha = 0.25,
    fill = "grey40"
  ) +
  geom_line(
    data = cond_line_ar,
    aes(x = ar_wth, y = epred),
    linewidth = 1.2,
    colour = "grey20"
  ) +
  geom_point(
    data = df_cues,
    aes(x = ar_wth, y = brain_mean, colour = cue),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_label_repel(
    data = df_cues,
    aes(x = ar_wth, y = brain_mean, label = participant, colour = cue),
    fill = alpha("white", 0.8),
    size = 5,
    fontface = "bold",
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    segment.size = 0.7,
    segment.alpha = 0.6,
    max.overlaps = 30,
    force = 3,
    seed = 812,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ factor(block, levels = names(block_labels_named)),
    labeller = as_labeller(block_labels_named),
    nrow = 2
  ) +
  scale_colour_manual(values = cue_colors_ar) +
  scale_x_continuous(breaks = seq(-5, 7, by = 1)) +
  guides(alpha = "none", linewidth = "none") +
  labs(
    x = "Valence rating (within-person centered)",
    y = "Anterior insula \u03b2",
    colour = "Cue",
    title = "Valence \u2192 anterior insula: conditional effect by block",
    subtitle = "Population-level line \u00b1 95% CI (fit_ai_ar_by_block) \u00b7 labels = participant ID"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 14, colour = "grey50")
  )

plot_a_ar
ggsave(
  file.path(PLOT_DIR, "ai_ar_by_block_conditional_effect.png"),
  plot = plot_a_ar,
  width = 20,
  height = 24,
  dpi = 150
)


# Fading vlines on 10-unit grid aligned to the 0-100 expectancy scale.
# Thickest line sits at the grand mean; decay outward in 10-unit steps.
grand_mean_exp <- mean(df$exp_btw)
exp_label_vals <- seq(0, 100, by = 10)
exp_break_vals <- exp_label_vals - grand_mean_exp

vlines_exp <- tibble(
  label = exp_label_vals,
  x = exp_break_vals,
  dist = abs(label - grand_mean_exp) / 10
) |>
  mutate(
    line_alpha = pmax(0.05, 0.55 * (0.72^dist)),
    line_lwd = pmax(0.10, 0.70 * (0.72^dist))
  )

plot_a_exp <- ggplot() +
  geom_vline(
    data = vlines_exp,
    aes(xintercept = x, linewidth = line_lwd, alpha = line_alpha),
    colour = "grey30"
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(
    data = cond_line_exp,
    aes(x = exp_wth, ymin = lo, ymax = hi),
    alpha = 0.25,
    fill = "grey40"
  ) +
  geom_line(
    data = cond_line_exp,
    aes(x = exp_wth, y = epred),
    linewidth = 1.2,
    colour = "grey20"
  ) +
  geom_point(
    data = df_cues,
    aes(x = exp_wth, y = brain_mean, colour = cue),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_label_repel(
    data = df_cues,
    aes(x = exp_wth, y = brain_mean, label = participant, colour = cue),
    fill = alpha("white", 0.8),
    size = 5,
    fontface = "bold",
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    segment.size = 0.7,
    segment.alpha = 0.6,
    max.overlaps = 30,
    force = 3,
    seed = 812,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ factor(block, levels = names(block_labels_named)),
    labeller = as_labeller(block_labels_named),
    nrow = 2
  ) +
  scale_colour_manual(values = cue_colors_exp) +
  # x-axis labels re-referenced to 0-100 by adding back the grand mean
  scale_x_continuous(breaks = exp_break_vals, labels = exp_label_vals) +
  guides(alpha = "none", linewidth = "none") +
  labs(
    x = "Expectancy rating (0\u2013100, grand-mean re-referenced)",
    y = "Anterior insula \u03b2",
    colour = "Cue",
    title = "Expectancy \u2192 anterior insula: conditional effect by block",
    subtitle = "Population-level line \u00b1 95% CI (fit_ai_exp_by_block) \u00b7 labels = participant ID"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 14, colour = "grey50")
  )

plot_a_exp
ggsave(
  file.path(PLOT_DIR, "ai_exp_by_block_conditional_effect.png"),
  plot = plot_a_exp,
  width = 20,
  height = 14,
  dpi = 150
)


# Plot B: same population line replicated per cue row so each cue × block
# panel shows how well the model line describes that cue's observations
cond_line_by_cue <- cond_line_val |>
  tidyr::crossing(
    cue = factor(
      c("CS+", "GS1", "GS2", "GS3"),
      levels = c("CS+", "GS1", "GS2", "GS3")
    )
  )

# -- geom_text version (simpler, more cluttered) --
# ggplot() +
#   geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
#   geom_ribbon(
#     data = cond_line_by_cue,
#     aes(x = val_wth, ymin = lo, ymax = hi),
#     alpha = 0.25,
#     fill = "grey40"
#   ) +
#   geom_line(
#     data = cond_line_by_cue,
#     aes(x = val_wth, y = epred),
#     linewidth = 0.9,
#     colour = "grey20"
#   ) +
#   geom_text(
#     data = df_cues,
#     aes(x = val_wth, y = brain_mean, label = participant, colour = cue),
#     size = 4,
#     alpha = 0.85
#   ) +
#   facet_grid(
#     cue ~ factor(block, levels = names(block_labels_named)),
#     labeller = labeller(block = as_labeller(block_labels_named))
#   ) +
#   scale_colour_manual(values = cue_colors_val, guide = "none") +
#   labs(
#     x = "Valence rating (within-person centered)",
#     y = "Anterior insula \u03b2",
#     title = "Valence \u2192 anterior insula: conditional effect by cue \u00d7 block",
#     subtitle = "Population-level line \u00b1 95% CI (fit_ai_val_by_block) \u00b7 labels = participant ID"
#   ) +
#   theme_minimal(base_size = 20) +
#   theme(
#     strip.text = element_text(face = "bold"),
#     legend.position = "none",
#     plot.subtitle = element_text(size = 14, colour = "grey50")
#   )

# -- geom_label_repel version with fading vlines + saved at large size --
plot_b <- ggplot() +
  geom_vline(
    data = vlines_df,
    aes(xintercept = x, linewidth = line_lwd, alpha = line_alpha),
    colour = "grey30"
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(
    data = cond_line_by_cue,
    aes(x = val_wth, ymin = lo, ymax = hi),
    alpha = 0.25,
    fill = "grey40"
  ) +
  geom_line(
    data = cond_line_by_cue,
    aes(x = val_wth, y = epred),
    linewidth = 0.9,
    colour = "grey20"
  ) +
  geom_point(
    data = df_cues,
    aes(x = val_wth, y = brain_mean, colour = cue),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_label_repel(
    data = df_cues,
    aes(x = val_wth, y = brain_mean, label = participant, colour = cue),
    fill = alpha("white", 0.8),
    size = 5,
    fontface = "bold",
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    segment.size = 0.7,
    segment.alpha = 0.6,
    max.overlaps = Inf,
    force = 3,
    seed = 812,
    show.legend = FALSE
  ) +
  facet_grid(
    cue ~ factor(block, levels = names(block_labels_named)),
    labeller = labeller(block = as_labeller(block_labels_named))
  ) +
  scale_colour_manual(values = cue_colors_val, guide = "none") +
  scale_x_continuous(breaks = seq(-5, 7, by = 1)) +
  guides(alpha = "none", linewidth = "none") +
  labs(
    x = "Valence rating (within-person centered)",
    y = "Anterior insula \u03b2",
    title = "Valence \u2192 anterior insula: conditional effect by cue \u00d7 block",
    subtitle = "Population-level line \u00b1 95% CI (fit_ai_val_by_block) \u00b7 labels = participant ID"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.subtitle = element_text(size = 14, colour = "grey50")
  )

plot_b
ggsave(
  file.path(PLOT_DIR, "ai_val_by_cue_block_conditional_effect.png"),
  plot = plot_b,
  width = 20,
  height = 14,
  dpi = 150
)

# Plot C: participant trajectory across blocks for each cue
# x = block (ordered), y = brain_mean, lines connect same participant,
# ID labeled at extinction — lets you track each person's arc per cue
df_traj <- df_cues |>
  mutate(block = factor(block, levels = names(block_labels_named)))

df_traj_label <- df_traj |>
  filter(block == "ext")

ggplot(df_traj, aes(x = block, y = brain_mean, group = participant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(alpha = 0.4, colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.5, colour = "grey50") +
  geom_text(
    data = df_traj_label,
    aes(label = participant),
    hjust = -0.2,
    size = 3.5,
    colour = "grey20"
  ) +
  scale_x_discrete(labels = block_labels_named) +
  facet_wrap(~cue, nrow = 2) +
  labs(
    x = "Block",
    y = "Anterior insula \u03b2",
    title = "Participant trajectories: anterior insula by cue across blocks",
    subtitle = "Each line = one participant \u00b7 ID labeled at extinction"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.subtitle = element_text(size = 14, colour = "grey50"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# testing using below

# GP by trial
v4_38_log_lik <- read_chains("model038", "V4", variables = "log_lik")
v5_38_log_lik <- read_chains("model038", "V5", variables = "log_lik")
v6_38_log_lik <- read_chains("model038", "V6", variables = "log_lik")
te_38_log_lik <- read_chains("model038", "TE", variables = "log_lik")
tpj_38_log_lik <- read_chains("model038", "TPJ", variables = "log_lik")
acc_38_log_lik <- read_chains("model038", "ACC", variables = "log_lik")
na_38_log_lik <- read_chains("model038", "NA", variables = "log_lik")
amy_38_log_lik <- read_chains("model038", "AMY", variables = "log_lik")
hip_38_log_lik <- read_chains("model038", "HIP", variables = "log_lik")
ai_38_log_lik <- read_chains("model038", "ant_insula", variables = "log_lik")
ofc_38_log_lik <- read_chains("model038", "OFC", variables = "log_lik")


# by cue by block by participant
v4_42_log_lik <- read_chains("model042", "V4", variables = "log_lik")
v5_42_log_lik <- read_chains("model042", "V5", variables = "log_lik")
v6_42_log_lik <- read_chains("model042", "V6", variables = "log_lik")
te_42_log_lik <- read_chains("model042", "TE", variables = "log_lik")
tpj_42_log_lik <- read_chains("model042", "TPJ", variables = "log_lik")
acc_42_log_lik <- read_chains("model042", "ACC", variables = "log_lik")
na_42_log_lik <- read_chains("model042", "NA", variables = "log_lik")
amy_42_log_lik <- read_chains("model042", "AMY", variables = "log_lik")
hip_42_log_lik <- read_chains("model042", "HIP", variables = "log_lik")
ai_42_log_lik <- read_chains("model042", "ant_insula", variables = "log_lik")
ofc_42_log_lik <- read_chains("model042", "OFC", variables = "log_lik")


#loo
v4_38_loo <- compute_loo(v4_38_log_lik)
v5_38_loo <- compute_loo(v5_38_log_lik)
v6_38_loo <- compute_loo(v6_38_log_lik)
te_38_loo <- compute_loo(te_38_log_lik)
tpj_38_loo <- compute_loo(tpj_38_log_lik)
acc_38_loo <- compute_loo(acc_38_log_lik)
na_38_loo <- compute_loo(na_38_log_lik)
amy_38_loo <- compute_loo(amy_38_log_lik)
hip_38_loo <- compute_loo(hip_38_log_lik)
ai_38_loo <- compute_loo(ai_38_log_lik)
ofc_38_loo <- compute_loo(ofc_38_log_lik)

v4_42_loo <- compute_loo(v4_42_log_lik)
v5_42_loo <- compute_loo(v5_42_log_lik)
v6_42_loo <- compute_loo(v6_42_log_lik)
te_42_loo <- compute_loo(te_42_log_lik)
tpj_42_loo <- compute_loo(tpj_42_log_lik)
acc_42_loo <- compute_loo(acc_42_log_lik)
na_42_loo <- compute_loo(na_42_log_lik)
amy_42_loo <- compute_loo(amy_42_log_lik)
hip_42_loo <- compute_loo(hip_42_log_lik)
ai_42_loo <- compute_loo(ai_42_log_lik)
ofc_42_loo <- compute_loo(ofc_42_log_lik)


# roi_label <- "V4"
# roi_label <- "V5"
# roi_label <- "V6"
roi_label <- "TE"
# roi_label <- "TPJ"
# roi_label <- "ACC"
# roi_label <- "Amygdala"
# roi_label <- "Hippocampus"
# roi_label <- "Anterior Insula"
# roi_label <- "OFC"

# mu <- read_chains("model042", "V4")
# trial <- read_trial_betas_chains("model038", "V4")
# mu <- read_chains("model042", "V5")
# trial <- read_trial_betas_chains("model038", "V5")
# mu <- read_chains("model042", "V6")
# trial <- read_trial_betas_chains("model038", "V6")
mu <- read_chains("model042", "TE")
trial <- read_trial_betas_chains("model038", "TE")
# mu <- read_chains("model042", "TPJ")
# trial <- read_trial_betas_chains("model038", "TPJ")
# mu <- read_chains("model042", "ACC")
# trial <- read_trial_betas_chains("model038", "ACC")
# mu <- read_chains("model042", "AMY")
# trial <- read_trial_betas_chains("model038", "Amygdala")
# mu <- read_chains("model042", "HIP")
# trial <- read_trial_betas_chains("model038", "Hippocampus")
# mu <- read_chains("model042", "ant_insula")
# trial <- read_trial_betas_chains("model038", "ant_insula")
# mu <- read_chains("model042", "OFC")
# trial <- read_trial_betas_chains("model038", "Orbital_Frontal_Complex")

# Headline summary plot
combined <- plot_roi_summary(trial, mu, roi_label)
print(combined)
save_roi_summary(combined, roi_label)


######
model042_ai_all <- read_chains(
  model_id = "model042",
  roi_id = "ant_insula",
  variables = "betas"
)

model042_ai_all_summary <- summarise_draws(
  model042_ai_all
)


roi_label <- "Anterior Insula"
# roi_label <- "Amygdala"

# Load the two model families for this ROI
# ai_mu <- read_chains("model039", "ant_insula")
ai_mu <- read_chains("model039", "ant_insula")
ai_trial <- read_trial_betas_chains("model038", "ant_insula")

# ai_mu <- read_chains("model039", "ant_insula")
# ai_mu <- read_chains("model042", "AMY")
# ai_trial <- read_trial_betas_chains("model038", "Amygdala")

# Headline summary plot
combined <- plot_roi_summary(ai_trial, ai_mu, roi_label)
print(combined)


ai_38_log_lik_mat <- as_draws_matrix(ai_38_log_lik)

# by cue by block by participant
ai_39_log_lik <- read_chains("model039", "ant_insula", variables = "log_lik")
ai_39_log_lik_mat <- as_draws_matrix(ai_39_log_lik)

# model 39 with t-distribution regularization of mu_betas
ai_42_log_lik <- read_chains("model042", "ant_insula", variables = "log_lik")
ai_42_log_lik_mat <- as_draws_matrix(ai_42_log_lik)

ai_38_loo <- compute_loo(ai_38_log_lik)
ai_39_loo <- compute_loo(ai_39_log_lik)
ai_42_loo <- compute_loo(ai_42_log_lik)

ai_38_loo
ai_39_loo
ai_42_loo

loo::loo_compare(
  ai_38_loo,
  ai_39_loo,
  ai_42_loo
)

amy_38_log_lik <- read_chains("model038", "Amygdala", variables = "log_lik")

amy_39_log_lik <- read_chains("model039", "AMY", variables = "log_lik")

amy_42_log_lik <- read_chains("model042", "AMY", variables = "log_lik")

amy_38_loo <- compute_loo(amy_38_log_lik)
amy_39_loo <- compute_loo(amy_39_log_lik)
amy_42_loo <- compute_loo(amy_42_log_lik)

amy_38_loo
amy_39_loo
amy_42_loo

loo::loo_compare(
  amy_38_loo,
  amy_39_loo,
  amy_42_loo
)

# Scalars come in cleanly:
stan_list <- jsonlite::fromJSON(
  "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/fmri_stan_list_single_trial_Ant_Ins_by_cue_pb03_res.json",
  simplifyVector = FALSE
)
censor_vec <- stan_list$usable_bold_indices_one_is_true %>% unlist()
par_vec <- stan_list$par %>% unlist()

par_vec <- par_vec[censor_vec == 1]

# log_lik is n_samples × n_obs (the volume-level matrix Stan returns)
ai_38_log_lik_par <- matrix(
  NA,
  nrow = nrow(ai_38_log_lik_mat),
  ncol = stan_list$n_par
)

for (p in 1:stan_list$n_par) {
  cols <- which(par_vec == p) # which volumes belong to p
  ai_38_log_lik_par[, p] <- rowSums(ai_38_log_lik_mat[, cols])
}

ai_38_loo_lopo <- loo::loo(ai_38_log_lik_par) # treats each column as one "observation"


ai_39_log_lik_par <- matrix(
  NA,
  nrow = nrow(ai_39_log_lik_mat),
  ncol = stan_list$n_par
)

for (p in 1:stan_list$n_par) {
  cols <- which(par_vec == p) # which volumes belong to p
  ai_39_log_lik_par[, p] <- rowSums(ai_39_log_lik_mat[, cols])
}


ai_39_loo_lopo <- loo::loo(ai_39_log_lik_par) # treats each column as one "observation"

ai_38_loo_lopo
ai_39_loo_lopo

loo_compare(
  ai_38_loo_lopo,
  ai_39_loo_lopo
)


model042_ai_all <- read_chains(
  model_id = "model042",
  roi_id = "ant_insula",
  variables = beta
)

model042_ai_all_summary <- summarise_draws(
  model042_ai_all
)

model042_amy_all <- read_chains(
  model_id = "model042",
  roi_id = "AMY",
  variables = NULL
)

model042_amy_all_summary <- summarise_draws(
  model042_amy_all
)

model042_all <- read_chains(
  model_id = "model042",
  roi_id = "HIP",
  variables = NULL
)

model042_all_summary_2 <- summarise_draws(
  model042_all
)
# ============================================================================
#  9. Example driver — anterior insula end-to-end
# ============================================================================
# Comment out / delete this section if you want to source this file purely
# as a library.

if (FALSE) {
  roi_label <- "V1F"
  v1f_mu <- read_chains("model042", "V1F")
  v1f_mu_df <- tidy_mu_betas(v1f_mu)
  posterior_plot <- plot_phase_cue_posterior(v1f_mu_df, roi_label)
  v1f_mu_full <- read_chains(
    "model042",
    "V1F",
    variables = c("mu_betas", "betas")
  )
  v1f_pp_cue_betas <- tidy_participant_betas(v1f_mu_full) %>%
    filter(j < 18)
  print(plot_participant_cue_betas(v1f_pp_cue_betas))
  v1f_pp <- compute_participant_metrics(v1f_mu_full)
  print(plot_participant_metrics(v1f_pp, free_y = T))

  roi_label <- "Anterior Insula"
  # roi_label <- "Amygdala"

  # Load the two model families for this ROI
  # ai_mu <- read_chains("model039", "ant_insula")
  ai_mu <- read_chains("model042", "ant_insula")
  ai_trial <- read_trial_betas_chains("model038", "ant_insula")

  # ai_mu <- read_chains("model039", "ant_insula")
  # ai_mu <- read_chains("model042", "AMY")
  # ai_trial <- read_trial_betas_chains("model038", "Amygdala")

  # Headline summary plot
  combined <- plot_roi_summary(ai_trial, ai_mu, roi_label)
  print(combined)
  # save_roi_summary(combined, roi_label)

  roi_label <- "Anterior Insula"
  # roi_label <- "Amygdala"

  # Load the two model families for this ROI
  # ai_mu <- read_chains("model039", "ant_insula")
  ai_mu <- read_chains("model039", "ant_insula")
  ai_trial <- read_trial_betas_chains("model038", "ant_insula")

  # ai_mu <- read_chains("model039", "ant_insula")
  # ai_mu <- read_chains("model042", "AMY")
  # ai_trial <- read_trial_betas_chains("model038", "Amygdala")

  # Headline summary plot
  combined <- plot_roi_summary(ai_trial, ai_mu, roi_label)
  print(combined)
  # save_roi_summary(combined, roi_label)

  # Metrics panel (generalization / selectivity / sharpening)
  ai_metrics <- compute_cue_metrics(tidy_trial_betas(ai_trial))
  metrics_panel <- plot_metrics_timeseries(ai_metrics)
  ts_panel <- plot_trial_timeseries(
    tidy_trial_betas(ai_trial),
    roi_label,
    ylim = roi_ylim[[roi_label]]
  )
  print(ts_panel / metrics_panel)

  # Per-participant metrics from model039's per-subject betas. Need to
  # reload the draws with `betas` included — the default loader only
  # pulls `mu_betas` for speed.
  ai_mu_full <- read_chains(
    "model042",
    "ant_insula",
    variables = c("mu_betas", "betas")
  )
  ai_pp <- compute_participant_metrics(ai_mu_full)
  print(plot_participant_metrics(ai_pp, free_y = T))

  # Loop over every model039 ROI and save the combined plot
  rois_to_plot <- chain_registry %>%
    # filter(model_id == "model039") %>%
    filter(model_id == "model042") %>%
    pull(roi_id)

  ai_mu_full <- read_chains(
    "model042",
    "ant_insula",
    variables = c("mu_betas", "betas")
  )
  ai_pp_metrics <- compute_participant_metrics(ai_mu_full)
  ai_pp_cue_betas <- tidy_participant_betas(ai_mu_full) %>%
    filter(j < 18)

  print(plot_participant_metrics(ai_pp_metrics, show_labels = T, free_y = T))

  ggsave(
    filename = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/anterior_insula_par_cue_metrics.png",
    width = 20,
    height = 7,
    dpi = 100
  )

  print(plot_participant_cue_betas(ai_pp_cue_betas))

  ggsave(
    filename = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/anterior_insula_par_cue.png",
    width = 20,
    height = 7,
    dpi = 100
  )

  amy_mu_full <- read_chains(
    "model042",
    "AMY",
    variables = c("mu_betas", "betas")
  )
  amy_pp_metrics <- compute_participant_metrics(amy_mu_full)
  amy_pp_cue_betas <- tidy_participant_betas(amy_mu_full) %>%
    filter(j < 18)

  print(plot_participant_metrics(amy_pp_metrics))

  ggsave(
    filename = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/amygdala_par_cue_metrics.png",
    width = 20,
    height = 7,
    dpi = 100
  )

  print(plot_participant_cue_betas(amy_pp_cue_betas))

  ggsave(
    filename = "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/amygdala_par_cue.png",
    width = 20,
    height = 7,
    dpi = 100
  )

  for (r in rois_to_plot) {
    label <- chain_registry %>%
      # filter(model_id == "model039", roi_id == r) %>%
      filter(model_id == "model042", roi_id == r) %>%
      pull(roi_label)
    has_m038 <- any(
      chain_registry$model_id == "model038" &
        chain_registry$roi_id == r
    )
    if (!has_m038) {
      next
    }

    message("Plotting ", label, "...")
    # mu_d <- read_chains("model039", r)
    mu_d <- read_chains("model042", r)
    trial_d <- read_trial_betas_chains("model038", r)
    save_roi_summary(plot_roi_summary(trial_d, mu_d, label), label)
  }
}
