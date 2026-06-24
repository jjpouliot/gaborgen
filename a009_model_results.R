# ============================================================================
#  Gaborgen24 fMRI — Model results: infrastructure + behavioral model
# ----------------------------------------------------------------------------
#  Consolidated from:
#    a008_model_loading_plots.R  (fMRI chain loading + Stan plotting functions)
#    a008_separate_plots_5.R     (m_final_less visualisations)
#
#  Layout
#  ------
#    1.  Libraries
#    2.  Configuration         (paths, palettes, intervals, y-limits)
#    3.  Beta metadata         (mu_betas[j] -> phase/cue;
#                               betas[trial_idx] -> cue/block/trial_per_cue)
#    4.  Chain registry        (one row per model x ROI)
#    5.  Path / loader funs    (chain_paths, chain_paths_lookup,
#                               read_chains, read_trial_betas_chains)
#    6.  Tidiers               (tidy_mu_betas, tidy_trial_betas,
#                               compute_cue_metrics, tidy_participant_betas,
#                               compute_participant_metrics)
#    7.  Plotters              (plot_phase_cue_posterior,
#                               plot_trial_timeseries,
#                               plot_metrics_timeseries,
#                               plot_participant_metrics,
#                               plot_participant_cue_betas,
#                               plot_roi_summary, save_roi_summary)
#    8.  LOO helper            (compute_loo)
#    9.  Behavioral data       (ratings loading + ratings2 — run once)
#   10.  build_roi_df()        generalized data-prep for any model042 ROI
#   11.  fit_m_final_less()    generalized brm fit for any df_fit
#   12.  Behavioral plot fns   (plot_partial_effects, plot_surface,
#                               plot_residuals, plot_loo_r2,
#                               make_m_final_less_plots, save_m_final_less_plots)
#   13.  Driver                (example calls; easy to loop over all ROIs)
# ============================================================================

# ---- 1. Libraries ----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(patchwork)
library(posterior)
library(loo)
library(tidybayes)
library(brms)


# ---- 2. Configuration ------------------------------------------------------

DATA_DIR <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/stan_chains"
PLOT_DIR <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/misc"
parent_folder <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI"

# Cue palettes — pretty labels for mu_betas plots, lowercase for
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
  "V5" = c(-0.1, 0.3),
  "TE" = c(-0.11, 0.2),
  "TPJ" = c(-0.05, 0.2),
  "ACC" = c(-0.1, 0.25),
  "Amygdala" = c(-0.075, 0.33),
  "Hippocampus" = c(-0.05, 0.15),
  "Anterior Insula" = c(-0.05, 0.33)
)


# ---- 3. Beta metadata ------------------------------------------------------

#  mu_betas[roi, j]  for j = 1..17  (model039 / model042)
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
    )
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
#  177-191  shock           : 15 trials
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
  "101",
  "102",
  "103",
  "106",
  "107",
  "108",
  "109",
  "113",
  "114",
  "115",
  "116",
  "117",
  "119",
  "120",
  "121",
  "122",
  "123",
  "125",
  "126",
  "127",
  "128",
  "129",
  "131",
  "132",
  "133",
  "134",
  "135",
  "137",
  "138",
  "140",
  "141",
  "143",
  "145",
  "149",
  "150",
  "151",
  "152",
  "153",
  "154",
  "155",
  "158",
  "159",
  "160",
  "161"
)


# ---- 4. Chain registry -----------------------------------------------------
# One row per (model_id, roi_id). `chain_nums` is a comma-separated string
# of chain numbers to load.

chain_registry <- tibble::tribble(
  ~model_id  , ~roi_id                   , ~roi_label          , ~suffix    , ~job_id    , ~chain_nums ,
  # --- model042 (hierarchical: mu_betas + per-subject betas, t-distribution regularization) ---
  "model042" , "V1F"                     , "V1F"               , "pb03_res" , "33726894" , "1,2"       ,
  "model042" , "V4"                      , "V4"                , "pb03_res" , "33218289" , "1,2"       ,
  "model042" , "V5"                      , "V5"                , "pb03_res" , "33233099" , "1,2"       ,
  "model042" , "V6"                      , "V6"                , "pb03_res" , "33233105" , "1,2"       ,
  "model042" , "TE"                      , "TE"                , "pb03_res" , "33233115" , "1,2"       ,
  "model042" , "TPJ"                     , "TPJ"               , "pb03_res" , "33233116" , "1,2"       ,
  "model042" , "ACC"                     , "ACC"               , "pb03_res" , "33233118" , "1,2"       ,
  "model042" , "NA"                      , "Nucleus Accumbens" , "pb03_res" , "33233130" , "1,2"       ,
  "model042" , "AMY"                     , "Amygdala"          , "pb03_res" , "33099488" , "1,2"       ,
  # "model042" , "AMY"                     , "Amygdala"          , "pb03_res_first_to_hab" , "34820009" , "1,2"       ,
  "model042" , "HIP"                     , "Hippocampus"       , "pb03_res" , "33233131" , "1,2"       ,
  "model042" , "ant_insula"              , "Anterior Insula"   , "pb03_res" , "33066991" , "1,2"       ,
  "model042" , "OFC"                     , "OFC"               , "pb03_res" , "33349421" , "1,2"       ,
  "model042" , "VMPFC"                   , "VMPFC"             , "pb03_res" , "34906917" , "1,2"       ,

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
  "model038" , "VMPFC"                   , "VMPFC"             , "pb03_res" , "34907532" , "1,2"       ,
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
  row <- registry |>
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

# Loader for hierarchical fits (model039 / model042). Defaults to reading
# only `mu_betas`. Pass variables = c("mu_betas", "betas") to get
# per-subject betas (needed for build_roi_df).
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

# Loader for single-trial fits (model038).
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
  draws |>
    spread_draws(mu_betas[i, j]) |>
    group_by(.draw, j) |>
    summarise(mu_beta = mean(mu_betas), .groups = "drop") |>
    left_join(meta, by = "j") |>
    filter(!is.na(phase))
}

# Long form for single-trial betas: one row per (draw, roi, trial), with
# cue / block / trial_per_cue attached.
tidy_trial_betas <- function(draws, meta = trial_beta_meta) {
  posterior::as_draws_df(draws) |>
    pivot_longer(
      cols = starts_with("betas["),
      names_to = c(".unused", "roi", "beta_index"),
      names_sep = "\\[|,|\\]",
      values_to = "beta_value"
    ) |>
    select(-.unused) |>
    mutate(
      roi = as.integer(roi),
      beta_index = as.integer(beta_index)
    ) |>
    left_join(meta, by = "beta_index")
}

# Per-draw generalization / selectivity / sharpening over the cue trials.
compute_cue_metrics <- function(trial_long) {
  trial_long |>
    filter(cue != "shock") |>
    group_by(.draw, trial_per_cue) |>
    summarise(
      csp = mean(beta_value[cue == "csp"]),
      gs1 = mean(beta_value[cue == "gs1"]),
      gs2 = mean(beta_value[cue == "gs2"]),
      gs3 = mean(beta_value[cue == "gs3"]),
      .groups = "drop"
    ) |>
    mutate(
      generalization = (csp + gs1 + gs2 + gs3) / 4,
      selectivity = csp - (gs1 + gs2 + gs3) / 3,
      sharpening = csp - gs1 - gs2 + gs3
    ) |>
    pivot_longer(
      c(generalization, selectivity, sharpening),
      names_to = "metric",
      values_to = "value"
    )
}

# Per-participant tidied betas from model039/042: averages over ROI within
# each draw and joins phase/cue labels.
tidy_participant_betas <- function(draws, meta = mu_beta_meta) {
  draws |>
    spread_draws(betas[participant, roi, j]) |>
    group_by(.draw, participant, j) |>
    summarise(beta = mean(betas), .groups = "drop") |>
    left_join(meta, by = "j")
}

# Per-participant, per-phase generalization / selectivity / sharpening
# from model039/042's per-subject `betas[participant, roi, j]`.
compute_participant_metrics <- function(draws, meta = mu_beta_meta) {
  draws |>
    spread_draws(betas[participant, roi, j]) |>
    group_by(.draw, participant, j) |>
    summarise(beta = mean(betas), .groups = "drop") |>
    left_join(meta, by = "j") |>
    filter(cue != "Shock") |>
    select(.draw, participant, phase, cue, beta) |>
    pivot_wider(names_from = cue, values_from = beta) |>
    mutate(
      avg_gen = (`CS+` + GS1 + GS2 + GS3) / 4,
      linear_gen = -1 * ((-3 * `CS+` - 1 * GS1 + 1 * GS2 + 3 * GS3) / 10),
      selectivity = `CS+` - (GS1 + GS2 + GS3) / 3,
      sharpening = `CS+` - GS1 - GS2 + GS3
    ) |>
    select(
      .draw,
      participant,
      phase,
      avg_gen,
      linear_gen,
      selectivity,
      sharpening
    ) |>
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
  p <- mu_long |>
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
  summarised <- trial_long |>
    group_by(.draw, cue, trial_per_cue, beta_index) |>
    summarise(avg_roi_draw = mean(beta_value), .groups = "drop") |>
    group_by(cue, trial_per_cue, beta_index) |>
    summarise(
      median_posterior = median(avg_roi_draw),
      lower = quantile(avg_roi_draw, ci_lower),
      upper = quantile(avg_roi_draw, ci_upper),
      .groups = "drop"
    )

  p <- summarised |>
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

# Generalization / selectivity / sharpening metrics over trials.
plot_metrics_timeseries <- function(metrics_long, breaks = trial_block_breaks) {
  metrics_long |>
    ggplot(aes(x = trial_per_cue, y = value, fill = metric, color = metric)) +
    geom_vline(xintercept = breaks) +
    tidybayes::stat_lineribbon(.width = .341, alpha = 0.25) +
    geom_hline(yintercept = 0) +
    theme_bw(base_size = 20)
}

# Per-participant posteriors of generalization / selectivity / sharpening.
plot_participant_metrics <- function(
  participant_metrics,
  .width = 0.68,
  point_size = 0.3,
  fatten = 1.8,
  alpha = 0.75,
  free_y = FALSE,
  show_labels = FALSE
) {
  half <- .width / 2

  per_participant <- participant_metrics |>
    mutate(participant = factor(participant)) |>
    group_by(participant, phase, metric) |>
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
      scales = if (free_y) "free_y" else NULL,
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
      axis.ticks.x = if (show_labels) element_line() else element_blank(),
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

# Per-participant beta per cue x block.
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

  per_participant <- participant_beta_long |>
    mutate(participant = factor(participant)) |>
    group_by(participant, phase, cue) |>
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
      axis.ticks.x = if (show_labels) element_line() else element_blank(),
      panel.spacing.x = unit(0.4, "lines"),
      legend.position = "right"
    ) +
    labs(
      x = "Participant",
      y = paste0("Per-participant β (median +/- ", round(.width * 100), "% CI)")
    )
}

# Headline figure: time-series on the left, phase x cue posterior on the right.
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

  posterior_plot + ts_plot + theme(legend.position = "right")
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
  fname <- file.path(dir, paste0(prefix, gsub(" ", "_", roi_label), ".png"))
  ggsave(
    filename = fname,
    plot = combined_plot,
    width = width,
    height = height,
    dpi = dpi
  )
  invisible(fname)
}


# ---- 8. LOO helper ---------------------------------------------------------

# Approximate leave-one-out CV. Requires log_lik to have been read in with
# the draws (pass variables = c("mu_betas", "log_lik") to the loader).
compute_loo <- function(draws) {
  ll <- draws |>
    posterior::subset_draws(variable = "log_lik") |>
    posterior::as_draws_array()
  r_eff <- loo::relative_eff(exp(ll), cores = 10)
  loo::loo(ll, r_eff = r_eff, cores = 10)
}


# ============================================================================
# ---- 9. Behavioral data (run once at source time) --------------------------
# ============================================================================
# Loads all Day 1 ratings, pivots to long form, and creates `ratings2` with
# within/between person splits. `ratings2` is used by build_roi_df() below.

Day1_ratings_paths <- list.files(
  path = parent_folder,
  pattern = "Day1.*ratings.dat$",
  recursive = TRUE,
  full.names = TRUE
)

for (i in seq_along(Day1_ratings_paths)) {
  current_ratings <- read.csv(Day1_ratings_paths[i])
  if (i == 1) {
    ratings_df <- current_ratings
  } else {
    ratings_df <- rbind.data.frame(ratings_df, current_ratings)
  }
}

ratings_df$partInd <- as.character(ratings_df$partInd)

# Experiment mislabels a rating column
names(ratings_df)[9] <- "ar_gs2"

ratings_df <- ratings_df |>
  select(!contains("Dur")) |>
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

ratings2 <- ratings_df |>
  filter(partInd %in% useable_participants) |>
  mutate(
    block = recode(
      ratInd,
      `1` = "hab",
      `2` = "acq_1",
      `3` = "acq_2",
      `4` = "ext"
    )
  ) |>
  group_by(partInd) |>
  mutate(
    val_btw = mean(val),
    val_wth = val - mean(val),
    ar_btw = mean(ar),
    ar_wth = ar - mean(ar),
    exp_btw = mean(exp),
    exp_wth = exp - mean(exp)
  ) |>
  ungroup() |>
  mutate(participant = as.integer(as.factor(partInd)))


# ============================================================================
# ---- 10. build_roi_df() — generalized data-prep for any ROI ---------------
# ============================================================================
#
# Loads per-participant betas for `roi_id` from `model_id`, joins them to
# `ratings2`, and returns a `df_fit` ready for fit_m_final_less().
#
# Arguments
#   roi_id    : string matching a roi_id in chain_registry (e.g. "ant_insula")
#   model_id  : string matching a model_id in chain_registry (default "model042")
#   j_max     : maximum j index to keep; 16 excludes shock (j=17), 17 includes it
#   registry  : chain_registry tibble (default uses the global)
#
# Returns a tibble with one row per (participant, block, condition/cue), ready
# to pass directly to fit_m_final_less().

build_roi_df <- function(
  roi_id,
  model_id = "model042",
  j_max = 16,
  registry = chain_registry
) {
  # -- label maps (same for every ROI) --------------------------------------
  cue_map <- c("CS+" = "csp", "GS1" = "gs1", "GS2" = "gs2", "GS3" = "gs3")
  phase_map <- c(
    "Habituation" = "hab",
    "Acquisition #1" = "acq_1",
    "Acquisition #2" = "acq_2",
    "Extinction" = "ext"
  )

  # -- 1. load per-subject betas --------------------------------------------
  draws <- read_chains(
    model_id = model_id,
    roi_id = roi_id,
    variables = c("mu_betas", "betas"),
    registry = registry
  )

  # -- 2. tidy: average over ROI dimension, attach phase/cue labels ---------
  pp_betas <- tidy_participant_betas(draws) |>
    filter(j <= j_max)

  # -- 3. summarise to brain_mean / brain_sd per participant x phase x cue --
  brain <- pp_betas |>
    group_by(participant, phase, cue) |>
    summarise(
      brain_mean = mean(beta),
      brain_sd = sd(beta),
      .groups = "drop"
    ) |>
    mutate(
      condition = cue_map[as.character(cue)],
      block = phase_map[as.character(phase)]
    )

  # -- 4. join ratings -------------------------------------------------------
  df <- brain |>
    inner_join(
      ratings2 |>
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
          exp_btw
        ),
      by = c("participant", "block", "condition")
    )

  n_dropped <- nrow(brain) - nrow(df)
  if (n_dropped > 0) {
    message(sprintf(
      "build_roi_df(%s): %d rows dropped on join (label mismatch or missing ratings)",
      roi_id,
      n_dropped
    ))
  }

  # -- 5. compute derived predictors ----------------------------------------
  gm <- mean(df$exp)

  df <- df |>
    mutate(
      exp_abs = (exp - gm) / 10, # absolute, grand-mean-centered
      exp_wth10 = exp_wth / 10, # within-person deviation (10 pp units)
      exp_btw10 = (exp_btw - gm) / 10 # person mean, centered
    )

  zc <- function(x) as.numeric(scale(x)) # z-score helper

  df <- df |>
    mutate(
      val_btw10 = (val_btw - mean(val)) / 10,
      ar_btw10 = (ar_btw - mean(ar)) / 10,
      aff_wth = rowMeans(cbind(zc(val_wth), zc(ar_wth))), # within-person aversiveness
      aff_btw = rowMeans(cbind(zc(val_btw), zc(ar_btw)))
    ) |>
    group_by(participant) |>
    mutate(
      aff_wth_pm = mean(aff_wth),
      expXaff_pm = mean(exp_wth10 * aff_wth) # person-mean of the interaction
    ) |>
    ungroup() |>
    mutate(
      threat_wth = rowMeans(cbind(zc(exp_wth), zc(val_wth), zc(ar_wth))),
      threat_btw = rowMeans(cbind(zc(exp_btw), zc(val_btw), zc(ar_btw)))
    )

  # -- 6. drop rows with any missing model predictors -----------------------
  df_fit <- df |>
    tidyr::drop_na(
      brain_mean,
      brain_sd,
      exp_wth10,
      aff_wth,
      exp_btw10,
      aff_btw,
      block,
      cue,
      participant
    )

  df_fit
}


# ============================================================================
# ---- 11. fit_m_final_less() — generalized brm fit -------------------------
# ============================================================================
#
# Fits the primary behavioral model to any `df_fit` returned by build_roi_df().
#
# Formula (m_final_less):
#   brain_mean | se(brain_sd, sigma = TRUE) ~
#     exp_wth10 * aff_wth * block +
#     block:cue +
#     exp_btw10 +
#     aff_btw +
#     (1 + exp_wth10 + aff_wth + exp_wth10:aff_wth | participant)
#
# Arguments
#   df_fit  : data frame from build_roi_df()
#   chains  : number of MCMC chains (default 4)
#   cores   : number of cores (default 4)
#   iter    : total iterations per chain including warmup (default 2000)
#   ...     : additional arguments passed to brm() (e.g. control, prior)
#
# Returns a brmsfit object.

fit_m_final_less <- function(
  df_fit,
  chains = 4,
  cores = 4,
  iter = 2000,
  ...
) {
  brm(
    brain_mean | se(brain_sd, sigma = TRUE) ~
      exp_wth10 *
      aff_wth *
      block +
      block:cue +
      exp_btw10 +
      aff_btw +
      (1 + exp_wth10 + aff_wth + exp_wth10:aff_wth | participant),
    data = df_fit,
    family = gaussian(),
    chains = chains,
    cores = cores,
    iter = iter,
    save_pars = save_pars(all = TRUE),
    backend = "cmdstanr",
    ...
  )
}


# ============================================================================
# ---- 12. Behavioral plot functions ----------------------------------------
# ============================================================================
# Generalized from a008_separate_plots_5.R. Every function accepts:
#   model     : a brmsfit from fit_m_final_less()
#   df_fit    : the data frame from build_roi_df() used to fit that model
#   roi_label : a string used in plot titles and y-axis labels
#
# Module-level constants shared across all plot functions.

.beh_block_labels <- c(
  hab = "Habituation",
  acq_1 = "Acquisition 1",
  acq_2 = "Acquisition 2",
  ext = "Extinction"
)
.beh_block_order <- c("hab", "acq_1", "acq_2", "ext")
.beh_cue_levels <- c("CS+", "GS1", "GS2", "GS3")
.beh_cue_colors <- c(
  "CS+" = "#E63946",
  "GS1" = "#2DC653",
  "GS2" = "#9B5DE5",
  "GS3" = "#457B9D"
)
.beh_affect_colors <- c(
  "Low affect (\u22121 SD)" = "#2166AC",
  "Avg affect (0)" = "#666666",
  "High affect (+1 SD)" = "#B2182B"
)
.beh_shared_theme <- theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 11, colour = "grey50")
  )

# ── Low-level helpers ─────────────────────────────────────────────────────────

# Posterior mean of epred for each row of newdata.
epred_vec <- function(model, newdata, re_form = NA) {
  posterior_epred(model, newdata = newdata, re_formula = re_form) |>
    colMeans()
}

# Add epred draws then summarise to median + 95% CrI per grid row.
epred_ci <- function(grid, model, n_draws = 1500, ...) {
  add_epred_draws(grid, model, re_formula = NA, ndraws = n_draws, ...) |>
    group_by(
      block,
      cue,
      across(any_of(c("exp_wth10", "aff_wth", "affect_level")))
    ) |>
    summarise(
      epred = median(.epred),
      lo = quantile(.epred, 0.025),
      hi = quantile(.epred, 0.975),
      .groups = "drop"
    )
}

# Partial residual for a given focal predictor.
#   focal = "exp"  -> expectancy (aff + between zeroed; cue kept)
#   focal = "aff"  -> affect     (exp + between zeroed; cue kept)
#   focal = "syn"  -> exp × aff surface (between zeroed; cue marginalised)
# re_form = NA  -> population partial residual (cloud includes between spread)
# re_form = NULL -> conditional partial residual (pure residual)
partial_resid <- function(
  model,
  data,
  focal = c("exp", "aff", "syn"),
  re_form = NA
) {
  focal <- match.arg(focal)
  mu_full <- epred_vec(model, data, re_form)
  ref <- data

  if (focal == "exp") {
    ref$aff_wth <- 0
    ref$exp_btw10 <- 0
    ref$aff_btw <- 0
    mu_focal <- epred_vec(model, ref, re_form)
  } else if (focal == "aff") {
    ref$exp_wth10 <- 0
    ref$exp_btw10 <- 0
    ref$aff_btw <- 0
    mu_focal <- epred_vec(model, ref, re_form)
  } else {
    ref$exp_btw10 <- 0
    ref$aff_btw <- 0
    mu_focal <- Reduce(
      `+`,
      lapply(.beh_cue_levels, function(cc) {
        r <- ref
        r$cue <- factor(cc, levels = .beh_cue_levels)
        epred_vec(model, r, re_form)
      })
    ) /
      length(.beh_cue_levels)
  }

  data |>
    mutate(
      .mu_full = mu_full,
      .mu_focal = mu_focal,
      .presid = brain_mean - mu_full + mu_focal
    )
}

# ── plot_partial_effects() ────────────────────────────────────────────────────
#
# Returns a named list with elements $exp, $aff, and $syn, each a ggplot of
# partial residuals for that predictor.

plot_partial_effects <- function(
  model,
  df_fit,
  roi_label,
  n_draws = 1500
) {
  aff_sd <- sd(df_fit$aff_wth)
  mean_bsd <- mean(df_fit$brain_sd)
  y_lab <- paste0(roi_label, " \u03b2  (partial residual)")

  df_obs <- df_fit |>
    mutate(
      cue = factor(cue, levels = .beh_cue_levels),
      block = factor(block, levels = .beh_block_order)
    )

  # --- grids ----------------------------------------------------------------
  exp_grid <- crossing(
    block = factor(.beh_block_order, levels = .beh_block_order),
    cue = factor(.beh_cue_levels, levels = .beh_cue_levels),
    exp_wth10 = seq(
      min(df_fit$exp_wth10),
      max(df_fit$exp_wth10),
      length.out = 35
    )
  ) |>
    mutate(aff_wth = 0, exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)

  aff_grid <- crossing(
    block = factor(.beh_block_order, levels = .beh_block_order),
    cue = factor(.beh_cue_levels, levels = .beh_cue_levels),
    aff_wth = seq(min(df_fit$aff_wth), max(df_fit$aff_wth), length.out = 35)
  ) |>
    mutate(exp_wth10 = 0, exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)

  syn_grid_all <- crossing(
    block = factor(.beh_block_order, levels = .beh_block_order),
    cue = factor(.beh_cue_levels, levels = .beh_cue_levels),
    exp_wth10 = seq(
      min(df_fit$exp_wth10),
      max(df_fit$exp_wth10),
      length.out = 35
    ),
    aff_wth = c(-aff_sd, 0, aff_sd)
  ) |>
    mutate(exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)

  # --- partial residuals ----------------------------------------------------
  pr_exp <- partial_resid(model, df_obs, "exp")
  pr_aff <- partial_resid(model, df_obs, "aff")
  pr_syn <- partial_resid(model, df_obs, "syn") |>
    mutate(aff_sd_units = aff_wth / aff_sd)

  # --- prediction lines -----------------------------------------------------
  exp_line <- epred_ci(exp_grid, model, n_draws = n_draws)
  aff_line <- epred_ci(aff_grid, model, n_draws = n_draws)

  syn_line <- add_epred_draws(
    syn_grid_all,
    model,
    re_formula = NA,
    ndraws = n_draws
  ) |>
    group_by(.draw, block, exp_wth10, aff_wth) |>
    summarise(.epred_avg = mean(.epred), .groups = "drop") |>
    group_by(block, exp_wth10, aff_wth) |>
    summarise(
      epred = median(.epred_avg),
      lo = quantile(.epred_avg, 0.025),
      hi = quantile(.epred_avg, 0.975),
      .groups = "drop"
    ) |>
    mutate(
      affect_level = factor(
        round(aff_wth, 6),
        levels = round(c(-aff_sd, 0, aff_sd), 6),
        labels = c(
          "Low affect (\u22121 SD)",
          "Avg affect (0)",
          "High affect (+1 SD)"
        )
      ),
      block = factor(block, levels = .beh_block_order)
    )

  # --- plots ----------------------------------------------------------------
  p_exp <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_ribbon(
      data = exp_line,
      aes(x = exp_wth10, ymin = lo, ymax = hi, fill = cue),
      alpha = 0.15
    ) +
    geom_line(
      data = exp_line,
      aes(x = exp_wth10, y = epred, colour = cue),
      linewidth = 0.9
    ) +
    geom_point(
      data = pr_exp,
      aes(x = exp_wth10, y = .presid, colour = cue),
      size = 1.3,
      alpha = 0.5
    ) +
    facet_wrap(~block, labeller = as_labeller(.beh_block_labels), nrow = 1) +
    scale_colour_manual(values = .beh_cue_colors, name = "Cue") +
    scale_fill_manual(values = .beh_cue_colors, name = "Cue") +
    labs(
      x = "Within-person expectancy (exp_wth10; 10 pp units)",
      y = y_lab,
      title = paste0("Expectancy \u2192 ", roi_label, "  |  m_final_less"),
      subtitle = "Lines: population partial effect (aff_wth = 0, between = 0). Dots: partial residuals"
    ) +
    .beh_shared_theme

  p_aff <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_ribbon(
      data = aff_line,
      aes(x = aff_wth, ymin = lo, ymax = hi, fill = cue),
      alpha = 0.15
    ) +
    geom_line(
      data = aff_line,
      aes(x = aff_wth, y = epred, colour = cue),
      linewidth = 0.9
    ) +
    geom_point(
      data = pr_aff,
      aes(x = aff_wth, y = .presid, colour = cue),
      size = 1.3,
      alpha = 0.5
    ) +
    facet_wrap(~block, labeller = as_labeller(.beh_block_labels), nrow = 1) +
    scale_colour_manual(values = .beh_cue_colors, name = "Cue") +
    scale_fill_manual(values = .beh_cue_colors, name = "Cue") +
    labs(
      x = "Within-person affect composite (aff_wth; z-scored avg. of val & ar)",
      y = y_lab,
      title = paste0("Affect \u2192 ", roi_label, "  |  m_final_less"),
      subtitle = "Lines: population partial effect (exp_wth10 = 0, between = 0). Dots: partial residuals"
    ) +
    .beh_shared_theme

  p_syn <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_point(
      data = pr_syn,
      aes(x = exp_wth10, y = .presid, fill = aff_sd_units),
      shape = 21,
      colour = NA,
      size = 1.7,
      alpha = 0.7
    ) +
    geom_ribbon(
      data = filter(syn_line, affect_level == "Low affect (\u22121 SD)"),
      aes(x = exp_wth10, ymin = lo, ymax = hi),
      fill = .beh_affect_colors[["Low affect (\u22121 SD)"]],
      alpha = 0.14
    ) +
    geom_ribbon(
      data = filter(syn_line, affect_level == "Avg affect (0)"),
      aes(x = exp_wth10, ymin = lo, ymax = hi),
      fill = .beh_affect_colors[["Avg affect (0)"]],
      alpha = 0.14
    ) +
    geom_ribbon(
      data = filter(syn_line, affect_level == "High affect (+1 SD)"),
      aes(x = exp_wth10, ymin = lo, ymax = hi),
      fill = .beh_affect_colors[["High affect (+1 SD)"]],
      alpha = 0.14
    ) +
    geom_line(
      data = syn_line,
      aes(x = exp_wth10, y = epred, colour = affect_level),
      linewidth = 1.1
    ) +
    facet_wrap(~block, labeller = as_labeller(.beh_block_labels), nrow = 1) +
    scale_colour_manual(values = .beh_affect_colors, name = NULL) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "grey85",
      high = "#B2182B",
      midpoint = 0,
      name = "Dot affect (SD)",
      guide = guide_colourbar(
        barheight = unit(0.4, "cm"),
        barwidth = unit(4, "cm")
      )
    ) +
    labs(
      x = "Within-person expectancy (exp_wth10; 10 pp units)",
      y = y_lab,
      title = paste0(
        "Expectancy \u00d7 Affect synergy \u2192 ",
        roi_label,
        "  |  m_final_less"
      ),
      subtitle = "Lines: cue-marginal partial effect at \u00b11 SD affect. Dots: partial residuals, shaded by affect"
    ) +
    .beh_shared_theme +
    theme(legend.key.width = unit(1.4, "cm"))

  list(exp = p_exp, aff = p_aff, syn = p_syn)
}

# ── plot_surface() ────────────────────────────────────────────────────────────
#
# 2-D expectancy × affect surface (median posterior prediction) with
# observed data overlaid. Dots sized by |residual|.
#
# Arguments
#   by_cue      : FALSE = single CS+ surface (1×4 block panels, default)
#                 TRUE  = 4-row cue × block grid — makes between-cue intercept
#                         differences directly visible as color differences
#   fill_limits : e.g. c(-0.01, 0.10) to fix the color range; NULL = auto

plot_surface <- function(
  model,
  df_fit,
  roi_label,
  n_draws = 400,
  by_cue = FALSE,
  fill_limits = NULL
) {
  mean_bsd <- mean(df_fit$brain_sd)
  y_lab <- paste0(roi_label, " \u03b2")

  df_obs <- df_fit |>
    mutate(
      cue = factor(cue, levels = .beh_cue_levels),
      block = factor(
        block,
        levels = .beh_block_order,
        labels = .beh_block_labels
      )
    )

  # ── y-axis in raw rating points instead of SD ────────────────────────────────
  # aff_wth is a z-composite, so it has no natural unit. Relabel its axis in the
  # average within-person change of valence & arousal (0–10 points) using the
  # linear map between the z-composite and the raw-averaged composite. The printed
  # R² says how exact the relabel is (≈1 when valence & arousal have similar
  # within-person SDs; if it drops, the two composites diverge and the raw-point
  # labels become approximate).
  aff_lin <- lm((val_wth + ar_wth) / 2 ~ aff_wth, data = df_fit)
  b0 <- coef(aff_lin)[[1]]
  b1 <- coef(aff_lin)[[2]]
  raw_brks <- seq(-4, 5, by = 1)
  z_pos <- (raw_brks - b0) / b1
  keepy <- z_pos >= min(df_fit$aff_wth) & z_pos <= max(df_fit$aff_wth)
  raw_labs <- ifelse(
    raw_brks > 0,
    paste0("+", raw_brks),
    as.character(raw_brks)
  )

  # Build prediction grid — include all cues when by_cue = TRUE
  cues_in_grid <- if (by_cue) .beh_cue_levels else "CS+"

  surf_grid <- crossing(
    block = factor(.beh_block_order, levels = .beh_block_order),
    cue = factor(cues_in_grid, levels = .beh_cue_levels),
    exp_wth10 = seq(
      min(df_fit$exp_wth10),
      max(df_fit$exp_wth10),
      length.out = 40
    ),
    aff_wth = seq(min(df_fit$aff_wth), max(df_fit$aff_wth), length.out = 40)
  ) |>
    mutate(exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)

  surf <- add_epred_draws(
    surf_grid,
    model,
    re_formula = NA,
    ndraws = n_draws
  ) |>
    group_by(block, cue, exp_wth10, aff_wth) |>
    summarise(epred = median(.epred), .groups = "drop") |>
    mutate(
      block = factor(
        block,
        levels = .beh_block_order,
        labels = .beh_block_labels
      )
    )

  res_obs <- df_obs |>
    mutate(
      .fitted = fitted(model, re_formula = NA)[, "Estimate"],
      .resid = brain_mean - .fitted
    )

  r_max <- max(abs(res_obs$.resid))
  halo_width <- 1.5
  vsize <- function(x) 1.8 + pmin(x^2 / r_max^2, 1) * 5.7

  obs_pts <- res_obs |>
    transmute(
      block,
      cue,
      exp_wth10,
      aff_wth,
      abs_resid = abs(.resid),
      size_color = vsize(abs(.resid)),
      size_black = vsize(abs(.resid)) + halo_width
    )

  # pretty() gives clean round-number contour breaks robust to small signal range
  contour_breaks <- pretty(surf$epred, n = 6)

  p <- ggplot(surf, aes(exp_wth10, aff_wth)) +
    geom_raster(aes(fill = epred), interpolate = TRUE) +
    geom_contour(
      aes(z = epred),
      breaks = contour_breaks,
      colour = "black",
      linewidth = 0.5,
      alpha = 0.6
    ) +
    geom_point(
      data = obs_pts,
      aes(exp_wth10, aff_wth, size = size_black),
      shape = 16,
      colour = "grey10",
      show.legend = FALSE
    ) +
    geom_point(
      data = obs_pts,
      aes(exp_wth10, aff_wth, colour = cue, size = size_color),
      shape = 16,
      alpha = 0.95
    ) +
    scale_fill_viridis_c(
      option = "turbo",
      name = y_lab,
      breaks = contour_breaks,
      labels = sprintf("%.2f", contour_breaks),
      guide = guide_colourbar(
        order = 1,
        barheight = unit(6, "cm"),
        barwidth = unit(0.5, "cm")
      )
    ) +
    scale_colour_manual(
      values = .beh_cue_colors,
      name = "Cue",
      guide = guide_legend(override.aes = list(shape = 16, size = 4))
    ) +
    scale_size_identity(guide = "none") +
    scale_x_continuous(
      breaks = c(-4, 0, 4, 8),
      labels = c("\u221240", "0", "+40", "+80"),
      name = "Within-person \u0394 expectancy (pp)",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = z_pos[keepy],
      labels = raw_labs[keepy],
      name = "Within-person \u0394 affect (rating pts)",
      expand = c(0, 0)
    ) +
    .beh_shared_theme +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.box = "vertical"
    )

  if (by_cue) {
    p +
      facet_grid(
        cue ~ block,
        labeller = labeller(block = as_labeller(.beh_block_labels))
      ) +
      labs(
        title = paste0(
          roi_label,
          ": expectancy \u00d7 affect surface by cue \u00d7 block"
        ),
        subtitle = "Rows = cues; warm-vs-cool differences between rows show cue discrimination"
      )
  } else {
    p +
      facet_wrap(~block, nrow = 1) +
      labs(
        title = paste0(
          "Expectancy \u00d7 affect surface (CS+ level) \u2014 ",
          roi_label
        ),
        subtitle = "Model prediction for CS+; dots = all cues, coloured; size \u221d |residual|\u00b2"
      )
  }
}

# ── plot_residuals() ──────────────────────────────────────────────────────────
#
# Returns list(obs, cell): individual residuals ordered by cue × block, and
# mean residuals per cue × block cell.

plot_residuals <- function(model, df_fit, roi_label) {
  df_obs <- df_fit |>
    mutate(
      cue = factor(cue, levels = .beh_cue_levels),
      block = factor(
        block,
        levels = .beh_block_order,
        labels = .beh_block_labels
      )
    ) |>
    mutate(
      .fitted = fitted(model, re_formula = NA)[, "Estimate"],
      .resid = brain_mean - .fitted
    )

  order_key <- df_obs |>
    arrange(cue, block, participant) |>
    mutate(obs_idx = row_number()) |>
    select(participant, block, cue, obs_idx)

  res_indexed <- df_obs |>
    left_join(order_key, by = c("participant", "block", "cue"))

  cue_seps <- order_key |>
    group_by(cue) |>
    summarise(max_idx = max(obs_idx), .groups = "drop") |>
    filter(as.integer(cue) < length(.beh_cue_levels))

  block_seps <- order_key |>
    group_by(cue, block) |>
    summarise(max_idx = max(obs_idx), .groups = "drop") |>
    filter(!max_idx %in% cue_seps$max_idx)

  cue_mids <- order_key |>
    group_by(cue) |>
    summarise(mid = mean(obs_idx), .groups = "drop")

  p_obs <- ggplot(res_indexed, aes(obs_idx, .resid, color = cue)) +
    geom_vline(
      data = cue_seps,
      aes(xintercept = max_idx + 0.5),
      color = "black",
      linewidth = 0.7
    ) +
    geom_vline(
      data = block_seps,
      aes(xintercept = max_idx + 0.5),
      color = "grey55",
      linewidth = 0.3
    ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    geom_point(size = 1.1, alpha = 0.5) +
    scale_color_manual(values = .beh_cue_colors, name = "Cue") +
    scale_x_continuous(breaks = cue_mids$mid, labels = cue_mids$cue) +
    labs(
      x = NULL,
      y = "Residual (obs \u2212 fitted)",
      title = paste0(
        "Individual residuals (cue \u00d7 block order) \u2014 ",
        roi_label
      ),
      subtitle = "Thick lines = cue boundaries; thin lines = block boundaries within cue"
    ) +
    .beh_shared_theme

  res_cell <- df_obs |>
    group_by(block, cue) |>
    summarise(mr = mean(.resid), se = sd(.resid) / sqrt(n()), .groups = "drop")

  p_cell <- ggplot(res_cell, aes(cue, mr, color = cue)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    geom_pointrange(aes(ymin = mr - 2 * se, ymax = mr + 2 * se)) +
    facet_wrap(~block, nrow = 1) +
    scale_color_manual(values = .beh_cue_colors, guide = "none") +
    labs(
      x = NULL,
      y = "Mean residual (obs \u2212 fitted)",
      title = paste0("Mean residual by cue \u00d7 block \u2014 ", roi_label)
    ) +
    .beh_shared_theme

  list(obs = p_obs, cell = p_cell)
}

# ── plot_loo_r2() ─────────────────────────────────────────────────────────────
#
# Posterior of LOO-adjusted R² (full and population-level predictions).

plot_loo_r2 <- function(model, roi_label, n_boot = 4000) {
  .ll <- log_lik(model)
  .nch <- model$fit@sim$chains
  .r_eff <- relative_eff(
    exp(.ll),
    chain_id = rep(seq_len(.nch), each = nrow(.ll) / .nch)
  )
  .psis <- psis(log_ratios = -.ll, r_eff = .r_eff)
  .y <- model$data$brain_mean

  .err_full <- E_loo(posterior_epred(model), .psis, log_ratios = -.ll)$value -
    .y
  .err_pop <- E_loo(
    posterior_epred(model, re_formula = NA),
    .psis,
    log_ratios = -.ll
  )$value -
    .y

  loo_R2_post <- function(err, y = .y, seed = 7314) {
    set.seed(seed)
    N <- length(y)
    w <- matrix(rexp(n_boot * N), nrow = n_boot)
    w <- w / rowSums(w)
    vy <- (N / (N - 1)) *
      (rowSums(sweep(w, 2, y^2, `*`)) - rowSums(sweep(w, 2, y, `*`))^2)
    ve <- (N / (N - 1)) *
      (rowSums(sweep(w, 2, err^2, `*`)) - rowSums(sweep(w, 2, err, `*`))^2)
    pmin(pmax(1 - ve / vy, -1), 1)
  }

  r2_long <- bind_rows(
    tibble(
      R2 = loo_R2_post(.err_full),
      type = "Full (incl. participant effects)"
    ),
    tibble(R2 = loo_R2_post(.err_pop), type = "Population only (fixed effects)")
  )

  ggplot(r2_long, aes(R2, fill = type)) +
    geom_density(alpha = 0.5, colour = NA) +
    scale_fill_manual(
      values = c(
        "Full (incl. participant effects)" = "#457B9D",
        "Population only (fixed effects)" = "#E63946"
      ),
      name = NULL
    ) +
    labs(
      x = expression("LOO-adjusted " * R^2),
      y = "Posterior density",
      title = paste0("Posterior of LOO-adjusted R\u00b2 \u2014 ", roi_label),
      subtitle = "Out-of-sample variance explained (Bayesian bootstrap)"
    ) +
    .beh_shared_theme
}

# ── make_m_final_less_plots() ─────────────────────────────────────────────────
#
# Master convenience function: runs all plot functions and returns a named list.
#
# Arguments
#   model     : brmsfit from fit_m_final_less()
#   df_fit    : data frame from build_roi_df()
#   roi_label : string for axis / title labels (e.g. "Amygdala")
#   n_draws   : posterior draws to use for partial effect lines (default 1500)
#   surf_draws: posterior draws for the surface plot (default 400; slow)
#
# Returns list(exp, aff, syn, surf, resid_obs, resid_cell, loo_r2)

make_m_final_less_plots <- function(
  model,
  df_fit,
  roi_label,
  n_draws = 1500,
  surf_draws = 400
) {
  message("Computing partial effect lines...")
  partial <- plot_partial_effects(model, df_fit, roi_label, n_draws = n_draws)

  message("Computing surface...")
  surf <- plot_surface(model, df_fit, roi_label, n_draws = surf_draws)

  message("Computing residual plots...")
  resids <- plot_residuals(model, df_fit, roi_label)

  message("Computing LOO R\u00b2...")
  r2 <- plot_loo_r2(model, roi_label)

  list(
    exp = partial$exp,
    aff = partial$aff,
    syn = partial$syn,
    surf = surf,
    resid_obs = resids$obs,
    resid_cell = resids$cell,
    loo_r2 = r2
  )
}

# ── save_m_final_less_plots() ─────────────────────────────────────────────────
#
# Save all plots from make_m_final_less_plots() to PLOT_DIR.

save_m_final_less_plots <- function(
  plots,
  roi_label,
  dir = PLOT_DIR,
  prefix = NULL,
  width = 15,
  height = 8,
  dpi = 150
) {
  if (is.null(prefix)) {
    prefix <- paste0(gsub(" ", "_", tolower(roi_label)), "_")
  }

  sv <- function(name, plot, w = width, h = height) {
    ggsave(file.path(dir, name), plot, width = w, height = h, dpi = dpi)
  }

  sv(paste0(prefix, "exp_m_final_less.png"), plots$exp)
  sv(paste0(prefix, "aff_m_final_less.png"), plots$aff)
  sv(paste0(prefix, "syn_m_final_less.png"), plots$syn)
  sv(paste0(prefix, "surf.png"), plots$surf)
  sv(paste0(prefix, "resid_obs_wc.png"), plots$resid_obs)
  sv(paste0(prefix, "resid_cell_wc.png"), plots$resid_cell)
  ggsave(
    file.path(dir, paste0(prefix, "loo_r2.png")),
    plots$loo_r2,
    width = 8,
    height = 6,
    dpi = dpi
  )

  invisible(dir)
}


# ============================================================================
# ---- 13. Driver — example calls for one or many ROIs ----------------------
# ============================================================================
# Uncomment and run interactively; do not leave these executing at source time.

if (FALSE) {
  # ── Single ROI (anterior insula) ────────────────────────────────────────

  df_ai <- build_roi_df("ant_insula", model_id = "model042")
  nrow(df_ai) # sanity check

  m_ai <- fit_m_final_less(df_ai)
  summary(m_ai)

  # LOO for the fit
  m_ai <- add_criterion(m_ai, "loo", moment_match = TRUE)
  m_ai$criteria$loo

  # Quick hypothesis test: is the exp × aff synergy non-zero in acquisition?
  hypothesis(
    m_ai,
    c(
      acq1 = "exp_wth10:aff_wth + exp_wth10:aff_wth:blockacq_1 = 0",
      acq2 = "exp_wth10:aff_wth + exp_wth10:aff_wth:blockacq_2 = 0"
    )
  )

  # ── All model042 ROIs ────────────────────────────────────────────────────

  roi_ids <- chain_registry |>
    filter(model_id == "model042") |>
    pull(roi_id)

  # Build all data frames first (fast: just loads posteriors + joins)
  roi_dfs <- setNames(
    lapply(roi_ids, build_roi_df),
    roi_ids
  )

  # Fit each ROI (slow: one full brm run per ROI)
  roi_fits <- setNames(
    lapply(roi_ids, function(r) fit_m_final_less(roi_dfs[[r]])),
    roi_ids
  )

  # ── fMRI summary plots for any ROI (from sections 5-7) ──────────────────

  roi_label <- "Anterior Insula"
  mu_draws <- read_chains("model042", "ant_insula")
  trial_draws <- read_trial_betas_chains("model038", "ant_insula")
  combined <- plot_roi_summary(trial_draws, mu_draws, roi_label)
  print(combined)
}


# Anterior Insula
roi_label <- "Anterior Insula"

df_ai <- build_roi_df("ant_insula", model_id = "model042")
cat("Rows:", nrow(df_ai), "\n")


m_ai <- fit_m_final_less(df_ai)

summary(m_ai)

# Generate partial effects plots (exp + aff + syn)
message("Computing partial effects for Anterior Insula...")
ai_partial <- plot_partial_effects(
  m_ai,
  df_ai,
  roi_label = roi_label,
  n_draws = 1000
)
print(ai_partial$exp)

print(ai_partial$aff)


print(ai_partial$syn)


# Residual plots
ai_resids <- plot_residuals(m_ai, df_ai, roi_label = roi_label)
print(ai_resids$cell)


# Surface plot
print(plot_surface(m_ai, df_ai, roi_label = roi_label, n_draws = 400))

print(plot_surface(m_ai, df_ai, roi_label, n_draws = 400, by_cue = TRUE))


# Amygdala
df_amy <- build_roi_df("AMY", model_id = "model042")
cat("Rows:", nrow(df_amy), "\n")


m_amy <- fit_m_final_less(df_amy)

summary(m_amy)


# Check the random effect SDs and brain_sd scale
cat("--- Amygdala brain_sd summary ---\n")
summary(df_amy$brain_sd)
cat("\n--- Anterior insula brain_sd (for comparison) ---\n")
summary(df_ai$brain_sd)
cat("\n--- Amygdala brain_mean summary ---\n")
summary(df_amy$brain_mean)


# Generate partial effects plots (exp + aff + syn)
message("Computing partial effects for Amygdala...")
amy_partial <- plot_partial_effects(
  m_amy,
  df_amy,
  roi_label = "Amygdala",
  n_draws = 1000
)
print(amy_partial$exp)

print(amy_partial$aff)


print(amy_partial$syn)


# Residual plots
amy_resids <- plot_residuals(m_amy, df_amy, roi_label = "Amygdala")
print(amy_resids$cell)

# Random effect collapse: all participant-level slope SDs are estimated at ~0. With amygdala brain_mean values roughly 5× smaller than anterior insula
# (median ≈ −0.005 vs a larger range for insula), the measurement error term absorbs essentially all the residual variance and there's nothing left
# for individual differences in slopes. The model is effectively a fixed-effects model for this ROI.
#
# Substantive effects: there's no evidence of an expectancy–affect synergy in the amygdala (all exp_wth10:aff_wth CIs span zero cleanly). The
# clearest signal is in the block:cue terms — the GS cues show reliably lower amygdala activation than CS+ across acquisition and habituation blocks,
# which is the cue discrimination pattern.

# Surface plot
print(plot_surface(m_amy, df_amy, roi_label = "Amygdala", n_draws = 400))

print(plot_surface(m_amy, df_amy, "Amygdala", n_draws = 400, by_cue = TRUE))
