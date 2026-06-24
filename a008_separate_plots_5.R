# =============================================================================
# a008_separate_plots.R  —  m_final_less visualisations (anterior insula)
#
# What changed vs. the previous version
# -------------------------------------
# 1. Expectancy / Affect / Synergy panels now plot PARTIAL RESIDUALS, not raw
#    brain_mean, so the dots live on the same scale as the partial line:
#
#       partial_resid_i = y_i - mu_full_i + mu_focal_i
#
#    mu_full  = population fit at the row's ACTUAL predictors (re_formula = NA)
#    mu_focal = population fit with every NON-focal predictor set to its
#               reference, focal predictor(s) + block (+ cue) kept at actual.
#    Leftover (y - mu_full) carries the residual AND the between-person spread
#    the population line ignores -> an honest cloud. Flip re_form to NULL in
#    partial_resid() for the tighter conditional (per-participant) version.
# 2. Synergy affect levels use +/- sd(aff_wth), so "+/-1 SD" is actually true
#    (aff_wth is a mean of two z-scores, so its SD < 1).
# 3. Synergy dots are coloured by aff_wth (the moderator), not cue, so the
#    moderation is visible; cue intercept is residualised out.
# 4. Surface: clipped to the per-block convex hull of observed (exp, aff) to
#    kill corner extrapolation, observed cells overlaid, y-axis in true SD,
#    x-axis relabelled as percentage points (not "Delta%").
# 5. Only m_final_less is kept (the _2 / _3 comparison plots were dropped).
# =============================================================================

library(patchwork)
library(tidyverse)
library(tidybayes)
library(brms)

# ── Labels / colours / theme ─────────────────────────────────────────────────
block_labels <- c(
  hab = "Habituation",
  acq_1 = "Acquisition 1",
  acq_2 = "Acquisition 2",
  ext = "Extinction"
)
block_order <- c("hab", "acq_1", "acq_2", "ext")
cue_levels <- c("CS+", "GS1", "GS2", "GS3")
cue_colors <- c(
  "CS+" = "#E63946",
  "GS1" = "#2DC653",
  "GS2" = "#9B5DE5",
  "GS3" = "#457B9D"
)

mean_bsd <- mean(df_fit$brain_sd)
aff_sd <- sd(df_fit$aff_wth) # honest +/-1 SD for the synergy levels
n_draws <- 1500

shared_theme <- theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 11, colour = "grey50")
  )

# Observed data (re-levelled for consistent facet order). Keep CODE labels for
# block/cue so posterior_epred() matches the levels the model was fit on.
df_obs <- df_fit |>
  mutate(
    cue = factor(cue, levels = cue_levels),
    block = factor(block, levels = block_order)
  )

# ── helper: add_epred → median + 95% CrI (for the partial LINES) ─────────────
epred_ci <- function(grid, model, ...) {
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

# =============================================================================
# PARTIAL-RESIDUAL MACHINERY
# =============================================================================
# posterior mean of epred for each row of `newdata`
epred_vec <- function(model, newdata, re_form = NA) {
  posterior_epred(model, newdata = newdata, re_formula = re_form) |>
    colMeans()
}

# partial_resid: returns df_obs with .presid for the chosen focal component.
#   focal = "exp"  -> expectancy line (aff, between zeroed; cue kept)
#   focal = "aff"  -> affect line     (exp, between zeroed; cue kept)
#   focal = "syn"  -> within x within surface, marginalised over cue
# re_form = NA   -> population partial residual (cloud = resid + between spread)
# re_form = NULL -> conditional partial residual (cloud = pure residual)
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
    # keep exp, aff, block; zero between; average focal over cues so the dots
    # match the cue-marginal synergy lines (cue intercept residualised out)
    ref$exp_btw10 <- 0
    ref$aff_btw <- 0
    mu_focal <- Reduce(
      `+`,
      lapply(cue_levels, function(cc) {
        r <- ref
        r$cue <- factor(cc, levels = cue_levels)
        epred_vec(model, r, re_form)
      })
    ) /
      length(cue_levels)
  }

  data |>
    mutate(
      .mu_full = mu_full,
      .mu_focal = mu_focal,
      .presid = brain_mean - mu_full + mu_focal
    )
}

pr_exp <- partial_resid(m_final_less, df_obs, "exp")
pr_aff <- partial_resid(m_final_less, df_obs, "aff")
pr_syn <- partial_resid(m_final_less, df_obs, "syn") |>
  mutate(aff_sd_units = aff_wth / aff_sd) # colour dots in SD units

# =============================================================================
# PARTIAL LINES (grids)
# =============================================================================
# ── 1. Expectancy grid (aff_wth = 0) ─────────────────────────────────────────
exp_grid <- crossing(
  block = factor(block_order, levels = block_order),
  cue = factor(cue_levels, levels = cue_levels),
  exp_wth10 = seq(min(df_fit$exp_wth10), max(df_fit$exp_wth10), length.out = 35)
) |>
  mutate(aff_wth = 0, exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)
exp_line <- epred_ci(exp_grid, m_final_less)

# ── 2. Affect grid (exp_wth10 = 0) ───────────────────────────────────────────
aff_grid <- crossing(
  block = factor(block_order, levels = block_order),
  cue = factor(cue_levels, levels = cue_levels),
  aff_wth = seq(min(df_fit$aff_wth), max(df_fit$aff_wth), length.out = 35)
) |>
  mutate(exp_wth10 = 0, exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)
aff_line <- epred_ci(aff_grid, m_final_less)

# ── 3. Synergy grid: marginal over cue, 3 affect levels at TRUE +/-1 SD ──────
syn_grid_all <- crossing(
  block = factor(block_order, levels = block_order),
  cue = factor(cue_levels, levels = cue_levels),
  exp_wth10 = seq(
    min(df_fit$exp_wth10),
    max(df_fit$exp_wth10),
    length.out = 35
  ),
  aff_wth = c(-aff_sd, 0, aff_sd)
) |>
  mutate(exp_btw10 = 0, aff_btw = 0, brain_sd = mean_bsd)

syn_line <- add_epred_draws(
  syn_grid_all,
  m_final_less,
  re_formula = NA,
  ndraws = n_draws
) |>
  group_by(.draw, block, exp_wth10, aff_wth) |>
  summarise(.epred_avg = mean(.epred), .groups = "drop") |> # marginal over cue
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
      labels = c("Low affect (−1 SD)", "Avg affect (0)", "High affect (+1 SD)")
    ),
    block = factor(block, levels = block_order)
  )

# =============================================================================
# PLOT 1: EXPECTANCY  (partial residuals)
# =============================================================================
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
  facet_wrap(~block, labeller = as_labeller(block_labels), nrow = 1) +
  scale_colour_manual(values = cue_colors, name = "Cue") +
  scale_fill_manual(values = cue_colors, name = "Cue") +
  labs(
    x = "Within-person expectancy (exp_wth10; 10 pp units)",
    y = "Anterior insula β  (partial residual)",
    title = "Expectancy → anterior insula  |  m_final_less",
    subtitle = "Lines: population partial effect (aff_wth = 0, between = 0). Dots: partial residuals (affect + between + cue removed)"
  ) +
  shared_theme

# =============================================================================
# PLOT 2: AFFECT  (partial residuals)
# =============================================================================
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
  facet_wrap(~block, labeller = as_labeller(block_labels), nrow = 1) +
  scale_colour_manual(values = cue_colors, name = "Cue") +
  scale_fill_manual(values = cue_colors, name = "Cue") +
  labs(
    x = "Within-person affect composite (aff_wth; z-scored avg. of val & ar)",
    y = "Anterior insula β  (partial residual)",
    title = "Affect → anterior insula  |  m_final_less",
    subtitle = "Lines: population partial effect (exp_wth10 = 0, between = 0). Dots: partial residuals (expectancy + between removed)"
  ) +
  shared_theme

# =============================================================================
# PLOT 3: SYNERGY  (partial residuals; dots coloured by affect, not cue)
# =============================================================================
affect_colors <- c(
  "Low affect (−1 SD)" = "#2166AC",
  "Avg affect (0)" = "#666666",
  "High affect (+1 SD)" = "#B2182B"
)

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
  # constant-fill ribbons (one layer per level) so the dot gradient can own
  # the single `fill` scale ggplot allows
  geom_ribbon(
    data = filter(syn_line, affect_level == "Low affect (−1 SD)"),
    aes(x = exp_wth10, ymin = lo, ymax = hi),
    fill = affect_colors[["Low affect (−1 SD)"]],
    alpha = 0.14
  ) +
  geom_ribbon(
    data = filter(syn_line, affect_level == "Avg affect (0)"),
    aes(x = exp_wth10, ymin = lo, ymax = hi),
    fill = affect_colors[["Avg affect (0)"]],
    alpha = 0.14
  ) +
  geom_ribbon(
    data = filter(syn_line, affect_level == "High affect (+1 SD)"),
    aes(x = exp_wth10, ymin = lo, ymax = hi),
    fill = affect_colors[["High affect (+1 SD)"]],
    alpha = 0.14
  ) +
  geom_line(
    data = syn_line,
    aes(x = exp_wth10, y = epred, colour = affect_level),
    linewidth = 1.1
  ) +
  facet_wrap(~block, labeller = as_labeller(block_labels), nrow = 1) +
  scale_colour_manual(values = affect_colors, name = NULL) +
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
    y = "Anterior insula β  (partial residual)",
    title = "Expectancy × Affect synergy → anterior insula  |  m_final_less",
    subtitle = "Lines: cue-marginal partial effect at ±1 SD affect. Dots: partial residuals (between + cue removed), shaded by affect"
  ) +
  shared_theme +
  theme(legend.key.width = unit(1.4, "cm"))

# =============================================================================
# PLOT 4: COMBINED SURFACE  (full surface; support + fit shown via dots)
# =============================================================================
# Per-observation population residual, reused by the surface dots and PLOT 5.
res_obs <- df_fit |>
  mutate(
    .fitted = fitted(m_final_less, re_formula = NA)[, "Estimate"],
    .resid = brain_mean - .fitted,
    block = factor(block, levels = block_order, labels = block_labels),
    cue = factor(cue, levels = cue_levels)
  )

surf_grid <- crossing(
  block = factor(block_order, levels = block_order),
  exp_wth10 = seq(
    min(df_fit$exp_wth10),
    max(df_fit$exp_wth10),
    length.out = 60
  ),
  aff_wth = seq(min(df_fit$aff_wth), max(df_fit$aff_wth), length.out = 60)
) |>
  mutate(
    cue = factor("CS+", levels = cue_levels),
    exp_btw10 = 0,
    aff_btw = 0,
    brain_sd = mean_bsd
  )

surf <- add_epred_draws(
  surf_grid,
  m_final_less,
  re_formula = NA,
  ndraws = 400
) |>
  group_by(block, exp_wth10, aff_wth) |>
  summarise(epred = median(.epred), .groups = "drop") |>
  mutate(block = factor(block, levels = block_order, labels = block_labels))

# dots: position = (exp, aff); colour = cue; size = closeness to the fitted value
obs_pts <- res_obs |>
  transmute(
    block,
    cue,
    exp_wth10,
    aff_wth,
    abs_resid = abs(.resid)
  )

# ── y-axis in raw rating points instead of SD ────────────────────────────────
# aff_wth is a z-composite, so it has no natural unit. Relabel its axis in the
# average within-person change of valence & arousal (0–10 points) using the
# linear map between the z-composite and the raw-averaged composite. The printed
# R² says how exact the relabel is (≈1 when valence & arousal have similar
# within-person SDs; if it drops, the two composites diverge and the raw-point
# labels become approximate).
aff_lin <- lm((val_wth + ar_wth) / 2 ~ aff_wth, data = df_fit)
message(sprintf(
  "aff_wth → rating-point relabel R² = %.4f (slope = %.3f pts per z-unit)",
  summary(aff_lin)$r.squared,
  coef(aff_lin)[[2]]
))
b0 <- coef(aff_lin)[[1]]
b1 <- coef(aff_lin)[[2]]
# raw_brks <- seq(-4, 4, by = 1)
raw_brks <- seq(-4, 5, by = 1)
z_pos <- (raw_brks - b0) / b1 # aff_wth position for each raw-point tick
keepy <- z_pos >= min(df_fit$aff_wth) & z_pos <= max(df_fit$aff_wth)
raw_labs <- ifelse(raw_brks > 0, paste0("+", raw_brks), as.character(raw_brks))

# ── dot size = |residual|: bigger dot = worse fit (poorly predicted stand out) ─
r_max <- max(obs_pts$abs_resid)
size_breaks <- c(0.05, 0.1, 0.2, 0.3, 0.4)
size_breaks <- size_breaks[size_breaks <= r_max]

# Map |residual| -> visual size with a power-2 transform so the dense low end is
# compressed and the rare large residuals visually jump out. Precompute the
# mapped size (rather than letting scale_size_continuous do it) so we can draw a
# slightly larger dark dot underneath each colored dot — that underlay reads as a
# thick outline that makes points pop against the turbo raster. The floor (1.8)
# keeps the smallest dots big enough to show their fill.
vsize <- function(x) 1.8 + pmin(x^2 / r_max^2, 1) * 5.7 # range 1.8..7.5
halo_width <- 1.5 # constant additive -> uniform ring thickness across dot sizes
obs_pts <- obs_pts |>
  mutate(
    size_color = vsize(abs_resid),
    # dark halo: colored dot is drawn ON TOP, so its interior is never eaten
    size_black = size_color + halo_width
  )
size_legend_breaks <- vsize(size_breaks) # identity positions for legend keys

contour_breaks <- seq(-0.05, 0.35, by = 0.05)

p_surf <- ggplot(surf, aes(exp_wth10, aff_wth)) +
  geom_raster(aes(fill = epred), interpolate = TRUE) +
  geom_contour(
    aes(z = epred),
    breaks = contour_breaks,
    colour = "black",
    linewidth = 0.5,
    alpha = 0.6
  ) +
  # dark underlay (slightly larger) -> reads as a thick proportional outline
  geom_point(
    data = obs_pts,
    aes(exp_wth10, aff_wth, size = size_black),
    shape = 16,
    colour = "grey10",
    show.legend = FALSE
  ) +
  # solid cue-coloured dot on top (uses the free `colour` scale, not `fill`)
  geom_point(
    data = obs_pts,
    aes(exp_wth10, aff_wth, colour = cue, size = size_color),
    shape = 16,
    alpha = 0.95
  ) +
  facet_wrap(~block, nrow = 1) +
  scale_fill_viridis_c(
    option = "turbo",
    name = "Insula β",
    breaks = contour_breaks,
    labels = sprintf("%.2f", contour_breaks),
    guide = guide_colourbar(
      order = 1,
      barheight = unit(6, "cm"),
      barwidth = unit(0.5, "cm")
    )
  ) +
  scale_colour_manual(
    values = cue_colors,
    name = "Cue",
    guide = guide_legend(order = 2, override.aes = list(shape = 16, size = 4))
  ) +
  scale_size_identity(
    name = "|obs − fitted|",
    breaks = size_legend_breaks,
    labels = size_breaks,
    guide = guide_legend(
      order = 3,
      override.aes = list(shape = 16, colour = "grey40")
    )
  ) +
  scale_x_continuous(
    breaks = c(-4, 0, 4, 8),
    labels = c("−40", "0", "+40", "+80"),
    name = "Within-person Δ expectancy (percentage points)",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = z_pos[keepy],
    labels = raw_labs[keepy],
    name = "Within-person Δ affect\n(rating points, avg. of unpleasantness & arousal)",
    expand = c(0, 0)
  ) +
  labs(
    title = "Combined expectancy × affect surface (CS+ level)",
    subtitle = "Shape identical across cues (no cue×exp / cue×aff term); only the level shifts. Dots = observations coloured by cue; size ∝ |residual|² (bigger = worse fit)."
  ) +
  shared_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(0.3, "cm"),
    legend.key.height = unit(0.65, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text = element_text(size = 9)
  )

# =============================================================================
# PLOT 5: RESIDUALS in cue × block order  (m_final_less only — already honest)
# =============================================================================
# (res_obs is computed above, before the surface plot)
order_key <- res_obs |>
  arrange(cue, block, participant) |>
  mutate(obs_idx = row_number()) |>
  select(participant, block, cue, obs_idx)

res_obs_indexed <- res_obs |>
  left_join(order_key, by = c("participant", "block", "cue"))

cue_seps <- order_key |>
  group_by(cue) |>
  summarise(max_idx = max(obs_idx), .groups = "drop") |>
  filter(as.integer(cue) < length(cue_levels))

block_seps <- order_key |>
  group_by(cue, block) |>
  summarise(max_idx = max(obs_idx), .groups = "drop") |>
  filter(!max_idx %in% cue_seps$max_idx)

cue_mids <- order_key |>
  group_by(cue) |>
  summarise(mid = mean(obs_idx), .groups = "drop")

p_resid_obs_wc <- ggplot(res_obs_indexed, aes(obs_idx, .resid, color = cue)) +
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
  scale_color_manual(values = cue_colors, name = "Cue") +
  scale_x_continuous(breaks = cue_mids$mid, labels = cue_mids$cue) +
  labs(
    x = NULL,
    y = "Residual (obs − fitted)",
    title = "Individual residuals in cue × block order  —  m_final_less",
    subtitle = "Thick lines = cue boundaries; thin lines = block boundaries within cue"
  ) +
  shared_theme

# cell-mean residuals (companion view)
res_cell_wc <- res_obs |>
  group_by(block, cue) |>
  summarise(mr = mean(.resid), se = sd(.resid) / sqrt(n()), .groups = "drop")

p_resid_cell_wc <- ggplot(res_cell_wc, aes(cue, mr, color = cue)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_pointrange(aes(ymin = mr - 2 * se, ymax = mr + 2 * se)) +
  facet_wrap(~block, nrow = 1) +
  scale_color_manual(values = cue_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Mean residual (obs − fitted)",
    title = "Mean residual by cue × block  —  m_final_less"
  ) +
  shared_theme

# =============================================================================
# LOO-ADJUSTED R²  (out-of-sample variance explained)
# =============================================================================
# Bayesian R² is in-sample and optimistic; LOO R² evaluates predictions on
# held-out observations (PSIS-LOO). We report two versions, matching the partial
# plots' conditioning:
#   full        = predictions include participant random effects (re_formula = NULL)
#   population  = population-level (fixed-effect) predictions (re_formula = NA),
#                 but the LOO weights still come from the actual (full) model.
# The posterior comes from a Bayesian bootstrap over observations, replicating
# brms:::.loo_R2. (Subgroup R² over small cells can legitimately be negative —
# that just means LOO predicts worse than the cell mean there — but these global
# values are well-positive.)
library(loo)

.ll <- log_lik(m_final_less) # S x N, for the PSIS-LOO weights
.nchains <- m_final_less$fit@sim$chains
.r_eff <- relative_eff(
  exp(.ll),
  chain_id = rep(seq_len(.nchains), each = nrow(.ll) / .nchains)
)
.psis <- psis(log_ratios = -.ll, r_eff = .r_eff)
.y <- m_final_less$data$brain_mean

# LOO point expectations -> LOO residuals, for full and population-level preds
.err_loo_full <- E_loo(
  posterior_epred(m_final_less),
  .psis,
  log_ratios = -.ll
)$value -
  .y
.err_loo_pop <- E_loo(
  posterior_epred(m_final_less, re_formula = NA),
  .psis,
  log_ratios = -.ll
)$value -
  .y

# Bayesian-bootstrap posterior of LOO R² from a vector of LOO residuals
loo_R2_post <- function(err, y = .y, n_boot = 4000, seed = 7314) {
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

r2_loo_full <- loo_R2_post(.err_loo_full)
r2_loo_pop <- loo_R2_post(.err_loo_pop)

r2_loo_long <- bind_rows(
  tibble(R2 = r2_loo_full, type = "Full (incl. participant effects)"),
  tibble(R2 = r2_loo_pop, type = "Population only (fixed effects)")
)

p_loo_r2 <- ggplot(r2_loo_long, aes(R2, fill = type)) +
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
    title = "Posterior of LOO-adjusted R² — m_final_less",
    subtitle = "Out-of-sample variance explained (Bayesian bootstrap over observations)"
  ) +
  shared_theme

# =============================================================================
# SAVE
# =============================================================================
width_val <- 15
height_val <- 8
sv <- function(name, plot) {
  ggsave(
    file.path(PLOT_DIR, name),
    plot,
    width = width_val,
    height = height_val,
    dpi = 150
  )
}

sv("ai_exp_m_final_less.png", p_exp)
sv("ai_aff_m_final_less.png", p_aff)
sv("ai_syn_m_final_less.png", p_syn)
sv("ai_surf.png", p_surf)
sv("ai_resid_obs_wc.png", p_resid_obs_wc)
sv("ai_resid_cell_wc.png", p_resid_cell_wc)
ggsave(
  file.path(PLOT_DIR, "ai_loo_r2.png"),
  p_loo_r2,
  width = 8,
  height = 6,
  dpi = 150
)
