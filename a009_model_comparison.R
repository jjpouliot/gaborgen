# ============================================================
# a009_model_comparison.R
# Compare m_final_less / m_final_less_2 / m_final_more across
# every figure type from a008_separate_plots.R.
#
# Run AFTER sourcing a008 so these are already in scope:
#   df_fit, df_obs, block_labels, block_order, cue_levels,
#   cue_colors, affect_colors, mean_bsd, n_draws, shared_theme,
#   epred_ci(), exp_grid, aff_grid, syn_grid_all, surf_grid,
#   contour_breaks
# and the fitted models (with loo added via add_criterion).
# ============================================================

library(tidyverse)
library(tidybayes)
library(patchwork)

# ---- the models to compare, in display (facet-row) order ----
models <- list(
  "less (slopes by block)" = m_final_less,
  "less_2 (no cue)" = m_final_less_2,
  "more (slopes by block x cue)" = m_final_more
)
# add m_final_less_3 here if you want the additive no-cue row too
model_levels <- names(models)

# pretty block factor, used everywhere so facet_grid needs no labeller
to_block <- function(x) factor(x, levels = block_order, labels = block_labels)
df_obs_p <- df_obs |> mutate(block = to_block(block))

# ============================================================
# 1. Tagged prediction tibbles (one row-set per model)
# ============================================================

exp_all <- imap(models, function(m, lab) {
  epred_ci(exp_grid, m) |> mutate(model = lab)
}) |>
  list_rbind() |>
  mutate(block = to_block(block), model = factor(model, levels = model_levels))

aff_all <- imap(models, function(m, lab) {
  epred_ci(aff_grid, m) |> mutate(model = lab)
}) |>
  list_rbind() |>
  mutate(block = to_block(block), model = factor(model, levels = model_levels))

# synergy: marginalise over cue within each draw, then summarise
syn_one <- function(model, label) {
  add_epred_draws(syn_grid_all, model, re_formula = NA, ndraws = n_draws) |>
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
        aff_wth,
        levels = c(-1, 0, 1),
        labels = c(
          "Low affect (-1 SD)",
          "Avg affect (0)",
          "High affect (+1 SD)"
        )
      ),
      model = label
    )
}
syn_all <- imap(models, ~ syn_one(.x, .y)) |>
  list_rbind() |>
  mutate(block = to_block(block), model = factor(model, levels = model_levels))

# ============================================================
# 2. Comparison panels: facet_grid(model ~ block)
# ============================================================

p_exp_cmp <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(
    data = df_obs_p,
    aes(exp_wth10, brain_mean, colour = cue),
    size = 0.9,
    alpha = 0.30
  ) +
  geom_ribbon(
    data = exp_all,
    aes(exp_wth10, ymin = lo, ymax = hi, fill = cue),
    alpha = 0.12
  ) +
  geom_line(
    data = exp_all,
    aes(exp_wth10, epred, colour = cue),
    linewidth = 0.9
  ) +
  facet_grid(model ~ block) +
  scale_colour_manual(values = cue_colors, name = "Cue") +
  scale_fill_manual(values = cue_colors, name = "Cue") +
  labs(
    x = "Within-person expectancy (10 pp units)",
    y = "Anterior insula \u03b2",
    title = "Expectancy \u2192 insula, by model",
    subtitle = "parallel cue lines = shared slope; fanning = cue-specific slope (only m_final_more can fan)"
  ) +
  shared_theme

p_aff_cmp <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(
    data = df_obs_p,
    aes(aff_wth, brain_mean, colour = cue),
    size = 0.9,
    alpha = 0.30
  ) +
  geom_ribbon(
    data = aff_all,
    aes(aff_wth, ymin = lo, ymax = hi, fill = cue),
    alpha = 0.12
  ) +
  geom_line(
    data = aff_all,
    aes(aff_wth, epred, colour = cue),
    linewidth = 0.9
  ) +
  facet_grid(model ~ block) +
  scale_colour_manual(values = cue_colors, name = "Cue") +
  scale_fill_manual(values = cue_colors, name = "Cue") +
  labs(
    x = "Within-person affect composite (z-scored avg. of val & ar)",
    y = "Anterior insula \u03b2",
    title = "Affect \u2192 insula, by model",
    subtitle = "same logic: parallel = shared affect slope, fanning = cue-specific affect slope"
  ) +
  shared_theme

p_syn_cmp <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(
    data = df_obs_p,
    aes(exp_wth10, brain_mean, colour = cue),
    size = 0.9,
    alpha = 0.25
  ) +
  geom_ribbon(
    data = syn_all,
    aes(exp_wth10, ymin = lo, ymax = hi, fill = affect_level),
    alpha = 0.16
  ) +
  geom_line(
    data = syn_all,
    aes(exp_wth10, epred, colour = affect_level),
    linewidth = 1.0
  ) +
  facet_grid(model ~ block) +
  scale_colour_manual(
    values = c(cue_colors, affect_colors),
    name = NULL,
    breaks = names(affect_colors)
  ) +
  scale_fill_manual(values = affect_colors, guide = "none") +
  labs(
    x = "Within-person expectancy (10 pp units)",
    y = "Anterior insula \u03b2",
    title = "Expectancy \u00d7 affect synergy, by model (marginalised over cue)",
    subtitle = "spread between affect levels = synergy strength"
  ) +
  shared_theme

p_exp_cmp
p_aff_cmp
p_syn_cmp

# ============================================================
# 3. The focused payoff: shared vs cue-specific slopes
#    (less = parallel, more = fanned) on one set of panels
# ============================================================

exp_two <- exp_all |>
  filter(model %in% c("less (slopes by block)", "more (slopes by block x cue)"))

p_slope_overlay <- ggplot(
  exp_two,
  aes(exp_wth10, epred, colour = cue, linetype = model)
) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~block, nrow = 1) +
  scale_colour_manual(values = cue_colors, name = "Cue") +
  scale_linetype_manual(
    values = c(
      "less (slopes by block)" = "solid",
      "more (slopes by block x cue)" = "21"
    ),
    name = "Model"
  ) +
  labs(
    x = "Within-person expectancy (10 pp units)",
    y = "Anterior insula \u03b2",
    title = "Shared vs cue-specific slopes",
    subtitle = "solid = parallel (less); dashed = cue-specific (more) \u2014 does the dashed fanning exceed noise?"
  ) +
  shared_theme

p_slope_overlay

# ============================================================
# 4. Residuals by cue x block, now including m_final_more
#    (more has full cell intercepts too, so its cell-MEAN
#     residual will sit ~0 like less; the discriminating
#     diagnostic for more is RMSE within cell, plotted below)
# ============================================================

get_resid <- function(model, name) {
  df_fit |>
    mutate(
      .fitted = fitted(model, re_formula = NA)[, "Estimate"],
      .resid = brain_mean - .fitted,
      model = name
    )
}

res_long <- imap(models, ~ get_resid(.x, .y)) |>
  list_rbind() |>
  mutate(
    block = to_block(block),
    cue = factor(cue, levels = cue_levels),
    model = factor(model, levels = model_levels)
  )

# within-cell RMSE: does letting slopes vary by cue tighten the fit?
res_rmse <- res_long |>
  group_by(model, block, cue) |>
  summarise(rmse = sqrt(mean(.resid^2)), .groups = "drop")

p_rmse <- ggplot(res_rmse, aes(cue, rmse, fill = cue)) +
  geom_col() +
  facet_grid(model ~ block) +
  scale_fill_manual(values = cue_colors, guide = "none") +
  labs(
    x = NULL,
    y = "within-cell RMSE (population-level)",
    title = "Fit tightness by cue x block \u2014 lower = better within-cell fit"
  ) +
  shared_theme

p_rmse

# ============================================================
# 5. ELPD by cell: where does m_final_more win/lose vs less?
#    (requires loo added to both via add_criterion)
# ============================================================

elpd_more_vs_less <- tibble(
  block = to_block(df_fit$block),
  cue = factor(df_fit$cue, levels = cue_levels),
  diff = m_final_more$criteria$loo$pointwise[, "elpd_loo"] -
    m_final_less$criteria$loo$pointwise[, "elpd_loo"]
) |>
  group_by(block, cue) |>
  summarise(
    elpd_gain = sum(diff),
    elpd_se = sqrt(n() * var(diff)),
    .groups = "drop"
  )

p_elpd_more <- ggplot(elpd_more_vs_less, aes(cue, elpd_gain, fill = cue)) +
  geom_col() +
  geom_linerange(
    aes(ymin = elpd_gain - elpd_se, ymax = elpd_gain + elpd_se),
    linewidth = 1.2
  ) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  facet_wrap(~block, nrow = 1) +
  scale_fill_manual(values = cue_colors, guide = "none") +
  labs(
    x = NULL,
    y = "ELPD gained by cue-specific slopes (more - less)",
    title = "Where letting slopes vary by cue helps (positive) or just overfits (~0 / negative)"
  ) +
  shared_theme

p_elpd_more

# overall verdict
print(loo_compare(m_final_less, m_final_less_2, m_final_more))

# ============================================================
# 6. Surfaces by model (CS+ cell). For m_final_more the surface
#    genuinely differs across cues, so swap cue in surf_grid or
#    facet by cue to see the cue-specific surfaces.
# ============================================================

surf_all <- imap(models, function(m, lab) {
  add_epred_draws(surf_grid, m, re_formula = NA, ndraws = 300) |>
    group_by(block, exp_wth10, aff_wth) |>
    summarise(epred = median(.epred), .groups = "drop") |>
    mutate(block = to_block(block), model = lab)
}) |>
  list_rbind() |>
  mutate(model = factor(model, levels = model_levels))

p_surf_cmp <- ggplot(surf_all, aes(exp_wth10, aff_wth, fill = epred)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(
    aes(z = epred),
    breaks = contour_breaks,
    colour = "black",
    linewidth = 0.4,
    alpha = 0.7
  ) +
  facet_grid(model ~ block) +
  scale_fill_viridis_c(option = "turbo", name = "Insula \u03b2") +
  labs(
    x = "Within-person expectancy",
    y = "Within-person affect",
    title = "Expectancy \u00d7 affect surface by model (CS+ cell)"
  ) +
  shared_theme

p_surf_cmp

# ── Save all comparison plots ─────────────────────────────────────────────────
ggsave(
  file.path(PLOT_DIR, "ai_exp_model_cmp.png"),
  p_exp_cmp,
  width = 22,
  height = 14,
  dpi = 150
)
ggsave(
  file.path(PLOT_DIR, "ai_aff_model_cmp.png"),
  p_aff_cmp,
  width = 22,
  height = 14,
  dpi = 150
)
ggsave(
  file.path(PLOT_DIR, "ai_syn_model_cmp.png"),
  p_syn_cmp,
  width = 22,
  height = 14,
  dpi = 150
)
ggsave(
  file.path(PLOT_DIR, "ai_slope_overlay.png"),
  p_slope_overlay,
  width = 22,
  height = 7,
  dpi = 150
)
ggsave(
  file.path(PLOT_DIR, "ai_rmse_model_cmp.png"),
  p_rmse,
  width = 22,
  height = 12,
  dpi = 150
)
ggsave(
  file.path(PLOT_DIR, "ai_elpd_more_vs_less.png"),
  p_elpd_more,
  width = 22,
  height = 7,
  dpi = 150
)
ggsave(
  file.path(PLOT_DIR, "ai_surf_model_cmp.png"),
  p_surf_cmp,
  width = 22,
  height = 14,
  dpi = 150
)
