library(tidyverse)

roi_ts_dir <- "/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/roi_timeseries"

# Pick one participant to inspect — change as needed
participant <- "149"
highpass <- signal::butter(n = 5, W = 0.026, type = "high")

# Load main ROI stats file
bold_file <- file.path(roi_ts_dir, paste0(participant, "_roi_stats_pb03.txt"))
bold_raw <- read_delim(bold_file, trim_ws = TRUE, show_col_types = FALSE)

# Pull NZMean columns and apply AFNI-style percent signal change centering
bold_scaled <- bold_raw |>
  select(starts_with("NZMean")) |>
  mutate(across(everything(), ~ 100 * . / mean(., na.rm = TRUE) - 100)) |>
  mutate(across(everything(), ~ signal::filtfilt(highpass, .))) |>
  mutate(time = row_number())

# Pivot to long for ggplot
bold_long <- bold_scaled |>
  pivot_longer(-time, names_to = "roi", values_to = "pct_signal_change")

# Plot a random subset of ROIs to keep the plot readable
set.seed(5)
set.seed(43)
rois_to_plot <- bold_long |>
  distinct(roi) |>
  slice_sample(n = min(6, n_distinct(bold_long$roi))) |>
  pull(roi)

bold_long |>
  filter(roi %in% rois_to_plot) |>
  ggplot(aes(x = time, y = pct_signal_change, color = roi)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  coord_cartesian(xlim = c(0, 100)) +
  facet_wrap(~roi, ncol = 2, scales = "free_y") +
  labs(
    title = paste("pb03 ROI Timeseries — Participant", participant),
    subtitle = "AFNI-style percent signal change (mean = 0)",
    x = "TR",
    y = "% Signal Change"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
