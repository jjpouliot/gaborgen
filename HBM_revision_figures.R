Oz_fft_df %>%
  group_by(participant, block, cue) %>%
  reframe(
    mean_ar = mean(ar),
    mean_val = mean(val),
    mean_exp = mean(exp),
  ) %>%
  ggplot(aes(x = cue, y = mean_ar)) +
  geom_quasirandom() +
  geom_line(
    aes(
      x = cue,
      y = mean_ar,
      group = participant,
      position = ggbeeswarm::position_quasirandom()
    ),
  ) +
  facet_wrap(~block, ncol = 1)

Oz_fft_df %>%
  group_by(participant, block, cue) %>%
  reframe(
    mean_ar = mean(ar),
    mean_val = mean(val),
    mean_exp = mean(exp),
  ) %>%
    ggplot() +
  geom_line(
    aes(
      x = cue,
      y = mean_ar,
      group = participant
    ), # still needed to connect the right rows
    position = ggbeeswarm::position_quasirandom(
      width = 0.1,
      groupOnX = T
    )
  ) +
  geom_quasirandom(
    aes(
      x = cue,
      y = mean_ar,
      color = participant
    ),
    width = 0.1,
    groupOnX = T
  ) +
  facet_wrap(~block, ncol = 1)


library(ggbeeswarm)

Oz_fft_df %>%
  group_by(participant, block, cue) %>%
  reframe(
    mean_ar = mean(ar),
    mean_val = mean(val),
    mean_exp = mean(exp)
  ) %>%
  ggplot() +
  # points: groupOnX = TRUE lives *outside* aes()
geom_line(aes(x = cue, 
  y = mean_ar,
  group=participant,
  position=position_quasirandom()
)
  ) 
  geom_quasirandom(
    aes(color = participant),
    width    = 0.1,
    groupOnX = TRUE
  ) +
  # lines: same position_jitter spec, with groupOnX=TRUE
  ) +
  facet_wrap(~block, ncol = 1) +
  theme_minimal()



kl_normal <- function(mu_prior, sd_prior, mu_post, sd_post) {
  log(sd_prior / sd_post) +
    (sd_post^2 + (mu_post - mu_prior)^2) / (2 * sd_prior^2) -
    0.5
}


Oz_fft_df %>%
  filter(participant %in% c(113, 118, 129)) %>%
  ggplot() +
  geom_point(aes(x = trial, y = zamp, colour = cue), size = .2) +
  facet_wrap(~participant, ncol = 1) +
  theme_bw() +
  theme(
    # remove the little strip “header” above each panel
    # strip.background = element_blank(),
    # strip.text        = element_blank(),
    # collapse spacing between panels
    panel.spacing = unit(0.1, "lines")
  )

Oz_fft_df$trial %>% max()

model003_fit_meta_data <- model003_fit$metadata()

CSP_parameters <- model003_fit_meta_data$model_params[
  str_detect(model003_fit_meta_data$model_params, "CSP_ass|scaling\\[")
]

CSP_df <- posterior::as_draws_df(
  model003_fit$draws(variables = CSP_parameters)
)

trial_posterior <- data.frame(
  trial = numeric(),
  CSP = numeric(),
  GS1 = numeric(),
  GS2 = numeric(),
  GS3 = numeric()
)

trial_posterior_par_5 <- data.frame(
  trial = numeric(),
  CSP = numeric(),
  GS1 = numeric(),
  GS2 = numeric(),
  GS3 = numeric()
)

trial_posterior_par_14 <- data.frame(
  trial = numeric(),
  CSP = numeric(),
  GS1 = numeric(),
  GS2 = numeric(),
  GS3 = numeric()
)

for (t in 1:176) {
  trial_posterior[t, 1] <- t

  trial_posterior[t, 2] <- mean(
    model003_df$`intercept[1]` +
      model003_df$`fatigue[1]` * t +
      unlist(CSP_df$`scaling[1]` * CSP_df[, t + 4])
  )

  trial_posterior[t, 3] <- mean(
    model003_df$`intercept[1]` +
      model003_df$`fatigue[1]` * t +
      unlist(CSP_df$`scaling[2]` * CSP_df[, t + 4])
  )

  trial_posterior[t, 4] <- mean(
    model003_df$`intercept[1]` +
      model003_df$`fatigue[1]` * t +
      unlist(CSP_df$`scaling[3]` * CSP_df[, t + 4])
  )

  trial_posterior[t, 5] <- mean(
    model003_df$`intercept[1]` +
      model003_df$`fatigue[1]` * t +
      unlist(CSP_df$`scaling[4]` * CSP_df[, t + 4])
  )

  trial_posterior_par_5[t, 1] <- t

  trial_posterior_par_5[t, 2] <- mean(
    model003_df$`intercept[5]` +
      model003_df$`fatigue[5]` * t +
      unlist(CSP_df$`scaling[1]` * CSP_df[, t + 4 + 176 * 4])
  )

  trial_posterior_par_5[t, 3] <- mean(
    model003_df$`intercept[5]` +
      model003_df$`fatigue[5]` * t +
      unlist(CSP_df$`scaling[2]` * CSP_df[, t + 4 + 176 * 4])
  )

  trial_posterior_par_5[t, 4] <- mean(
    model003_df$`intercept[5]` +
      model003_df$`fatigue[5]` * t +
      unlist(CSP_df$`scaling[3]` * CSP_df[, t + 4 + 176 * 4])
  )

  trial_posterior_par_5[t, 5] <- mean(
    model003_df$`intercept[5]` +
      model003_df$`fatigue[5]` * t +
      unlist(CSP_df$`scaling[4]` * CSP_df[, t + 4 + 176 * 4])
  )

  trial_posterior_par_14[t, 1] <- t

  trial_posterior_par_14[t, 2] <- mean(
    model003_df$`intercept[14]` +
      model003_df$`fatigue[14]` * t +
      unlist(CSP_df$`scaling[1]` * CSP_df[, t + 4 + 176 * 13])
  )

  trial_posterior_par_14[t, 3] <- mean(
    model003_df$`intercept[14]` +
      model003_df$`fatigue[14]` * t +
      unlist(CSP_df$`scaling[2]` * CSP_df[, t + 4 + 176 * 13])
  )

  trial_posterior_par_14[t, 4] <- mean(
    model003_df$`intercept[14]` +
      model003_df$`fatigue[14]` * t +
      unlist(CSP_df$`scaling[3]` * CSP_df[, t + 4 + 176 * 13])
  )

  trial_posterior_par_14[t, 5] <- mean(
    model003_df$`intercept[14]` +
      model003_df$`fatigue[14]` * t +
      unlist(CSP_df$`scaling[4]` * CSP_df[, t + 4 + 176 * 13])
  )
}

cue_color <- c("red1", "green1", "purple1", "blue1")

# don't freak out if it looks like the CS+ effect happens before a shock, those trials are probably missing because of movement
(trial_posterior %>%
  ggplot() +
  geom_point(
    data = Oz_fft_df[Oz_fft_df$participant == "113", ],
    aes(x = trial, y = zamp, color = cue, shape = paired)
  ) +
  geom_line(aes(x = trial, y = CSP), color = cue_color[1]) +
  geom_line(aes(x = trial, y = GS1), color = cue_color[2]) +
  geom_line(aes(x = trial, y = GS2), color = cue_color[3]) +
  geom_line(aes(x = trial, y = GS3), color = cue_color[4]) +
  scale_color_manual(values = cue_color) +
  scale_shape_manual(values = c(16, 15)) +
  ggtitle("Participant 1") +
  theme_classic()) /
  (trial_posterior_par_5 %>%
    ggplot() +
    geom_point(
      data = Oz_fft_df[Oz_fft_df$participant == "118", ],
      aes(x = trial, y = zamp, color = cue, shape = paired)
    ) +
    geom_line(aes(x = trial, y = CSP), color = cue_color[1]) +
    geom_line(aes(x = trial, y = GS1), color = cue_color[2]) +
    geom_line(aes(x = trial, y = GS2), color = cue_color[3]) +
    geom_line(aes(x = trial, y = GS3), color = cue_color[4]) +
    scale_color_manual(values = cue_color) +
    scale_shape_manual(values = c(16, 15)) +
    ggtitle("Participant 5") +
    theme_classic()) /
  (trial_posterior_par_14 %>%
    ggplot() +
    #shocks
    # geom_vline(xintercept = 34) +
    # geom_vline(xintercept = 37) +
    geom_point(
      data = Oz_fft_df[Oz_fft_df$participant == "129", ],
      aes(x = trial, y = zamp, color = cue, shape = paired)
    ) +
    geom_line(aes(x = trial, y = CSP), color = cue_color[1]) +
    geom_line(aes(x = trial, y = GS1), color = cue_color[2]) +
    geom_line(aes(x = trial, y = GS2), color = cue_color[3]) +
    geom_line(aes(x = trial, y = GS3), color = cue_color[4]) +
    scale_color_manual(values = cue_color) +
    scale_shape_manual(values = c(16, 15)) +
    ggtitle("Participant 14") +
    theme_classic()) &
  theme(text = element_text(size = 20))


which(gaborgen_stan_list$paired[gaborgen_stan_list$participant == 14] == 1)
