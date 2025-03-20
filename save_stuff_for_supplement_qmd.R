model002_fit_meta_data <- model002_fit$metadata()

model002_fit_relevant_parameters <- model002_fit_meta_data$model_params[
  str_detect(model002_fit_meta_data$model_params, "learning_rate\\[")]

model002_fit_learning_rate_posteriors <- model002_fit$draws(variables = model002_fit_relevant_parameters,format = "df")


CSP_parameters <- model002_fit_meta_data$model_params[
  str_detect(model002_fit_meta_data$model_params, "CSP_ass|scaling\\[")] 

CSP_df <- posterior::as_draws_df(
  model002_fit$draws(variables = CSP_parameters))

model002_fit_relevant_parameters <- model002_fit_meta_data$model_params[
  str_detect(model002_fit_meta_data$model_params, "scaling\\[")]

model002_fit_scaling_posteriors <- model002_fit$draws(variables = model002_fit_relevant_parameters,format = "df")


model004_fit_meta_data <- model004_fit$metadata()

model004_fit_relevant_parameters <- model004_fit_meta_data$model_params[
  str_detect(model004_fit_meta_data$model_params, "learning_unpaired\\[|learning_paired\\[|scaling\\[")]

model004_fit_learning_rate_posteriors <- model004_fit$draws(variables = model004_fit_relevant_parameters,format = "df")



cue_color <- c("red1","green1", "purple1", "blue1")




save(cue_color, model002_fit_loo, model003_fit_loo, Oz_fft_df, 
     model002_fit_learning_rate_posteriors, model002_fit_scaling_posteriors,
     model004_fit_learning_rate_posteriors,
     gaborgen_stan_list, 
     CSP_df,
     file = "/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/supplemental_data.RData")


Oz_fft_df %>% 
  ggplot() +
  geom_density(aes(x = zamp, color = participant)) +
  theme_bw()

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot() +
  geom_density(aes(x = amplitude_15Hz_fft, 
                   color = stan_par_id,
                   ),linewidth = 1,
               alpha = .4) +
  theme_bw() +
  theme(text = element_text(size = 20))

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0)) +
  geom_density(aes(x = zamp, 
                   color = stan_par_id,
                   ),linewidth = 1,
               alpha = .4) +
  theme_bw() +
  theme(text = element_text(size = 20))

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot() +
  geom_density_ridges(aes(x = amplitude_15Hz_fft, 
                   color = cue,
                   y = block
                   ),linewidth = 1,
               alpha = .4) +
  scale_color_manual(values = cue_color) +
  facet_wrap(~stan_par_id) +
  theme_bw()

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot(aes(x = amplitude_15Hz_fft, y = block, color = cue)) +
  stat_density_ridges(scale = 1,
                      rel_min_height = 0.01,
                      panel_scaling = FALSE,
                      from = 0, 
                      alpha = 0) +
  scale_color_manual(values = cue_color) +
  coord_cartesian(xlim = c(0,.45)) +
  facet_wrap(~stan_par_id) +
  theme_bw() +
  theme(text = element_text(size = 20))

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot(aes(x = amplitude_15Hz_fft, y = block, color = cue)) +
  geom_jitter(height = .1,size = .65, alpha = 1) +
  scale_color_manual(values = cue_color) +
  coord_cartesian(xlim = c(0,.45)) +
  facet_wrap(~stan_par_id) +
  theme_bw() +
  theme(text = element_text(size = 20))

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot(aes(x = zamp, y = block, color = cue)) +
  stat_density_ridges(scale = 1,
                      rel_min_height = 0.01,
                      panel_scaling = FALSE,
                      alpha = 0) +
  scale_color_manual(values = cue_color) +
  facet_wrap(~stan_par_id) +
  theme_bw()

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id),
         block = factor(block, levels = rev(unique(block)))) %>% 
  ggplot(aes(x = zamp, y = block, color = cue)) +
  geom_jitter(height = .1,size = .65, alpha = 1) +
  scale_color_manual(values = cue_color) +
  facet_wrap(~stan_par_id) +
  theme_bw()

Oz_fft_df %>% 
  mutate(stan_par_id = factor(stan_par_id)) %>% 
  group_by(stan_par_id) %>% 
  reframe(average_amp = mean(amplitude_15Hz_fft)) %>% 
  arrange(desc(average_amp)) %>% 
  print(n = 999)
