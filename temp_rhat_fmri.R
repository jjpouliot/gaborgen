library(tidyverse)
library(posterior)

d <- as_draws_rvars(example_draws("multi_normal"))
rhat(d$Sigma)

warmup_draws <- read.csv("/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/misc/temp_fmri_mod3_warmup_chains.csv")

warmup_draws %>% 
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = delta))

warmup_draws %>% 
  select(!c("Betas.13",
            "Betas.14",
            "Betas.15",
            "Betas.16",
            "Betas.17",
            "Betas.18"
            )) %>% 
  pivot_longer(starts_with("Beta")) %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_density(aes(x = value, color = name))

hold <- warmup_draws %>%
  group_by(chain) %>%
  mutate(.iteration = row_number(),
         .chain = chain) 

hold2 <- as_draws_df(hold)

hold3 <- hold2 %>%
  filter(.iteration < 201)

hold4 <- hold3 %>%
  ungroup() %>% 
  select(.chain,.iteration,Betas.5) %>% 
  pivot_wider(names_from = .chain, values_from = Betas.5) 



m <- as.matrix(hold4[, 2:5])  # convert selected columns to plain matrix
rhat(m)

rhat(data.frame(hold4[2:5]))
