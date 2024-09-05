library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)

source("code/custom_functions.R") # function to extract posteriors and get mean and sd of predictors

# load model
brm_taxon_mac_fine_vel = readRDS(file = "models/brm_taxon_mac_fine_vel.rds")

# rerun with weaker priors
brm_taxon_mac_fine_vel_weak = update(brm_taxon_mac_fine_vel,
                                     prior = c(prior(normal(3, 4), class = "Intercept"),
                                               prior(exponential(2), class = "sd"),
                                               prior(normal(0, 2), class = "b")))

saveRDS(brm_taxon_mac_fine_vel_weak, file = "models/brm_taxon_mac_fine_vel_weak.rds")



# extract posteriors ------------------------------------------------------
brm_taxon_mac_fine_vel_weak = readRDS(file = "models/brm_taxon_mac_fine_vel_weak.rds")

mod_d = brm_taxon_mac_fine_vel_weak$data

# get posteriors
velocity_preds_weak = get_snail_posts(model = brm_taxon_mac_fine_vel_weak, 
                                 predictor1 = "velocity_s", 
                                 predictor2 = "macrophyte_s",
                                 predictor3 = "fines_s") %>% 
  mutate(predictor = "velocity_s", 
         pred_values = velocity_s,
         model = "weak priors")

macrophyte_preds_weak = get_snail_posts(model = brm_taxon_mac_fine_vel_weak, 
                                   predictor2 = "velocity_s", 
                                   predictor1 = "macrophyte_s",
                                   predictor3 = "fines_s") %>% 
  mutate(predictor = "macrophyte_s", 
         pred_values = macrophyte_s,
         model = "weak priors")

fines_preds_weak = get_snail_posts(model = brm_taxon_mac_fine_vel_weak, 
                              predictor3 = "velocity_s", 
                              predictor2 = "macrophyte_s",
                              predictor1 = "fines_s") %>% 
  mutate(predictor = "fines_s",
         pred_values = fines_s,
         model = "weak priors")

post_preds_taxa = readRDS(file = "posteriors/post_preds_taxa.rds") %>% mutate(model = "main model priors")

all_weak_preds = bind_rows(velocity_preds_weak, 
                           macrophyte_preds_weak, 
                           fines_preds_weak) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (pred_values*sd) + mean,
         raw_predictor = case_when(predictor == "fines_s" ~ "a) % Fine Sediment",
                                   predictor == "macrophyte_s" ~ "b) % Macrophytes",
                                   TRUE ~ "c) Stream Velocity"))
  

all_preds = bind_rows(all_weak_preds, post_preds_taxa)


# plot --------------------------------------------------------------------

plot_prior_sensitivity = all_preds %>% 
  ggplot(aes(x = pred_values, y = .epred + 1)) + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = model),
                  alpha = 0.3, linewidth = 0.1) + 
  ggh4x::facet_grid2(model~raw_predictor, scales = "free_x") +
  # geom_point(data = pred_data, aes(y = total_snail + 1),
             # size = 0.1, shape = ".") +
  scale_y_log10() +
  labs(x = "Predictor Value",
       y = "Total Snails per Quadrat + 1",
       fill = "",
       color = "") +
  guides(fill = "none",
         color = "none") +
  scale_fill_brewer(type = "qual") +
  theme_default() +
  theme(strip.text = element_text(hjust = 0)) +
  NULL

ggsave(plot_prior_sensitivity, file = "plots/plot_prior_sensitivity.jpg", width = 6.5, height = 3.5)
