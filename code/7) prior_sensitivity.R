library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)

source("code/custom_functions.R") # function to extract posteriors and get mean and sd of predictors


# prior vs posterior ------------------------------------------------------
post_preds_regression = readRDS(file = "posteriors/post_preds_regression.rds") %>% 
  mutate(model = "Posterior")

brm_taxon_mac_fine_vel = readRDS(file = "models/brm_taxon_mac_fine_vel.rds")

snail_density = readRDS("data/snail_density.rds") 

velocity_mean_sd = tibble(mean = attributes(snail_density$velocity_s)$`scaled:center`,
                          sd = attributes(snail_density$velocity_s)$`scaled:scale`,
                          predictor = "velocity_s")

macrophyte_mean_sd = tibble(mean = attributes(snail_density$macrophyte_s)$`scaled:center`,
                            sd = attributes(snail_density$macrophyte_s)$`scaled:scale`,
                            predictor = "macrophyte_s")

fines_mean_sd = tibble(mean = attributes(snail_density$fines_s)$`scaled:center`,
                       sd = attributes(snail_density$fines_s)$`scaled:scale`,
                       predictor = "fines_s")

predictors_mean_sd = bind_rows(velocity_mean_sd, macrophyte_mean_sd, fines_mean_sd)


# rerun with priors
# brm_taxon_mac_fine_vel_prior = update(brm_taxon_mac_fine_vel, sample_prior = "only")
# saveRDS(brm_taxon_mac_fine_vel_prior, file = "models/brm_taxon_mac_fine_vel_prior.rds")
brm_taxon_mac_fine_vel_prior = readRDS(file = "models/brm_taxon_mac_fine_vel_prior.rds")


# get samples of prior and wrangle
# get model data
mod_d = brm_taxon_mac_fine_vel_prior$data

# get raw data
snail_density = readRDS(file = "data/snail_density.rds")

# make predictors
taxon_predictor_grid_prior = get_predictor_grid(brm_taxon_mac_fine_vel_prior, velocity_s, macrophyte_s, fines_s) %>% 
  expand_grid(taxon = unique(mod_d$taxon))

total_predictor_grid_prior = get_predictor_grid(brm_taxon_mac_fine_vel_prior, velocity_s, macrophyte_s, fines_s)

# get priors

taxon_prior = taxon_predictor_grid_prior %>%
  add_epred_draws(brm_taxon_mac_fine_vel_prior, re_formula = ~ (1 + macrophyte_s + fines_s + velocity_s|taxon)) %>% 
  pivot_longer(cols = ends_with("_s")) %>% 
  filter(value != 0) %>% 
  mutate(taxon_order = case_when(taxon == "Pyrgulopsis" ~ "b) *Pyrgulopsis*",
                                 taxon == "Fossaria" ~ "c) *Fossaria*",
                                 TRUE ~ "d) *Physa*")) %>% 
  mutate(taxon = case_when(taxon == "Pyrgulopsis" ~ "*Pyrgulopsis*",
                           taxon == "Fossaria" ~ "*Fossaria*",
                           taxon == "Physa" ~ "*Physa*")) 


# bind priors 
prior_preds_regression = taxon_prior %>% 
  mutate(predictor = name,
         formatted_predictor = case_when(name == "fines_s" ~ "% Fine Sediment",
                                         name == "macrophyte_s" ~ "% Macrophytes",
                                         TRUE ~ "Stream Velocity"),
         numbered_predictor = case_when(name == "fines_s" ~ "a) % Fine Sediment",
                                        name == "macrophyte_s" ~ "b) % Macrophytes",
                                        TRUE ~ "c) Stream Velocity")) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (value*sd) + mean) %>% 
  mutate(model = "Prior")


# combine
prior_post_preds_regression = bind_rows(prior_preds_regression, post_preds_regression %>% 
                                          filter(taxon != "All Snails"))


# plot
library(ggh4x)
library(ggtext)

prior_post_plot = prior_post_preds_regression %>% 
  filter(.draw <= 100) %>% 
  ggplot(aes(x = pred_values, y = .epred + 1)) + 
  stat_lineribbon(.width = c(0.95), aes(fill = model),
                  alpha = 0.6, linewidth = 0.1) +
  facet_grid2(taxon~numbered_predictor, scales = "free_x") +
  scale_y_log10() +
  labs(x = "Predictor Value",
       y = "Total Snails per Quadrat + 1",
       fill = "",
       color = "") +
  ggthemes::scale_fill_colorblind() +
  theme_default() +
  theme(strip.text.y = element_markdown(),
        strip.text.x = element_text(hjust = 0),
        legend.position = "top") + 
  guides(colour = guide_legend(override.aes = list(alpha = 0.9)))+
  NULL

ggsave(prior_post_plot, file = "plots/prior_post_plot.tif", width = 6.5, height = 8, bg = "white")
ggsave(prior_post_plot, file = "plots/prior_post_plot.jpg", width = 6.5, height = 8, bg = "white")


# prior sensitivity -------------------------------------------------------

# load model
brm_taxon_mac_fine_vel = readRDS(file = "models/brm_taxon_mac_fine_vel.rds")

# rerun with weaker priors
brm_taxon_mac_fine_vel_weak = update(brm_taxon_mac_fine_vel,
                                     prior = c(prior(normal(3, 4), class = "Intercept"),
                                               prior(exponential(2), class = "sd"),
                                               prior(normal(0, 2), class = "b")))

saveRDS(brm_taxon_mac_fine_vel_weak, file = "models/brm_taxon_mac_fine_vel_weak.rds")


# extract posteriors 
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
  theme(strip.text = element_markdown(hjust = 0)) +
  NULL

ggsave(plot_prior_sensitivity, file = "plots/plot_prior_sensitivity.jpg", width = 6.5, height = 3.5)
ggsave(plot_prior_sensitivity, file = "plots/plot_prior_sensitivity.tif", width = 6.5, height = 3.5)
