library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)
library(ggh4x)
library(ggtext)


source("code/custom_functions.R")

# load data ---------------------------------------------------------------
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

# load models
brm_taxon_mac_fine_vel = readRDS("models/brm_taxon_mac_fine_vel.rds")
brm_total_snails = readRDS("models/brm_total_snails.rds")

# density -------------------------------------------------------------

# get model data
mod_d = brm_taxon_mac_fine_vel$data

# get raw data
snail_density = readRDS(file = "data/snail_density.rds")

# make predictors
taxon_predictor_grid = get_predictor_grid(brm_total_snails, velocity_s, macrophyte_s, fines_s) %>% 
  expand_grid(taxon = unique(mod_d$taxon))

total_predictor_grid = get_predictor_grid(brm_total_snails, velocity_s, macrophyte_s, fines_s)

# get posteriors

taxon_posts = taxon_predictor_grid %>%
  add_epred_draws(brm_taxon_mac_fine_vel, re_formula = ~ (1 + macrophyte_s + fines_s + velocity_s|taxon)) %>% 
  pivot_longer(cols = ends_with("_s")) %>% 
  filter(value != 0) %>% 
  mutate(taxon_order = case_when(taxon == "Pyrgulopsis" ~ "b) *Pyrgulopsis*",
                                 taxon == "Fossaria" ~ "c) *Fossaria*",
                                 TRUE ~ "d) *Physa*")) %>% 
  mutate(taxon = case_when(taxon == "Pyrgulopsis" ~ "*Pyrgulopsis*",
                           taxon == "Fossaria" ~ "*Fossaria*",
                           taxon == "Physa" ~ "*Physa*",
                           TRUE ~ "All Snails")) 

total_posts = total_predictor_grid %>%
  add_epred_draws(brm_total_snails, re_formula = NA) %>% 
  pivot_longer(cols = ends_with("_s")) %>% 
  filter(value != 0) %>% 
  mutate(taxon = "All Snails",
         taxon_order = "a) All Snails") 


# bind posteriors 
post_preds_regression = bind_rows(total_posts, taxon_posts) %>% 
         mutate(predictor = name,
                formatted_predictor = case_when(name == "fines_s" ~ "% Fine Sediment",
                                     name == "macrophyte_s" ~ "% Macrophytes",
                                     TRUE ~ "Stream Velocity"),
                numbered_predictor = case_when(name == "fines_s" ~ "a) % Fine Sediment",
                                   name == "macrophyte_s" ~ "b) % Macrophytes",
                                   TRUE ~ "c) Stream Velocity")) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (value*sd) + mean) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "All Snails"))

saveRDS(post_preds_regression, file = "posteriors/post_preds_regression.rds")

# make raw data similar to posteriors

pred_data_taxon = mod_d %>% 
  pivot_longer(cols = c(velocity_s, macrophyte_s, fines_s),
               values_to = "pred_values",
               names_to = "predictor")  %>% 
  mutate(name = predictor,
         formatted_predictor = case_when(name == "fines_s" ~ "% Fine Sediment",
                                         name == "macrophyte_s" ~ "% Macrophytes",
                                         TRUE ~ "Stream Velocity"),
         numbered_predictor = case_when(name == "fines_s" ~ "a) % Fine Sediment",
                                        name == "macrophyte_s" ~ "b) % Macrophytes",
                                        TRUE ~ "c) Stream Velocity")) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (pred_values*sd) + mean, 
         .epred = taxon_snail) %>% 
  mutate(taxon = case_when(taxon == "Pyrgulopsis" ~ "*Pyrgulopsis*",
                           taxon == "Fossaria" ~ "*Fossaria*",
                           taxon == "Physa" ~ "*Physa*",
                           TRUE ~ "All Snails")) 


pred_data_total = pred_data_taxon %>% group_by(.epred, pred_values, site_f, quadrat_f, transect_f, 
                                               name, predictor, formatted_predictor, numbered_predictor) %>% 
  reframe(.epred = sum(.epred)) %>% 
  mutate(taxon = "All Snails")

pred_data_all = bind_rows(pred_data_taxon, pred_data_total) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "All Snails"))

# plot_posteriors 

plot_density_regression = post_preds_regression %>% 
  filter(.draw <= 1000) %>% 
  ggplot(aes(x = pred_values, y = .epred + 1)) + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = predictor),
                  alpha = 0.3, linewidth = 0.1) +
  facet_grid2(taxon~numbered_predictor, scales = "free_x") +
  geom_point(data = pred_data_all,
             size = 0.1, shape = ".") +
  scale_y_log10() +
  labs(x = "Predictor Value",
       y = "Total Snails per Quadrat + 1",
       fill = "",
       color = "") +
  guides(fill = "none",
         color = "none") +
  ggthemes::scale_fill_colorblind() +
  theme_default() +
  theme(strip.text.y = element_markdown(),
        strip.text.x = element_text(hjust = 0)) +
  NULL

ggsave(plot_density_regression, file = "plots/plot_density_regression.tif", width = 6.5, height = 8, bg = "white")
ggsave(plot_density_regression, file = "plots/plot_density_regression.jpg", width = 6.5, height = 8, bg = "white")

# capture probabilities -----------------------------------------------
capture_prob_taxon = mod_d %>% 
  select(taxon, site_f, taxon_snail) %>% 
  mutate(zi = 1 - case_when(taxon_snail == 0 ~ 0, TRUE ~ 1)) %>% 
  mutate(taxon_order = case_when(taxon == "Pyrgulopsis" ~ "c) *Pyrgulopsis*",
                                 taxon == "Fossaria" ~ "b) *Fossaria*",
                                 TRUE ~ "d) *Physa*"))

capture_prob_data = bind_rows(capture_prob_taxon, capture_prob_taxon %>% mutate(taxon_order = "a) All Snails"))

capture_probs_taxa = mod_d %>% 
  select(taxon, site_f) %>% 
  distinct() %>% 
  mutate(velocity_s = 0, macrophyte_s = 0, fines_s = 0) %>% 
  add_epred_draws(brm_taxon_mac_fine_vel, allow_new_levels = T, re_formula = ~ (1 + taxon|site_f), dpar = "zi") %>% 
  left_join(capture_prob_taxon %>% ungroup %>% distinct(taxon, taxon_order))

capture_probs_overall = tibble(site_f = unique(snail_density$site_f)) %>%
  mutate(macrophyte_s = 0, fines_s = 0, velocity_s = 0) %>% 
  add_epred_draws(brm_total_snails, re_formula = ~(1|site_f), dpar = "zi") %>% 
  mutate(taxon = "All Snail Taxa",
         taxon_order = "a) All Snails")

capture_probs = bind_rows(capture_probs_taxa, capture_probs_overall)

saveRDS(capture_probs, file = "posteriors/capture_probs.rds")

plot_capture_probs = capture_probs %>% 
  ggplot(aes(x = site_f, y = 1 - zi, fill = taxon_order)) +
  geom_line(data = . %>% filter(.draw <= 500),
            aes(group = .draw, color = taxon_order), alpha = 0.02) + 
  # geom_violin(linewidth = 0.2) + 
  geom_boxplot(aes(group = site_f), outlier.shape = NA, width = 0.5) +
  geom_jitter(data = capture_prob_data, width = 0.2, height = 0,
              alpha = 0.6, shape = "|") +
  labs(x = "(upstream)                Site                (downstream)",
       y = "Capture Probability") +
  guides(color = "none",
         fill = "none") +
  facet_wrap(~taxon_order, ncol = 1) +
  scale_color_brewer(type = 'qual') + 
  scale_fill_brewer(type = 'qual') +
  theme_default() +
  theme(strip.text = element_markdown(hjust = 0)) +
  NULL

ggsave(plot_capture_probs, file = "plots/plot_capture_probs.jpg", width = 5, height = 9)
ggsave(plot_capture_probs, file = "plots/plot_capture_probs.tif", width = 5, height = 9, bg = "white")

