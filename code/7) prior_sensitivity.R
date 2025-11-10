library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)
library(ggh4x)
library(ggtext)

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
         formatted_predictor = case_when(name == "fines_s" ~ "% Fine sediment",
                                         name == "macrophyte_s" ~ "% Macrophytes",
                                         TRUE ~ "Stream Velocity"),
         numbered_predictor = case_when(name == "fines_s" ~ "A % Fine sediment",
                                        name == "macrophyte_s" ~ "B % Macrophytes",
                                        TRUE ~ "C Stream velocity")) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (value*sd) + mean) %>% 
  mutate(model = "Prior")


# combine
prior_post_preds_regression = bind_rows(prior_preds_regression, post_preds_regression) %>% 
  filter(taxon != "All snails") %>% 
  mutate(taxon = case_when(taxon == "*Fossaria*" ~ "*Galba*", 
                           T ~ taxon)) %>% 
  mutate(taxon_order = case_when(taxon == "All snails" ~ "A All snails",
                                 taxon == "Pyrgulopsis" ~ "B *Pyrgulopsis*",
                                 taxon == "Galba" ~ "C *Galba*",
                                 taxon == "Fossaria" ~ "C *Galba*",
                                 TRUE ~ "D *Physa*"))



make_prior_post = function(posteriors = NA){
  
  plot_data = posteriors 
  
  plot_data %>% 
    filter(.draw <= 1000) %>% 
    ggplot(aes(x = pred_values, y = .epred*4 + 1)) + 
    stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = model),
                    alpha = 0.6, linewidth = 0.1) +
    facet_grid2(taxon ~ numbered_predictor, scales = "free_x") +
    scale_y_log10(labels = scales::comma) +
    labs(x = "Predictor value",
         y = bquote("Snail density (no./m"^2 ~ "+ 1)"),
         fill = "",
         color = "") +
    guides(color = "none") +
    ggthemes::scale_fill_colorblind() +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.y = element_markdown(),
          strip.text.x = element_text(hjust = 0)) +
    NULL
}

prior_post_preds_regression_list = prior_post_preds_regression %>% group_by(taxon) %>% group_split()

figure_s1_list = list()

for(i in 1:length(prior_post_preds_regression_list)){
  figure_s1_list[[i]] = make_prior_post(posteriors = prior_post_preds_regression_list[[i]])
}

library(patchwork)

(plot_density_regression_prior = (figure_s1_list[[3]] + theme(axis.text.x = element_blank()))/
    # (figure_s1_list[[4]] + theme(axis.text.x = element_blank(),
    #                             strip.text.x = element_blank()))/
    (figure_s1_list[[1]] + theme(axis.text.x = element_blank(),
                                 strip.text.x = element_blank()))/
    (figure_s1_list[[2]] + theme(strip.text.x = element_blank()))/
    plot_layout(axis_titles = "collect",
                guides = "collect")
)


ggsave(plot_density_regression_prior, file = "plots/fig_s1_plot_density_regression.tif", width = 6.5, height = 7, bg = "white")
ggsave(plot_density_regression_prior, file = "plots/fig_s1_plot_density_regression.jpg", width = 6.5, height = 7, bg = "white")


# prior sensitivity -------------------------------------------------------

# load model
brm_total_snails = readRDS(file = "models/brm_total_snails.rds")

# rerun with weaker priors
# brm_total_snails_weak = update(brm_total_snails,
#                                      prior = c(prior(normal(3, 4), class = "Intercept"),
#                                                prior(exponential(2), class = "sd"),
#                                                prior(normal(0, 2), class = "b")))
# 
# saveRDS(brm_total_snails_weak, file = "models/brm_total_snails_weak.rds")
brm_total_snails_weak = readRDS(file = "models/brm_total_snails_weak.rds")


# extract posteriors 
brm_total_snails_weak = readRDS(file = "models/brm_total_snails_weak.rds")

mod_d = brm_total_snails_weak$data

# get posteriors
velocity_preds_weak = tibble(velocity_s = seq(min(mod_d$velocity_s),
                                              max(mod_d$velocity_s),
                                              length.out = 30)) %>% 
  mutate(macrophyte_s = 0,
         fines_s = 0) %>% 
  add_epred_draws(brm_total_snails_weak, re_formula = NA) %>% 
  mutate(predictor = "velocity_s",
         model = "weak priors")
  
  
macrophyte_preds_weak = tibble(macrophyte_s = seq(min(mod_d$macrophyte_s),
                                              max(mod_d$macrophyte_s),
                                              length.out = 30)) %>% 
  mutate(velocity_s = 0,
         fines_s = 0) %>% 
  add_epred_draws(brm_total_snails_weak, re_formula = NA) %>% 
  mutate(predictor = "macrophyte_s",
         model = "weak priors")

fines_preds_weak = tibble(fines_s = seq(min(mod_d$fines_s),
                                                 max(mod_d$fines_s),
                                                 length.out = 30)) %>% 
  mutate(velocity_s = 0,
         macrophyte_s = 0) %>% 
  add_epred_draws(brm_total_snails_weak, re_formula = NA) %>% 
  mutate(predictor = "fines_s",
         model = "weak priors")

velocity_preds_main = tibble(velocity_s = seq(min(mod_d$velocity_s),
                                              max(mod_d$velocity_s),
                                              length.out = 30)) %>% 
  mutate(macrophyte_s = 0,
         fines_s = 0) %>% 
  add_epred_draws(brm_total_snails, re_formula = NA) %>% 
  mutate(predictor = "velocity_s",
         model = "main priors")


macrophyte_preds_main = tibble(macrophyte_s = seq(min(mod_d$macrophyte_s),
                                                   max(mod_d$macrophyte_s),
                                                   length.out = 30)) %>% 
  mutate(velocity_s = 0,
         fines_s = 0) %>% 
  add_epred_draws(brm_total_snails, re_formula = NA) %>% 
  mutate(predictor = "macrophyte_s",
         model = "main priors")

fines_preds_main = tibble(fines_s = seq(min(mod_d$fines_s),
                                        max(mod_d$fines_s),
                                        length.out = 30)) %>% 
  mutate(velocity_s = 0,
         macrophyte_s = 0) %>% 
  add_epred_draws(brm_total_snails, re_formula = NA) %>% 
  mutate(predictor = "fines_s",
         model = "main priors")



all_main_preds = bind_rows(velocity_preds_weak, 
                           macrophyte_preds_weak, 
                           fines_preds_weak,
                           velocity_preds_main, 
                           macrophyte_preds_main, 
                           fines_preds_main) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = case_when(predictor == "velocity_s" ~ velocity_s,
                                 predictor == "macrophyte_s" ~ macrophyte_s,
                                 TRUE ~ fines_s)) %>% 
  mutate(pred_values = (pred_values*sd) + mean,
         raw_predictor = case_when(predictor == "fines_s" ~ "A % Fine sediment",
                                   predictor == "macrophyte_s" ~ "B % Macrophytes",
                                   TRUE ~ "C Stream velocity"))
  


# plot --------------------------------------------------------------------

plot_prior_sensitivity1 = all_main_preds %>% 
  mutate(model = case_when(model == "main priors" ~ "Main model priors",
                           model == "weak priors" ~ "Weak priors")) %>%
  filter(model == "Main model priors") %>% 
  ggplot(aes(x = pred_values, y = .epred*4 + 1)) + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = model),
                  alpha = 0.3, linewidth = 0.1) + 
  ggh4x::facet_grid2(model~raw_predictor, scales = "free_x") +
  # geom_point(data = pred_data, aes(y = total_snail + 1),
             # size = 0.1, shape = ".") +
  scale_y_log10() +
  labs(x = "Predictor Value",
       y = bquote("Snail density (no./m"^2 ~ "+ 1)"),
       fill = "",
       color = "") +
  guides(fill = "none",
         color = "none") +
  scale_fill_brewer(type = "qual") +
  theme_default() +
  theme(strip.text = element_markdown(hjust = 0)) +
  NULL

plot_prior_sensitivity2 = all_main_preds %>% 
  mutate(model = case_when(model == "main priors" ~ "Main model priors",
                           model == "weak priors" ~ "Weak priors")) %>%
  filter(model != "Main model priors") %>% 
  ggplot(aes(x = pred_values, y = .epred*4 + 1)) + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), fill = "purple3",
                  alpha = 0.3, linewidth = 0.1) + 
  ggh4x::facet_grid2(model~raw_predictor, scales = "free_x") +
  # geom_point(data = pred_data, aes(y = total_snail + 1),
  # size = 0.1, shape = ".") +
  scale_y_log10() +
  labs(x = "Predictor Value",
       y = bquote("Snail density (no./m"^2 ~ "+ 1)"),
       fill = "",
       color = "") +
  guides(fill = "none",
         color = "none") +
  scale_fill_brewer(type = "qual") +
  theme_default() +
  theme(strip.text = element_markdown(hjust = 0)) +
  NULL

library(patchwork)

plot_prior_sensitivity = (plot_prior_sensitivity1 + theme(axis.text.x = element_blank()))/(plot_prior_sensitivity2 + theme(strip.text.x = element_blank())) +
  plot_layout(axis_titles = "collect")

ggsave(plot_prior_sensitivity, file = "plots/plot_prior_sensitivity.jpg", width = 6.5, height = 3.5)
ggsave(plot_prior_sensitivity, file = "plots/plot_prior_sensitivity.tif", width = 6.5, height = 3.5)
