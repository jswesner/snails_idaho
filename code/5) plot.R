library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)
library(ggh4x)
library(ggtext)
library(patchwork)


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

# make predictor
taxon_predictor_grid = get_predictor_grid(brm_total_snails, velocity_s, macrophyte_s, fines_s) %>% 
  expand_grid(taxon = unique(mod_d$taxon))

total_predictor_grid = get_predictor_grid(brm_total_snails, velocity_s, macrophyte_s, fines_s)

# get posteriors

taxon_posts = taxon_predictor_grid %>%
  add_epred_draws(brm_taxon_mac_fine_vel, re_formula = ~ (1 + macrophyte_s + fines_s + velocity_s|taxon)) %>% 
  pivot_longer(cols = ends_with("_s")) %>% 
  filter(value != 0) %>% 
  mutate(taxon_order = case_when(taxon == "Pyrgulopsis" ~ "B *Pyrgulopsis*",
                                 taxon == "Galba" ~ "C *Galba*",
                                 taxon == "Fossaria" ~ "C *Galba*",
                                 TRUE ~ "D *Physa*")) %>% 
  mutate(taxon = case_when(taxon == "Pyrgulopsis" ~ "*Pyrgulopsis*",
                           taxon == "Galba" ~ "*Galba*",
                           taxon == "Fossaria" ~ "*Galba*",
                           taxon == "Physa" ~ "*Physa*",
                           TRUE ~ "All snails")) 

total_posts = total_predictor_grid %>%
  add_epred_draws(brm_total_snails, re_formula = NA) %>% 
  pivot_longer(cols = ends_with("_s")) %>% 
  filter(value != 0) %>% 
  mutate(taxon = "All snails",
         taxon_order = "A All snails") 


# bind posteriors 
post_preds_regression = bind_rows(total_posts, taxon_posts) %>% 
         mutate(predictor = name,
                formatted_predictor = case_when(name == "fines_s" ~ "% Fine sediment",
                                     name == "macrophyte_s" ~ "% Macrophytes",
                                     TRUE ~ "Stream velocity"),
                numbered_predictor = case_when(name == "fines_s" ~ "A % Fine sediment",
                                   name == "macrophyte_s" ~ "B % Macrophytes",
                                   TRUE ~ "C Stream velocity")) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (value*sd) + mean) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "All snails"))

saveRDS(post_preds_regression, file = "posteriors/post_preds_regression.rds")

# make raw data similar to posteriors

pred_data_taxon = mod_d %>% 
  pivot_longer(cols = c(velocity_s, macrophyte_s, fines_s),
               values_to = "pred_values",
               names_to = "predictor")  %>% 
  mutate(name = predictor,
         formatted_predictor = case_when(name == "fines_s" ~ "% Fine sediment",
                                         name == "macrophyte_s" ~ "% Macrophytes",
                                         TRUE ~ "Stream velocity"),
         numbered_predictor = case_when(name == "fines_s" ~ "A % Fine sediment",
                                        name == "macrophyte_s" ~ "B % Macrophytes",
                                        TRUE ~ "C Stream velocity")) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (pred_values*sd) + mean, 
         .epred = taxon_snail) %>% 
  mutate(taxon = case_when(taxon == "Pyrgulopsis" ~ "*Pyrgulopsis*",
                           taxon == "Galba" ~ "*Galba*",
                           taxon == "Fossaria" ~ "*Galba*",
                           taxon == "Physa" ~ "*Physa*",
                           TRUE ~ "All snails")) 


pred_data_total = pred_data_taxon %>% group_by(.epred, pred_values, site_f, quadrat_f, transect_f, 
                                               name, predictor, formatted_predictor, numbered_predictor) %>% 
  reframe(.epred = sum(.epred)) %>% 
  mutate(taxon = "All snails")

pred_data_all = bind_rows(pred_data_taxon, pred_data_total) %>% 
  mutate(taxon = as.factor(taxon),
         taxon = fct_relevel(taxon, "All snails"))

# plot_posteriors 

make_fig2 = function(posteriors = NA, 
         raw_data = NA){
  
  plot_data = posteriors
  
  plot_data %>% 
    filter(.draw <= 1000) %>% 
    ggplot(aes(x = pred_values, y = .epred*4 + 1)) + 
    stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = predictor),
                    alpha = 0.3, linewidth = 0.1) +
    facet_grid2(taxon~numbered_predictor, scales = "free_x") +
    geom_point(data = pred_data_all %>% 
                 filter(taxon == unique(plot_data$taxon)),
               size = 0.1, shape = ".") +
    scale_y_log10() +
    labs(x = "Predictor value",
         y = bquote("Total snails/m"^2 ~ "+ 1"),
         fill = "",
         color = "") +
    guides(fill = "none",
           color = "none") +
    ggthemes::scale_fill_colorblind() +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.y = element_markdown(),
          strip.text.x = element_text(hjust = 0)) +
    NULL
}

post_preds_regression_list = post_preds_regression %>% group_by(taxon) %>% group_split()

figure_2_list = list()

for(i in 1:length(post_preds_regression_list)){
  figure_2_list[[i]] = make_fig2(posteriors = post_preds_regression_list[[i]])
}

plot_density_regression = (figure_2_list[[1]] + theme(axis.text.x = element_blank()))/
  (figure_2_list[[4]] + theme(axis.text.x = element_blank(),
                              strip.text.x = element_blank()))/
  (figure_2_list[[2]] + theme(axis.text.x= element_blank(),
                              strip.text.x = element_blank()))/
  (figure_2_list[[3]] + theme(axis.text.x = element_blank(),
                              strip.text.x = element_blank()))/
  plot_layout(axis_titles = "collect")


ggsave(plot_density_regression, file = "plots/fig_2_plot_density_regression.tif", width = 6.5, height = 8, bg = "white")
ggsave(plot_density_regression, file = "plots/fig_2_plot_density_regression.jpg", width = 6.5, height = 8, bg = "white")

# capture probabilities -----------------------------------------------
capture_prob_taxon = mod_d %>% 
  select(taxon, site_f, taxon_snail) %>% 
  mutate(zi = 1 - case_when(taxon_snail == 0 ~ 0, TRUE ~ 1)) %>% 
  mutate(taxon_order = case_when(taxon == "Pyrgulopsis" ~ "B *Pyrgulopsis*",
                                 taxon == "Galba" ~ "C *Galba*",
                                 taxon == "Fossaria" ~ "C *Galba*",
                                 TRUE ~ "D *Physa*"))

capture_prob_data = bind_rows(capture_prob_taxon, capture_prob_taxon %>% mutate(taxon_order = "A All snails",
                                                                                taxon = "All snails")) 

capture_probs_taxa = mod_d %>% 
  select(taxon, site_f) %>% 
  distinct() %>% 
  mutate(velocity_s = 0, macrophyte_s = 0, fines_s = 0) %>% 
  add_epred_draws(brm_taxon_mac_fine_vel, allow_new_levels = T, re_formula = ~ (1 + taxon|site_f), dpar = "zi") %>% 
  left_join(capture_prob_taxon %>% ungroup %>% distinct(taxon, taxon_order))

capture_probs_overall = tibble(site_f = unique(snail_density$site_f)) %>%
  mutate(macrophyte_s = 0, fines_s = 0, velocity_s = 0) %>% 
  add_epred_draws(brm_total_snails, re_formula = ~(1|site_f), dpar = "zi") %>% 
  mutate(taxon = "All snails",
         taxon_order = "A All snails")

capture_probs = bind_rows(capture_probs_taxa, capture_probs_overall)

saveRDS(capture_probs, file = "posteriors/capture_probs.rds")



make_fig3 = function(posteriors = NA, 
                     raw_data = NA){
  
  plot_data = posteriors
  colors = unique(plot_data$fill)
  
  plot_data %>% 
    ggplot(aes(x = site_f, y = 1 - zi), fill = colors) +
    geom_line(data = . %>% filter(.draw <= 500),
              aes(group = .draw),
              color = colors,
              alpha = 0.02) + 
    # geom_violin(linewidth = 0.2) + 
    geom_boxplot(aes(group = site_f), outlier.shape = NA, width = 0.5) +
    geom_jitter(data = capture_prob_data %>% 
                  filter(taxon == unique(plot_data$taxon)), width = 0.2, height = 0,
                alpha = 0.6, shape = "|") +
    labs(x = "Upstream    \u27F5  Reach  \u27F6      Downstream",
         y = "Occurrence probability") +
    guides(color = "none",
           fill = "none") +
    facet_wrap(~taxon_order, ncol = 1) +
    scale_color_brewer(type = 'qual') + 
    scale_fill_brewer(type = 'qual') +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_markdown(hjust = 0)) +
    NULL
  
}

capture_probs_list = capture_probs %>% 
  group_by(taxon_order) %>% 
  mutate(PANEL = case_when(taxon == "Pyrgulopsis" ~ "2",
                           taxon == "Fossaria" ~ "3",
                           taxon == "Physa" ~ "4",
                           TRUE ~ "1")) %>% 
  left_join(color_list) %>% 
  group_split()

figure_3_list = list()

for(i in 1:length(capture_probs_list)){
  figure_3_list[[i]] = make_fig3(posteriors = capture_probs_list[[i]])
}

plot_capture_probs = (figure_3_list[[1]] + theme(axis.text.x = element_blank()))/
  (figure_3_list[[2]] + theme(axis.text.x = element_blank()))/
  (figure_3_list[[3]] + theme(axis.text.x= element_blank()))/
  (figure_3_list[[4]] + theme(axis.text.x = element_blank()))/
  plot_layout(axis_titles = "collect")

ggsave(plot_capture_probs, file = "plots/fig_3_plot_capture_probs.jpg", width = 5, height = 9)
ggsave(plot_capture_probs, file = "plots/fig_3_plot_capture_probs.tif", width = 5, height = 9, bg = "white")


# mean density by site ----------------------------------------------------

site_posts = brm_total_snails$data %>% 
  select(-total_snail) %>% 
  distinct() %>% 
  add_epred_draws(brm_total_snails) %>% 
  group_by(.draw, site_f) %>% 
  reframe(.epred = mean(.epred)) %>% 
  mutate(taxon = "All snails")

site_posts_taxon = brm_taxon_mac_fine_vel$data %>% 
  select(-taxon_snail) %>% 
  distinct() %>% 
  add_epred_draws(brm_taxon_mac_fine_vel) %>% 
  group_by(.draw, site_f, taxon) %>% 
  reframe(.epred = mean(.epred))

site_posts_all = bind_rows(site_posts, site_posts_taxon) %>% 
  mutate(taxon = case_when(taxon == "Pyrgulopsis" ~ "*Pyrgulopsis*",
                           taxon == "Galba" ~ "*Galba*",
                           taxon == "Fossaria" ~ "*Galba*",
                           taxon == "Physa" ~ "*Physa*",
                           TRUE ~ "All snails")) %>% 
  left_join(post_preds_regression  %>% ungroup %>% distinct(taxon, taxon_order))

saveRDS(site_posts_all, file = "posteriors/site_posts_all.rds")

temp = pred_data_all %>% 
  left_join(post_preds_regression  %>% ungroup %>% distinct(taxon, taxon_order)) %>% 
  ggplot(aes(x = site_f, y = .epred + 1)) + 
  geom_point(aes(fill = taxon_order), shape = 21, color = "black",
             size = 1) +
  facet_wrap(~taxon_order) +
  tidybayes::stat_pointinterval(data = site_posts_all) + 
  scale_y_log10() +
  labs(x = "(upstream)  <---------  Reach --------->       (downstream)",
       y = "Total Snails per Quadrat + 1",
       fill = "",
       color = "") +
  guides(fill = "none",
         color = "none") +
  scale_fill_brewer(type = 'qual') + 
  theme_classic() +
  theme(strip.text = element_markdown(hjust = 0),
        strip.background = element_blank()) +
  NULL


temp_data = ggplot_build(temp)$data

color_list = temp_data[[1]] %>% distinct(fill, PANEL)

site_posts_all = readRDS(file = "posteriors/site_posts_all.rds")



make_fig_1 = function(posteriors = NA, raw_data = NA){ 
  plot_data = posteriors
  colors = unique(plot_data$fill)
  
  plot_data %>% 
  ggplot(aes(x = site_f, y = .epred*4 + 1*4)) +
  geom_line(data = . %>% filter(.draw <= 500),
            aes(group = .draw),
            color = colors,
            alpha = 0.02) + 
  # geom_violin(linewidth = 0.2) +
    geom_boxplot(aes(group = site_f), outlier.shape = NA, width = 0.5,
                 fill = colors) +
  geom_jitter(data = raw_data %>% 
                left_join(post_preds_regression %>% ungroup %>% distinct(taxon, taxon_order)) %>%
                filter(taxon_order == unique(plot_data$taxon_order)),
              size = 0.4,
              width = 0.05, height = 0,
              alpha = 0.6) +
  labs(x = "",
       y = bquote("Total snails/m"^2 ~ "+ 1")) +
  guides(color = "none",
         fill = "none") +
  facet_wrap(~taxon_order) +
  scale_color_brewer(type = 'qual') + 
  scale_fill_brewer(type = 'qual') +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000"),
                limits = c(1, NA)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_markdown(hjust = 0)) +
  NULL
}

site_posts_all_list = site_posts_all %>% arrange(taxon_order) %>% 
  group_by(taxon_order) %>% 
  mutate(PANEL = case_when(taxon == "*Pyrgulopsis*" ~ "2",
                           taxon == "*Galba*" ~ "3",
                           taxon == "*Physa*" ~ "4",
                           TRUE ~ "1")) %>% 
  left_join(color_list) %>% 
  group_split()

fig_1_list = list()

for(i in 1:length(site_posts_all_list)){
  fig_1_list[[i]] = make_fig_1(raw_data = pred_data_all,
                               posteriors = site_posts_all_list[[i]])
}

library(patchwork)
plot_density_site = (fig_1_list[[1]] + theme(axis.text = element_blank()))/
 (fig_1_list[[2]] + theme(axis.text = element_blank()))/
  (fig_1_list[[3]] + theme(axis.text = element_blank()))/
  fig_1_list[[4]] +
  plot_layout(axis_titles = "collect")

plot_density_site = plot_density_site + labs(x = "Upstream    \u27F5  Reach  \u27F6      Downstream")

ggsave(plot_density_site, file = "plots/fig_1_plot_density_site.jpg", width = 5, height = 9)
ggsave(plot_density_site, file = "plots/fig_1_plot_density_site.tif", width = 5, height = 9, bg = "white")



