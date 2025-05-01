library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)

# load data -------------------------------------------------------
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

# load models -------------------------------------------------------
# brm_intonly = readRDS(file = "models/brm_intonly.rds")
# brm_int1 = readRDS(file = "models/brm_int.rds")
# brm_intonly_randzi = readRDS(file = "models/brm_intonly_randzi.rds")
# brm_velocity = readRDS(file = "models/brm_velocity.rds")
# brm_velocity_elev = readRDS(file = "models/brm_velocity_elev.rds")
# brm_velocity_fines = readRDS(file = "models/brm_velocity_fines.rds")
# brm_fines = readRDS(file = "models/brm_fines.rds")
# brm_macrophytes = readRDS("models/brm_macrophytes.rds")
# brm_mac_vel_fines = readRDS("models/brm_mac_vel_fines.rds")
brm_taxon_mac_fine_vel = readRDS("models/brm_taxon_mac_fine_vel.rds")

# calculate effects  -------------------------------------------------------
fixed_effects = as_draws_df(brm_taxon_mac_fine_vel) %>% 
  select(.draw, b_Intercept, ends_with("_s"), -starts_with("r_taxon"),
         -starts_with("sd_"), -starts_with("cor_")) %>% 
  pivot_longer(cols = ends_with("_s")) %>% 
  mutate(name = str_remove(name, "b_"),
         name = str_remove(name, "_s")) %>% 
  rename(b_slope = value)


random_effects = as_draws_df(brm_taxon_mac_fine_vel) %>% 
  select(.draw, -starts_with("b_"), starts_with("r_taxon"),
         -starts_with("sd_"), -starts_with("cor_")) %>% 
  pivot_longer(cols = -.draw) %>% 
  mutate(name = str_sub(name, 9, 30),
         name = str_remove(name, "]")) %>% 
  separate(name, into = c("taxon", "parameter")) %>% 
  pivot_wider(names_from = parameter, values_from = value) %>% 
  pivot_longer(cols = c(-.draw, -taxon, -Intercept)) %>% 
  rename(r_int = Intercept,
         r_slope = value)

post_int_slopes = fixed_effects %>% 
  right_join(random_effects) %>% 
  mutate(intercept = b_Intercept + r_int,
         slope = b_slope + r_slope)  


# summarize effects -------------------------------------------------------

slope_taxa_raw = post_int_slopes %>% 
  pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "parameter_value") %>% 
  group_by(taxon, name, parameter) %>% 
  mutate(exp_parameter_value = exp(parameter_value))

slope_overall_raw = post_int_slopes %>% 
  pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "parameter_value") %>% 
  group_by(.draw, name, parameter) %>% 
  reframe(parameter_value = mean(parameter_value)) %>% 
  mutate(exp_parameter_value = exp(parameter_value)) %>% 
  mutate(taxon = "total_snails")

slopes_raw = bind_rows(slope_taxa_raw, slope_overall_raw)

(slope_table_raw = slopes_raw %>% 
  group_by(name, parameter, taxon) %>% 
  median_qi(parameter_value) %>% 
  filter(parameter == "slope") %>% 
  select(name, taxon, parameter, parameter_value, .lower, .upper) %>% 
  mutate(across(where(is.numeric), round, 1)) %>% 
  mutate(parameters = paste0(parameter_value, " (", .lower, ",", .upper, ")")) %>% 
  select(-parameter_value, -.lower, -.upper) %>% 
  pivot_wider(names_from = name, values_from = parameters)
)

write_csv(slope_table_raw, file = "tables/slope_table_raw.csv")

slope_table_percentchange = post_int_slopes %>% 
  pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "parameter_value") %>% 
  group_by(taxon, name, parameter) %>% 
  mutate(exp_parameter_value = exp(parameter_value)) %>% 
  filter(parameter == "slope") %>% 
  mutate(percent_change = exp_parameter_value - 1) %>% 
  select(.draw, name, taxon, percent_change) %>% 
  median_qi(parameter_value = percent_change) %>% 
  select(taxon, name, parameter, parameter_value, .lower, .upper) %>% 
  mutate(across(where(is.numeric), round, 1)) %>% 
  mutate(parameters = paste0(parameter_value, " (", .lower, ",", .upper, ")")) %>% 
  select(-parameter_value, -.lower, -.upper) %>% 
  pivot_wider(names_from = name, values_from = parameters)

slope_table_percentchange_perunit = post_int_slopes %>% 
  pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "parameter_value") %>% 
  group_by(taxon, name, parameter) %>% 
  mutate(exp_parameter_value = exp(parameter_value)) %>% 
  filter(parameter == "slope") %>% 
  mutate(percent_change = exp_parameter_value) %>% 
  select(.draw, name, taxon, percent_change) %>% 
  left_join(predictors_mean_sd %>% mutate(name = str_remove(predictor, "_s")) %>% select(-predictor)) %>% 
  mutate(percent_change = 1 - (percent_change/sd)) %>% 
  median_qi(parameter_value = percent_change) %>% 
  select(taxon, name, parameter, parameter_value, .lower, .upper) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  mutate(parameters = paste0(parameter_value, " (", .lower, ",", .upper, ")")) %>% 
  select(-parameter_value, -.lower, -.upper) %>% 
  pivot_wider(names_from = name, values_from = parameters)

# summarize probabilities -------------------------------------------------

slopes_raw %>% 
  group_by(taxon, name, parameter) %>% 
  filter(parameter == "slope") %>% 
  reframe(prob_negative = sum(parameter_value < 0)/max(.draw)) %>% 
  arrange(name)


slopes_raw %>%
  select(.draw, taxon, parameter_value) %>% 
  pivot_wider(names_from = taxon, values_from = parameter_value) %>% 
  mutate(p_f = Pyrgulopsis - Fossaria,
         pyr_ph = Pyrgulopsis - Physa) %>% 
  group_by(parameter, name) %>% 
  reframe(prob_higher = sum(pyr_ph > 0)/max(.draw))

# check plotting ----------------------------------------------------------

post_int_slopes %>% 
  expand_grid(x = seq(-2, 2, length.out = 30)) %>% 
  mutate(.epred = exp(intercept + slope*x)) %>% 
  ggplot(aes(x = x, y = .epred, fill = name)) + 
  stat_lineribbon(.width = 0.95) + 
  facet_wrap(taxon ~ name) +
  scale_y_log10()
  



# capture probabilities ---------------------------------------------------

capture_probs = readRDS(file = "posteriors/capture_probs.rds")

capture_probs %>% group_by(site_f, taxon) %>%
  median_qi(1-zi) %>% 
  filter(site_f == "1" | site_f == "10")


# site densities ----------------------------------------------------------

site_posts_all = readRDS(file = "posteriors/site_posts_all.rds")

table_density_site = site_posts_all %>% 
  group_by(taxon, site_f) %>% 
  median_qi(.epred) %>% 
  rename(median = .epred,
         low95 = .lower,
         high95 = .upper) %>% 
  select(-.width, -.point, -.interval) %>% 
  mutate(taxon = str_remove(taxon, "\\*"),
         taxon = str_remove(taxon, "\\*"))

write_csv(table_density_site, file = "tables/table_density_site.csv")
