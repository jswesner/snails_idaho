library(tidyverse)
library(brms)

# load data ---------------------------------------------------------------
snail_density_long = readRDS("data/snail_density.rds") %>% 
  pivot_longer(cols = c(Fossaria, Physa, Pyrgulopsis), names_to = "taxon", values_to = "taxon_snail")


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

# load pre-fit models
brm_intonly = readRDS(file = "models/brm_intonly.rds")
brm_int1 = readRDS(file = "models/brm_int.rds")
brm_intonly_randzi = readRDS(file = "models/brm_intonly_randzi.rds")
brm_velocity = readRDS(file = "models/brm_velocity.rds")
brm_velocity_elev = readRDS(file = "models/brm_velocity_elev.rds")
brm_velocity_fines = readRDS(file = "models/brm_velocity_fines.rds")
brm_fines = readRDS(file = "models/brm_fines.rds")



library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)

# load data
snail_density = readRDS("data/snail_density.rds")

# prior predictive --------------------------------------------------------
# Hall et al. 2001 found an 'extreme' abundance of snails in Yellowstone of 483,000/m2
# So our prior model should be skeptical of a number that big. Our sample area was 0.25m2
# To get 483000m2 of snail in our sample area, we would need 0.25*483000 = 120,750 snails per sample
# So the prior should be skeptical of numbers in that range (i.e., intercepts > 11 after natural log)
# Hall Jr, R., Tank, J., & Dybdahl, M. (2001). Exotic snails dominate nitrogen and carbon cycling in a highly productive stream. UW-National Park Service Research Station Annual Reports, 25, 72-77.

# get a typical snail density using NEON data
invertebrate_taxonomy <- read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/aqua-sync_data/invertebrate-taxonomy.csv")

snail_size_data <- read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/aqua-sync_data/invertebrate-size-data.csv") %>% 
  left_join(invertebrate_taxonomy %>% distinct(taxon, class)) %>% 
  filter(class == "Gastropoda") %>% 
  group_by(sample, sampling_area) %>% 
  reframe(count = sum(count)) %>% 
  mutate(density = count/sampling_area,
         density_0.25m2 = density*0.25)

mean(snail_size_data$density_0.25m2)
median(snail_size_data$density_0.25m2)
sd(snail_size_data$density_0.25m2)

log(50)

prior_pred = brm(total_snail ~ 1 + (1|site_f/transect_f/quadrat_f),
                 data = snail_density,
                 family = zero_inflated_negbinomial(),
                 iter = 1000, chains = 1,
                 prior = c(prior(normal(3, 2), class = "Intercept"),
                           prior(exponential(4), class = "sd")),
                 sample_prior = "only")

saveRDS(prior_pred, file = "models/prior_pred.rds")

as_draws_df(prior_pred) %>% as_tibble() %>% 
  ggplot(aes(x = exp(b_Intercept))) +
  geom_histogram() +
  scale_x_log10() +
  geom_vline(xintercept = 50) +
  geom_vline(xintercept = 120750)

# fit prior model of taxa ------------------------------------------------------

brm_prior = brm(bf(taxon_snail ~ 1 + velocity_s + fines_s +
                     (1 + velocity_s + fines_s|taxon) + (1|site_f/transect_f/quadrat_f),
                   zi ~ 1 + (1 + taxon|site_f)),
                data = snail_density_long,
                family = zero_inflated_negbinomial(),
                prior = c(prior(normal(3, 2), class = "Intercept"),
                          prior(exponential(4), class = "sd"),
                          prior(normal(0, 1), class = "b")),
                chains = 1, iter = 1000,
                sample_prior = "only")

saveRDS(brm_prior, file = "models/brm_prior.rds")

mod_d = brm_prior$data

# plot
get_snail_posts <- function(model = brm_taxon_mac_fine_vel, 
                            model_data = mod_d, 
                            predictor1 = "velocity_s", 
                            predictor2 = "macrophyte_s",
                            predictor3 = "fines_s") {
  
  # Check if predictor columns exist and contain finite values
  if (!predictor1 %in% names(model_data)) {
    stop("Predictor1 column does not exist in model_data.")
  }
  
  predictor1_data <- model_data[[predictor1]]
  
  if (all(is.na(predictor1_data)) || !any(is.finite(predictor1_data))) {
    stop("Predictor1 column contains no finite values.")
  }
  
  # Dynamically unquote column names
  predictor1_sym <- sym(predictor1)
  predictor2_sym <- sym(predictor2)
  predictor3_sym <- sym(predictor3)
  
  pred_grid <- model_data %>%
    select(-taxon_snail) %>% 
    ungroup() %>% 
    distinct(taxon) %>%  
    expand_grid(!!predictor1_sym := seq(min(predictor1_data, na.rm = TRUE),
                                        max(predictor1_data, na.rm = TRUE),
                                        length.out = 30)) %>% 
    mutate(!!predictor2_sym := 0,
           !!predictor3_sym := 0, 
           site_f = "new",
           quadrat_f = "new",
           transect_f = "new")
  
  pred_grid %>% 
    add_epred_draws(model, allow_new_levels = TRUE, dpar = T,
                    re_formula =  ~ (1 + macrophyte_s + fines_s + velocity_s | taxon))
}


# get posteriors
velocity_priors = get_snail_posts(model = brm_prior, model_data = brm_prior$data, 
                                 predictor1 = "velocity_s", 
                                 predictor2 = "macrophyte_s",
                                 predictor3 = "fines_s") %>% 
  mutate(predictor = "velocity_s", 
         pred_values = velocity_s)

macrophyte_priors = get_snail_posts(model = brm_prior, model_data = brm_prior$data, 
                                   predictor2 = "velocity_s", 
                                   predictor1 = "macrophyte_s",
                                   predictor3 = "fines_s") %>% 
  mutate(predictor = "macrophyte_s", 
         pred_values = macrophyte_s)

fines_priors = get_snail_posts(model = brm_prior, model_data = brm_prior$data, 
                              predictor3 = "velocity_s", 
                              predictor2 = "macrophyte_s",
                              predictor1 = "fines_s") %>% 
  mutate(predictor = "fines_s",
         pred_values = fines_s)
  
post_priors_taxa = bind_rows(velocity_priors, macrophyte_priors, fines_priors) %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (pred_values*sd) + mean,
         raw_predictor = case_when(predictor == "fines_s" ~ "a) % Fine Sediment",
                                   predictor == "macrophyte_s" ~ "b) % Macrophytes",
                                   TRUE ~ "c) Stream Velocity"))


pred_data_taxon = mod_d %>% 
  pivot_longer(cols = c(velocity_s, macrophyte_s, fines_s),
               values_to = "pred_values",
               names_to = "predictor") %>% 
  left_join(predictors_mean_sd) %>% 
  mutate(pred_values = (pred_values*sd) + mean,
         raw_predictor = case_when(predictor == "fines_s" ~ "a) % Fine Sediment",
                                   predictor == "macrophyte_s" ~ "b) % Macrophytes",
                                   TRUE ~ "c) Stream Velocity")) 

# load posteriors to compare
post_preds_taxa = readRDS(file = "posteriors/post_preds_taxa.rds")

all_posts_taxa = bind_rows(post_preds_taxa %>% mutate(model = "posterior"), post_priors_taxa %>% mutate(model = "prior"))

saveRDS(all_posts_taxa, file = "posteriors/all_posts_taxa.rds")

# plot --------------------------------------------------------------------
all_posts_taxa = readRDS(file = "posteriors/all_posts_taxa.rds")


priors_density_taxa = all_posts_taxa %>% 
  ggplot(aes(x = pred_values, y = .epred + 1)) + 
  stat_lineribbon(.width = c(0.95), aes(fill = model),
                  alpha = 0.8, linewidth = 0.1) + 
  facet_grid2(taxon~raw_predictor, scales = "free") +
  # geom_point(data = pred_data_taxon, aes(y = taxon_snail + 1),
  #            size = 0.1, shape = ".") +
  scale_y_log10() +
  labs(x = "Predictor Value",
       y = "Total Snails per Quadrat",
       fill = "",
       color = "") +
  guides(color = "none") +
  ggthemes::scale_fill_colorblind() +
  theme_default() +
  theme(strip.text = element_text(hjust = 0)) +
  NULL

ggsave(priors_density_taxa, file = "plots/priors_density_taxa.jpg", width = 6.5, height = 6.5)
