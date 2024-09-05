library(tidyverse)
library(brms)

# load data
snail_density = readRDS("data/snail_density.rds")
snail_density_long = snail_density %>% 
  pivot_longer(cols = c(Fossaria, Physa, Pyrgulopsis), names_to = "taxon", values_to = "taxon_snail")

# fit prior model ------------------------------------------------------

brm_prior = brm(bf(taxon_snail ~ 1 + macrophyte_s + fines_s + velocity_s + 
                     (1 + macrophyte_s + fines_s + velocity_s|taxon) + (1|site_f/transect_f/quadrat_f),
                   zi ~ 1 + (1 + taxon|site_f)),
                data = snail_density_long,
                family = zero_inflated_negbinomial(),
                prior = c(prior(normal(3, 2), class = "Intercept"),
                          prior(exponential(4), class = "sd"),
                          prior(normal(0, 1), class = "b")),
                chains = 1, iter = 1000,
                sample_prior = "only")

saveRDS(brm_prior, file = "models/brm_prior.rds")


# fit models by taxon --------------------------------------------------------------

brm_velocity_fines_taxon = brm(bf(taxon_snail ~ 1 + velocity_s + fines_s +
                                   (1 + velocity_s + fines_s|taxon) + (1|site_f/transect_f/quadrat_f),
                           zi ~ 1 + (1 + taxon|site_f)),
                        data = snail_density_long,
                        family = zero_inflated_negbinomial(),
                        prior = c(prior(normal(3, 2), class = "Intercept"),
                                  prior(exponential(4), class = "sd"),
                                  prior(normal(0, 1), class = "b")),
                        chains = 4, iter = 2000)

saveRDS(brm_velocity_fines_taxon, file = "models/brm_velocity_fines_taxon.rds")


brm_taxon_int = update(brm_velocity_fines_taxon, 
                       formula = bf(taxon_snail ~ 1 +
                                      (1 |taxon) + (1|site_f/transect_f/quadrat_f),
                                    zi ~ 1 + (1 + taxon|site_f)),
                       chains = 4, iter = 2000)

saveRDS(brm_taxon_int, file = "models/brm_taxon_int.rds")

brm_taxon_macrophyte = update(brm_velocity_fines_taxon, 
                              formula = bf(taxon_snail ~ 1 + macrophyte_s +
                                             (1 + macrophyte_s|taxon) + (1|site_f/transect_f/quadrat_f),
                                           zi ~ 1 + (1 + taxon|site_f)),
                              iter = 2000, chains = 4, newdata = snail_density_long)

saveRDS(brm_taxon_macrophyte, file = "models/brm_taxon_macrophyte.rds")

brm_taxon_velocity = update(brm_velocity_fines_taxon, 
                              formula = bf(taxon_snail ~ 1 + velocity_s +
                                             (1 + velocity_s|taxon) + (1|site_f/transect_f/quadrat_f),
                                           zi ~ 1 + (1 + taxon|site_f)),
                              iter = 2000, chains = 4, newdata = snail_density_long)

saveRDS(brm_taxon_velocity, file = "models/brm_taxon_velocity.rds")

brm_taxon_fines = update(brm_velocity_fines_taxon, 
                            formula = bf(taxon_snail ~ 1 + fines_s +
                                           (1 + fines_s|taxon) + (1|site_f/transect_f/quadrat_f),
                                         zi ~ 1 + (1 + taxon|site_f)),
                            iter = 2000, chains = 4, newdata = snail_density_long)

saveRDS(brm_taxon_fines, file = "models/brm_taxon_fines.rds")

brm_taxon_velocity_mac = brm(bf(taxon_snail ~ 1 + velocity_s + macrophyte_s +
                                    (1 + velocity_s + macrophyte_s|taxon) + (1|site_f/transect_f/quadrat_f),
                                  zi ~ 1 + (1 + taxon|site_f)),
                               data = snail_density_long,
                               family = zero_inflated_negbinomial(),
                               prior = c(prior(normal(3, 2), class = "Intercept"),
                                         prior(exponential(4), class = "sd"),
                                         prior(normal(0, 1), class = "b")),
                               chains = 4, iter = 2000)

saveRDS(brm_taxon_velocity_mac, file = "models/brm_taxon_velocity_mac.rds")

brm_taxon_fines_mac = brm(bf(taxon_snail ~ 1 + macrophyte_s + fines_s +
                                    (1 + macrophyte_s + fines_s|taxon) + (1|site_f/transect_f/quadrat_f),
                                  zi ~ 1 + (1 + taxon|site_f)),
                               data = snail_density_long,
                               family = zero_inflated_negbinomial(),
                               prior = c(prior(normal(3, 2), class = "Intercept"),
                                         prior(exponential(4), class = "sd"),
                                         prior(normal(0, 1), class = "b")),
                               chains = 4, iter = 2000)

saveRDS(brm_taxon_fines_mac, file = "models/brm_taxon_fines_mac.rds")


brm_taxon_mac_fine_vel = update(brm_velocity_fines_taxon, 
                         formula = bf(taxon_snail ~ 1 + macrophyte_s + fines_s + velocity_s +
                                        (1 + macrophyte_s + fines_s + velocity_s|taxon) + (1|site_f/transect_f/quadrat_f),
                                      zi ~ 1 + (1 + taxon|site_f)),
                         iter = 2000, chains = 4, newdata = snail_density_long)

saveRDS(brm_taxon_mac_fine_vel, file = "models/brm_taxon_mac_fine_vel.rds")





# fit models for all snails -----------------------------------------------

brm_total_snails = brm(bf(total_snail ~ 1 + macrophyte_s + fines_s + velocity_s + (1|site_f/transect_f/quadrat_f),
       zi ~ 1 + (1|site_f)),
    data = snail_density,
    family = zero_inflated_negbinomial(),
    prior = c(prior(normal(3, 2), class = "Intercept"),
              prior(exponential(4), class = "sd"),
              prior(normal(0, 1), class = "b")),
    chains = 4, iter = 2000)

saveRDS(brm_total_snails, file = "models/brm_total_snails.rds")
