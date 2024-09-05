library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)

# compare models ----------------------------------------------------
brm_velocity_fines_taxon = readRDS(file = "models/brm_velocity_fines_taxon.rds")
brm_taxon_int = readRDS(file = "models/brm_taxon_int.rds")
brm_taxon_macrophyte = readRDS(file = "models/brm_taxon_macrophyte.rds")
brm_taxon_velocity = readRDS(file = "models/brm_taxon_velocity.rds")
brm_taxon_fines = readRDS(file = "models/brm_taxon_fines.rds")
brm_taxon_fines_mac = readRDS(file = "models/brm_taxon_fines_mac.rds")
brm_taxon_velocity_mac = readRDS(file = "models/brm_taxon_velocity_mac.rds")
brm_taxon_mac_fine_vel = readRDS(file = "models/brm_taxon_mac_fine_vel.rds")

# waic_a_taxon = waic(brm_velocity_fines_taxon)
# waic_b_taxon = waic(brm_taxon_int)
# waic_c_taxon = waic(brm_taxon_macrophyte)
# waic_d_taxon = waic(brm_taxon_velocity)
# waic_e_taxon = waic(brm_taxon_fines)
# waic_f_taxon = waic(brm_taxon_mac_fine_vel)
# waic_g_taxon = waic(brm_taxon_fines_mac)
# waic_h_taxon = waic(brm_taxon_velocity_mac)
# 
# waic_list = list(waic_a_taxon,
#                  waic_b_taxon, 
#                  waic_c_taxon, 
#                  waic_d_taxon, 
#                  waic_e_taxon, 
                 # waic_f_taxon,
                 # waic_g_taxon,
                 # waig_f_taxon)
# 
# saveRDS(waic_list, file = "models/waic_list.rds")

waic_list = readRDS(file = "models/waic_list.rds")

loo_table = loo_compare(waic_list)

loo_table_csv = as_tibble(loo_table) %>% mutate(model = rownames(loo_table))

write_csv(loo_table_csv, file = "tables/loo_table.csv")

# quick calculation of confint of diff from 2nd model
-0.3 + c(-2, 2)*1.4

