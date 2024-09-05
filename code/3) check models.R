library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)

brm_velocity_fines_taxon = readRDS(file = "models/brm_velocity_fines_taxon.rds")
brm_taxon_int = readRDS(file = "models/brm_taxon_int.rds")
brm_taxon_macrophyte = readRDS(file = "models/brm_taxon_macrophyte.rds")
brm_taxon_velocity = readRDS(file = "models/brm_taxon_velocity.rds")
brm_taxon_fines = readRDS(file = "models/brm_taxon_fines.rds")
brm_taxon_mac_fine_vel = readRDS(file = "models/brm_taxon_mac_fine_vel.rds")


# pp checks
pp_check(brm_velocity_fines_taxon, type = "hist") 
pp_check(brm_taxon_int, type = "hist") 
pp_check(brm_taxon_macrophyte, type = "hist") 
pp_check(brm_taxon_velocity, type = "hist") 
pp_check(brm_taxon_fines, type = "hist") 
pp_check(brm_taxon_mac_fine_vel, type = "hist") 

# pp checks
pp_check(brm_velocity_fines_taxon, type = "stat") 
pp_check(brm_taxon_int, type = "stat") 
pp_check(brm_taxon_macrophyte, type = "stat") 
pp_check(brm_taxon_velocity, type = "stat") 
pp_check(brm_taxon_fines, type = "stat") 
pp_check(brm_taxon_mac_fine_vel, type = "stat")