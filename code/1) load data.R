library(tidyverse)
library(janitor)
library(lubridate)
library(readxl)

# load data
macrophytes = read_csv("data/macrophytes_snails.csv") %>% 
  pivot_longer(cols = starts_with("site_"), names_to = "site", values_to = "perc_macrophyte") %>% 
  mutate(site = parse_number(site))

lat_long_elev <- read_excel("data/lat_long_elev.xlsx")

snail_density = read_excel("data/2024.Birch.snail.density.final.xlsx", sheet = "Data") %>% 
  mutate(month = month(date),
         season = case_when(month <= 6 ~ "spring",
                             TRUE ~ "fall"),
         season = as.factor(season),
         season = fct_relevel(season, "spring")) %>% 
  left_join(lat_long_elev) %>% 
  left_join(macrophytes) %>% 
  mutate(site_f = as.factor(site),
         transect_f = as.factor(transect),
         quadrat_f = as.factor(quadrat),
         zero_one = case_when(total_snail == 0 ~ 0, TRUE ~ 1),
         velocity_s = scale(velocity),
         depth_s = scale(depth),
         fines_s = scale(fines),
         elevation_s = scale(elevation),
         distance_s = scale(distance),
         macrophyte_s = scale(perc_macrophyte),
         sample_area_m2 = 0.5*0.5) 

saveRDS(snail_density, file = "data/snail_density.rds")

snail_definitions = read_excel("data/2024.Birch.snail.density.final.xlsx", sheet = "Data description")

# exploratory plots

snail_density %>% 
  ggplot(aes(x = site, y = velocity, color = season)) + 
  geom_point() +
  facet_wrap(~season) +
  # scale_y_log10() +
  NULL

snail_density %>% 
  ggplot(aes(x = total_snail, color = season)) + 
  geom_histogram(bins = 50) +
  facet_wrap(~season) 


snail_density %>% 
  ggplot(aes(x = elevation_s, y = total_snail, color = season)) + 
  geom_jitter() +
  facet_wrap(~season) +
  scale_y_log10()

snail_density %>% 
  ggplot(aes(x = macrophyte_s, y = total_snail, color = season)) + 
  geom_jitter() +
  facet_wrap(~season) +
  scale_y_log10()

plot(elevation ~ velocity, data = snail_density)

snail_density %>% 
  pivot_longer(cols = c(smallgrav, largegrav, cobble, boulder, bedrock)) %>% 
  ggplot(aes(x = fines, y = value)) + 
  geom_jitter() + 
  facet_wrap(~name, scales = "free")
