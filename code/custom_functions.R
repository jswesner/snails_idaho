# function to get factorial combinations of predictor values
get_predictor_grid <- function(model, ...) {
  predictors <- quos(...)
  
  model$data %>%
    select(!!!predictors) %>%
    distinct() %>%
    summarise(across(everything(), list(min = min, max = max))) %>%
    pivot_longer(cols = everything()) %>%
    mutate(min_max = str_sub(name, -3, -1),
           name = str_remove(name, "_max"),
           name = str_remove(name, "_min")) %>%
    pivot_wider(names_from = min_max, values_from = value) %>%
    rowwise() %>%
    mutate(sequence = list(seq(min, max, length.out = 50))) %>%
    unnest(sequence) %>%
    select(name, sequence) %>%
    rename(x = sequence) %>%
    rownames_to_column() %>%
    pivot_wider(names_from = "name", values_from = "x") %>%
    replace_na(list(velocity_s = 0, macrophyte_s = 0, fines_s = 0))
}


# function to get snail posteriors for different predictors
get_snail_posts <- function(model = brm_taxon_mac_fine_vel, 
                            model_data = mod_d, 
                            response = "taxon_snail",
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
  response_sym <- sym(response)
  
  pred_grid <- model_data %>%
    select(-!!response_sym) %>% 
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


# function to get snail posteriors for different predictors no taxon
get_snail_posts_no_taxon <- function(model = brm_taxon_mac_fine_vel, 
                            model_data = mod_d, 
                            response = "taxon_snail",
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
  response_sym <- sym(response)
  
  pred_grid <- model_data %>%
    select(-!!response_sym) %>% 
    ungroup() %>%  
    select(-taxon) %>% 
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
                    re_formula =  NA)
}


# get parameters for scaled variables -------------------------------------

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

rm(velocity_mean_sd)
rm(macrophyte_mean_sd)
rm(fines_mean_sd)


