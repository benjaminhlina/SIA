# ---- load packages ---- 
{
  library(dplyr)
  library(here)
  library(openxlsx)
  library(purrr)
  library(readr)
  library(stringr)
}

# ---- bring in model fits and effects ----

# distance models 
dis_fit <- read_rds(here("Results", 
                         "distance_isotope_model_fit.rds"))
dis_effect <- read_rds(here("Results", 
                            "distance_isotope_model_effects.rds"))

# temp models 

temp_fit <- read_rds(here("Results", 
                          "temp_isotope_model_fit.rds"))
temp_effect <- read_rds(here("Results", 
                             "temp_isotope_model_effects.rds"))

# depth models 

depth_fit <- read_rds(here("Results", 
                           "depth_isotope_model_fit.rds"))
depth_effect <- read_rds(here("Results", 
                              "depth_isotope_model_effects.rds"))


# ---- view all fits and merge them together ----

glimpse(dis_fit)
glimpse(temp_fit)
glimpse(depth_fit)


sia_metric_fit <- bind_rows(dis_fit, temp_fit, depth_fit)

sia_metric_fit



# ---- view all effects and merge them together -----

glimpse(dis_effect)
glimpse(temp_effect)
glimpse(depth_effect)


dis_effect

temp_effect <- temp_effect %>% 
  mutate(
    df = NA, 
    statistic = NA, 
    chi_sq = NA
  ) %>% 
  rename(
    terms = term, 
    p_value = p.value
  ) %>% 
  select(id, metric, terms, df, statistic, chi_sq, p_value)

depth_effect <- depth_effect %>% 
  mutate(
    df = NA, 
    statistic = NA, 
    chi_sq = NA
  ) %>% 
  rename(
    terms = term, 
    p_value = p.value
  ) %>% 
  select(id, metric, terms, df, statistic, chi_sq, p_value)


sia_metric_effects <- bind_rows(dis_effect, temp_effect, depth_effect)



# ---- write excel files ----


write.xlsx(sia_metric_fit, here("Results", 
                                "sia_metric_model_selection.xlsx"))


write.xlsx(sia_metric_effects, here("Results", 
                                    "sia_metric_effect_selection.xlsx"))

