# ---- bring in packages -----
library(dplyr)
library(ggplot2)
library(here)
library(lubridate)
library(MixSIAR)
library(nicheROVER)
library(readr)
library(rjags)
library(SIBER)
library(SIBERG)
library(stringr)
library(sf)
library(tidyr)
source(here::here("R", 
                  "Data Cleaning", 
                  "julian_date_reorder.r"
                  
))

# ---- bring in acoustic data -----

at <- read_rds(here("Saved Data", 
                    "kenauk lake trout 2017 - 2021.rds"))

# ---- bring in stable isotope data -----
df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))

# view structures 


glimpse(at)
glimpse(df)


df_short <- df %>% 
  select(sample, c_13, n_15, tl_mm, fl_mm, girth_mm, wt_g)

df_short
glimpse(at)

unique(at$sensor_unit)

at_temp <- at %>% 
  filter(sensor_unit %in% "°C") %>% 
  mutate(
    fish_basin = str_replace(fish_basin, " Basin", ""),
    doy = yday(date)
  ) %>% 
  group_by(floy_tag, fish_basin, date, doy, month, 
           month_number, season, year) %>% 
  summarise(
    mean_temp = mean(sensor_value), 
    sem_temp = sd(sensor_value) / sqrt(n())
  ) %>% 
  ungroup()
at_temp


ati <- at_temp %>% 
  left_join(df_short, by = c("floy_tag" = "sample")) 



ati <- ati %>% 
  mutate(
    ydoy = days(date)
  )



length(unique(ati$floy_tag))



ggplot(data = ati, aes(x = date, y = mean_temp)) + 
  geom_point(size = 3, aes(colour = c_13)) + 
  geom_errorbar(aes(ymin = mean_temp - sem_temp, 
                    ymax = mean_temp + sem_temp), width = 0.15) + 
  facet_wrap(. ~ floy_tag, ncol = 5) 

ati_doy <- ati %>% 
  group_by(floy_tag, fish_basin, ydoy, c_13, n_15) %>% 
  summarise(
    temp_mean = mean(mean_temp), 
    temp_sem = sd(mean_temp) / sqrt(n())
  ) %>% 
  ungroup()


ati_doy <- ati_doy %>% 
  mutate(
    date = as.Date(ydoy, origin = "2021-04-30"), 
    month_abb = month(date, label = TRUE)
  )


# ati_doy %>% 
#   filter(ydoy %in% seq(25, 350, 65)) %>% 
#   group_by(month_abb) %>% 
#   summarise() %>% 
#   ungroup() %>% 
#   mutate(
#     month_abb = forcats::fct_relevel(month_abb, "May", "Jul", "Oct",
#                                      "Dec", "Feb", "Apr")
#   ) %>% 
# 
#   .$month_abb -> month_label 
# month_label


ggplot(data = ati_doy, aes(x = ydoy, y = temp_mean)) + 
  geom_point(size = 3, aes(colour = c_13), alpha = 0.5) + 
  # geom_errorbar(aes(ymin = temp_mean - temp_sem, 
  #                   ymax = temp_mean + temp_sem), width = 0.15) + 
  facet_wrap(. ~ fish_basin, ncol = 5) + 
  # scale_x_continuous(breaks = seq(25, 350, 65), 
  #                    label = month_label) +
  theme_bw(base_size = 15) + 
  scale_y_continuous(breaks = seq(0, 17.5, 2.5)) + 
  scale_colour_viridis_c(begin = 0.25, end = 0.85, 
                         option = "D", 
                         name = expression(paste(delta ^ 13, "C"))) + 
  labs(x = "Date", 
       y = "Daily Temperature (°C)") -> p
p



ggsave(filename = here("Plots", 
                       "temp_d13c_basin.png"), 
       plot = p, 
       height = 4.37, width = 8.34 * 3)












ggplot(data = ati_doy, aes(x = ydoy, y = temp_mean)) + 
  geom_point(size = 3, aes(colour = n_15), alpha = 0.5) + 
  # geom_errorbar(aes(ymin = temp_mean - temp_sem, 
  #                   ymax = temp_mean + temp_sem), width = 0.15) + 
  facet_wrap(. ~ fish_basin, ncol = 5) + 
  scale_colour_viridis_c(begin = 0.25, end = 0.85, 
                         option = "D") -> p1

p
p1