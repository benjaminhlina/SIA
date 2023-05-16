# ---- bring in packages -----
library(dplyr)
library(ggplot2)
library(here)
library(lubridate)
library(readr)
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

at_depth <- at %>% 
  filter(sensor_unit %in% "m") %>% 
  mutate(
    fish_basin = str_replace(fish_basin, " Basin", "") %>% 
      factor(., level = c("East", "West", "North"))
  ) %>% 
  group_by(floy_tag, fish_basin) %>% 
  summarise(
    mean_depth = mean(sensor_value), 
    var_depth = var(sensor_value),
    sd_depth = sd(sensor_value), 
    sem_depth = sd(sensor_value) / sqrt(n())
  ) %>% 
  ungroup()
at_depth


ati <- at_depth %>% 
  left_join(df_short, by = c("floy_tag" = "sample")) 

ati



length(unique(ati$floy_tag))



ggplot(data = ati, aes(x = c_13, y = mean_depth)) + 
  geom_point(size = 3, aes(colour = fish_basin)) + 
  stat_ellipse(aes(colour = fish_basin), 
               type = "norm", linetype = 1,
               linewidth = 1) +
  geom_errorbar(aes(ymin = mean_depth - sem_depth, 
                    ymax = mean_depth + sem_depth), width = 0.05) +
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", 
                         name = "Basin") + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = "Mean deptherature Use (*degree*C*)")
ggplot(data = ati, aes(x = n_15, y = mean_depth)) + 
  geom_point(size = 3, aes(colour = fish_basin)) + 
  stat_ellipse(aes(colour = fish_basin), 
               type = "norm", linetype = 1,
               linewidth = 1) +
  geom_errorbar(aes(ymin = mean_depth - sem_depth, 
                    ymax = mean_depth + sem_depth), width = 0.05) +
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", 
                         name = "Basin") + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = "Mean Depth Use (m)")
ggplot(data = ati, aes(x = c_13, y = n_15)) + 
  geom_point(size = 4, aes(fill = mean_depth, 
                           shape = fish_basin), 
             stroke = 0.7, colour = "black") + 
  scale_shape_manual(name = "Basin", 
                     values = c(21, 22, 23)) + 
  # geom_errorbar(aes(ymin = mean_depth - sem_depth, 
  #                   ymax = mean_depth + sem_depth), width = 0.05) +
  scale_fill_viridis_c(
    begin = 0.2,
    # end = 0.9, 
    option = "D", 
    breaks = rev(seq(5, 20, 5)), 
    direction = -1,
    
    guide = guide_colorbar(reverse = TRUE), 
    
    name = "Depth Use (m)") +
  theme_bw(base_size = 15) + 
  theme(
    # panel.grid = element_line(colour = alpha("black", alpha = 0.35)),
    axis.text = element_text(colour = "black"), 
    # legend.background = element_blank(),
    legend.position = c(0.17, 0.885), 
    legend.box = "horizontal"
  ) +
  labs(
    x = expression(paste(delta ^ 13, "C")),
    y = expression(paste(delta ^ 15, "N")),
  ) -> p 

# p
ggsave(filename = here("Plots", 
                       "temp use and isotopes", 
                       "d13C_vs_d15N_mean_depth.png"), 
       width = 11, height = 8.5, plot = p)
# ggplot(data = ati, aes(z= c_13, y = mean_depth)) + 
#   geom_boxplot(size = 3, aes(colour = fish_basin)) + 
#   geom_errorbar(aes(ymin = mean_depth - sem_depth, 
#                     ymax = mean_depth + sem_depth), width = 0.05) +
#   scale_colour_viridis_d(begin = 0.25, end = 0.85, 
#                          option = "D", 
#                          name = "Basin") + 
#   labs(x = expression(paste(delta ^ 13, "C")), 
#        y = "Mean deptherature Use (*degree*C*)")
# ggplot(data = ati, aes(x = c_13, y = mean_depth)) + 
#   geom_point(size = 3, aes(colour = fish_basin)) + 
#   geom_errorbar(aes(ymin = mean_depth - sd_depth, 
#                     ymax = mean_depth + sd_depth), width = 0.05) +
#   scale_colour_viridis_d(begin = 0.25, end = 0.85, 
#                          option = "D", 
#                          name = "Basin") + 
#   labs(x = expression(paste(delta ^ 13, "C")), 
#        y = "Mean deptherature Use (*degree*C*)")
# facet_wrap(. ~ floy_tag, ncol = 5) 

ati_doy <- ati %>% 
  group_by(floy_tag, fish_basin, ydoy, c_13, n_15) %>% 
  summarise(
    depth_mean = mean(mean_depth), 
    depth_sem = sd(mean_depth) / sqrt(n())
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


ggplot(data = ati_doy, aes(x = ydoy, y = depth_mean)) + 
  geom_point(size = 3, aes(colour = c_13), alpha = 0.5) + 
  # geom_errorbar(aes(ymin = depth_mean - depth_sem, 
  #                   ymax = depth_mean + depth_sem), width = 0.15) + 
  facet_wrap(. ~ fish_basin, ncol = 5) + 
  # scale_x_continuous(breaks = seq(25, 350, 65), 
  #                    label = month_label) +
  theme_bw(base_size = 15) + 
  scale_y_continuous(breaks = seq(0, 17.5, 2.5)) + 
  scale_colour_viridis_c(begin = 0.25, end = 0.85, 
                         option = "D", 
                         name = expression(paste(delta ^ 13, "C"))) + 
  labs(x = "Date", 
       y = "Daily deptherature (°C)") -> p
p



ggsave(filename = here("Plots", 
                       "depth_d13c_basin.png"), 
       plot = p, 
       height = 4.37, width = 8.34 * 3)












ggplot(data = ati_doy, aes(x = ydoy, y = depth_mean)) + 
  geom_point(size = 3, aes(colour = n_15), alpha = 0.5) + 
  # geom_errorbar(aes(ymin = depth_mean - depth_sem, 
  #                   ymax = depth_mean + depth_sem), width = 0.15) + 
  facet_wrap(. ~ fish_basin, ncol = 5) + 
  scale_colour_viridis_c(begin = 0.25, end = 0.85, 
                         option = "D") -> p1

p
p1