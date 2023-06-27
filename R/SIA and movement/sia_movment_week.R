# ---- bring in packages -----
library(broom.mixed)
library(car)
library(dplyr)
library(DHARMa)
library(emmeans)
library(fitdistrplus)
library(ggplot2)
library(glmmTMB)
library(GLMMadaptive)
library(gratia)
library(here)
library(lme4)
library(mgcv)
library(purrr)
library(readr)
library(rstanarm)
library(sf)
library(tidyr)


# ---- bring in dataframe ------

df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))


edges_sf <- read_rds(here("Saved Data", 
                          "lkt_edges_sf_weeks.rds"))

glimpse(edges_sf)
# ---- Bring in shapefile ------
p_lake <- st_read(dsn = here::here("Shapefiles",
                                   "."),
                  layer = "plake_edit_wo_link") %>% 
  
  dplyr::select(AREA_GEO, geometry)
p_lake_utm <- p_lake %>% 
  st_transform(., crs = 32618)


edges_sf %>% 
  group_by(from_to, cost_dist, geometry) %>% 
  summarise() %>% 
  ggplot() + 
  geom_sf(data = p_lake_utm, fill = NA) + 
  geom_sf(aes(colour = cost_dist), linewidth = 1) + 
  scale_colour_viridis_c(option = "B", end = 0.75) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(
    x = "Latitude", 
    y = "Longitude"
  )

glimpse(df)
df <- df %>% 
  filter(c_13 != is.na(c_13))

edges_sf <- edges_sf %>% 
  filter(from != to)

glimpse(edges_sf)

summary(edges_sf$weight)

# ---- calculate distance for each month for each year ---- 

# first caclcute the distance traveled between places for each 
# momth for each year 
dis_m <- edges_sf %>% 
  st_drop_geometry() %>% 
  dplyr::select(floy_tag, fish_basin, weeks, month, from_nam, to_nam, 
                season, year, 
                weight, cost_dist) %>% 
  mutate(
    dis = weight * cost_dist, 
    dis_km = dis / 1000, 
    dis_cubed = dis_km ^ (1 / 3)
  )
# determine the total distance traveled per month per year sum up all movement
dis_sum  <- dis_m %>% 
  group_by(floy_tag, fish_basin, weeks, month, year) %>% 
  summarise(
    dis_sum_km = sum(dis_km),
    dis_sum = sum(dis_cubed)
  ) %>% 
  ungroup()

dis_sum

# summarized this across each month so we have the average for each month 
dis_m_sum <- dis_sum %>% 
  group_by(floy_tag, fish_basin, weeks) %>% 
  summarise(
    dis_km = mean(dis_sum),
    dis_sd = sd(dis_sum),
    dis_sem = sd(dis_sum) / sqrt(n()),
    dis_var = var(dis_sum)
  ) %>% 
  ungroup()


dis_m_sum 
# look at distrubtion of monthly mean movement fore each month 
ggplot(data = dis_m_sum, aes(x = dis_km)) + 
  geom_histogram() + 
  facet_wrap(~ weeks) + 
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  labs(
    x = expression(paste("Movement (", sqrt("km / month", "" ^ 3), 
                         ")")),
    y = "Frequency")

# caclaute the monthly mean movment sd and sem 
dis_mean_overall <- dis_m_sum %>% 
  group_by(floy_tag) %>% 
  summarise(
    mean_dis = mean(dis_km),
    sd_dis = sd(dis_km),
    sem_dis = sd(dis_km) / sqrt(n()), 
    var_dis = var(dis_km)
  ) %>% 
  ungroup()


# average movement regardless of time period 
# dis_mean_overall <- dis_m %>% 
#   group_by(floy_tag) %>% 
#   summarise(
#     mean_dis = mean(dis_cubed),
#     sd_dis = sd(dis_cubed),
#     sem_dis = sd(dis_cubed) / sqrt(n()), 
#     var_dis = var(dis_cubed)
#   ) %>% 
#   ungroup()


# look at distribution of monthly mean movement fore each month 
ggplot(data = dis_mean_overall, aes(x = mean_dis)) + 
  geom_histogram() +  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  labs(
    x = expression(paste("Average Monthly Movement (", sqrt("km", "" ^ 3),")")),
    y = "Frequency")

# ---- join this data to stable isotope data -----


df_movment_overall <- df %>% 
  left_join(dis_mean_overall, by = c("sample" = "floy_tag"))


df_movment_overall <- df_movment_overall %>% 
  filter(mean_dis != is.na(mean_dis) & 
           c_13 != is.na(c_13))

dis_mean_overall_s <- dis_mean_overall %>%
  filter(floy_tag != "07478")

