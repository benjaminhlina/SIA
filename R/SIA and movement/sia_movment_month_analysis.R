# ---- bring in packages -----
library(broom)
library(car)
library(dplyr)
library(DHARMa)
library(emmeans)
library(fitdistrplus)
library(ggplot2)
library(glmmTMB)
library(here)
library(lme4)
library(readr)
library(sf)

# ---- bring in dataframe ------

df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))


edges_sf <- read_rds(here("Saved Data", 
                          "lkt_edges_sf.rds"))

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
  theme(panel.grid = element_blank())



glimpse(df)
glimpse(edges_sf)


df <- df %>% 
  filter(c_13 != is.na(c_13))

glimpse(edges_sf)

edges_sf <- edges_sf %>% 
  filter(from != to)

glimpse(edges_sf)

# ---- calculate distance for each month for each year ---- 

dis_my <- edges_sf %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  dplyr::select(floy_tag, fish_basin, month, from_nam, to_nam, 
                season, year, 
                weight, cost_dist) %>% 
  mutate(
    dis = weight * cost_dist, 
    dis_km = dis / 1000, 
    dis_cubed = dis_km ^ (1 / 3)
  )



dis_my 


dis_my %>% 
  group_by(floy_tag, month) %>% 
  summarise(
    dis_mean = mean(dis_km), 
    dis_cubed = dis_mean ^ (1 / 3)
  ) %>% 
  ungroup()


# dis_m <- dis_my %>% 
#   group_by(
#     floy_tag, fish_basin, month
#   ) %>% 
#   summarise(
#     sum_dis = sum(dis_cubed),
#     sum_dis_raw = sum(dis_km), 
#   ) %>% 
#   ungroup()
# 
# dis_m
# 
# dis_mean <- dis_m %>% 
#   group_by(floy_tag, fish_basin) %>% 
#   summarise(
#     mean_dis_raw = mean(sum_dis_raw), 
#     sem_dis_raw = sd(sum_dis_raw) / sqrt(n()), 
#     mean_dis = mean(sum_dis), 
#     sem_dis = sd(sum_dis) / sqrt(n())
#   ) %>% 
#   ungroup()
# 
# dis_mean

df_movment <- df %>% 
  left_join(dis_mean, by = c("sample" = "floy_tag"))




df_movment

glimpse(df_movment)
df_movment_heard <- df_movment %>% 
  filter(mean_dis != is.na(mean_dis))

glimpse(df_movment_heard)

# we will use mean_dis 
# ---- look at distrubtions ----
descdist(df_movment_heard$mean_dis)
descdist(df_movment_heard$c_13)
descdist(df_movment_heard$n_15)
ggplot(data = df_movment_heard, aes(x = mean_dis)) +
  geom_histogram()
ggplot(data = df_movment_heard, aes(x = c_13)) +
  geom_histogram()
ggplot(data = df_movment_heard, aes(x = n_15)) +
  geom_histogram()

# --- create model for mean_dis ~ c_13 ----- 
m <- glmmTMB(mean_dis ~ c_13 * fish_basin + (1|sample),
             data = df_movment_heard, family = Gamma(link = "log"),
             REML = TRUE)


res <- simulateResiduals(m)
plot(res)

Anova(m, type = "III")
m1 <- glmmTMB(c_13 ~ mean_dis * fish_basin + (1|sample),
             data = df_movment_heard, family = gaussian(link = "identity"),
             REML = TRUE)


res1 <- simulateResiduals(m1)
plot(res1)

Anova(m1)

# --- create model for mean_dis ~ n_13 ----- 
m2 <- glmmTMB(mean_dis ~ n_15 * fish_basin + (1|sample),
             data = df_movment_heard, family = Gamma(link = "log"),
             REML = TRUE)


res2 <- simulateResiduals(m2)
plot(res2)

Anova(m2, type = "II")

m3 <- glmmTMB(n_15 ~ mean_dis * fish_basin + (1|sample),
             data = df_movment_heard, family = Gamma(link = "log"),
             REML = TRUE)


res3 <- simulateResiduals(m3)
plot(res3)

Anova(m1)




# 
# m1 <- glmmTMB(mean_dis ~ c_13 * fish_basin + (1|sample),
#              data = df_movment_heard, family = gaussian(link = "log"))
# res1 <- simulateResiduals(m1)
# plot(res1)
# 
# Anova(m1)
