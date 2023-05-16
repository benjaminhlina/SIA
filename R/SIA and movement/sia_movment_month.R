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


dis_mean_month <- dis_my %>% 
  group_by(floy_tag, fish_basin, month) %>% 
  summarise(
    dis_mean_km = mean(dis_km),
    dis_mean_km_sd = sd(dis_km),
    dis_mean_km_sem = sd(dis_km) / sqrt(n()),
    dis_mean_km_var = var(dis_km),
    dis_mean = mean(dis_cubed),
    dis_sd = sd(dis_cubed),
    dis_sem = sd(dis_cubed) / sqrt(n()), 
    dis_var = var(dis_cubed)
  ) %>% 
  ungroup()

dis_mean_month

ggplot(data = dis_mean_month, aes(x = month, y = dis_mean)) + 
  geom_point(aes(fill = fish_basin), shape = 21, colour = "black", 
             stroke = 0.8) + 
  scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                       option = "D", 
                       # name = expression(
                       #    scale_fill_viridis_c(begin = 0.25, end = 0.95, 
                       # option = "D", 
                       # name = expression(
                       #   paste("Movement (", sqrt("km / month", "" ^ 3), 
                       #              ")")),  
                       # # alpha = 0.65
                       # ) + 
                       # alpha = 0.65
  ) 


ggplot(data = dis_mean_month, aes(x = dis_mean)) + 
  geom_histogram() + 
  facet_wrap(~ month)

dis_mean <- dis_mean_month %>% 
  filter(month %in% c("November", "December", 
                      "January", "February", 
                      "March", "April")) %>% 
  group_by(floy_tag) %>% 
  summarise(
    dis_mean_km_o = mean(dis_mean_km),
    dis_mean_km_sd_o = sd(dis_mean_km),
    dis_mean_km_sem_o = sd(dis_mean_km) / sqrt(n()),
    dis_mean_km_var_o = var(dis_mean_km),
    dis_mean_o = mean(dis_mean),
    dis_sd_o = sd(dis_mean),
    dis_sem_o = sd(dis_mean) / sqrt(n()), 
    dis_var_o = var(dis_mean)
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


glimpse(df_movment)
df_movment_heard <- df_movment %>% 
  filter(dis_mean_o != is.na(dis_mean_o))


# ------------------ ANALYSIS -------------------------------------------
df_movment_heard <- df_movment_heard %>% 
  filter(c_13 != is.na(c_13))
glimpse(df_movment_heard)
# ---- look at distrubtions ----
descdist(df_movment_heard$dis_mean_o)
descdist(df_movment_heard$c_13)
descdist(df_movment_heard$n_15)
ggplot(data = df_movment_heard, aes(x = dis_mean_o)) +
  geom_histogram()
ggplot(data = df_movment_heard, aes(x = c_13)) +
  geom_histogram()
ggplot(data = df_movment_heard, aes(x = n_15)) +
  geom_histogram()

# --- create model for mean_dis ~ c_13 ----- 
m <- glmmTMB(dis_mean_o ~ c_13 * basin + (1|sample),
             data = df_movment_heard, 
             family = gaussian(link = "identity"),
             REML = TRUE)


res <- simulateResiduals(m)
plot(res)
Anova(mod = m)

m1 <- glmmTMB(c_13 ~ dis_mean_o * basin + (1|sample),
              data = df_movment_heard, 
              family = gaussian(link = "identity"),
              REML = TRUE)


res1 <- simulateResiduals(m1)
plot(res1)

Anova(m1)

# --- create model for mean_dis ~ n_13 ----- 
m2 <- glmmTMB(dis_mean_o ~ n_15 * basin + (1|sample),
              data = df_movment_heard, 
              family = gaussian(link = "identity"),
              REML = TRUE)



res2 <- simulateResiduals(m2)
plot(res2)

Anova(m2, type = "II")

m3 <- glmmTMB(log(n_15) ~ dis_mean_o * basin + (1|sample),
              data = df_movment_heard, 
              family = gaussian(link = "identity"),
              REML = TRUE)


res3 <- simulateResiduals(m3)
plot(res3)

Anova(m3)




# 
# m1 <- glmmTMB(mean_dis ~ c_13 * fish_basin + (1|sample),
#              data = df_movment_heard, family = gaussian(link = "log"))
# res1 <- simulateResiduals(m1)
# plot(res1)
# 
# Anova(m1)

# ---- isotopes by basin with basic ellipses plot -----
ggplot(data = df_movment_heard, 
       aes(x = c_13, y = n_15)) + 
  geom_point(size = 5, aes(shape = basin, 
                           fill = dis_mean), 
             colour = "black") +
  scale_shape_manual(values = c(21:23), 
                     name = "Basin") + 
  scale_fill_viridis_c(begin = 0.25, end = 0.95, 
                       option = "D", 
                       name = expression(
                         paste("Movement (", sqrt("km / month", "" ^ 3), 
                                    ")")),  
                       # alpha = 0.65
                       ) + 
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.22, 0.89), 
        # legend.title.align = 0.5,
        legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p 
p


ggsave(filename = here("Plots",
                       "movement and isotopes",
                       "c_13_n_15_distance_traveled.png"),
       width = 11, height = 8.5, plot = p)


ggplot(data = df_movment_heard, 
       aes(x = n_15, y = dis_mean)) + 
  geom_point(size = 5, aes(fill = basin), 
             colour = "black", 
             shape = 21, stroke = 0.8) +
  geom_errorbar(aes(ymin = dis_mean - dis_sem, 
                    ymax = dis_mean + dis_sem), width = 0.05) +
  # scale_shape_manual(values = c(21:23), 
  #                    name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin",
                       alpha = 0.75
                       ) +
  scale_y_continuous(breaks = seq(0, 8, 2), 
                     limits = c(0, 8.5)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.92, 0.92), 
        # legend.title.align = 0.5,
        legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 15, "N")), 
       y = expression(paste("Movement (", sqrt("km / month", "" ^ 3), 
               ")"))) -> p1
p1

ggsave(filename = here("Plots", 
                       "movement and isotopes", 
                       "d15N_distance_traveled.png"), 
       width = 11, height = 8.5, plot = p1)


ggplot(data = df_movment_heard, 
       aes(x = c_13, y = dis_mean_o)) + 
  geom_point(size = 5, aes(fill = basin), 
             colour = "black", 
             shape = 21, stroke = 0.8) +
  geom_errorbar(aes(ymin = dis_mean_o - dis_sem_o, 
                    ymax = dis_mean_o + dis_sem_o), width = 0.05) +
  # scale_shape_manual(values = c(21:23), 
  #                    name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin",
                       alpha = 0.75
                       ) +
  scale_y_continuous(breaks = seq(0, 10, 2), 
                     limits = c(0, 9)) +
  theme_bw(base_size = 15) +
  
  geom_text(
    aes(x = (-29 + -2), y = 8.75), label = "Pelagic",
    size = 9, vjust = 0, hjust = 0, check_overlap = TRUE) +
  geom_text(
    aes(x = (-27.3 + 1), y = 8.75), label = "Littoral",
    size = 9, vjust = 0, hjust = 0, check_overlap = TRUE) +
  geom_segment(aes(x = -27.28392, y = 8.65, xend = -24.28392, yend = 8.65),
               arrow = arrow(length = unit(0.5, "cm"))) +
geom_segment(aes(x = -29, y = 8.65, xend = -32, yend = 8.65),
               arrow = arrow(length = unit(0.5, "cm"))) +
  
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.92, 0.92), 
        # legend.title.align = 0.5,
        legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste("Movement (", sqrt("km / month", "" ^ 3), 
               ")"))) -> p2
p2


ggsave(filename = here("Plots", 
                       "movement and isotopes", 
                       "d13C_distance_traveled.png"), 
       width = 11, height = 8.5, plot = p2)

ggplot(data = df_movment_heard, 
       aes(x = n_15, y = dis_mean_km)) + 
  geom_point(size = 5, aes(fill = basin), 
             colour = "black", 
             shape = 21, stroke = 0.8) +
  geom_errorbar(aes(ymin = dis_mean_km - dis_mean_km_sem, 
                    ymax = dis_mean_km + dis_mean_km_sem), width = 0.05) +
  # scale_shape_manual(values = c(21:23), 
  #                    name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin",
                       alpha = 0.75
  ) +
  # scale_y_continuous(breaks = seq(0, 175, 25)) + 
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.92, 0.92), 
        # legend.title.align = 0.5,
        legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste("Movement (km  ", month ^ -1, ")"))) -> p3


p3
glimpse(df_movment_heard)

ggplot(data = df_movment_heard, 
       aes(x = c_13, y = dis_mean_km)) + 
  geom_point(size = 5, aes(fill = basin), 
             colour = "black", 
             shape = 21, stroke = 0.8) +
  geom_errorbar(aes(ymin = dis_mean_km - dis_mean_km_sem, 
                    ymax = dis_mean_km + dis_mean_km_sem), width = 0.05) +
  # scale_shape_manual(values = c(21:23), 
  #                    name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin",
                       alpha = 0.75
  ) +
  # scale_y_continuous(breaks = seq(0, 175, 25)) + 
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.92, 0.92), 
        # legend.title.align = 0.5,
        legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste("Movement (km  ", month ^ -1, ")"))) -> p4

p1
p2
p3
p4