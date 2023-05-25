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
  dplyr::select(floy_tag, fish_basin, month, from_nam, to_nam, 
                season, year, 
                weight, cost_dist) %>% 
  mutate(
    dis = weight * cost_dist, 
    dis_km = dis / 1000, 
    dis_cubed = dis_km ^ (1 / 3)
  )
# determine the total distance traveled per month per year sum up all movement
dis_sum  <- dis_m %>% 
  group_by(floy_tag, fish_basin, month, year) %>% 
  summarise(
    dis_sum_km = sum(dis_km),
    dis_sum = sum(dis_cubed)
  ) %>% 
  ungroup()

dis_sum

# summarized this across each month so we have the average for each month 
dis_m_sum <- dis_sum %>% 
  group_by(floy_tag, fish_basin, month) %>% 
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
  facet_wrap(~ month) + 
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

# --- join this data to stable isotope data -----


df_movment_overall <- df %>% 
  left_join(dis_mean_overall, by = c("sample" = "floy_tag"))


df_movment_overall <- df_movment_overall %>% 
  filter(mean_dis != is.na(mean_dis) & 
           c_13 != is.na(c_13))

dis_mean_overall_s <- dis_mean_overall %>%
  filter(floy_tag != "07478")
## $\delta^{15}$N *versus* Movement ± sem with colour as basin for isothermal and stratification

# ---- Plots ----

ggplot(data = df_movment_overall, 
       aes(x = mean_dis, y = n_15)) + 
  geom_point(size = 4, aes(fill = basin), 
             colour = "black", 
             shape = 21, stroke = 0.8) +
  # geom_errorbar(aes(xmin = mean_dis - sem_dis, 
  #                   xmax = mean_dis + sem_dis), width = 0.05) +
  # scale_shape_manual(values = c(21:23), 
  #                    name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin",
                       alpha = 0.75
  ) +
  scale_y_continuous(breaks = seq(9, 15, 1)) +
  # facet_wrap(. ~ therm) + 
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.92, 0.85), 
        # legend.title.align = 0.5,
        legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(y = expression(paste(delta ^ 15, "N")), 
       x = expression(paste("Mean Monthly Movement (", sqrt("km", "" ^ 3),")"))
  ) -> p1
p1
ggsave(filename = here("Plots",
                       "movement and isotopes",
                       "n_15_distance_traveled_overall.png"),
       width = 11, height = 8.5, plot = p1)

write_rds(p1, here("Saved Plots",
                   "n_15_distance_traveled_overall.rds"))

## $\delta^{13}$C *versus* Movement ± sem with colour as basin for isothermal and stratification


ggplot(data = df_movment_overall, 
       aes(y = c_13, x = mean_dis)) + 
  geom_point(size = 4, aes(fill = basin), 
             colour = "black", 
             shape = 21, stroke = 0.8) +
  # geom_errorbar(aes(xmin = mean_dis - sem_dis, 
  #                   xmax = mean_dis + sem_dis), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin",
                       alpha = 0.75
  ) +
  scale_y_continuous(breaks = rev(seq(-23, -31, -1))) + 
  # scale_y_continuous(breaks = seq(0, 10, 2), 
  #                    limits = c(0, 9)) +
  theme_bw(base_size = 15) +
  # facet_wrap(.~ therm) + 
  #   geom_text(
  #     aes(x = (-29 + -2), y = 8.75), label = "Pelagic",
  #     size = 9, vjust = 0, hjust = 0, check_overlap = TRUE) +
  #   geom_text(
  #     aes(x = (-27.3 + 1), y = 8.75), label = "Littoral",
  #     size = 9, vjust = 0, hjust = 0, check_overlap = TRUE) +
  #   geom_segment(aes(x = -27.28392, y = 8.65, xend = -24.28392, yend = 8.65),
  #                arrow = arrow(length = unit(0.5, "cm"))) +
  # geom_segment(aes(x = -29, y = 8.65, xend = -32, yend = 8.65),
  #                arrow = arrow(length = unit(0.5, "cm"))) +
#   
theme(axis.text = element_text(colour = "black"),
      panel.grid = element_blank(), 
      legend.position = c(0.92, 0.85), 
      # legend.title.align = 0.5,
      legend.box = "horizontal",
      legend.background = element_blank()) + 
  labs(y = expression(paste(delta ^ 13, "C")), 
       x = expression(paste("Mean Monthly Movement (", sqrt("km", "" ^ 3), 
                            ")"))) -> p2

p2

ggsave(filename = here("Plots",
                       "movement and isotopes",
                       "c_13_distance_traveled_overall.png"),
       width = 11, height = 8.5, plot = p2)



write_rds(p2, here("Saved Plots",
                   "c_13_distance_traveled_overall.rds"))
## $\delta^{15}$N *versus* Movement ± sem with colour as basin

ggplot(data = df_movment_overall, 
       aes(x = c_13, y = n_15)) + 
  geom_point(size = 4, aes(fill = mean_dis, 
                           shape = basin), 
             colour = "black", stroke = 0.8) + 
  scale_shape_manual(values = c(21:23),
                     name = "Basin") +
  scale_fill_viridis_c(begin = 0.25, end = 0.95, 
                       option = "D",
                       alpha = 0.75, 
                       name = expression(paste("Movement (", sqrt("km", "" ^ 3), 
                                               ")"))) +
  # scale_y_continuous(breaks = seq(0, 8, 2), 
  #                    limits = c(0, 8.5)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.92, 0.80), 
        # legend.title.align = 0.5,
        # legend.box = "horizontal",
        legend.background = element_blank()) + 
  labs(
    x = expression(paste(delta ^ 13, "C")), 
    y = expression(paste(delta ^ 15, "N"))) -> p3
p3


ggsave(filename = here("Plots",
                       "movement and isotopes",
                       "c_13_n_15_distance_traveled.png"),
       width = 11, height = 8.5, plot = p3)

write_rds(p2, here("Saved Plots",
                   "c_13_n_15_distance_traveled.rds"))
# ---- ANALYSIS -----
# ---- look at distributions ----


glimpse(df_movment_overall)
descdist(df_movment_overall$c_13)
descdist(df_movment_overall$n_15)

ggplot(data = df_movment_overall, aes(x = c_13)) +
  geom_histogram()

ggplot(data = df_movment_overall, aes(x = n_15)) +
  geom_histogram()


glimpse(df_movment_overall)

# ---- create model for c_13 ~ mean_dis -----
options(contrasts = c("contr.sum", "contr.poly"))

# m <- gam(c_13 ~ mean_dis * basin,
#          family = scat,
#          method = "REML",
#          data = df_movment_overall)

m <- lm(c_13 ~ mean_dis * basin,
        contrasts = list(basin = "contr.sum"),
         # family = scat, 
         # method = "REML",
         data = df_movment_overall)


# appraise(m)
# anova.gam(m)
# ---- create specific stuff for model saving -----
hist(residuals(m))
res <- simulateResiduals(m)
plot(res)
# appraise(m)
par(mfrow = c(2, 2))
# gam.check(m)
plot(m)
par(mfrow = c(1, 1))

me_d13c <- Anova(m, type = "III")
# m_aov <- anova.gam(m)
# 
# 
# me_d13c <- m_aov$pTerms.table %>% 
#   as_tibble(rownames = "terms") %>% 
#   janitor::clean_names()
me_d13c

# ---- model section ----  
m1 <- update(m, ~ mean_dis)
m2 <- update(m, ~ basin)
m3 <- update(m, ~ mean_dis + basin)

# ---- create model list for model selection ------
model_list <- list(m, m1, m2)
# give the elements useful names
names(model_list) <- c("m", 
                       "m1", "m2")

glance(m)
# get the summaries using `lapply

summary_list <- lapply(model_list, function(x) tidy(x, parametric = TRUE))
glance_list <- lapply(model_list, glance)

glance_summary <- map_df(glance_list, ~as.data.frame(.x), .id = "id") %>% 
  mutate(model = lapply(model_list, formula) %>%
           as.character() 
  ) %>% 
  dplyr::select(model, id:df.residual) %>% 
  arrange(AIC)


glance_summary

glance_summary <- glance_summary %>% 
  mutate(
         delta_AIC = AIC - first(AIC), 
         AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)
glance_summary

glance_summary %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_movement_model_selection.xlsx"))

# ---- create specific stuff for model saving -----
car::Anova(m)
summary(m)

main_effects <- tidy(car::Anova(m))



ind_effects <- tidy(m)


# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_movment_lmer_main_effect.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_movment_lmer_ind_effect.xlsx"))

# ---- Nitrogen vs movement -----

m3 <- glmmTMB(n_15 ~ mean_dis * basin + (1|sample),
              data = df_movment_overall,
              family = Gamma(link = "log"),
              control = glmmTMBControl(optimizer=optim,
                                     optArgs = list(method = "BFGS")),
              REML = TRUE)

# model fit evaluation 
res3 <- simulateResiduals(m3)
plot(res3)
Anova(m3)
# model section 
m4 <- update(m3, ~ mean_dis + (1|sample))
m5 <- update(m3, ~ mean_dis)


model_list <- list(m3, m4, m5)
# give the elements useful names
names(model_list) <- c("m3", 
                       "m4", "m5")

glance(m)
# get the summaries using `lapply

summary_list <- lapply(model_list, function(x) tidy(x, parametric = TRUE))
glance_list <- lapply(model_list, glance)

glance_summary <- map_df(glance_list, ~as.data.frame(.x), .id = "id") %>% 
  mutate(model = lapply(model_list, formula) %>%
           as.character() 
  ) %>% 
  dplyr::select(model, id:df.residual) %>% 
  arrange(AIC)


glance_summary

glance_summary <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)
glance_summary

glance_summary %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_movement_model_selection.xlsx"))

# ---- create specific stuff for model saving -----
car::Anova(m3)
# summary(m3)

main_effects <- tidy(car::Anova(m3))



ind_effects <- tidy(m3)


# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_movment_lmer_main_effect.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_movment_lmer_ind_effect.xlsx"))

