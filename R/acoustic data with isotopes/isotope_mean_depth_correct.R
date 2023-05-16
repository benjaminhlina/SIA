# ---- bring in packages -----
library(broom.mixed)
library(car)
library(dplyr)
library(DHARMa)
library(fitdistrplus)
library(ggplot2)
library(glmmTMB)
library(here)
library(lubridate)
library(openxlsx)
library(purrr)
library(readr)
library(stringr)
library(sf)
library(tidyr)

# ---- bring in acoustic data -----

at <- read_rds(here("Saved Data", 
                    "kenauk lake trout 2017 - 2021.rds"))

# ---- bring in stable isotope data -----
df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))

# view structures 


glimpse(at)
glimpse(df)

# ---- grab only the columns you need -----
df_short <- df %>% 
  dplyr::select(sample, c_13, n_15, tl_mm, fl_mm, girth_mm, wt_g)

df_short
glimpse(at)

unique(at$sensor_unit)

# ---- mean depth for each month ==== 
at_my <- at %>% 
  filter(sensor_unit == "m") %>% 
  mutate(
    fish_basin = str_replace(fish_basin, " Basin", "") %>% 
      factor(., level = c("East", "West", "North"))
  ) %>% 
  group_by(floy_tag, fish_basin, 
           month, year 
  ) %>% 
  summarise(
    depth = mean(sensor_value), 
    var = var(sensor_value),
    sd = sd(sensor_value), 
    sem = sd(sensor_value) / sqrt(n())
  ) %>% 
  ungroup()

at_depth <- at_my %>% 
  group_by(floy_tag, fish_basin, 
           month
  ) %>% 
  summarise(
    md = mean(depth), 
    var = var(depth),
    sd = sd(depth), 
    sem = sd(depth) / sqrt(n())
  ) %>% 
  ungroup()
at_depth

# ---- mean depth across months ----- 
at_depth_sum <- at_depth %>% 
  group_by(floy_tag, fish_basin) %>% 
  summarise(
    mean_depth = mean(md), 
    sd = sd(md), 
    var_depth = var(md), 
    sem = sd(md) / sqrt(n())
  ) %>% 
  ungroup()
# print(n = 100, at_depth)

# ---- join mean depth across months to isotope data ----
ati <- at_depth_sum %>% 
  left_join(df_short, by = c("floy_tag" = "sample")) 


# ---- plot mean depth across months for d13C -----
ggplot(data = ati, aes(x = c_13, y = mean_depth)) + 
  geom_point(size = 4, aes(fill = fish_basin), 
             shape = 21, stroke = 0.8) + 
  # stat_ellipse(aes(colour = fish_basin), 
  #              type = "norm", linetype = 1,
  # linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1) +  
  geom_errorbar(aes(ymin = mean_depth - sem, 
                    ymax = mean_depth + sem), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin") + 
  scale_x_continuous(breaks = rev(seq(-26, -31, -1))) + 
  scale_y_reverse() + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.85, 0.9)
  ) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = "Mean Monthly Depth Use (m)") -> p 


p
ggsave(filename = here("Plots", 
                       "depth use and isotopes", 
                       "mean_depth_month_d13c.png"), plot = p, 
       width = 11, height = 8.5)


write_rds(p, here("Saved Plots",
                   "c_13_depth_overall.rds"))


# ---- plot mean depth across months for d13C -----
ggplot(data = ati, aes(x = n_15, y = mean_depth)) + 
  geom_point(size = 4, aes(fill = fish_basin), 
             shape = 21, stroke = 0.8) + 
  # stat_ellipse(aes(colour = fish_basin), 
  #              type = "norm", linetype = 1,
  # linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1) +  
  geom_errorbar(aes(ymin = mean_depth - sem, 
                    ymax = mean_depth + sem), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin") + 
  scale_x_continuous(breaks = seq(9, 15, 1)) + 
  scale_y_reverse() + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.85, 0.9)
  ) + 
  labs(x = expression(paste(delta ^ 15, "N")), 
       y = "Mean Monthly Depth Use (m)") -> p1 


p1

ggsave(filename = here("Plots", 
                       "depth use and isotopes", 
                       "mean_depth_month_d15n.png"), plot = p1, 
       width = 11, height = 8.5)


write_rds(p1, here("Saved Plots",
                  "n_15_depth_overall.rds"))


# ---- look at distributions ----
ati <- ati %>% 
  filter(c_13 != is.na(c_13))

glimpse(ati)
descdist(ati$mean_depth)
descdist(ati$c_13)
descdist(ati$n_15)

ggplot(data = ati, aes(x = mean_depth)) +
  geom_histogram()

ggplot(data = ati, aes(x = c_13)) +
  geom_histogram()

ggplot(data = ati, aes(x = n_15)) +
  geom_histogram()


glimpse(ati)

# ---- create model for c_13 ~ mean_depth -----
# m <- mixed_model(fixed = dis_mean_o ~ c_13 + basin, 
#                  random = ~ c_13 | sample, 
#                  data = df_movment_heard_o,
#                   family =)
m <- glmmTMB(c_13 ~ mean_depth * fish_basin + (1|floy_tag),
             family = gaussian(link = "identity"),
             data = ati,
             REML = TRUE)
# model fit evalauteion 
res <- simulateResiduals(m)

plot(res)

Anova(mod = m)
# model section 
m1 <- update(m, ~ mean_depth + (1|floy_tag), 
             # control = glmmTMBControl(optimizer=optim,
             #                          optArgs = list(method = "BFGS")
                                      )
)
m2 <- update(m, ~ mean_depth)

# m1 <- glmmTMB(mean_depth ~ c_13 * basin + (1|sample),
#               data = ati,
#               family = Gamma(link = "log"),
#               REML = TRUE)
# 
# res1 <- simulateResiduals(m1)
# 
# plot(res1)
# Anova(m1)
# 
# 
# res2 <- simulateResiduals(m2)

# make list of only the models you hav! 



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
                                  "d13C_depth_model_selection.xlsx"))

# ---- create specific stuff for model saving -----
car::Anova(m)
summary(m)

main_effects <- tidy(car::Anova(m))



ind_effects <- tidy(m)


# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_depth_lmer_main_effect.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_depth_lmer_ind_effect.xlsx"))

# ---- Nitrogen vs movement -----

m3 <- glmmTMB(n_15 ~ mean_depth * fish_basin + (1|floy_tag),
              data = ati,
              family = Gamma(link = "log"),
              # control = glmmTMBControl(optimizer=optim,
              #                          optArgs = list(method = "BFGS")),
              REML = TRUE)

# model fit evaluation 
res3 <- simulateResiduals(m3)
plot(res3)
Anova(m3)
# model section 
m4 <- update(m3, ~ mean_depth + (1|floy_tag), 
             # control = glmmTMBControl(optimizer=optim,
             #                          optArgs = list(method = "BFGS"))
             )
m5 <- update(m3, ~ mean_depth)


model_list <- list(m3, m4, m5)
# give the elements useful names
names(model_list) <- c("m3", 
                       "m4", "m5")

glance(m3)
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
                                  "d15N_depth_model_selection.xlsx"))

# ---- create specific stuff for model saving -----
car::Anova(m3)
# summary(m3)

main_effects <- tidy(car::Anova(m3))



ind_effects <- tidy(m3)


# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_depth_lmer_main_effect.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_depth_lmer_ind_effect.xlsx"))

main_effects <- tidy(car::Anova(m4))



ind_effects <- tidy(m4)


# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_depth_lmer_main_effect_no_basin.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_depth_lmer_ind_effect_no_basin.xlsx"))

