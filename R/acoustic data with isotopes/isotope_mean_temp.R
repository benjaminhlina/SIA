# ---- bring in packages -----
{
  library(broom.mixed)
  library(car)
  library(dplyr)
  library(DHARMa)
  library(emmeans)
  library(fitdistrplus)
  library(ggplot2)
  library(glmmTMB)
  library(here)
  library(lme4)
  library(MASS)
  library(openxlsx)
  library(purrr)
  library(readr)
  library(robustbase)
  library(stringr)
  library(tidyr)
}
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

# ---- mean temp for each month ==== 
at_t <- at %>% 
  filter(sensor_unit %in% "°C") %>% 
  mutate(
    fish_basin = str_replace(fish_basin, " Basin", "") %>% 
      factor(., level = c("East", "West", "North"))
  ) %>% 
  group_by(floy_tag, fish_basin, month, year 
  ) %>% 
  summarise(
    temp = mean(sensor_value), 
    var = var(sensor_value),
    sd = sd(sensor_value), 
    sem = sd(sensor_value) / sqrt(n())
  ) %>% 
  ungroup()


at_temp <- at_t %>% 
  group_by(floy_tag, fish_basin, month) %>% 
  summarise(
    te = mean(temp), 
    var = var(temp),
    sd = sd(temp), 
    sem = sd(temp) / sqrt(n())
  ) %>% 
  ungroup()
at_temp

# ---- mean temp across months ----- 
at_temp_sum <- at_temp %>% 
  group_by(floy_tag, fish_basin) %>% 
  summarise(
    mean_temp = mean(te), 
    sd = sd(te), 
    var_temp = var(te), 
    sem = sd(te) / sqrt(n()), 
    range = max(te) - min(te)
  ) %>% 
  ungroup()


# ---- join mean temp across months to isotope data ----
ati <- at_temp_sum %>% 
  left_join(df_short, by = c("floy_tag" = "sample")) 
# 
ati_s <- ati %>%
  filter(floy_tag != "07478")
# ---- plot mean temp across months for d13C -----
ggplot(data = ati, aes(y = c_13, x = mean_temp)) + 
  geom_point(size = 4, aes(fill = fish_basin), 
             shape = 21, stroke = 0.8) + 
  # stat_ellipse(aes(colour = fish_basin), 
  #              type = "norm", linetype = 1,
  # linewidth = 1) +
  # geom_errorbar(aes(xmin = mean_temp - sem, 
  #                   xmax = mean_temp + sem), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin") + 
  scale_y_continuous(breaks = rev(seq(-26, -31, -1))) + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.85, 0.9)
  ) + 
  labs(y = expression(paste(delta ^ 13, "C")), 
       x = "Mean Monthly Temperature Use (°C)") -> p 


p
ggsave(filename = here("Plots", 
                       "temp use and isotopes", 
                       "mean_temp_month_d13c.png"), plot = p, 
       width = 11, height = 8.5)


write_rds(p, here("Saved Plots",
                  "c_13_temp_overall.rds"))

# ---- nitrogren vs. temp -----
ggplot(data = ati, aes(y = n_15, x = mean_temp)) + 
  geom_point(size = 4, aes(fill = fish_basin), 
             shape = 21, stroke = 0.8) + 
  # stat_ellipse(aes(colour = fish_basin), 
  #              type = "norm", linetype = 1,
  # linewidth = 1) +
  # geom_errorbar(aes(xmin = mean_temp - sem, 
  #                   xmax = mean_temp + sem), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin") + 
  scale_y_continuous(breaks = seq(9, 15, 1)) + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.85, 0.9)
  ) + 
  labs(y = expression(paste(delta ^ 15, "N")), 
       x = "Mean Monthly Temperature Use (°C)") -> p1


p1
ggsave(filename = here("Plots", 
                       "temp use and isotopes", 
                       "mean_temp_month_d15n.png"), plot = p1, 
       width = 11, height = 8.5)


write_rds(p1, here("Saved Plots",
                   "n_15_temp_overall.rds"))



# ---- look at distributions ----
ati <- ati %>% 
  filter(c_13 != is.na(c_13))



glimpse(ati)
ati <- ati %>% 
  mutate(
    temp_p_sem = mean_temp + sem,
    temp_m_sem = mean_temp - sem, 
    sem_range = temp_p_sem - temp_m_sem,
    sd_range =  (mean_temp + sd) - (mean_temp - sd)
  )

descdist(ati$c_13)
descdist(ati$n_15)


ggplot(data = ati, aes(x = c_13)) +
  geom_histogram()

ggplot(data = ati, aes(x = n_15)) +
  geom_histogram()


glimpse(ati)

# ---- create model for c_13 ~ mean_temp -----

options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(c_13 ~ mean_temp * fish_basin,
        contrasts = list(fish_basin = "contr.sum"),
        data = ati)
drop1(m, .~., test = "F")
# model section 
m1 <- update(m, ~ mean_temp) 
m2 <- update(m, ~ fish_basin) 
m3 <- update(m, ~ mean_temp + fish_basin) 

# ---- create model list for model selection ------
model_list <- list(m, m1, m2, m3)
# give the elements useful names
names(model_list) <- c("m", 
                       "m1", "m2", "m3")

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
                                  "d13C_temp_model_selection.xlsx"))

# ---- create specific stuff for model saving -----
car::Anova(m)
summary(m)

main_effects <- tidy(car::Anova(m))



ind_effects <- tidy(m)


# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_temp_lmer_main_effect.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d13C Habitat",
                                  "d13C_temp_lmer_ind_effect.xlsx"))

# ---- Nitrogen vs movement -----
glimpse(ati)

m3 <- glmmTMB(n_15 ~ mean_temp * fish_basin + (1|floy_tag) +
                # (1|temp_p_sem) + (1| temp_m_sem),
              data = ati,
              family = Gamma(link = "identity"),
              # control = glmmTMBControl(optimizer=optim,
              #                          optArgs = list(method = "BFGS")),
              REML = TRUE)

# model fit evaluation 
res3 <- simulateResiduals(m3)
plot(res3)
Anova(m3, type = "II")
# model section 
m4 <- update(m3, ~ mean_temp * sem + (1|floy_tag), 
             # control = glmmTMBControl(optimizer=optim,
             #                          optArgs = list(method = "BFGS"))
)
m5 <- update(m3, ~ mean_temp)


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
                                  "d15N_temp_model_selection.xlsx"))

# ---- create specific stuff for model saving -----
car::Anova(m3)
# summary(m3)

main_effects <- tidy(car::Anova(m3))



ind_effects <- tidy(m3)
summary(m3)

ind_effects
# main_effects %>% 
main_effects %>% 
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_temp_lmer_main_effect.xlsx"))
ind_effects %>%
  openxlsx::write.xlsx(here::here("results",
                                  "d15N Habitat",
                                  "d15N_temp_lmer_ind_effect.xlsx"))


