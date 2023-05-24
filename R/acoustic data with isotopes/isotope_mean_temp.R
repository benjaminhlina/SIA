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

# ---- nitrogren vs. temp plot -----
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

gs_d13c <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)
gs_d13c

# ---- create specific stuff for model saving -----

drop1(m, .~., test = "F")
shapiro.test(residuals(m))

par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))

Anova(mod = m, type = 3)



me_d13c <- tidy(car::Anova(m,  type = 3))
me_d13c

# ---- Nitrogen vs movement with large fish  -----
glimpse(ati)


m3 <- lm(n_15 ~ mean_temp * fish_basin, 
         data = ati,
         contrasts = list(fish_basin = contr.sum))

Anova(m3, type = 3)

# model section 
m4 <- update(m3, ~ mean_temp)
m5 <- update(m3, ~ fish_basin)
m6 <- update(m3, ~ mean_temp + fish_basin)



model_list <- list(m3, m4, m5, m6)
# give the elements useful names
names(model_list) <- c("m3", 
                       "m4", "m5", "m6")

glance(m6)
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

gs_d15n_outlier <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)
gs_d15n


# ---- create specific stuff for model saving -----
car::Anova(m3, type = "III")

par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))
shapiro.test(residuals(m3))

me_d15n_outlier <- tidy(car::Anova(m3, type = "III"))
cooksD <- cooks.distance(m3)
influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential

# the one fish that we tagged that is at a higher trophic level really pulls 
# the results to accomedate this we will regress a second time without it 

# ---- Nitrogen vs movement without large fish  -----
glimpse(ati_s)


m7 <- lm(n_15 ~ mean_temp * fish_basin, 
         data = ati_s,
         contrasts = list(fish_basin = contr.sum),
)

Anova(m7, type = 3)

# model section 
m8 <- update(m7, ~ mean_temp)
m9 <- update(m7, ~ fish_basin)
m10 <- update(m7, ~ mean_temp + fish_basin)



model_list <- list(m7, m8, m9, m10)
# give the elements useful names
names(model_list) <- c("m7", 
                       "m8", "m9", "m10")

glance(m7)
# get the summaries using `lapply

summary_list <- lapply(model_list, function(x) tidy(x, parametric = TRUE))
glance_list <- lapply(model_list, glance)
# glance_list <- lapply(glance_list, function(x) mutate_at(x, 
#                                                          .vars = "logLik", 
#                                                          as.numeric))


glance_summary <- map_df(glance_list, ~as.data.frame(.x), .id = "id") %>% 
  mutate(model = lapply(model_list, formula) %>%
           as.character()
  ) %>% 
  dplyr::select(model, id:df.residual) %>% 
  arrange(AIC)


glance_summary

gs_d15n <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)
gs_d15n

# ---- create specific stuff for model saving -----
car::Anova(m7, type = "III")
par(mfrow = c(2, 2))
plot(m7)
par(mfrow = c(1, 1))
shapiro.test(residuals(m7))

me_d15n <- tidy(car::Anova(m7, type = "III"))


# ---- combine main effects and model fit together and save as an RDS ----
# we will combine the rds for all three metrics later in another script 

model_fit <- bind_rows(list(d13c = gs_d13c, 
                            d15n = gs_d15n,
                            d15n_outlier = gs_d15n_outlier), 
                       .id = "id") %>% 
  mutate(
    metric = "temp"
  ) %>% 
  dplyr::select(id, metric, model, r.squared:df.residual)


model_effects <- bind_rows(list(d13c = me_d13c,
                                d15n = me_d15n,
                                d15n_outlier = me_d15n_outlier), 
                           .id = "id") %>% 
  mutate(
    metric = "temp"
  ) %>% 
  dplyr::select(id, metric, term, p.value)

model_fit

model_effects

write_rds(model_fit, here("Results", 
                          "temp_isotope_model_fit.rds"))
write_rds(model_effects, here("Results", 
                              "temp_isotope_model_effects.rds"))

# ---- plot mean temp across months for d13C -----

descdist(ati$mean_temp)


m7 <- glm(mean_temp ~ c_13 * n_15,
          data = ati,
          family = Gamma()
          )

Anova(m7, type = 3)

res <- simulateResiduals(m7)

plot(res)
shapiro.test(residuals.glm(m7))



summary(ati$mean_temp)
#
nd <- expand_grid(
  c_13 = seq(-25.5, -31.5, length.out = 100),
  n_15 = seq(8.5, 15, length.out = 100),
  mean_temp = seq(4, 8, length.out = 100)

)
#
fits <- predict(m7, newdata = nd, type = "response")

preds <- bind_cols(nd, fits) %>%
  rename(
    fit = ...4
  )

preds

ggplot() +

  geom_raster(data = preds, aes(x = c_13,
                              y = n_15, fill = fit), ) +
  geom_point(data = ati, size = 4,
             aes(y = n_15, x = c_13,
                 fill = mean_temp),
             shape = 21, stroke = 0.8
             ) +
#   # stat_ellipse(aes(colour = fish_basin),
#   #              type = "norm", linetype = 1,
#   # linewidth = 1) +
#   # geom_errorbar(aes(xmin = mean_temp - sem,
#   #                   xmax = mean_temp + sem), width = 0.05) +
  scale_fill_viridis_c(name = "Temperature Use (°C)",
                       option = "D",alpha = 0.5,
                       breaks = seq(4, 8, 1),
                       limit = c(4, 8)
                       ) +
  # scale_colour_viridis_c(
  #   # begin = 0.25, end = 0.85,
  #                      option = "D",
  #                      name = "Observed\nTemperature Use (°C)",
  #                      ) +
  scale_x_continuous(breaks = rev(seq(-25, -31, -1))) + 
  scale_y_continuous(breaks = seq(8, 15, 1)) + 
  # scale_colour_viridis_c(
  # scale_y_continuous(breaks = rev(seq(-26, -31, -1))) +
  # facet_wrap(.~ fish_basin) +
  coord_cartesian(expand = FALSE) +
  theme_bw(base_size = 15) +
  theme(
    legend.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    # legend.position = c(0.85, 0.9)
  ) +
  labs(
    x = expression(paste(delta ^ 13, "C")),
    y = expression(paste(delta ^ 15, "N"))) -> p3


p3
ggsave(filename = here("Plots",
                       "temp use and isotopes",
                       "mean_temp_month_d13c_raster_pred.png"), plot = p3,
       width = 11, height = 8.5)

write_rds(p3, here("Saved Plots",
                   "d13c_d15n_temp_predicted.rds"))