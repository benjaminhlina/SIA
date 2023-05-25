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
    sem = sd(md) / sqrt(n()), 
    range = max(md) - min(md)
  ) %>% 
  ungroup()
# print(n = 100, at_depth)

# ---- join mean depth across months to isotope data ----
ati <- at_depth_sum %>% 
  left_join(df_short, by = c("floy_tag" = "sample")) 


ati_s <- ati %>%
  filter(floy_tag != "07478")
# ---- plot mean depth across months for d13C -----
ggplot(data = ati, aes(y = c_13, x = mean_depth)) + 
  geom_point(size = 4, aes(fill = fish_basin), 
             shape = 21, stroke = 0.8) + 
  # stat_ellipse(aes(colour = fish_basin), 
  #              type = "norm", linetype = 1,
  # linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +  
  # geom_errorbar(aes(xmin = mean_depth - sem, 
  #                   xmax = mean_depth + sem), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin") + 
  scale_y_continuous(breaks = rev(seq(-26, -31, -1))) + 
  # scale_x_reverse() + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.85, 0.9)
  ) + 
  labs(y = expression(paste(delta ^ 13, "C")), 
       x = "Mean Monthly Depth Use (m)") -> p 


p
ggsave(filename = here("Plots", 
                       "depth use and isotopes", 
                       "mean_depth_month_d13c.png"), plot = p, 
       width = 11, height = 8.5)


write_rds(p, here("Saved Plots",
                  "c_13_depth_overall.rds"))


# ---- plot mean depth across months for d13C -----
ggplot(data = ati, aes(y = n_15, x = mean_depth)) + 
  geom_point(size = 4, aes(fill = fish_basin), 
             shape = 21, stroke = 0.8) + 
  # stat_ellipse(aes(colour = fish_basin), 
  #              type = "norm", linetype = 1,
  # linewidth = 1) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +  
  # geom_errorbar(aes(xmin = mean_depth - sem, 
  #                   xmax = mean_depth + sem), width = 0.05) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", 
                       name = "Basin") + 
  scale_y_continuous(breaks = seq(9, 15, 1)) + 
  # scale_y_reverse() + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.85, 0.9)
  ) + 
  labs(y = expression(paste(delta ^ 15, "N")), 
       x = "Mean Monthly Depth Use (m)") -> p1 


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
descdist(ati$c_13)
descdist(ati$n_15)

ggplot(data = ati, aes(x = c_13)) +
  geom_histogram()

ggplot(data = ati, aes(x = n_15)) +
  geom_histogram()


glimpse(ati)

# ---- create model for c_13 ~ mean_depth -----

m <- lm(c_13 ~ mean_depth * fish_basin,
        data = ati, 
        contrasts = list(fish_basin = "contr.sum"))
# ---- create specific stuff for model saving -----

drop1(m, .~., test = "F")
shapiro.test(residuals(m))
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
cooksD <- cooks.distance(m)
influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential

# a few outliers but their affect is pretty maginal 
car::Anova(m, type = "III")


me_d13c_d <- tidy(car::Anova(m))


# model section 
m1 <- update(m, ~ mean_depth)
m2 <- update(m, ~ fish_basin)
m3 <- update(m, ~ mean_depth + fish_basin)


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

gs_d13c_d <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)

gs_d13c_d

# ---- Nitrogen vs movement -----

descdist(ati$n_15)

m3 <- glm(n_15 ~ mean_depth * fish_basin,
         data = ati, family = Gamma(link = "identity"), 
         contrasts = list(fish_basin = "contr.sum")
         )

res <- simulateResiduals(m3)
plot(res)
# ---- create specific stuff for model saving -----
drop1(m3, .~., test = "F")
shapiro.test(residuals(m3))

par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))

car::Anova(m3, type = "III")

me_d15n_outlier <- tidy(car::Anova(m3, type = "III"))

cooksD <- cooks.distance(m3)
influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential

# our large fish that is in at a higher trophic level is quite an outlier 
# we will run the regression again without it and report both 

# model section 
m4 <- update(m3, ~ mean_depth)
m5 <- update(m3, ~ fish_basin)
m6 <- update(m3, ~ mean_depth + fish_basin)


model_list <- list(m3, m4, m5, m6)
# give the elements useful names
names(model_list) <- c("m3", 
                       "m4", "m5", "m6")

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

gs_d15n_d_outlier <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)

gs_d15n_d_outlier



# ---- create model without large fish for d15n vs depth -----

ati_s <- ati_s  %>% 
  filter(n_15 != is.na(n_15))
glimpse(ati_s)
unique(is.na(ati_s$n_15))

descdist(ati_s$n_15)
ggplot(data = ati_s, aes(x = n_15)) + 
  geom_histogram()


m7 <- lm(n_15 ~ mean_depth * fish_basin, 
         data = ati_s,
         contrasts = list(fish_basin = contr.sum)
)

Anova(m7, type = 3)
# ---- create specific stuff for model saving -----
car::Anova(m7, type = "III")
par(mfrow = c(2, 2))
plot(m7)
par(mfrow = c(1, 1))
shapiro.test(residuals(m7))

me_d15n_d <- tidy(car::Anova(m7, type = "III"))


# we can see our single outlier pulls so much that it affects the signficance
# model section 
m8 <- update(m7, ~ mean_depth)
m9 <- update(m7, ~ fish_basin)
m10 <- update(m7, ~ mean_depth + fish_basin)



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

gs_d15n_d <- glance_summary %>% 
  mutate(
    delta_AIC = AIC - first(AIC), 
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) %>% 
  dplyr::select(model:AIC, delta_AIC, AIC_weight, BIC:df.residual)
gs_d15n_d

# ---- combine main effects and model fit together and save as an RDS ----
# we will combine the rds for all three metrics later in another script 

model_fit <- bind_rows(list(d13c = gs_d13c_d, 
                            d15n = gs_d15n_d,
                            d15n_outlier = gs_d15n_d_outlier), 
                       .id = "id") %>% 
  mutate(
    metric = "depth"
  ) %>% 
  dplyr::select(id, metric, model, r.squared:df.residual)

model_fit


model_effects <- bind_rows(list(d13c = me_d13c_d, 
                                d15n = me_d15n_d,
                                d15n_outlier = me_d15n_outlier), 
                           .id = "id") %>% 
  mutate(
    metric = "depth"
  ) %>% 
  dplyr::select(id, metric, term, p.value)

model_fit

model_effects

write_rds(model_fit, here("Results", 
                          "depth_isotope_model_fit.rds"))
write_rds(model_effects, here("Results", 
                              "depth_isotope_model_effects.rds"))

# ---- plot mean temp across months for d13C -----
descdist(ati$mean_depth)


m11 <- glm(mean_depth ~ c_13 * n_15,
          data = ati,
          family = Gamma(link = "log")
)

Anova(m7, type = 3)

res <- simulateResiduals(m7)

plot(res)
shapiro.test(residuals.glm(m7))
summary(ati$n_15)

# Can't predict as the data is not signgicant to predict 

# nd <- expand_grid(
#   c_13 = seq(-25.5, -31.5, length.out = 100),
#   n_15 = seq(8.5, 15, length.out = 100),
#   mean_depth = seq(0, 20, length.out = 100)
# )
# #
# fits <- predict(m7, newdata = nd, type = "response")
# 
# preds <- bind_cols(nd, fits) %>%
#   rename(
#     fit = ...4
#   )
# 
# preds
# #
# #
# #
# #
# #
# #
# #
# ggplot() +
#   geom_raster(data = preds, aes(x = c_13,
#                                 y = n_15, fill = fit), ) +
#   geom_point(data = ati, size = 4,
#              aes(y = n_15, x = c_13,
#                  fill = mean_depth),
#              shape = 21, stroke = 0.8
#   ) +
#   #   # stat_ellipse(aes(colour = fish_basin),
#   #   #              type = "norm", linetype = 1,
#   #   # linewidth = 1) +
#   #   # geom_errorbar(aes(xmin = mean_depth - sem,
#   #   #                   xmax = mean_depth + sem), width = 0.05) +
#   scale_fill_viridis_c(name = "Depth Use (m)",
#                        option = "D",alpha = 0.5 ,trans = "reverse", 
#                        direction = 1
#                        # breaks = seq(4.5, 7.5, 1),
#                        # limit = c(4, 8)
#   ) +
#   scale_x_continuous(breaks = rev(seq(-25, -31, -1))) + 
#   scale_y_continuous(breaks = seq(8, 15, 1)) + 
#   # scale_colour_viridis_c(
#   #   # begin = 0.25, end = 0.85,
#   #   option = "D",
#   #   trans = "reverse", 
#   #   direction = 1,
#   #   name = "Observed\nDepth Use (m)",
#   # ) +
#   # scale_y_continuous(breaks = rev(seq(-26, -31, -1))) +
#   
# 
#   # facet_wrap(.~ fish_basin) +
#   coord_cartesian(expand = FALSE) +
#   theme_bw(base_size = 15) +
#   theme(
#     legend.title = element_text(hjust = 0.5),
#     panel.grid = element_blank(),
#     # legend.position = c(0.85, 0.9)
#   ) +
#   labs(
#     x = expression(paste(delta ^ 13, "C")),
#     y = expression(paste(delta ^ 15, "N"))) -> p3
# 
# 
# p3
# ggsave(filename = here("Plots",
#                        "depth use and isotopes",
#                        "mean_depth_month_d13c_raster_pred.png"), plot = p3,
#        width = 11, height = 8.5)
# 
# write_rds(p3, here("Saved Plots",
#                   "d13c_d15n_depth_predicted.rds"))
# 
