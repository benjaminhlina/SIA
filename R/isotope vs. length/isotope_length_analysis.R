# ---- bring in packages -----
library(broom)
library(car)
library(dplyr)
library(emmeans)
library(ggplot2)
library(here)
library(lme4)
library(readr)

# ---- bring in dataframe ------

df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))


glimpse(df)

# ---- isotopes by basin with basic ellipses plot -----
ggplot(data = df, aes(x = c_13, y = n_15)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin)) +
  # stat_ellipse(aes(colour = basin), type = "norm", linetype = 2,  
               # size = 1) + 
  stat_ellipse(aes(colour = basin), type = "t",  linewidth = 1) + 
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Basin") + 
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.83), 
          legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p 
p

# ggsave(filename = here("Plots",
#                        "sai_test_plot.png"), plot = p,
#        height = 4.37, width = 8.34)



# ---- clean data for d15N vs length analyzis ----- 
glimpse(df)

df_1 <- df %>%
  filter(!n_15 == is.na(n_15))

# ---- plot regression of d15N vs lenght ---- 
ggplot(data = df, aes(x = tl_mm, y = n_15)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin)) + 
  stat_smooth(method = "lm", linewidth = 1, colour = "black") +
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", name = "Basin") + 
  theme_bw(base_size = 15) +
  
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "R2")),
                        label.x = 0.125,
                        label.y = 0.97) +
  scale_y_continuous(breaks = seq(9, 15, 1)) + 
  scale_x_continuous(breaks = seq(350, 750, 50)) + 
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.83),
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = "Total Length (mm)", 
       y = expression(paste(delta ^ 15, "N"))) -> p1
p1 

# ggsave(filename = here("Plots",
#                        "length_vs_N_15.png"), plot = p1,
#        height = 4.37, width = 8.34)
# ---- assess distribution and run analysis on all points ----- 

fitdistrplus::descdist(df_1$n_15)
norm <- fitdistrplus::fitdist(df_1$n_15, distr = "norm")

plot(norm)


ggplot(data = df_1, aes(x = n_15)) + 
  geom_histogram()


m <- lm(data = df, n_15 ~ tl_mm)
par(mfrow = c(2, 2))
plot(m)

car::Anova(m)

summary(m)

m_augment <- augment(m, df_1)

glimpse(m_augment)

# ggplot diagnosistic plots 
ggplot(data = m_augment, 
       aes(x = tl_mm, y = n_15)) + 
  geom_point() + 
  geom_smooth(method="lm", color="red",se=FALSE) +
  geom_segment(aes(xend = tl_mm, 
                   yend = .fitted), linetype="dashed") + 
  theme_classic(base_size = 15) + 
  scale_x_continuous(breaks = seq(350, 750, 50)) + 
  scale_y_continuous(breaks = seq(9, 15, 1)) + 
  labs(x = "Total Length (mm)", 
       y = expression(paste(delta ^ 15, "N")))
ggplot(data = m_augment, 
       aes(x = tl_mm, .resid)) + 
  geom_point() + 
  theme_classic(base_size = 15) + 
  theme(
    plot.title = element_text(hjust = 0.5)
  ) + 
  scale_x_continuous(breaks = seq(350, 750, 50)) + 
  # scale_y_continuous(breaks = seq(9, 15, 1)) + 
  geom_hline(yintercept = 0, color = "red", linetype='dashed') + 
  labs(y = "Residuals",
       x = "Total Length (mm)", 
       title = "Residual plot the lions regression model.")


# ---- remove outliers and redo the analysis ---- 
# remove the two larger fish 
df_2 <- df_1 %>% 
  dplyr::select(sample, n_15, tl_mm, basin) %>% 
  arrange(tl_mm) %>% 
  filter(tl_mm < 675)


# assess normality and variances 
fitdistrplus::descdist(df_2$n_15)
norm <- fitdistrplus::fitdist(df_2$n_15, distr = "norm")

plot(norm)


ggplot(data = df_2, aes(x = n_15)) + 
  geom_histogram()

moments::kurtosis(df_2$n_15)
moments::skewness(df_2$n_15)
# use linear model to assess differences 
m1 <- lm(data = df_2, n_15 ~ tl_mm)
par(mfrow = c(2, 2))
plot(m1)

Anova(m1)

summary(m1)

m_augment_2 <- augment(m1, df_2)

glimpse(m_augment_2)
# add in basin as a covariate 
m2 <- lm(data = df_2, n_15 ~ tl_mm * basin)
par(mfrow = c(2, 2))
plot(m2)

Anova(m2)

summary(m2)


m_augment_3 <- augment(m2, df_2)

glimpse(m_augment_2)


compares <- emmeans(m2, pairwise ~ basin)
compares


# ---- plot multiple regressions ---- 
ggplot(data = df_2, aes(x = tl_mm, y = n_15)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin)) + 
  stat_smooth(method = "lm", linewidth = 1, colour = "black") +
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", name = "Basin") + 
  theme_bw(base_size = 15) +
  
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "R2")),
                        label.x = 0.125,
                        label.y = 0.97) +
  scale_y_continuous(breaks = seq(9, 15, 1)) + 
  scale_x_continuous(breaks = seq(350, 750, 50)) + 
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.83),
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = "Total Length (mm)", 
       y = expression(paste(delta ^ 15, "N"))) -> p2
p2

ggplot(data = df_2, aes(x = tl_mm, y = n_15, group = basin)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin),) + 
  stat_smooth(method = "lm", linewidth = 1, aes(colour = basin), 
              # fullrange = TRUE
              ) +
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", name = "Basin") + 
  theme_bw(base_size = 15) +
  
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "R2")),
                        label.x = 0.115, 
                        label.y = c(0.935, 0.905, 0.875),
                        inherit.aes = TRUE) +
  scale_y_continuous(breaks = seq(9, 15, 1)) + 
  scale_x_continuous(breaks = seq(350, 750, 50)) + 
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.91),
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = "Total Length (mm)", 
       y = expression(paste(delta ^ 15, "N"))) -> p3
p3

# ggsave(filename = here("Plots",
#                        "length and isotopes",
#                        "length_vs_N_15_basin.png"), plot = p3,
#        height = 8.5, width = 11)

# ---- analyze d13C vs lenght ----- 

ggplot(data = df, aes(x = tl_mm, y = c_13)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin)) + 
  # stat_ellipse(aes(colour = basin), type = "norm", linetype = 1,
  # linewidth = 1) +
  stat_smooth(method = "lm", linewidth = 1, colour = "black") +
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", name = "Basin") + 
  theme_bw(base_size = 15) +
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "R2")),
                        label.x = 0.125,
                        label.y = 0.97) +
  scale_y_continuous(breaks = rev(seq(-23, -34, -1))) + 
  scale_x_continuous(breaks = seq(350, 750, 50)) + 
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.83),
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = "Total Length (mm)", 
       y = expression(paste(delta ^ 13, "C"))) -> p4
p4

# ggsave(filename = here("Plots",
#                        "length_vs_C_13.png"), plot = p2,
#        height = 4.37, width = 8.34)


# ---- run mixxed effects d15N vs tl ---- 

m3 <- lmer(n_15 ~ tl_mm + (1|basin), data = df_2)
Anova(m3)


summary(m3)
plot(m3)


