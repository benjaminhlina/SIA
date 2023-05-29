# ---- bring in packages -----
library(dplyr)
library(DHARMa)
library(fitdistrplus)
library(forcats)
library(emmeans)
library(ggplot2)
library(here)
library(nicheROVER) # possible look at using this vs SIBER
library(purrr)
library(patchwork)
library(readr)
library(tidyr)

# ---- bring in data ----


df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))

glimpse(df)


glimpse(df)
df <- df %>% 
  filter(c_13 != is.na(c_13) & n_15 != is.na(n_15))

descdist(df$c_13)
descdist(df$n_15)

ggplot(data = df, aes(x = c_13)) +
  geom_histogram()

ggplot(data = df, aes(x = n_15)) +
  geom_histogram()
mean(df$c_13)
sum <- df %>% 
  group_by(
    basin
  ) %>% 
  summarise(
    mean_c = mean(c_13), 
    sem_c = sd(c_13) / sqrt(n()),
    mean_n = mean(n_15), 
    sem_n = sd(n_15) / sqrt(n())
  ) %>% 
  ungroup()

sum
# model 

m <- lm(c_13 ~ basin, data = df, 
        contrasts = list(basin = "contr.sum"))

plot(residuals(m))
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))

car::Anova(m, type = 3)
summary(m)


em <- emmeans(m, ~ basin, type = "sidak")

contrast(em, method = "pairwise", type = "fdr")


m1 <- glm(n_15 ~ basin, data = df, family = Gamma(link = "inverse"), 
          contrasts = list(basin = "contr.sum"))

summary(m1)
car::Anova(m1, type = 3)
res <- simulateResiduals(m1)
plot(res)
