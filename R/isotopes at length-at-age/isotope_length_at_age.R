# ---- bring in packages -----
library(dplyr)
library(ggplot2)
library(here)
library(MixSIAR)
library(nicheROVER) # possible look at using this vs SIBER
library(readr)
library(rjags)
library(SIBER)
library(SIBERG)
library(sf)
library(tidyr)


# ---- bring in length @ age data ----- 

la <- read_rds(here("Saved Data", 
                    "cleaned_length_at_age_raw.rds"))

df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))

vb <- read_rds(here("Saved Data", 
                    "von_b_bays_coef.rds"))


# look at data 
glimpse(df)
glimpse(la)
glimpse(vb)


# ---- manipulate vb coefficents to long and wide ------


# vb_long <- vb %>% 
#   pivot_longer(cols = -metric, 
#                names_to = "coef", 
#                values_to = "est")
# 
# 
# vb_wide <- vb_long %>% 
#   group_by(metric) %>% 
#   mutate(id = row_number()) %>% 
#   pivot_wider(names_from = "metric", 
#               values_from = "est") 
# 
# vb_long
# vb_wide

vb_model_param <- vb %>% 
  filter(metric %in% "50%")


vb_model_param


glimpse(df)

df <- df %>% 
  mutate(
    linf = vb_model_param$linf, 
    k = vb_model_param$k, 
    t0 = vb_model_param$t0, 
    est_age = seq(0, 25, length.out = 65), 
    est_lgn = linf * (1 - exp(-k * (est_age - t0))), 
    age_est = (log(1 - (tl_mm / linf)) / -k ) + t0
)


write_rds(df, file = here("Saved data", 
                         "tag_iso_est_age.rds"))
glimpse(df)
min(df$tl_mm)

df_test <- tibble(
  tl = seq(350, 750, 10), 
  linf = vb_model_param$linf, 
  k = vb_model_param$k, 
  t0 = vb_model_param$t0
)


df_test %>% 
  mutate(
    age_est = (log(1 - (tl / linf)) / -k ) + t0
  ) %>% 
  ggplot() + 
  geom_line(aes(x = est_lgn, y = age_est),
            linewidth = 1) 

ggplot(data = df) + 
  geom_line(aes(x = tl_mm, y = age_est),
               linewidth = 1)  +
  geom_point(aes(x = tl_mm, y = age_est,
                 colour = c_13),
             size = 3, alpha = 0.5) +

  scale_colour_viridis_c(begin = 0.25, end = 0.85, 
                         option = "D", 
                         name = expression(paste(delta ^ 13, "C"))) +
  theme_bw(base_size = 15) + 
  labs(x = "Estimated Age (yr)", 
       y = "Total Length (mm)") -> p 

p
ggplot(data = df) + 
  geom_line(aes(x = est_age, y = est_lgn),
               linewidth = 1)  +
  geom_point(aes(x = age_est, y = tl_mm, 
                 colour = n_15), 
             size = 3, alpha = 0.5) + 
  
  scale_colour_viridis_c(begin = 0.25, end = 0.85, 
                         option = "D", 
                         name = expression(paste(delta ^ 15, "N"))) +
  theme_bw(base_size = 15) + 
  labs(x = "Estimated Age (yr)", 
       y = "Total Length (mm)") -> p1
p1
ggsave(filename = here("Plots", 
                       "est_age_lgt_c13.png"), plot = p, 
       height = 4.37, width = 8.34)
ggsave(filename = here("Plots", 
                       "est_age_lgt_n15.png"), plot = p1, 
       height = 4.37, width = 8.34)






# 
# est_lgn = linf * (1 - exp(-k * (est_age - t0))) 
# 
# (est_lgn / linf) + 1 = exp(-k * (est_age - t0))
# log((est_lgn / linf) + 1) = (-k) * (est_age - t0))
# 
# log((est_lgn / linf) + 1) + t0 / -k