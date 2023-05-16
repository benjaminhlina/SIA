# ---- bring in packages -----
library(dplyr)
library(ggplot2)
library(here)
library(MixSIAR)
library(nicheROVER)
library(readr)
library(SIBER)
library(SIBERG)
library(stringr)
library(sf)
library(tidyr)
# ---- bring in SIA data frame ------
df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))


# View(df)
glimpse(df)

# ---- select only isotopes and basin groups -----
df_SIBRE <- df %>% 
  select(c_13, n_15, basin)


# ---- change basin to numerical group, add community for SIBRE -----
df_SIBRE <- df_SIBRE %>% 
  mutate(
    group = case_when(
      basin %in% "East" ~ 1,
      basin %in% "West" ~ 2,
      basin %in% "North" ~ 3,
    ),
    community = 1
  ) %>% 
  select(c_13, n_15, group, community) %>% 
  rename(
    iso1 = c_13, 
    iso2 = n_15
  ) %>% 
  drop_na() %>% 
  as.data.frame() # working in SIBRE it cannnot handle tibbles 

# ---- create SIBRE object -----
df_SIBRE_obj <- createSiberObject(df_SIBRE)
df_SIBRE_obj

# ---- save SIBRE Object for Bayesian analysis -----
write_rds(df_SIBRE_obj, here("Saved Data", 
                             "lkt_pap_SIBRE_obj.rds"))

# ---- SIBRE object using base plot ------
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, 
                             lty = 1, lwd = 2)
group.hulls.args    <- list(lty = 2, col = "grey20")



par(mfrow=c(1,1))
plotSiberObject(df_SIBRE_obj,
                ax.pad = 2,
                hulls = F,
                # community.hulls.args = community.hulls.args,
                # ellipses = T,
                # group.ellipses.args = group.ellipses.args,
                group.hulls = T,
                group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)



# ---- plot ellipses with ggplot -----

ggplot(data = df, aes(x = c_13, y = n_15)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin)) +
  stat_ellipse(aes(group = basin, 
                   color = basin),
               fill = NA, 
               alpha = 0.70, 
               level = 0.95,
               type = "norm",
               geom = "polygon") +  
  scale_colour_viridis_d(begin = 0.25, end = 0.85, 
                         option = "D", name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Basin") + 
  scale_x_continuous(breaks = rev(seq(-24, -34, -1))) +
  scale_y_continuous(breaks = seq(7, 15, 1)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.83), 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p 
p


ggsave(filename = here("Plots",
                       "sai__elips.png"),
       plot = p,
       height = 4.37, width = 8.34)

# ---- create convex hull using {sf}, this was quite clever Ben -----

df_convex_hull <- df %>% 
  filter(c_13 != is.na(c_13)) %>% 
  st_as_sf(coords = c("c_13", "n_15")) %>% 
  group_by(basin) %>%
  summarize(geometry = st_union(geometry)) %>% 
  st_convex_hull() %>% 
  st_cast("POINT") %>%
  mutate(c_13 = st_coordinates(.)[,1],
         n_15 = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  as_tibble()


# ---- ggplot ellipses and convex hull ------
ggplot() + 
  geom_point(data = df, size = 3, alpha = 0.70, 
             aes(x = c_13, y = n_15, colour = basin)) +
  stat_ellipse(
    data = df, aes(x = c_13, y = n_15, 
                   group = interaction(basin, community), 
                   fill = basin, 
                   color = basin), 
    alpha = 0.10, 
    level = 0.95,
    type = "norm",
    geom = "polygon") +  
  geom_polygon(data = df_convex_hull, fill = NA, 
               colour = "black", 
               aes(y = n_15, x = c_13, linetype = basin), linewidth = 0.75) + 
  scale_colour_viridis_d(begin = 0.23, end = 0.85, 
                         option = "D", name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Basin") + 
  scale_linetype(name = "Basin") + 
  scale_x_continuous(breaks = rev(seq(-23, -34, -1))) +
  scale_y_continuous(breaks = seq(7, 15, 1)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.07, 0.86), 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p1
# p1

# ggsave(filename = here("Plots",
#                        "sai_convex_hull_elips.png"), 
#        plot = p1,
#        height = 4.37, width = 8.34)
