# ---- bring in packages -----
library(dplyr)
library(ellipse)
library(forcats)
library(factoextra)
library(ggplot2)
library(here)
library(purrr)
library(patchwork)
library(readr)
library(SIBER)
library(sf)
library(tidyr)

# ---- bring in data ----


df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))

glimpse(df)

# ---- select isotopes to analyze for cluster analysis ----- 
df_cn <- df %>% 
  select(c_13, n_15) %>% 
  drop_na()

#  scale the data
df_scale <- scale(df_cn)



### Elbow method (look at the knee)
# Elbow method for kmeans
fviz_nbclust(df_scale, kmeans) +
  geom_vline(xintercept = 3, linetype = 2)

km_cn <- kmeans(df_cn, 3, nstart = 25)


df_cn$groups <- km_cn$cluster


ggplot(data = df_cn, aes(x = c_13, y = n_15)) + 
  geom_point(aes(colour = as.factor(groups)), size = 3) + 
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "Groups")


# ---- save kcluster mean dataframe ----- 
# 
# write_rds(x = df_cn, file = here("Saved Data", 
#                                  "cleaned_lkt_tagged_sia_cluster.rds"))
# ---- create SIBRE object for Bayesian analysis of clusters -----
df_cn <- df_cn %>% 
  rename(iso1 = c_13, 
         iso2 = n_15, 
         group = groups) %>% 
  mutate(
    community = 1
  ) %>% 
  as.data.frame()

df_sibre <- createSiberObject(df_cn)
# ---- set up parameters to run JAGS -----

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10 ^ 6   # number of iterations to run the model for

parms$n.burnin <- 10^6 # discard the first set of values
parms$n.thin <- 100000    # thin the posterior by this many
parms$n.chains <- 4       # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

parms
priors

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(df_sibre, parms, priors)

glimpse(ellipses.posterior)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
sea.b <- siberEllipses(ellipses.posterior)


sea.b_tb <- as_tibble(sea.b) %>% 
  rename(group_1 = V1, 
         group_2 = V2, 
         group_3 = V3) %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(cols = -id,
               names_to = "group", 
               values_to = "sea_b"
  ) %>% 
  mutate(id = 1:nrow(.), 
         group = factor(group)
  )


group.ML <- groupMetricsML(df_sibre)




siberDensityPlot(sea.b, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

points(1:ncol(sea.b), group.ML[3,], col="red", pch = "x", lwd = 2)


# ---- extract ellispse for plotting ---- 
# how many of the posterior draws do you want?
n.posts <- 14 # max number is 80 



# decide how big an ellipse you want to draw
p.ell <- 0.95

# for a standard ellipse use
# p.ell <- pchisq(1, 2)

str(ellipses.posterior)


# a list to store the results
all_ellipses <- list()
# i <- 1
# j <- 1
# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}


ellipse_df <- bind_rows(all_ellipses, .id = "id") %>% 
  rename(
    c_13 = x,
    n_15 = y
  )
#   basin = rep
# ) %>%
# mutate(
#  basin = factor(
#    case_when(
#    basin == 1 ~ "East",
#    basin == 2 ~ "West",
#    basin == 3 ~ "North"
#  ), level = c("East", "West", "North")
# )
# )

# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)


ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

df_cn$group <- as.factor(df_cn$group)
ellipse_df$group <- as.factor(ellipse_df$group)


# ---- plot ellipse from Baysian estimates ----- 
ggplot(data = df_cn, aes(x = iso1, y = iso2)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = group)) +
  geom_polygon(data = ellipse_df,
               mapping = aes(x = c_13, y = n_15,
                             group = interaction(rep, group),
                             color = group,
                             fill = NULL),
               fill = NA,
               alpha = 0.2) + 
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "Cluster") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Cluster") + 
  scale_x_continuous(breaks = rev(seq(-20, -34, -1))) +
  scale_y_continuous(breaks = seq(6, 20, 2)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.10, 0.91), 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p 


# p
ggsave(filename = here("Plots", 
                       "Bayesian Niche cluster", 
                       "estiatmed_bays_ellisps_cluster.png"), plot = p, 
       height = 8.5, width = 11)


# ---- use sf to create sf objects to get ellipse size -----

ellipse_df_sf <- st_as_sf(ellipse_df, coords = c("c_13", "n_15")) %>% 
  dplyr::group_by(group, rep) %>%
  dplyr::summarize(do_union = FALSE) %>% 
  sf::st_cast("POLYGON") %>% 
  ungroup() %>% 
  mutate(
    area = st_area(.)
  ) %>% 
  select(group, rep, area, geometry)

ellipse_df_sf


ggplot() + 
  geom_sf(data = ellipse_df_sf, aes(colour = group), fill = NA, linewidth = 1)
ggplot(data = ellipse_df_sf, aes(x = area)) + 
  geom_histogram()
ggplot(data = ellipse_df_sf, aes(y = area, x = group)) + 
  geom_boxplot(aes(fill = group), alpha = 0.5, width = 0.1) + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Group")


tapply(ellipse_df_sf$area, ellipse_df_sf$group, mean)
