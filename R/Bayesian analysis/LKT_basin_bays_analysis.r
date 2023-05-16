# ---- bring in packages -----
library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(rjags)
library(SIBER)
library(sf)
library(tidyr)

# ---- bring SIBRE object & isotopic data -----
df_SIBRE_obj <- read_rds(here("Saved Data", 
                              "lkt_pap_SIBRE_obj.rds"))


glimpse(df_SIBRE_obj)
df_SIBRE_obj

df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia.rds"))

glimpse(df)
# ---- determine convex TA, SEA and SEAc -----
# hull total area (TA), standard Ellipse Area (SEA), and 
# its corresponding small sample size corrected version (SEAc)

group.ML <- groupMetricsML(df_SIBRE_obj)


group_ml <- bind_cols(metric = rownames(group.ML), 
                      as_tibble(group.ML)) %>% 
  pivot_longer(cols = -metric, 
               names_to = "group", 
               values_to = "value") %>% 
  mutate(
    basin = factor(
      case_when(
        group == 1.1 ~ "East",
        group == 1.2 ~ "West",
        group == 1.3 ~ "North"
      ), level = c("West", "North","East")
    )
  )

group_ml
# these metrics are informative but we are unable to compare niches 
# statistical sense as we lack a measure of the uncertainty 
# around each estimate.

# ---- Use Bayesian Inferences to compare basin Niches ------
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
ellipses.posterior <- siberMVN(df_SIBRE_obj, parms, priors)

glimpse(ellipses.posterior)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
sea.b <- siberEllipses(ellipses.posterior)


sea.b_tb <- as_tibble(sea.b) %>% 
  rename(West = V1, 
         North = V2, 
         East = V3) %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(cols = -id,
               names_to = "basin", 
               values_to = "sea_b"
  ) %>% 
  mutate(id = 1:nrow(.), 
         basin = factor(basin, 
                        level = c("West", "North", "East")
         )
  )







siberDensityPlot(sea.b, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

points(1:ncol(sea.b), group.ML[3,], col="red", pch = "x", lwd = 2)


ggplot(data = sea.b_tb, aes(x = basin, y = sea_b)) +
  
  geom_violin(width = 0.15, outlier.colour = NA) + 
  # stat_summary(
  #   geom = "point", fun = mean, 
  #   colour = "black", size = 3
  # )
  
  geom_point(data = group_ml %>% 
               filter(metric == "SEAc"), aes(x = basin, y = value), 
             colour = "red", shape = 4, size = 3) + 
  scale_y_continuous(breaks = seq(0, 15, 5), 
                     limits = c(0, 15)) + 
  theme_classic(base_size = 15) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title =  "SIBER ellipses on each group",
       x = "Basin", 
       y = expression("Standard Ellipse Area " ('\u2030' ^2))
  )


# 
# SEA.B.credibles <- lapply(
#   as.data.frame(SEA.B), 
#   function(x,...){tmp<-hdrcde::hdr(x)$hdr},
#   prob = cr.p)
# 
# print(SEA.B.credibles)
# 
# 
# 
# SEA.B.modes <- lapply(
#   as.data.frame(SEA.B), 
#   function(x,...){tmp<-hdrcde::hdr(x)$mode},
#   prob = cr.p, all.modes=T)
# 
# print(SEA.B.modes)
# 
# 
# Pg1.1_lt_g1.2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
# print(Pg1.1_lt_g1.2)
# Pg1.1_lt_g1.3 <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B)
# print(Pg1.1_lt_g1.3)
# 
# 
# 
# Pg1.1_lt_g2.1 <- sum( SEA.B[,2] < SEA.B[,3] ) / nrow(SEA.B)
# print(Pg1.1_lt_g2.1)
# 
# 
# overlap.G1.2.G1.3 <- maxLikOverlap("1.2", "1.3", df_SIBRE_obj, p = 0.95, n =)
# 
# prop.of.first <- as.numeric(overlap.G1.2.G1.3["overlap"] / overlap.G1.2.G1.3["area.1"])
# print(prop.of.first)
# 
# 
# bayes.overlap.G2.G3 <- bayesianOverlap("1.2", "1.3", ellipses.posterior, 
#                                        draws = 10, p.interval = 0.95,
#                                        n = 360)
# print(bayes.overlap.G2.G3)



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

ellipse_df <- ellipse_df %>% 
  mutate(
    basin = factor(
      case_when(
        group == 1 ~ "East",
        group == 2 ~ "West",
        group == 3 ~ "North"
      ), level = c("East", "West", "North")
    )
  )

ellipse_df

unique(ellipse_df$rep)

ggplot(data = df, aes(x = c_13, y = n_15)) + 
  geom_point(size = 3, alpha = 0.70, aes(colour = basin)) +
  geom_polygon(data = ellipse_df,
               mapping = aes(x = c_13, y = n_15,
                             group = interaction(rep, basin),
                             color = basin,
                             fill = NULL),
               fill = NA,
               alpha = 0.2) + 

scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                       option = "D", name = "Basin") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Basin") + 
  scale_x_continuous(breaks = rev(seq(-24, -34, -1))) +
  scale_y_continuous(breaks = seq(7, 15, 1)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.93, 0.15), 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p 
  

p
ggsave(filename = here("Plots", 
                       "Bayesian Niche cluster",
                       "estimated_bayes_ellipse.png"), plot = p, 
       height = 8.5, width = 11)




# ---- use sf to create sf objects to get ellipse size -----


ellipse_df_sf <- st_as_sf(ellipse_df, coords = c("c_13", "n_15")) %>% 
  dplyr::group_by(basin, rep) %>%
  dplyr::summarize(do_union = FALSE) %>% 
  sf::st_cast("POLYGON") %>% 
  ungroup() %>% 
  mutate(
    area = st_area(.)
  ) %>% 
  select(basin, rep, area, geometry)




ggplot() + 
  geom_sf(data = ellipse_df_sf, aes(colour = basin), fill = NA, linewidth = 1)
ggplot(data = ellipse_df_sf, aes(x = area)) + 
  geom_histogram()
ggplot(data = ellipse_df_sf, aes(y = area, x = basin)) + 
  geom_boxplot(aes(fill = basin), alpha = 0.5, width = 0.1) + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "Basin")


tapply(ellipse_df_sf$area, ellipse_df_sf$basin, mean)


