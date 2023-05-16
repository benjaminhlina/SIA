

# ---- bring in packages -----
library(dplyr)
library(ellipse)
library(forcats)
library(ggplot2)
library(here)
library(nicheROVER) # possible look at using this vs SIBER
library(purrr)
library(patchwork)
library(readr)
library(tidyr)

# ---- bring in data ----


df <- read_rds(here("Saved Data", 
                    "cleaned_lkt_tagged_sia_cluster.rds"))

glimpse(df)

# ---- calculated the isotopic means for each groups -----

aggregate(df[1:2], df[3], mean, na.rm = TRUE)



# ---- generate parameter draws from the "default" posteriors of each fish ----
nsample <- 1e3

fish_par <- tapply(1:nrow(df), df$groups,
                   function(ii) niw.post(nsamples = nsample, 
                                         X = df[ii, 1:2]))
# clrs <- c("black", "red", "blue")
# niche.par.plot(fish_par, col = clrs, plot.index = 2)
# str(fish_par)
# ---- separate niw.post object into sigma and mu tibbles ----- 
# sigma

df_sigma <- map(fish_par, pluck, 2) %>% 
  imap(~ as_tibble(.x) %>% 
         mutate( 
           metric = "sigma", 
           id = c("c_13", "n_15"),
           groups = .y
         )
  ) %>%
  bind_rows() %>% 
  pivot_longer(cols = -c("id", "groups", "metric"),
               names_to = "isotope", 
               values_to = "post_sample"
  )  %>% 
  separate(isotope, into = c("isotopes", "sample_number"), sep = "\\.")


df_sigma_cn <- df_sigma %>% 
  filter(id != isotopes)
df_sigma_wide <- df_sigma %>%
  select(id:post_sample) %>% 
  pivot_wider(names_from = id, 
              values_from = post_sample)

# df_sigma_n <- df_sigma %>% 
#   filter(id %in% "n_15" & isotopes %in% "n_15")

# mu 
df_mu <- bind_rows(
  fish_par$'1'$mu %>% 
    as_tibble() %>% 
    mutate(groups = factor("1")),
  fish_par$'2'$mu %>% 
    as_tibble() %>% 
    mutate(groups = factor("2")),
  fish_par$'3'$mu %>% 
    as_tibble() %>% 
    mutate(groups = factor("3"))
) %>% 
  group_by(groups) %>% 
  mutate(
    sample_number = 1:1000
  ) %>% 
  ungroup()





# ---- density plots ggplot ----

# ggplot(data = df_niw_post, aes(x = post_sample)) + 
#   geom_density(aes(fill = groups), alpha = 0.5) + 
#   facet_grid(isotopes ~ metric, scales = "free_x")

p <- ggplot(data = df_mu, aes(x = c_13)) +
  geom_density(aes(fill = groups), alpha = 0.35) + 
  scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                       option = "D", name = "groups") + 
  scale_y_continuous(breaks = seq(0, 4, 1), 
                     # limits = c(0, 3)
  ) + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    legend.position = c(0.10, 0.88)
  ) + 
  labs(x = expression(paste(mu[~delta^13],""[C])),
       y = expression(paste("p(", mu[~delta^13],""[C], " ", "| X)")))


p1 <- ggplot(data = df_mu, aes(x = n_15)) +
  geom_density(aes(fill = groups), alpha = 0.5) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                       option = "D", name = "groups") + 
  scale_y_continuous(breaks = seq(0, 6, 1), 
                     # limits = c(0, 3)
  ) +  
  theme_bw(base_size = 15) + 
  theme(legend.position = "none", 
        panel.grid = element_blank()) + 
  labs(x = expression(paste(mu[~delta^15],""[N])),
       y = expression(paste("p(", mu[~delta^15],""[N]," ", "| X)")),
  )

p2 <- ggplot(data = df_sigma_cn, aes(x = post_sample)) +
  geom_density(aes(fill = groups), alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                       option = "D", name = "groups") + 
  scale_y_continuous(breaks = seq(0, 8, 1), 
                     # limits = c(0, 3)
  ) +  
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid = element_blank()) + 
  labs(x = expression(paste(Sigma[~delta^13],""[C], " ", ""[~delta^15],""[N])),
       y = expression(paste("p(",Sigma[~delta^13],""[C], " ", 
                            ""[~delta^15],""[N]," ", "| X)")))


p3 <- p + p1 + p2

p3
# 
# ggsave(filename = here("Plots",
#                        "nicheROVER plots",
#                        "posterior_mean_density_isotope_kmean_groups.png"),
#        height = 8.5, width = 11 * 2, plot = p3)




# 
# clrs <- c("black", "blue", "orange")
# 
# df_dat <- as.data.frame(df)
# fish_data <- tapply(1:nrow(df_dat), df$groups, function(ii) X = df_dat[ii,4:5]) 
#
# np <- niche.plot(niche.par = fish_par, niche.data = fish_data, pfrac = .05,
#            iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
#            col = clrs, xlab = expression("Isotope Ratio (per mil)"))
# 
# np
# ---- create ellipse plot from posterior samples ----- 
# decide how big an ellipse you want to draw
p.ell <- 0.95

# for a standard ellipse use
# p.ell <- pchisq(1, 2)


groups_name <- unique(df_sigma_wide$groups)

all_ellipses <- list()

# for loop over groupss and sample number within groups to create each 
# sample numbers ellipi 
for (i in 1:3) {
  
  sigma_groups <- df_sigma_wide %>% 
    filter(groups %in% groups_name[i])
  
  mu_groups <- df_mu %>% 
    filter(groups %in% groups_name[i])
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  
  for(j in 1:length(unique(sigma_groups$sample_number))) {
    sigma_ind <- sigma_groups %>%
      filter(sample_number %in% sample_number[j]) %>% 
      select(c_13, n_15)
    Sigma <- as.matrix(sigma_ind, 2,2)
    row.names(Sigma) <- c("c_13", "n_15")
    Sigma
    mu <- mu_groups %>%
      filter(sample_number %in% sample_number[j]) %>% 
      select(sample_number, c_13, n_15) %>% 
      pivot_longer(cols = -sample_number, 
                   names_to = "isotope", 
                   values_to = "mu") %>% 
      .$mu
    
    
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}


# combine ellipose list into dataframe and add groups names back in 
ellipse_df <- bind_rows(all_ellipses, .id = "id") %>% 
  mutate(
    groups = factor(id)
  ) %>% 
  as_tibble()


# randomly sample 20 ellipise out of 1000 
ellipse_df %>% 
  group_by(groups, rep) %>% 
  nest() %>%
  group_by(groups) %>% 
  slice_sample(n = 10, replace = TRUE) %>% 
  ungroup() %>% 
  unnest(cols = c(data)) -> random_ellipse 


random_ellipse

# ---- create niche plots -----
ggplot() + 
  geom_polygon(data = random_ellipse,
               mapping = aes(x = c_13, y = n_15,
                             group = interaction(rep, groups),
                             color = groups,
                             fill = NULL),
               fill = NA,
               linewidth = 0.5) + 
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "groups", 
                         alpha = 0.2
                         ) + 
  # scale_fill_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", name = "groups") + 
  scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  scale_y_continuous(breaks = seq(2, 24, 2)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p4 
# p4

ggplot() + 
  geom_density(data = df, aes(x = n_15, 
                              fill = as.factor(groups),
                              # colour = groups
  ), 
  alpha = 0.2, 
  linewidth = 0.8) +
  # scale_colour_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", name = "groups") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85,
                       option = "D", name = "groups") +
  scale_x_continuous(breaks = seq(8, 16, 1), 
                     limits = c(8, 16)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 15, "N")), 
       y = "Density") -> p5
# p5

ggplot() + 
  geom_density(data = df, aes(x = c_13, 
                              fill = as.factor(groups)
                              # colour = groups
  ), 
  alpha = 0.2, 
  linewidth = 0.8) +
  # scale_colour_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", name = "groups") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85,
                       option = "D", name = "groups") +
  scale_x_continuous(breaks = rev(seq(-22, -34, -1)),
                     limits = rev(c(-22, -34)))  +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = c(0.1, 0.84), 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = "Density") -> p6

# p6

ggplot() + 
  geom_point(data = df, aes(x = c_13, y = n_15,
                            colour = as.factor(groups)),
             size = 3, alpha = 0.70) +
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "groups") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "groups") + 
  scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N"))) -> p7



p8 <- p5 + p4 + p7 + p6

p8




# ggsave(filename = here("Plots",
#                        "nicheROVER plots",
#                        "estimated_ellispes_nicheROVER_clusters.png"), plot = p8,
#        height = 8.5, width = 11)




# sigma_list <- df_sigma_wid %>% 
#   groups_split(groups) %>%  
#   map(~ split(., as.factor(.$sample_number)))
# 
# sigma_list_n <- df_sigma_wid %>%
#   select(-isotopes) %>% 
#   groups_by(groups, sample_number) %>%
#   nest() %>% 
#   ungroups() %>% 
#   map(.$data, as.matrix(2, 2)
#   )


# map(~ split(., as.factor(.$sample_number)))
# 
# fish_par
# 
# 
# ellipse_df <- ellipse(x = k[[1]][[1000]],
#                       centre = , level = p.ell)




# ---- NICHE overlap calculation and plotting ----- 

# ---- calculate % overlap from posterior sampling ----
# Calculate how much each groupsing overlaps 
over_stat <- overlap(fish_par, nreps = nsample, nprob = 1e3, 
                     alpha = c(0.95
                               ,0.99
                               )
                     )

overlap.plot(over_stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# convert to dataframe to be able to plot and look at staistically 
over_stat_df <- over_stat %>% 
  as_tibble(rownames = "groups_a") %>% 
  mutate(
    id = 1:nrow(.), 
    groups_a = factor(groups_a)
  ) %>% 
  pivot_longer(cols = -c(id, groups_a), 
               names_to = "groups_b", 
               values_to = "mc_nr")  %>% 
  separate(groups_b, into = c("groups_c", "sample_number", "percentage"), 
           sep = "\\.") %>% 
  select(-id) %>% 
  rename(groups_b = groups_c) %>% 
  mutate(
    groups_b =  factor(groups_b), 
    mc_nr_perc = mc_nr * 100
  )

# filter for just the 95% 
over_stat_df_95 <- over_stat_df %>% 
  filter(percentage %in% "95%")

# calculate the mean overlap for each groups 
over_mean <- apply(over_stat, c(1:2, 4), mean) * 100 %>% 
  round(2) 

# convert to dataframe to be able to plot and compare 
over_mean_df <- over_mean %>% 
  as_tibble(rownames = "groups_a") %>% 
  mutate(
    id = 1:nrow(.), 
    groups_a = factor(groups_a)
  ) %>% 
  pivot_longer(cols = -c(id, groups_a), 
               names_to = "groups_b", 
               values_to = "mean_mc_nr")  %>% 
  separate(groups_b, into = c("groups_c", "percentage"), 
           sep = "\\.") %>% 
  select(-id) %>% 
  rename(groups_b = groups_c) %>% 
  mutate(
    groups_b =  factor(groups_b)
  )

# filter so that you have just the 95% 
over_mean_df_95 <- over_mean_df %>% 
  filter(percentage %in% "95%")

# calculate confidence intervals around the mean 
over_cred <- apply(over_stat * 100, c(1:2, 4), quantile, 
                   prob = c(0.025, 0.975), na.rm = TRUE)
# convert to dataframe for plotting 
over_cred_df <- over_cred %>% 
  as_tibble(rownames = "qual") %>% 
  mutate(
    id = 1:nrow(.)
  ) %>% 
  pivot_longer(cols = -c(id, qual), 
               names_to = "groups_b", 
               values_to = "quantile_mc_nr") %>% 
  separate(groups_b, into = c("groups_a", "groups_b", "percentage"), 
           sep = "\\.") %>% 
  select(groups_a, groups_b, percentage, qual, quantile_mc_nr) %>% 
  mutate(across(starts_with("groups"), ~factor(.x)
  )
  )

# filter out just 95% 
over_cred_df_95 <- over_cred_df %>% 
  filter(percentage == "95%")


# ----- plot 95 % niche size  ----- 
p9 <- ggplot(data = over_stat_df_95, aes(x = mc_nr_perc)) + 
  geom_histogram(bins = 50, colour = "white", aes(fill = groups_a)) + 
  geom_vline(data = over_mean_df_95, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_cred_df_95, aes(xintercept = quantile_mc_nr), 
             colour = "black", linewidth = 1, linetype = 6) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "groups") + 
  ggh4x::facet_grid2(groups_a ~ groups_b, 
                     independent = "y",
                     scales = "free_y") + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(colour = "black"), 
    legend.position = c(0.75, 0.24)
  ) + 
  
  labs(x = paste("Overlap Probability (%)", "\u2013", 
                 "Niche Region Size: 95%"), 
       y = "Frequency")

p9

ggsave(filename = here("Plots",
                       "nicheROVER plots",
                       "niche_percent_overlap_hist_cluster.png"), plot = p9,
       height = 8.5, width = 11)
p10 <- ggplot(data = over_stat_df_95, aes(x = mc_nr_perc)) + 
  geom_density(aes(fill = groups_a)) + 
  geom_vline(data = over_mean_df_95, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_cred_df_95, aes(xintercept = quantile_mc_nr), 
             colour = "black", linewidth = 1, linetype = 6) +
  scale_fill_viridis_d(begin = 0.25, end = 0.85, 
                       option = "D", name = "groups", 
                       alpha = 0.5) + 
  ggh4x::facet_grid2(groups_a ~ groups_b, 
                     independent = "y",
                     scales = "free_y") + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(colour = "black"), 
    legend.position = c(0.76, 0.24)
  ) + 
  
  labs(x = paste("Overlap Probability (%)", "\u2013", 
                 "Niche Region Size: 95%"), 
       y = "Frequency")

p10

ggsave(filename = here("Plots",
                       "nicheROVER plots",
                       "niche_percent_overlap_density_cluster.png"), plot = p10,
       height = 8.5, width = 11)
# ---- determine the size of the niche for each posterior sample -----
fish_size <- sapply(fish_par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# convert to a dataframe for plotting 
fish_size_df <- fish_size %>% 
  as_tibble() %>% 
  mutate(
    id = 1:nrow(.)
  ) %>% 
  pivot_longer(
    cols = -id, 
    names_to = "groups", 
    values_to = "niche_size"
  ) %>% 
  mutate(
    id = 1:nrow(.), 
    groups = factor(groups)
  )

fish_size_df

# calculate mean niche size and the standard deviation of the mean 

niche_size_mean <- rbind(est = colMeans(fish_size),
                         se = apply(fish_size, 2, 
                                    function(x) sd(x) / sqrt(length(x)))) %>% 
  as_tibble(rownames = "metric") %>% 
  pivot_longer(cols = -metric, 
               names_to = "groups", 
               values_to = "values") %>% 
  mutate(
    groups = factor(groups)
  ) %>% 
  pivot_wider(id_cols = groups, names_from = metric, values_from = values)

# openxlsx::write.xlsx(x = niche_size_mean, here("results", 
#                                                "nicheROVER_niche_size_mean.xlsx"))

# ---- plot niche size with mean and sd ----- 

ggplot(data = fish_size_df) + 
  geom_violin(
    aes(x = groups, y = niche_size),
    width = 0.2
  ) + 
  geom_point(data = niche_size_mean, aes(x = groups, y = est), 
             # size = 3
  ) +
  geom_errorbar(data = niche_size_mean, aes(x = groups, 
                                            ymin = est - se, 
                                            ymax = est + se), 
                width = 0.03) +
  scale_y_continuous(breaks = seq(10, 80, 10)) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "groups", 
       y = "Niche Size") -> p11

ggsave(filename = here("Plots",
                       "niche size",
                       "niche_size_viloin_cluter.png"), plot = p11,
       height = 8.5, width = 11)



ggplot(data = fish_size_df) + 
  geom_boxplot(
    aes(x = groups, y = niche_size),
    width = 0.05
  ) + 
  geom_point(data = niche_size_mean, aes(x = groups, y = est), 
             # size = 2.5
  ) +
  geom_errorbar(data = niche_size_mean, aes(x = groups, 
                                            ymin = est - se, 
                                            ymax = est + se), 
                width = 0.03) +
  scale_y_continuous(breaks = seq(10, 80, 10)) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "groups", 
       y = "Niche Size") -> p12

ggsave(filename = here("Plots",
                       "niche size",
                       "niche_size_boxplot_cluster.png"), plot = p12,
       height = 8.5, width = 11)


# ---- look at distribution and test differences in niche size ---- 
ggplot(data = fish_size_df, aes(x = niche_size)) + 
  geom_histogram()


fitdistrplus::descdist(fish_size_df$niche_size)


gamma_size <- fitdistrplus::fitdist(fish_size_df$niche_size, distr = "gamma", 
                                    method = "mme") 

lnorm_size <- fitdistrplus::fitdist(fish_size_df$niche_size, distr = "lnorm", 
                                    method = "mme") 

plot(gamma_size)

car::leveneTest(data = fish_size_df, niche_size ~ groups)
tapply(fish_size_df$niche_size, fish_size_df$groups, 
       shapiro.test)


par(mfrow = c(2, 2))




kruskal.test(data = fish_size_df, niche_size ~ groups)

dunn.test::dunn.test(x = fish_size_df$niche_size, g = fish_size_df$groups)
