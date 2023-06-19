# ---- bring in packages -----
library(dplyr)
library(ellipse)
library(forcats)
library(ggplot2)
library(ggtext)
library(glue)
library(here)
library(lemon)
library(nicheROVER) # possible look at using this vs SIBER
library(purrr)
library(patchwork)
library(readr)
library(tidyr)

# ---- bring in data ----


df <- fish %>% 
  janitor::clean_names()

# ---- calculated the isotopic means for each species -----

aggregate(df[2:4], df[1], mean, na.rm = TRUE)

# remove the samples that didn't run 
df <- df %>% 
  filter(d13c != is.na(d13c)) 

# ---- generate parameter draws from the "default" posteriors of each fish ----
nsample <- 1e3
fish_par <- df %>% 
  split(.$species) %>% 
  map(~ select(., d15n, d13c, d34s)) %>%  
  map(~niw.post(nsample = nsample, X = .))


# ---- separate niw.post object into sigma and mu tibbles ----- 
df_sigma <- map(fish_par, pluck, 2) %>% 
  imap(~ as_tibble(.x) %>% 
         mutate( 
           metric = "sigma", 
           id = c("d15n", "d13c", "d34s"),
           species = .y
         )
  ) %>%
  bind_rows() %>% 
  pivot_longer(cols = -c("id", "species", "metric"),
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
#   filter(id %in% "d15n" & isotopes %in% "d15n")

# ---- mu -----

df_mu <- map(fish_par, pluck, 1) %>% 
  imap(~ as_tibble(.x) %>% 
         mutate( 
           metric = "mu", 
           species = .y
         )
  ) %>%
  bind_rows() %>% 
  mutate(
    species = factor(species, 
                     levels = c("ARCS", "BDWF", "LKWF", "LSCS"))
  ) %>% 
  group_by(species) %>% 
  mutate(
    sample_number = 1:1000
  ) %>% 
  ungroup()
unique(df_mu$species)


glimpse(df_mu)

df_mu_long <- df_mu %>% 
  pivot_longer(cols = -c(metric, species, sample_number), 
               names_to = "isotope", 
               values_to = "mu_est") %>% 
  mutate(
    element = case_when(
      isotope == "d15n" ~ "N",
      isotope == "d13c" ~ "C",
      isotope == "d34s" ~ "S",
    ), 
    neutron = case_when(
      isotope == "d15n" ~ 15,
      isotope == "d13c" ~ 13,
      isotope == "d34s" ~ 34,
    ) 
  )

df_sigma <- df_sigma %>%
  mutate(
    
    )
  )

# ---- density plots ggplot ----



# ggplot(data = df_mu_long, aes(x = mu_est)) +
#   geom_density(aes(fill = species), alpha = 0.5) +
#   scale_fill_viridis_d(begin = 0.25, end = 0.75, 
#                        option = "D", name = "Species") + 
#   facet_grid(y_label ~ isotope_label, 
#              # repeat.tick.labels = TRUE, 
#              scales = "free", 
#              # strip.position = "bottom",
#              labeller = label_parsed) +
#   
#   theme_bw(base_size = 15) + 
#   theme(
#     panel.grid = element_blank(), 
#     strip.background = element_blank(), 
#     strip.text = element_text(size = 17),
#     strip.placement = "outside") + 
#   labs(
#     x = "",
#     y = expression(paste("p(", mu[~delta]," "[isotope], " ", "| X)")))
# )




plots <- df_mu_long %>% 
  split(.$isotope) %>% 
  imap(
    ~ ggplot(data = ., aes(x = mu_est)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                           option = "D", name = "Species") + 
      theme_bw(base_size = 15) + 
      theme(panel.grid = element_blank(), 
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(), 
            legend.position = "none"
      ) + 
      labs(
        x = paste("\u00b5<sub>\U03B4</sub>", "<sub><sup>", 
                  unique(.$neutron), "</sup></sub>", 
                  "<sub>",unique(.$element), "</sub>", sep = ""),
        y = paste0("p(\u00b5 <sub>\U03B4</sub>","<sub><sup>", 
                   unique(.$neutron), "</sub></sup>", 
                   "<sub>",unique(.$element),"</sub>",
                   " | X)"), sep = "")
  )




p <- plots$d15n + 
  theme(legend.position = c(0.10, 0.88)) 

p1 <- plots$d13c
p2 <- plots$d34s


p3 <- p + p1 + p2 

p3


plots <- df_sigma %>% 
  group_split(id, isotopes) %>% 
  imap(
    ~ ggplot(data = ., aes(x = post_sample)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75, 
                           option = "D", name = "Species") +
      theme_bw(base_size = 15) + 
      theme(panel.grid = element_blank()) +
      labs( 
      )
  )

plots[[1]] +
  labs(
    y = "Density",
    x = expression(paste(Sigma[~delta^15],""[N])),
  )
# ---- niche.plot ------
# 
# clrs <- c("black", "blue", "orange")
# 
# df_dat <- as.data.frame(df)
# fish_data <- tapply(1:nrow(df_dat), df$species, function(ii) X = df_dat[ii,4:5]) 
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


species_name <- unique(df_sigma_wide$species)

all_ellipses <- list()

# for loop over speciess and sample number within species to create each 
# sample numbers ellipi 
for (i in 1:4) {
  
  sigma_species <- df_sigma_wide %>% 
    filter(species %in% species_name[i])
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      select(d15n, d13c, d34s)
    
    Sigma <- as.matrix(sigma_ind, 3, 3)
    row.names(Sigma) <- c("d15n", "d13c", "d34s")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      select(sample_number, d15n, d13c, d34s) %>% 
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

# combine ellipose list into dataframe and add species names back in 
ellipse_df <- bind_rows(all_ellipses, .id = "id") %>% 
  mutate(
    species = factor(
      case_when(
        id == "1" ~ "ARCS",
        id == "2" ~ "BDWF",
        id == "3" ~ "LKWF",
        id == "4" ~ "LSCS"
      ), level = c("ARCS", "BDWF", "LKWF", "LSCS")
    )
  ) %>% 
  as_tibble()


# randomly sample 20 ellipise out of 1000 
ellipse_df %>% 
  group_by(species, rep) %>% 
  nest() %>%
  group_by(species) %>% 
  slice_sample(n = 10, replace = TRUE) %>% 
  ungroup() %>% 
  unnest(cols = c(data)) -> random_ellipse 


random_ellipse

# ---- create niche plots -----
ggplot() + 
  geom_polygon(data = random_ellipse,
               mapping = aes(x = d13c, y = d15n,
                             group = interaction(rep, species),
                             color = species,
                             fill = NULL),
               fill = NA,
               linewidth = 0.5) + 
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species", 
                         # alpha = 0.35
  ) + 
  # scale_fill_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", name = "species") + 
  scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  scale_y_continuous(breaks = seq(6, 16, 2)) +
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
  geom_density(data = df, aes(x = d15n, 
                              fill = species,
                              # colour = species
  ), 
  alpha = 0.35, 
  linewidth = 0.8) +
  # scale_colour_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", name = "species") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
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
  geom_density(data = df, aes(x = d13c, 
                              fill = species
                              # colour = species
  ), 
  alpha = 0.35, 
  linewidth = 0.8) +
  # scale_colour_viridis_d(begin = 0.25, end = 0.85, 
  #                      option = "D", name = "species") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
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
  geom_point(data = df, aes(x = d13c, y = d15n,
                            # colour = species, 
                            fill = species, 
  ),
  shape = 21, colour = "black", 
  stroke = 0.8,
  size = 3, alpha = 0.70) +
  # scale_colour_viridis_d(begin = 0.25, end = 0.75, 
  #                        option = "D", name = "species") + 
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
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

# p7

p8 <- p5 + p4 + p7 + p6

p8





# sigma_list <- df_sigma_wid %>% 
#   group_split(species) %>%  
#   map(~ split(., as.factor(.$sample_number)))
# 
# sigma_list_n <- df_sigma_wid %>%
#   select(-isotopes) %>% 
#   group_by(species, sample_number) %>%
#   nest() %>% 
#   ungroup() %>% 
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
# Calculate how much each grouping overlaps 
over_stat <- overlap(fish_par, nreps = nsample, nprob = 1e3, 
                     alpha = c(0.95, 0.99))

# convert to dataframe to be able to plot and look at staistically 
over_stat_df <- over_stat %>% 
  as_tibble(rownames = "species_a") %>% 
  mutate(
    id = 1:nrow(.), 
    species_a = factor(species_a, 
                       level = c("East", 
                                 "West",
                                 "North"))
  ) %>% 
  pivot_longer(cols = -c(id, species_a), 
               names_to = "species_b", 
               values_to = "mc_nr")  %>% 
  separate(species_b, into = c("species_c", "sample_number", "percentage"), 
           sep = "\\.") %>% 
  select(-id) %>% 
  rename(species_b = species_c) %>% 
  mutate(
    species_b =  factor(species_b, 
                        level = c("East", 
                                  "West",
                                  "North")), 
    mc_nr_perc = mc_nr * 100
  )

# filter for just the 95% 
over_stat_df_95 <- over_stat_df %>% 
  filter(percentage %in% "95%")

# calculate the mean overlap for each species 
over_mean <- apply(over_stat, c(1:2, 4), mean) * 100 %>% 
  round(2) 

# convert to dataframe to be able to plot and compare 
over_mean_df <- over_mean %>% 
  as_tibble(rownames = "species_a") %>% 
  mutate(
    id = 1:nrow(.), 
    species_a = factor(species_a, 
                       level = c("East", 
                                 "West",
                                 "North"))
  ) %>% 
  pivot_longer(cols = -c(id, species_a), 
               names_to = "species_b", 
               values_to = "mean_mc_nr")  %>% 
  separate(species_b, into = c("species_c", "percentage"), 
           sep = "\\.") %>% 
  select(-id) %>% 
  rename(species_b = species_c) %>% 
  mutate(
    species_b =  factor(species_b, 
                        level = c("East", 
                                  "West",
                                  "North")))

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
               names_to = "species_b", 
               values_to = "quantile_mc_nr") %>% 
  separate(species_b, into = c("species_a", "species_b", "percentage"), 
           sep = "\\.") %>% 
  select(species_a, species_b, percentage, qual, quantile_mc_nr) %>% 
  mutate(across(starts_with("species"), ~factor(.x, 
                                                level = c("East", 
                                                          "West",
                                                          "North"))
  )
  )
over_cred_df
# filter out just 95% 
over_cred_df_95 <- over_cred_df %>% 
  filter(percentage == "95%")


# ---- plot 95 % niche size  ----- 
p9 <- ggplot(data = over_stat_df_95, aes(x = mc_nr_perc)) + 
  geom_histogram(bins = 50, colour = "white", aes(fill = species_a)) + 
  geom_vline(data = over_mean_df_95, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_cred_df_95, aes(xintercept = quantile_mc_nr), 
             colour = "black", linewidth = 1, linetype = 6) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       alpha = 0.8,
                       option = "D", name = "species") + 
  ggh4x::facet_grid2(species_a ~ species_b, 
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
                       "niche_percent_overlap_hist.png"), plot = p9,
       height = 8.5, width = 11)

p10 <- ggplot(data = over_stat_df_95, aes(x = mc_nr_perc)) + 
  geom_density(aes(fill = species_a)) + 
  geom_vline(data = over_mean_df_95, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_cred_df_95, aes(xintercept = quantile_mc_nr), 
             colour = "black", linewidth = 1, linetype = 6) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species", 
                       alpha = 0.35) + 
  ggh4x::facet_grid2(species_a ~ species_b, 
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

# p10

ggsave(filename = here("Plots",
                       "nicheROVER plots", 
                       "niche_percent_overlap_density.png"), plot = p10,
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
    names_to = "species", 
    values_to = "niche_size"
  ) %>% 
  mutate(
    id = 1:nrow(.), 
    species = factor(species, 
                     level = c("East", 
                               "West",
                               "North"))
  )

fish_size_df

# calculate mean niche size and the standard deviation of the mean 

niche_size_mean <- rbind(est = colMeans(fish_size),
                         sd = apply(fish_size, 2, 
                                    function(x) sd(x)),
                         se = apply(fish_size, 2, 
                                    function(x) sd(x) / sqrt(length(x)))) %>% 
  as_tibble(rownames = "metric") %>% 
  pivot_longer(cols = -metric, 
               names_to = "species", 
               values_to = "values") %>% 
  mutate(
    species = factor(species, 
                     level = c("East", 
                               "West",
                               "North"))
  ) %>% 
  pivot_wider(id_cols = species, names_from = metric, values_from = values)

niche_size_mean
# ---- plot niche size with mean and sd ----- 

ggplot(data = fish_size_df) + 
  geom_violin(
    aes(x = species, y = niche_size),
    width = 0.2
  ) + 
  geom_point(data = niche_size_mean, aes(x = species, y = est), 
             # size = 3
  ) +
  geom_errorbar(data = niche_size_mean, aes(x = species, 
                                            ymin = est - se, 
                                            ymax = est + se), 
                width = 0.03) +
  scale_y_continuous(breaks = seq(10, 80, 10)) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "species", 
       y = "Niche Size") -> p11

ggsave(filename = here("Plots",
                       "niche size",
                       "niche_size_viloin.png"), plot = p11,
       height = 8.5, width = 11)



ggplot(data = fish_size_df) + 
  geom_boxplot(
    aes(x = species, y = niche_size),
    width = 0.05
  ) + 
  geom_point(data = niche_size_mean, aes(x = species, y = est), 
             # size = 2.5
  ) +
  geom_errorbar(data = niche_size_mean, aes(x = species, 
                                            ymin = est - se, 
                                            ymax = est + se), 
                width = 0.03) +
  scale_y_continuous(breaks = seq(10, 80, 10)) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "species", 
       y = "Niche Size") -> p12

ggsave(filename = here("Plots",
                       "niche size",
                       "niche_size_boxplot.png"), plot = p12,
       height = 8.5, width = 11)


# ---- look at distribution and test differences in niche size ---- 
ggplot(data = fish_size_df, aes(x = niche_size)) + 
  geom_histogram()


fitdistrplus::descdist(fish_size_df$niche_size)
norm_size <- fitdistrplus::fitdist(fish_size_df$niche_size, distr = "norm") 


gamma_size <- fitdistrplus::fitdist(fish_size_df$niche_size, distr = "gamma", 
                                    method = "mme") 

lnorm_size <- fitdistrplus::fitdist(fish_size_df$niche_size, distr = "lnorm", 
                                    method = "mme") 

plot(norm_size)
plot(gamma_size)
plot(lnorm_size)

car::leveneTest(data = fish_size_df, niche_size ~ species)
tapply(fish_size_df$niche_size, fish_size_df$species, 
       shapiro.test)
p12

par(mfrow = c(2, 2))

m <- lm(data = fish_size_df, niche_size ~ species)

plot(m)

car::Anova(m)
summary(m)


kruskal.test(data = fish_size_df, niche_size ~ species)

dunn.test::dunn.test(x = fish_size_df$niche_size, g = fish_size_df$species)
