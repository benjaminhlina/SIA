# ---- bring in packages -----

library(ggplot2)
library(here)
library(patchwork)
library(readr)

# --- bring in plots to arrange ----
p <- read_rds(here("Saved Plots", 
                   "c_13_distance_traveled_overall.rds")) 
p1 <- read_rds(here("Saved Plots", 
                    "n_15_distance_traveled_overall.rds"))
p2 <- read_rds(here("Saved Plots",
                    "c_13_depth_overall.rds"))
p3 <- read_rds(here("Saved Plots",
                    "n_15_depth_overall.rds"))
p4 <- read_rds(here("Saved Plots",
                    "c_13_temp_overall.rds"))
p5 <- read_rds(here("Saved Plots",
                    "n_15_temp_overall.rds"))
# ---- remove lengends -----
p <- p + 
  theme(
    legend.position = "none", 
    plot.tag.position  = c(0.20, 0.96)
  )
p1 <- p1 + 
  theme(
    legend.position = "none",
    plot.tag.position  = c(0.20, 0.96)
  )
p2 <- p2 + 
  theme(
    legend.position = "none",
    plot.tag.position  = c(0.21, 0.96)
  )
p3 <- p3 + 
  theme(
    legend.position = "none",
    plot.tag.position  = c(0.21, 0.96)
  )
p4 <- p4 + 
  theme(
    legend.position = "none",
    plot.tag.position  = c(0.20, 0.96)
  )
p5 <- p5 + 
  theme(
    # legend.position = "none"
    # legend.position = c(0.85, 0.21)
    legend.position = c(0.15, 0.75), 
    plot.tag.position  = c(0.20, 0.96),
    legend.background = element_blank()
    
    )
# ---- combine all plots ----- 
p6 <- p + p2 + p4 + p1 + p3 + p5 + 
  plot_annotation(tag_levels = "a") 
  

# p6
# ---- save -----
ggsave(filename = here("Plots", 
                       "Habitat Use and Isotopes", 
                       "dis_temp_depth_isotope.png"), plot = p6,
       height = 8.5, width = 11 * 1.25)

