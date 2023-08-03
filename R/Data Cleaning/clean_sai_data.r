# bring in packages -----
library(dplyr)
library(ggplot2)
library(here)
library(nicheROVER)
library(SIBER)
library(SIBERG)
library(readr)
library(stringr)
library(tidyr)


# bring in test dataframe ------


df <- read_csv(here("Data", 
                    "sai_lt_test.csv")) %>% 
  janitor::clean_names()

glimpse(df)

df$sample <- as.character(df$sample)
# bring in fish metaddata 

fish_meta <- read_csv(here("Data", 
                           "sai_lt_metadata.csv")) %>% 
  janitor::clean_names()

glimpse(fish_meta)


# calcuate wieghts -----

fish_meta <- fish_meta %>% 
  mutate(wt_g = if_else(is.na(wt_g), true = 
                          round((10 ^ (-5.218126 +  3.038077 * log10(tl_mm))), 
                                digits = 0),
                        false = wt_g)
  )

View(fish_meta)

lt_meta_f <- fish_meta %>% 
  filter(vemco_tag != is.na(vemco_tag)) %>% 
  group_by(year, vemco_tag, vemco_type, min_delay, max_delay) %>% 
  summarise(
    n = n(), 
    tl = mean(tl_mm), 
    tl_sem = sd(tl_mm) / sqrt(n()), 
    wt = mean(wt_g), 
    wt_sem = sd(wt_g) / sqrt(n())
  ) %>% 
  ungroup()

lt_meta_f %>% 
  openxlsx::write.xlsx(., here("Results", 
                               "table_1_sai_a.xlsx"))





fish_meta <- fish_meta %>% 
  mutate(basin = str_remove(basin, pattern = " Basin") %>% 
           factor(levels = c("East", "West", "North"))
  )

glimpse(fish_meta)

# View(fish_meta)

# join test dataframe with fish_metadata 

df <- df %>% 
  left_join(fish_meta, by = c("sample" = "tag_id"))




# View(df)
glimpse(df)

write_rds(x = df, file = here("Saved Data", 
                              "cleaned_lkt_tagged_sia.rds"))
