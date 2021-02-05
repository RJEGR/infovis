rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)
library(tidyverse)
library(ggh4x)


obj <- read_xlsx(files[1])
mtd_file <- list.files(path = dir, pattern = "CIAD_MappingFile.txt$", full.names = TRUE)

mtd <- read.csv(mtd_file, sep = "\t") %>% 
  select(Sample_IDs, Tissue, Time) %>%
  rename(Index = Sample_IDs)


obj %>%
  select_if(is.double) %>%
  names() -> colNames

obj %>%
  pivot_longer(cols = all_of(colNames), 
               values_to = "RA", 
               names_to = "Index") %>% 
  inner_join(mtd) -> obj_longer

# arrange by factor


obj_longer %>%
  arrange(desc(Category)) %>%
  mutate(SuperPathway = factor(SuperPathway, 
                               levels = unique(SuperPathway))) -> dataViz

# view(dataViz)

# Which noised data we've? (look at x-axis)

dataViz %>%
  group_by(Index) %>%
  mutate(Z = RA - mean(RA) / sd(RA)) %>%
  ggplot(aes(RA, fill = cut(RA, 100))) +
  geom_histogram(show.legend = F, bins = 100) +
  facet_wrap(~ Tissue)


obj_longer %>%
  # filter(RA > 12000) %>%
  group_by(SuperPathway, Tissue, Category) %>%
  summarise(RA = sum(RA)) %>%
  group_by(Tissue) %>%
  mutate(Rank = rank(RA)) -> dataViz

LCategory <- dataViz %>% group_by(Category) %>% summarise(t = sum(RA)) %>% arrange(desc(t)) %>% pull(Category)

library(ggsci)

dataViz %>%
  mutate(Category = factor(Category, 
                               levels = LCategory)) %>%
  mutate(SuperPathway = forcats::fct_reorder(SuperPathway, Rank)) %>%
  ggplot() +
  geom_col(aes(x = SuperPathway ,
               y = RA, fill = Category)) +
  labs(y = "Relative Abundance of predicted gene (%)") +
  ggsci::scale_fill_ucscgb() +
  facet_grid(Category ~ Tissue , scales = "free", space = "free_y") +
  coord_flip() +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank()) -> barPlot

# barPlot

ggsave(barPlot, filename = "barKegg.png", path = dir, 
       width = 18, height = 10)
