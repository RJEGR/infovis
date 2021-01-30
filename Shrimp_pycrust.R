rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)

obj <- read_xlsx(files[1])


library(tidyverse)

mtd <- read.csv(list.files(path = dir, pattern = "txt$", full.names = TRUE), sep = "\t") %>% select(Sample_IDs, Tissue, Time) %>%
  rename(Index = Sample_IDs)


obj %>%
  select_if(is.double) %>%
  names() -> colNames

obj %>%
  pivot_longer(cols = colNames, values_to = "RA", 
               names_to = "Index") %>% 
  inner_join(mtd) -> obj_longer

# arrange by factor

reorder_within <- function(x, by, within, fun = mean, sep = "_",...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = mean)
}

scale_x_reordered <- function(..., sep = "_") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

obj_longer %>%
  arrange(desc(Category)) %>%
  mutate(SuperPathway = factor(SuperPathway, 
                               levels = unique(SuperPathway))) -> dataViz

# view(dataViz)

# Which noised data we've? (look at x-axis)

dataViz %>%
  group_by(Index) %>%
  mutate(Z = RA - mean(RA) / sd(RA)) %>%
  # filter(abs(Z) < 3) %>%
  ggplot(aes(RA, fill = cut(RA, 100))) +
  geom_histogram(show.legend = F, bins = 100) +
  facet_wrap(~ Tissue)


dataViz %>%
  ungroup() %>%
  filter(RA > 12000) %>%
  # ggplot(aes(x = reorder_within(SuperPathway, RA, Tissue) , 
  #            y = RA, fill = Category)) +
  ggplot() +
  geom_col(aes(x = reorder_within(SuperPathway, RA, Tissue) , 
               y = RA, fill = Category)) +
  labs(y = "Relative Abundance of predicted gene (%)") +
  coord_flip() +
  facet_grid(~ Tissue , scales = "free", space = "free_x") +
  theme_bw(base_size = 16, base_family = "GillSans") -> barPlot

barPlot

ggsave(barPlot, filename = "barKegg.png", path = dir, 
       width = 18, height = 10)





