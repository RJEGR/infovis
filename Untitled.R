# Nanodrop

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/nanodrop/'
pattern_f <- 'csv'
  
file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

df <- lapply(file, read.csv)

head(df <- do.call(rbind, df))

m <- round(min(df$Total_RNA_ug))-1
M <- round(max(df$Total_RNA_ug)) + m

getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

df %>%
  ggplot() +
  geom_col(aes(x = id, y = Total_RNA_ug, fill = as.factor(pH))) +
  # geom_col(aes(x = id, y = Total_RNA_ug), fill = 'white') +
  facet_grid(~ hpf,scales = 'free_x') +
  scale_y_continuous(expression("Total RNA" ~ (µg)), breaks = seq(m,M, by = 15)) +
  scale_x_discrete("",labels = df$P)+
  scale_fill_manual(values = getPalette[-3]) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) +
  guides(fill = guide_legend(""))

# test ggdensity

df %>%
  ggplot(aes(A260_A280)) +
  geom_histogram()

    