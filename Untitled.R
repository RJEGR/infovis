# Nanodrop

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/nanodrop/'

pattern_f <- 'csv'


qqfun <- function(x) {
  x <- x[x > 0]
  qq <- qqnorm(x, plot.it = F) 
  qq %>% as_tibble()
}

zfun <- function(x) {
  x <- x[x > 0]
  z <- c((x - mean(x)) / sd(x))
  # z %>% as_tibble()
  return(z)
}


file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

df <- lapply(file, read.csv)

head(df <- do.call(rbind, df))

m <- round(min(df$Total_RNA_ug))-1
M <- round(max(df$Total_RNA_ug)) + m

getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

my_theme <-
  ggplot2::theme_bw(base_size = 10, base_family = "GillSans") +
  ggplot2::theme(
    legend.position = 'top',
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank())

theme_set(my_theme)

df %>%
  ggplot() +
  geom_col(aes(x = id, y = Total_RNA_ug, fill = as.factor(pH))) +
  facet_grid(~ hpf,scales = 'free_x') +
  scale_y_continuous(expression("Total RNA" ~ (µg)), breaks = seq(m,M, by = 15)) +
  scale_x_discrete("",labels = df$P)+
  scale_fill_manual("", values = getPalette[-3]) -> p

df %>%
  ggplot() +
  geom_col(aes(x = id, y = Total_RNA_ug, fill = as.factor(pH))) +
  facet_grid(~ rspnW_uL,scales = 'free_x') +
  scale_y_continuous(expression("Total RNA" ~ (µg)), breaks = seq(m,M, by = 15)) +
  scale_x_discrete("",labels = df$P)+
  scale_fill_manual("", values = getPalette[-3])

df %>%
  mutate(A280 = A260/A260_A280) %>%
  mutate(col = ifelse(A260_A280 > 5, 'red', 'black')) %>%
  ggplot(aes(A280, A260)) +
  geom_smooth(method = "lm", linetype="dashed", 
    size = 0.5, alpha=0.5, se = F) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001, label.y = 70) +
  geom_point(aes(color = col), alpha = 0.7) +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = expression(A[280]), y = expression(A[260])) +
  theme_bw() +
  theme(legend.position = 'none') -> p1

df %>%
  mutate(A260_A280 = ifelse(A260_A280 > 5, NA, A260_A280)) %>%
  ggplot() +
  geom_col(aes(x = id, y = A260_A280, fill = as.factor(pH))) +
  facet_grid(~ hpf,scales = 'free_x') +
  labs(y = expression(A[260]/A[280])) +
  scale_x_discrete("",labels = df$P) +
  scale_fill_manual("", values = getPalette[-3]) -> p2

library(patchwork)
p1+p2+plot_layout(widths = c(1,4))
p/p2

df %>%
  mutate(A280 = A260/A260_A280) %>%
  mutate(col = ifelse(A260_A280 > 5, 'red', 'black')) %>%
  ggplot(aes(A260_A280, Total_RNA_ug)) +
  facet_grid(pH ~ hpf,scales = 'free') +
  theme_bw() +
  geom_point()

# test ggdensity
# Not good enough <---
# remotes::install_github("jamesotto852/ggdensity")

library("ggdensity")

df %>%
  mutate(A280 = A260/A260_A280) %>%
  ggplot(aes(A280, A260)) + 
  # coord_equal() +
  # geom_density_2d_filled() +
  geom_hdr_lines() +
  geom_point(shape = 21) +
  theme_bw() +
  labs(x = expression(A[280]), y = expression(A[260])) 
  # facet_grid(hpf ~ pH, scales = 'free_y') 

# as boxplot ----

df %>%
  ggplot() +
  geom_boxplot(aes(x = pH, 
    y = Total_RNA_ug, fill = as.factor(pH))) +
  facet_grid(~ hpf,scales = 'free_x') +
  scale_y_continuous(expression("Total RNA" ~ (µg)), breaks = seq(m,M, by = 15)) +
  scale_x_discrete("",labels = df$P)+
  scale_fill_manual("", values = getPalette[-3]) +
  theme_bw()
