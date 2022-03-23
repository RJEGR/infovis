rm(list = ls())
# 
# test https://www.danieldsjoberg.com/gtsummary/index.html

options(stringsAsFactors = FALSE)

library(tidyverse)

Lcanal <- c('Canal-1', 'Canal-3', 'Canal-2', 'Canal-4')

# RColorBrewer::display.brewer.pal(12, 'Paired')

getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

getPalette <- structure(getPalette, names = Lcanal) 


path <- '~/Documents/DOCTORADO/pH_measures/2022_/asentamiento'
# 

file <- 'asentamiento_canal_1_2_agua_marina.dat$'

file <- list.files(path, pattern = file, full.names = TRUE)

df <- lapply(file, read.table)

head(df1 <- do.call(rbind, df))

namesL <- c('Canal-1', 'Canal-2', 'Canal-3', 'Canal-4')


head(times <- df1[1:4])
head(mVolts <- df1[5:8])
head(tbl1 <- df1[10:13])

names(tbl1) <- namesL


tbl <- tbl1

names(mVolts) <- namesL

head(dff <- cbind(times, tbl))

dff %>% distinct(V1)

x_breaks <- nrow(dff)*5 / 3600

dff %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = all_of(namesL), values_to = 'Obs') %>%
  drop_na(Obs) %>%
  filter(name %in% c('Canal-1', 'Canal-2')) -> df_longer

saveRDS(df_longer, file = paste0(path, '/df_longer_control.rds'))

df_longer %>%
  group_by(name, V1) %>%
  mutate(z = (Obs - mean(Obs)) / sd(Obs),
    Obs = ifelse(abs(z) >= 2, NA, Obs)) %>% 
  mutate(name = factor(name, levels = namesL)) %>%
  mutate(Day = substr(V1, 1,2), id = (id*5)/ 3600) %>%
  mutate(Day = as.numeric(Day)) %>%
  drop_na(Obs) %>%
  # sample_n(1000) %>%
  ungroup() %>%
  # mutate(name = ifelse(name %in% c('Canal-3')), NA, name)
  ggplot(aes(id, Obs, color = name, group = name)) +
  # geom_line() +
  # Span Controls the amount of smoothing for the default loess smoother.
  # geom_smooth(span = 0.7,
  # linetype="dashed", size = 0.5, alpha=0.5,
  # se = TRUE, na.rm = TRUE) +
  geom_point(alpha = 0.5, size = 0.1) +
  scale_color_manual(values = getPalette) +
  scale_y_continuous(breaks = seq(7,8.15, by = 0.1)) +
  scale_x_continuous(breaks = seq(0, x_breaks, by = 12)) +
  labs(y = 'pH', x = 'Hour') +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank()) +
  guides(colour = guide_legend("")) -> psave

psave

library(rstatix)

IC <- function(x, conf = 0.99) {
  a <- mean(x)
  sd <- sd(x)
  n <- length(x)
  err <- qnorm(p = conf, mean = a, sd = sd, lower.tail = F) * sd / sqrt(n)
  
  return(err)
}

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

source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

df_longer %>%
  drop_na(Obs) %>%
  group_by(V1, name) %>%
  mutate(value = Obs) %>%
  summarise(
    a = mean(value), sd = sd(value), IC = IC(value),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) -> out_stats

out_stats %>%
  filter(name %in% 'Canal-2') %>%
  ggplot(aes(y = a, x = V1, fill = name)) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  # facet_grid(.~Buffer, scales = 'free_y') +
  geom_point() +
  geom_line(aes(group = name), orientation = "x") +
  # geom_bar(stat="identity", fill = 'white', color = 'grey67') + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  theme(legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1))


