rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

Lcanal <- c('Canal-1', 'Canal-3', 'Canal-2', 'Canal-4')

# RColorBrewer::display.brewer.pal(12, 'Paired')

getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

getPalette <- structure(getPalette, names = Lcanal) 

path <- '~/Documents/DOCTORADO/pH_measures/2022_/'

file <- '.dat$'

file <- list.files(path, pattern = file, full.names = TRUE)

r.tb <- function(file) {
  df <- read.table(file)
  samtm <- basename(file)
  samtm <- sapply(strsplit(samtm, "[.]"), `[`, 1)
  
  df <- mutate(df, samtm = samtm)
 
  }

df <- lapply(file, r.tb)

head(df1 <- do.call(rbind, df))

# head(df1 <- rbind(df[[7]])) # , df[[2]]

head(times <- df1[1:4])
head(mVolts <- df1[5:8])
head(tbl1 <- df1[10:13])

namesL <- c('Canal-1', 'Canal-2', 'Canal-3', 'Canal-4')


into <- c("date", "exp", "Stage", 'hpf')

df1 %>% as_tibble() %>%
  separate(col = samtm, 
    into = into, sep = "_") %>%
  select(all_of(c(names(df1[10:13]), into))) -> tbl


names(mVolts) <- namesL


# mVolts %>%
#   pivot_longer(cols = names(tbl)) %>%
#   mutate(g = ifelse(name %in% c('Canal-1', 'Canal-3'), 'pH-7.6', 'pH-7.8')) %>%
#   # filter(name %in% c('Canal-3', 'Canal-4')) %>%
#   ggplot(aes(value, fill = name)) +
#   scale_fill_brewer('') +
#   geom_histogram() + 
#   facet_grid(~ g) +
#   theme_bw(base_size = 14, base_family = "GillSans")
# 

names(tbl)[1:4] <- namesL

tbl %>%
  pivot_longer(cols = namesL) %>%
  mutate(g = ifelse(name %in% c('Canal-1', 'Canal-3'), 
    'pH-7.6', 'pH-7.8')) %>%
  # filter(name %in% c('Canal-3', 'Canal-4')) %>%
  ggplot(aes(value, fill = name, group = name), alpha = 0.7) +
  scale_fill_manual(values = getPalette) +
  geom_histogram(bins = 100) + 
  geom_rug(aes(color = name)) +
  scale_color_manual(values = getPalette) + 
  facet_grid(~ g) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  scale_x_continuous(name = 'pH', 
    breaks = seq(6.5,8.15, by = 0.25), limits = c(6.5,8.15)
  ) +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_blank())


head(dff <- cbind(times, tbl))


library(gtsummary)

dff %>%
  pivot_longer(cols = namesL, values_to = 'Obs') -> dff_longer

dff_longer

dff_longer %>% 
  select(V1, name, Obs) %>% # V3, Stage, hpf,
  filter(name %in% c('Canal-3', 'Canal-4')) %>%
  drop_na(Obs) %>%
  # sample_n(100) %>% 
  tbl_summary(by = name,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)")) %>%
  add_ci() %>%
  add_p() -> md_table


dff %>%
  pivot_longer(cols = namesL, values_to = 'Obs') %>%
  mutate(name = factor(name, levels = namesL)) %>%
  mutate(Exp = ifelse(name %in% c('Canal-1', 'Canal-3'), 7.6, 7.8)) %>%
  mutate(dif = abs(Exp - Obs)) %>%
  mutate(name = factor(name, levels = Lcanal)) %>%
  ggplot(aes(y = Obs, x = V1, color = name)) +
  scale_color_manual('', values = getPalette) +
  stat_boxplot(geom ='errorbar', width = 0.3, 
    position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
    outlier.alpha = 0) +
  stat_summary(aes(group = name), fun=mean, geom="point", shape=23, 
    size=1, position = position_dodge(0.6), color = 'black') +
  scale_y_continuous(name = 'pH', 
    breaks = seq(7.0,8.15, by = 0.1), limits = c(7.0,8.15)
  ) +
  labs(x = 'Fecha') +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank()) +
  guides(colour = guide_legend("")) -> psave


ggsave(psave, filename = paste0(path, "raw_pH_boxplot.png"), width = 7, height = 4)

IC <- function(x, conf = 0.95) {
  a <- mean(x)
  sd <- sd(x)
  # n <- length(x)
  n <- sum(x > 0)
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


dff %>%
  pivot_longer(cols = namesL, values_to = 'Obs') %>%
  mutate(name = factor(name, levels = namesL)) %>%
  drop_na(Obs) %>%
  group_by(name) %>% 
  summarise(min = min(Obs), max = max(Obs), 
    a = mean(Obs), sd = sd(Obs), IC = IC(Obs),
    upper = a+IC, lower = a-IC, n = n(),
    rCal = Rcalculate_coeff(Obs),
    rCri = Rcritical_coeff(n),
    gaussian = ifelse(rCal > rCri, TRUE, FALSE)) %>% view()

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

dff %>%
  pivot_longer(cols = namesL, values_to = 'Obs') %>%
  mutate(name = factor(name, levels = namesL)) %>%
  drop_na(Obs) %>%
  group_by(name) %>%
  summarise(qqfun(Obs)) %>%
  mutate(z = zfun(y)) %>%
  mutate(outlier = ifelse(abs(z)>3, TRUE, FALSE)) -> outliersdf


outliersdf %>%
  ggplot(aes(y, fill = name), alpha = 0.5) +
  geom_histogram(bins = 100) +
  geom_rug(aes(color = name)) +
  labs(y = 'pH') +
  scale_fill_manual('', values = getPalette) +
  scale_color_manual('', values = getPalette) +
  facet_grid( outlier ~., scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 4),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank()) +
  guides(colour = guide_legend("")) -> psave

ggsave(psave, filename = paste0(path, "raw_outliers_hist.png"), width = 7, height = 4)

outliersdf %>%
  # filter(outlier == FALSE) %>%
  ggplot(aes(x, y, group = name, color = outlier)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  geom_point(size = 2.5, alpha = 0.5) + # , label.y = 2.5
  geom_rug(aes(color = outlier)) +
  facet_wrap(~name) +
  ggpubr::stat_cor(aes(group = name), 
    method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "Expected", 
    y = expression("Observed "~ pH),
    caption = "") +
  scale_color_manual(name = expression("Outlier-"~sigma), 
    values = c('black', 'red')) +
  scale_y_continuous(name = 'pH', 
    breaks = seq(6.5,8.15, by = 0.1), limits = c(6.5,8.15)) +
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 4),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank()) +
  guides(colour = guide_legend("")) -> psave

ggsave(psave, filename = paste0(path, "raw_qqplot.png"), width = 6, height = 6)


# mean pH after clean outliers ---

outliersdf %>%
  filter(outlier == FALSE) %>%
  mutate(g = ifelse(name %in% c('Canal-1', 'Canal-3'), 
    'pH-7.6', 'pH-7.8')) %>%
  drop_na(y) %>%
  group_by(name) %>% 
  summarise(min = min(y), max = max(y), 
    a = mean(y), sd = sd(y), IC = IC(y),
    upper = a+sd, lower = a-sd, n = n(),
    rCal = Rcalculate_coeff(y),
    rCri = Rcritical_coeff(n),
    gaussian = ifelse(rCal > rCri, TRUE, FALSE)) %>% view()


# outliersdf %>%
#   filter(outlier == FALSE) %>%
#   mutate(g = ifelse(name %in% c('Canal-1', 'Canal-3'), 
#     'pH-7.6', 'pH-7.8')) %>%
#   ggplot(aes(y = y, x = name, color = g, group = name)) +
#   # facet_grid(outlier ~.) +
#   scale_color_brewer('') +
#   stat_boxplot(geom ='errorbar', width = 0.3, 
#     position = position_dodge(0.6)) +
#   geom_boxplot(width = 0.3, position = position_dodge(0.6), 
#     outlier.alpha = 0) +
#   stat_summary(fun=mean, geom="point", shape=23, 
#     size=1, position = position_dodge(0.6)) +
#   scale_y_continuous(name = 'pH', 
#     breaks = seq(6.5,8.15, by = 0.25), limits = c(6.5,8.15)) +
#   scale_x_discrete(limits = Lcanal) +
#   labs(x = 'Fecha') +
#   theme_classic(base_size = 12, base_family = "GillSans")

outliersdf %>%
  filter(outlier == FALSE) %>%
  mutate(g = ifelse(name %in% c('Canal-1', 'Canal-3'), 
    'pH-7.6', 'pH-7.8')) %>%
  filter(name %in% c('Canal-3', 'Canal-4')) %>%
  ggplot(aes(y, fill = name, group = name), alpha = 0.7) +
  geom_histogram(bins = 100) + 
  geom_rug(aes(color = name)) +
  # scale_fill_brewer('') +
  # scale_color_brewer('') + 
  scale_fill_manual('', values = getPalette) +
  scale_color_manual('', values = getPalette) +
  # facet_grid(~ g) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  scale_x_continuous(name = 'pH', 
    breaks = seq(6.5,8.15, by = 0.25), limits = c(6.5,8.15)
  ) +
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave

ggsave(psave, filename = paste0(path, "clean_hist.png"), width = 5, height = 3)

# 
# 
# 
outliersdf %>%
  filter(outlier == FALSE) %>%
  drop_na(y) %>%
  group_by(name) %>%
  mutate(n = n())


outliersdf %>% 
  ungroup() %>%
  mutate(id = rep(1:nrow(dff), 4)) %>%
  select(id, names(outliersdf)) %>%
  filter(outlier == F) %>% 
  pivot_wider(id_cols = id, names_from = name, 
    values_from = y, values_fill = NA) -> clean_dff
  
x_breaks <- nrow(dff)*5 / 3600

# library(ggsci)

# gcolor <- structure(c('#a50f15', '#de2d26'), 
#   names = c('Canal-1', 'Canal-2'))

# Redoing plots w/ clean outliers ----

dff %>% mutate(id = 1:nrow(.), hpf = (id*5)/3600) %>% 
  mutate(g = rep(1:(nrow(dff)/12)+9, each = 12))

dff %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = namesL, values_to = 'Obs') %>%
  drop_na(Obs) %>%
  group_by(name) %>%
  mutate(z = (Obs - mean(Obs)) / sd(Obs),
    Obs = ifelse(abs(z) > 3, NA, Obs)) %>%
  drop_na(Obs) %>%
  mutate(name = factor(name, levels = namesL)) %>%
  mutate(Canal = ifelse(name %in% c('Canal-1','Canal-2'), 
    'Up', 'Down')) %>%
  filter(name %in% c('Canal-3', 'Canal-4')) %>%
  ggplot(aes((id*5)/3600, Obs, 
    color = name, group = name, shape = name)) +
  geom_line() +
  scale_color_manual(values = getPalette) +
  scale_y_continuous(breaks = seq(7,8.15, by = 0.1)) +
  scale_x_continuous(breaks = seq(0, x_breaks, by = 2)) +
  # facet_grid( Canal ~., scales = "free_y") +
  labs(y = 'pH', x = 'Hpf') +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank()) +
  guides(colour = guide_legend("")) -> psave

ggsave(psave, 
  filename = paste0(path, "/pH_time_series_data_cleaned.png"), 
  width = 12, height = 4)

psave + facet_grid(~ Stage, scales = 'free_x', space = 'free_x') -> psave

ggsave(psave, 
  filename = paste0(path, "/pH_time_series_data_cleaned_grid.png"), 
  width = 12, height = 4)

library(clock)
# install.packages('clock')

dff <- as_tibble(dff)

dff %>%
  separate(col = V1, into = c('day', 'month', 'year'), sep = '/') 


# volts ----
# c('V5','V6','V7','V8','V9')


path <- '~/Documents/DOCTORADO/pH_measures/2022_/'

file <- '*_150222.csv$'

file <- list.files(path, pattern = file, full.names = TRUE)

df <- lapply(file, read.table)


df1 <- do.call(rbind, df)

namesL <- c('Canal-1', 'Canal-2', 'Canal-3', 'Canal-4')

names(df1)[5:9] <- c(namesL, 'Canal-x')

# mVolts

