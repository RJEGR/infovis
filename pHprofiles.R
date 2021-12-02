rm(list = ls())

options(stringsAsFactors = FALSE)


path <- '~/Documents/DOCTORADO/pH_measures/'

file <- '.dat$'

file <- list.files(path, pattern = file, full.names = TRUE)

# read_tbl <- function(x) {
#   df <- read.table(x, sep = "")
#   return(df)
# }

df <- lapply(file, read.table)

df1 <- rbind(df[[1]], df[[2]])

df2 <- rbind(df[[3]], df[[4]], df[[5]])

times <- rbind(df1[1:4], df2[1:4])

tbl1 <- df1[11:14]
tbl2 <- df2[10:13]

namesL <- c('Canal-1', 'Canal-2', 'Canal-4', 'Canal-3')
# c('Control', 'Experimental', 'Pecera-1', 'Pecera-2')

names(tbl1) <- namesL
names(tbl2) <- namesL


tbl <- rbind(tbl1, tbl2)

library(tidyverse)

names(tbl) <-namesL

# tbl %>%
#   pivot_longer(cols = names(tbl)) %>%
#   ggplot(aes(value)) +
#   geom_histogram() + facet_grid(~ name)

# namesL <- c('R1', 'Pecera-R1', 'R2', 'Pecera-R2')

head(dff <- cbind(times, tbl))

dff %>%
  pivot_longer(cols = names(tbl)) %>%
  mutate(name = factor(name, levels = namesL)) %>%
  # filter(name != 'R1') %>%
  mutate(group = ifelse(!name %in% 'Control', 'Experimental', 'Control')) %>%
  ggplot(aes(value+0.0375, name, color = V1)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
    outlier.alpha = 0) +
  stat_summary(fun=mean, geom="point", shape=23, 
    size=1, position = position_dodge(0.6)) +
  scale_x_continuous(breaks = seq(7,8.15, by = 0.25), limits = c(7,8.15)) +
  labs(y = '', x = 'pH') +
  theme_bw(base_size = 14, base_family = "GillSans") +
  guides(colour = guide_legend("Fecha")) -> psave
  # facet_grid(group ~., space = 'free_y', scales = 'free')

ggsave(psave, filename = paste0(path, "pH_boxplot.png"), width = 7, height = 4)

dff %>%
  pivot_longer(cols = names(tbl)) %>%
  mutate(name = factor(name, levels = namesL)) %>%
  mutate(V3 = ifelse(grepl('^a.',V3), 'a.m', 'p.m')) %>%
  # mutate(V1 = ifelse(V1 %in% '19/11/2021' & V3 %in% 'a.m', NA, V1)) %>%
  drop_na(V1) %>%
  group_by(V1, name) %>%
  summarise(a = mean(value+0.0375)) %>%
  pivot_wider(names_from = name, values_from = a) %>%
  select(namesL)



# dff %>% as_tibble() %>%
  # mutate(V1 = gsub(c('/'), '-', V1)) %>%
  # mutate(V1 = lubridate::as_datetime(V1))
  
# as.POSIXct(mydataframe$Time, format = "%m-%d-%Y %H:%M:%S")
# 
# dff %>% pivot_longer(cols = names(tbl)) %>%
#   mutate(name = factor(name, levels = namesL)) %>%
#   mutate(V2 = as.POSIXct(V2,  format = "%H:%M:%S")) %>%
#   ggplot(aes(V2, value+0.2, color = V1)) +
#   # geom_line() +
#   geom_point(alpha = 0.5) +
#   # scale_x_time(breaks = scales::date_breaks('30 mins')) +
#   facet_grid( name ~.) +
#   labs(y = 'pH', x = 'Hora')

#

# mutate(Time = as.POSIXct(Time, origin = "2018-01-01", tz = "GMT")) %>%
#   ggplot(aes(Time, Value)) +
#   geom_point() +
#   scale_x_datetime(date_breaks = "30 min", date_labels = "%H:%M")

dff %>% pivot_longer(cols = names(tbl)) %>%
  # sample_n(1000) %>%
  mutate(name = factor(name, levels = namesL)) %>%
  mutate(V2 = as.POSIXct(V2,  format = "%H:%M:%S")) %>% # %H:%M:%S
  mutate(V3 = ifelse(grepl('^a.',V3), 'a.m', 'p.m')) %>%
  mutate(g = ifelse(name %in% 'Canal-1', 'Control', 'Experimental')) %>%
  mutate(V1 = ifelse(V1 %in% '19/11/2021' & V3 %in% 'a.m', NA, V1)) %>%
  drop_na(V1) %>%
  ggplot(aes(V2, value+0.0375, color = name)) +
  # geom_line() +
  geom_point(alpha = 0.5) +
  scale_x_datetime(breaks = scales::date_breaks("3 hours"), date_labels = "%H") +
  # scale_y_continuous(breaks = seq(7,8.15, by = 0.2), limits = seq(7,8.15, by = 0.2)) +
  ggh4x::facet_nested(~ V1 + V3, scales = "free") +
  labs(y = 'pH', x = 'Hora') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(panel.spacing = unit(0, "lines"), legend.position = "bottom") +
  guides(colour = guide_legend("")) -> psave

ggsave(psave, filename = paste0(path, "pH_time_series.png"), width = 10, height = 7)

  