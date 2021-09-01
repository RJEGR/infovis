rm(list = ls())

options(stringsAsFactors = FALSE)


path <- '~/Documents/DOCTORADO'
file <- 'meta_analisis_OA_essays_mollusk.csv'

file <- list.files(path, pattern = file, full.names = TRUE)


df <- read.csv(file)
# 'Tiempo' # (horas)

summary(df)

library(tidyverse)
library(ggplot2)

head(df)

df %>% group_by(Stage, Nombre) %>% tally()

df %>% ggplot(aes(x = Stage, y = pH, color = Nombre)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
               outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
               size=1, position = position_dodge(0.6)) +
  # geom_point(position = position_dodge(0.6)) +
  ggsci::scale_color_d3(name = '') +
  labs(y = '', x = '') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top") -> psave

psave +
  geom_boxplot(aes(y = pCO2/8)) +
  scale_y_continuous(name = 'pH', sec.axis = sec_axis(~.*8, 
                                         name = expression(~pCO[2])))

ggsave(psave, path = path, filename = 'laboratory_pH_essays.png', width = 4, height = 2)
  
# only_cgigas

df %>% 
  filter(grepl('Crassostrea', sp)) %>%
  ggplot(aes(x = sp, y = pH, color = Stage)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
               outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
               size=1, position = position_dodge(0.6)) +
  # geom_point(position = position_dodge(0.6)) +
  ggsci::scale_color_d3(name = '') +
  labs(y = '', x = '') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top") 

# by time of dev

df %>% 
  filter(grepl('Crassostrea', sp)) %>%
  ggplot(aes(x = dias, y = pH, color = sp)) +
  # geom_point(alpha = 0.5, color = 'white') +
  ggsci::scale_color_d3(name = '') +
  scale_x_continuous(breaks = seq(0, 55, by = 2)) -> p 


p + 
  geom_point(aes(x = horas/24), alpha = 0.7) +
  scale_x_continuous(name = 'Dias', sec.axis = sec_axis(~.*24, 
                                                      name = 'Horas')) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top")

# segmentar por etapas del desarrollo el grafico?



