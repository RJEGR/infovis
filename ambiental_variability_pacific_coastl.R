rm(list = ls())

options(stringsAsFactors = FALSE)


path <- '~/Documents/DOCTORADO/'
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

# psave +
#   geom_boxplot(aes(y = pCO2/2500)) +
#   scale_y_continuous(name = 'pH', sec.axis = sec_axis(~.*2500, 
#                                          name = expression(~pCO[2])))

ggsave(psave, path = path, filename = 'laboratory_pH_essays.png', width = 4, height = 2)
  
# only_cgigas

df %>% 
  filter(grepl('Crassostrea', sp)) %>%
  mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  ggplot(aes(x = sp, y = pH, color = Stage)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
               outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
               size=1, position = position_dodge(0.6)) +
  # geom_point(position = position_dodge(0.6)) +
  ggsci::scale_color_aaas(name = '') +
  labs(y = '', x = '') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, 
                                   vjust = 1, size = 12)) -> psave

ggsave(psave, path = path, 
       filename = 'laboratory_pH_essays_oyster.png', 
       width = 4, height = 3)

df %>% 
  filter(grepl('Crassostrea', sp)) %>%
  mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  group_by(sp, Stage) %>% tally()

# by time of dev

df %>% 
  filter(grepl('Crassostrea', sp)) %>%
  mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  ggplot(aes(x = dias, y = pH, color = sp)) +
  # geom_point(alpha = 0.5, color = 'white') +
  ggsci::scale_color_d3(name = '') +
  scale_x_continuous(breaks = seq(0, 55, by = 2)) -> p 


p + 
  geom_point(aes(x = horas/24), alpha = 0.7) +
  scale_x_continuous(name = 'Dias', sec.axis = sec_axis(~.*24, 
                                                      name = 'Horas')) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top")

# segmentar por etapas del desarrollo el grafico?
# Efecto ----
# test Hedge's d https://github.com/mtorchiano/effsize/

file <- 'meta_analisis_OA_essays_mollusk_efecto.csv'

file <- list.files(path, pattern = file, full.names = TRUE)


df1 <- read.csv(file)

df1 %>% 
  as_tibble() %>%
  filter(grepl('Crassostrea', sp)) %>%
  group_by(Group, Efecto) %>% 
  mutate(Efecto = ifelse(Efecto == 0, 'Neutral',
                                     ifelse(Efecto == 1, 'Pos', 'Neg'))) %>%
  tally() %>% 
  # mutate(Efecto = ifelse(grepl('*concha$', Efecto), 'Calcificación', Efecto)) %>%
  pivot_wider(names_from = Efecto, values_from = n, values_fill = 0) %>%
  view()

df1 %>% distinct(Medicion)


df1 %>% 
  as_tibble() %>%
  filter(grepl('Crassostrea', sp)) %>%
  mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  mutate(Efecto = ifelse(Efecto == 0, 'Neutral',
                         ifelse(Efecto == 1, 'Pos', 'Neg'))) -> df2
df2 %>%  
  ggplot(aes(x = Tiempo, y = deltapH, color = Group)) + # color = as.factor(Efecto
  geom_point(alpha = 1, size = 4, shape = 4) + # , shape = 3
  labs(y = 'pH (Delta)', x = 'Horas') +
  theme_bw(base_size = 14, base_family = "GillSans") +
  ggsci::scale_color_d3(name = '') +
  # scale_y_continuous(limits = c()) +
  facet_wrap(~ Efecto)
  # scale_x_continuous(breaks = seq(0, 55, by = 2)) 

library(ggalluvial)

ggplot(data = df2,
       aes(axis1 = Group, axis2 = as.factor(Tiempo), axis3 = as.factor(Efecto))) +
  scale_x_discrete(limits = c("Parametro", "Horas", "Efecto"), expand = c(.2, .05))
  stat_alluvium(aes(fill = as.factor(Efecto)), lode.guidance = "backfront") +
  geom_stratum() +
  stat_stratum(alpha = .25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  ggsci::scale_fill_d3(name = '') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top") 


