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
  ggplot(aes(x = sp, y = pH, color = Stage)) + # 
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
               outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
               size=1, position = position_dodge(0.6)) +
  # geom_point(position = position_dodge(0.6)) +
  # ggsci::scale_color_d3(name = '') +
  scale_color_manual(name = '', values = c('grey56', 'black')) +
  labs(y = '', x = '') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, 
                                   vjust = 1, size = 12)) -> psave

ggsave(psave, path = path, 
       filename = 'laboratory_pH_essays_oyster.png', 
       width = 4, height = 3)


# only abulon

df %>% 
  filter(grepl('Haliotis', sp)) %>%
  mutate(sp = str_replace(sp, "Haliotis", "H.")) %>%
  ggplot(aes(x = sp, y = pH)) + # 
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
    outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
    size=1, position = position_dodge(0.6)) +
  # scale_color_manual(name = '', values = c('grey56', 'black')) +
  labs(y = '', x = '') +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    axis.text.x = element_text(angle = 45,
      hjust = 1, 
      vjust = 1, size = 12)) -> psave

ggsave(psave, path = path, 
  filename = 'laboratory_pH_essays_abulon.png', 
  width = 4, height = 3)


df %>% 
  filter(grepl('Haliotis', sp)) %>%
  mutate(sp = str_replace(sp, "Haliotis", "H.")) %>%
  group_by(sp, Stage) %>% tally()

# by time of dev

df %>%
  # filter(grepl('Crassostrea', sp)) %>%
  # mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  ggplot(aes(pH)) +
  geom_histogram(bins = 30)


df %>% 
  filter(grepl('Haliotis', sp)) %>%
  mutate(sp = str_replace(sp, "Haliotis", "H.")) %>%
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


df %>% 
  filter(grepl('Crassostrea', sp)) %>%
  mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  ggplot(aes(pH)) +
  geom_histogram(bins = 30)


file <- 'meta_analisis_OA_essays_mollusk_efecto.csv'

file <- list.files(path, pattern = file, full.names = TRUE)


df1 <- read.csv(file)

df1 %>% 
  as_tibble() %>%
  filter(grepl('Crassostrea', sp)) %>%
  filter(Multiestresor == FALSE) %>% 
  group_by(Medicion, Efecto) %>% 
  mutate(Efecto = ifelse(Efecto == 0, 'Neutral',
                                     ifelse(Efecto == 1, 'Pos', 'Neg'))) %>%
  tally() %>% 
  # mutate(Efecto = ifelse(grepl('*concha$', Efecto), 'Calcificación', Efecto)) %>%
  pivot_wider(names_from = Efecto, values_from = n, values_fill = 0) 

# circoz

df1 %>% 
  as_tibble() %>%
  filter(grepl('Crassos', sp)) %>%
  # filter(grepl('Haliotis', sp)) %>%
  filter(Multiestresor == FALSE) -> tbl

n <- length(unique(tbl$Group))
grid.col <- ggsci::pal_d3(alpha = 0.5)(n)
names(grid.col) <- sort(unique(tbl$Group))

effPalette <- c('blue', 'red', 'grey68')
names(effPalette) <- c('Positivo', 'Negativo', 'Neutral')

TimeL <- c('1', '2', '3', '6', '20', '30')

tbl %>%
  group_by(Group, Efecto, Tiempo) %>%
  mutate(Efecto = ifelse(Efecto == 0, 'Neutral',
    ifelse(Efecto == 1, 'Positivo', 'Negativo'))) %>%
  tally() %>%
  mutate(Tiempo = factor(Tiempo, levels = TimeL)) %>%
  mutate(n = ifelse(n > 1, '*', '')) %>%
  ggplot(aes(y = Group, x = as.character(Tiempo))) +
  geom_tile(aes(fill = Efecto), color = "black", size = 0.5) +
  geom_text(aes(label = n)) +
  scale_fill_manual(values = effPalette) +
  scale_x_discrete(limits = TimeL) +
  theme_classic(base_family = "GillSans", base_size = 16) +
  labs(x = 'Tiempo (Días despues de la fertilización)', y = 'Medición') +
  facet_wrap(~Efecto) +
  theme(legend.position = 'none', strip.background = element_blank()) -> pheat

ggsave(pheat, path = path, filename = 'effect_heatmap.png', width = 7, height = 3)

library(circlize)

circos.clear()

circos.par(start.degree = 0, gap.degree = 4, 
  track.margin = c(-0.01, 0.01), 
  points.overflow.warning = FALSE)


tbl %>% distinct(Tiempo)

tbl %>%
  select(Group, Efecto, Tiempo) %>%
  mutate(Efecto = ifelse(Efecto == 0, 'Neutral',
    ifelse(Efecto == 1, 'Positivo', 'Negativo'))) %>%
  data.frame() %>%
  with(., table(Group, Efecto)) %>% 
  chordDiagram(
    grid.col = c(grid.col, effPalette),
    directional = -1,
    diffHeight = mm_h(5), target.prop.height = mm_h(4),
    # annotationTrack = "grid",
    # direction.type = c("arrows", "diffHeight"),
    preAllocateTracks = 1,
    small.gap = 5, big.gap = 10)

df1 %>% ggplot(aes(x = deltapH, y = Group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
    outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
    size=1, position = position_dodge(0.6)) +
  # scale_y_discrete(position = "right") +
  theme_classic(base_family = "GillSans", base_size = 16) +
  labs(x = 'pH', y = '') +
  theme(legend.position = 'none', strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y  = element_blank()) -> pbox
# 
# ggsave(pbox, path = path, filename = 'effect_boxplot.png', width = 4, height = 4)


pbox + pheat + theme(axis.title.y = element_blank()) + plot_layout(widths = c(1, 2))

df1 %>% 
  as_tibble() %>%
  filter(grepl('Crassostrea', sp)) %>%
  mutate(sp = str_replace(sp, "Crassostrea", "C.")) %>%
  mutate(Efecto = ifelse(Efecto == 0, 'Neutral',
                         ifelse(Efecto == 1, 'Pos', 'Neg'))) -> df2
df2 %>%  
  ggplot(aes(x = Tiempo, y = deltapH, color = Group)) + 
  geom_segment(aes(x = 0, xend = Tiempo, y= 0, yend = deltapH)) +
  geom_point(alpha = 1, size = 4, aes(shape = Efecto)) + # , shape = 3
  labs(y = 'pH (Delta)', x = 'Horas') +
  theme_bw(base_size = 14, base_family = "GillSans") +
  ggsci::scale_color_d3(name = '') 

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


