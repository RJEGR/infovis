# Ricardo Gomez Reyes, 2021

# Visualising silicate profile
# Graficar perfiles de H4SiO4 para Atlántico N y S, Indico y Pacífico N y S

# A16N2013: 30oN y 6oS
# I09N_2016: 26 oS
# P16N2015 : 30oN y 17oSt
# (notar que los datos de 2015
#   vienen en dos archivos)

# Graficar los 5 perfiles en la misma gráfica
# Graficar datos entre 0 y 5000 m de profundidad

# DEPTH METERS,
# SILCAT UMOL/KG

rm(list=ls())

library(tidyverse)
library(patchwork)
library(scales)

cols <- c('SECT_ID','GROUP', 'LATITUDE', 'LONGITUDE', 'CTDPRS', 'SILCAT', 'SALNTY', 'PHSPHT', 'OXYGEN', 'NITRAT',
          'ALKALI', 'PH') # 'CO2',
files <- list.files(pattern = 'csv', 
                    path = '~/Documents/DOCTORADO/CSV/', full.names = T)

x <- lapply(files, function(x) {
  y <- read.csv(x, skip = 1, comment.char = "#", stringsAsFactors = F)
  y <- data.frame(y, GROUP = basename(x))
  y <- y %>% select_at(vars(-contains(c('FLAG', 'TMP'))))
  # y <- y %>% select_at(vars(contains(cols))) # -contains('FLAG'), vars(contains(cols)
  y <- y %>% rename_with(~ toupper(gsub("_SWS", "_TOT", .x, fixed = TRUE)))
  return(y)})
 
x <- do.call(rbind, x) %>% 
  as_tibble() %>% drop_na() 

x %>%
  mutate(SECT_ID = as.factor(SECT_ID)) %>%
  mutate(GROUP = str_replace_all(GROUP, '.csv$', '')) %>%
  mutate(SECT_ID = str_replace_all(SECT_ID, '.[1-9]$', '')) %>%
  mutate_if(is.character, function(x) {str_replace_all(x, ' ', '')}) %>%
  mutate_at(vars(contains(cols[-c(1:2)])), as.double) -> x
  



# con respecto a la latitud
# A16N2013: 30N y 6S
# I09N_2016: 26 S
# P16N2015 : 30N y 17S

# eliminar 0.5 unidades de distancia

x %>%
  filter(SECT_ID %in% 'A16N') %>%
  filter(round(LATITUDE) == 30 | round(LATITUDE) == -6) -> A16N

A16N %>%
  mutate(SECT_ID = ifelse(round(LATITUDE) == 30, 'A16N_30N', 'A16N_6S')) -> A16N


x %>%
  filter(round(LATITUDE) == -26 & 
           SECT_ID %in% 'I09N' )  -> I09N

x %>%
  filter(SECT_ID %in% 'P16N') %>%
  filter(round(LATITUDE) == 30 | round(LATITUDE) == -17) -> P16N

P16N %>%
  mutate(SECT_ID = ifelse(round(LATITUDE) == 30, 'P16N_30N', 'P16N_17S')) -> P16N

y <- rbind(A16N, I09N, P16N)

y %>%
  mutate(ocean = '') %>%
  mutate(ocean = ifelse(grepl('^A', SECT_ID), 'Atlantic',ocean )) %>%
  mutate(ocean = ifelse(grepl('^I', SECT_ID), 'Indian',ocean )) %>%
  mutate(ocean = ifelse(grepl('^P', SECT_ID), 'Pacific',ocean )) -> y

y %>%
  mutate(CTDPRS = as.numeric(CTDPRS),
         SILCAT = as.numeric(SILCAT)) %>%
  filter(SILCAT > 0 & CTDPRS <= 5000) %>%
  ggplot(aes(y = CTDPRS, x = SILCAT, color = ocean, 
             shape = SECT_ID)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_line(orientation = 'y') +
  # geom_path() +
  scale_y_reverse() +  
  theme_bw(base_family = "GillSans", base_size = 18) +
  scale_x_continuous(position = "top") +
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1)) +
  labs(x = expression(~H[4]~SiO[4]~(µmol~Kg^{-1})), y = expression(~Pressure~(dbar)))  +
  scale_shape_discrete(name = 'Coordinates') +
  theme(legend.position = "none") -> p1

y %>%
  mutate(CTDPRS = as.numeric(CTDPRS),
         PHSPHT = as.numeric(PHSPHT)) %>%
  filter(PHSPHT > 0 & CTDPRS <= 5000) %>%
  ggplot(aes(y = CTDPRS, x = PHSPHT, color = ocean, 
             shape = SECT_ID)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_line(orientation = 'y') +
  # geom_path() +
  scale_y_reverse() +  
  theme_bw(base_family = "GillSans", base_size = 18) +
  scale_x_continuous(position = "top") +
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1))+
  labs(x = expression(~PO[4]~(µmol~Kg^{-1})), y = expression(~Pressure~(dbar)))  +
  scale_shape_discrete(name = 'Coordinates') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "right")  -> p2



saveP <- p1+p2
saveP
# Pacific 32N profiling SILICATE and PHSPTH

library(scales)

x %>%
  filter(SECT_ID %in% 'P16N') %>%
  filter(round(LATITUDE) == 32) %>%
  mutate_at(vars(c('SILCAT', 'PHSPHT', 'CTDPRS', 'SALNTY')), as.numeric) %>% 
  filter(PHSPHT > 0 & CTDPRS <= 5000) %>%
  # mutate()
  ggplot(data = ., aes(y = CTDPRS)) +
  scale_y_reverse() + 
  geom_line(aes(x = PHSPHT), colour = 'red', orientation = 'y') +
  geom_line(aes(x = SILCAT/50), color = 'blue', orientation = 'y') +
  scale_x_continuous(expression(~PO[4]~(µmol~Kg^{-1})),
                     sec.axis = sec_axis(~.*50, 
                                         name = expression(~H[4]~SiO[4]~(µmol~Kg^{-1})))) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  labs(y = expression(~Pressure~(dbar))) +
  theme(axis.text.x.top = element_text(color = "blue"),
        axis.text.x.bottom = element_text(color = "red"),
        axis.ticks.x.top = element_line(color = "blue"),
        axis.ticks.x.bottom = element_line(color = 'red'),
        axis.title.x.top = element_text(color = "blue"),
        axis.title.x.bottom = element_text(color = 'red')) -> p3 


# facet nutrients
cols <- c('SILCAT', 'PHSPHT', 'NITRAT', 'OXYGEN', 'ALKALI', 'PH_TOT', 'SALNTY') # , 
# expression(~H[4]~SiO[4]~(µmol~Kg^{-1}))

nitrL <- quote(NO[3]~(µmol~L^{-1})) 
fosL <- quote(PO[4]~(µmol~L^{-1}))
silL <- quote(H[4]~SiO[4]~(µmol~Kg^{-1}))
OxL <- quote(O[2]~(µmol~L^{-1}))

labels <- c(silL, fosL, nitrL, OxL,  'Alcalinidad', 'pH', 'Salinidad') # , 

y %>%
  filter(CTDPRS <= 5000) %>%
  pivot_longer(cols = cols, names_to = 'var') %>%
  filter(value > 0) %>%
  mutate(var = factor(var, levels = cols, labels = labels)) -> y_longer


y_longer %>% group_by(SECT_ID, var) %>% filter(value > mean(value))

# CTDPRS > 200 & CTDPRS <= 1000
y_longer %>%
  mutate(depth = '') %>%
  mutate(depth = ifelse(CTDPRS <= 200, '0-200', 
                        ifelse(between(CTDPRS, 200, 1000), '200-1000', 
                               ifelse(between(CTDPRS, 1000, 4000), '1000-4000',
                                      ifelse(CTDPRS > 4000, '> 4000', depth))))) %>%
  mutate(depth = factor(depth, levels = c('0-200', '200-1000', '1000-4000', '> 4000'))) %>%
  # group_by(depth) %>% count() %>% pull(n) %>% sum() # sanity check
  ggplot(aes(y = value, x = depth)) +
  facet_grid(var~ocean,  scales = 'free',
             labeller = labeller(var = label_parsed)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.4, position = position_dodge(0.6)) +
  stat_summary(fun=mean, geom="point", shape=23, size=1, color = 'red', position = position_dodge(0.6)) +
  labs(y = '') +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) -> psave

ggsave(psave, path = '~/Documents/DOCTORADO/', filename = 'facet_vars_box.png', width = 10, height = 10)

y_longer %>% group_by(GROUP, var, ocean) %>% 
  mutate(z = (value - mean(value))/sd(value) ) %>%
  filter(!abs(z) > 3) -> y_longer
  # group_by(var) %>%
  # summarise(q = quantile(z)) %>%
  # ggplot(aes(x=var, y = abs(q))) + geom_boxplot()


y_longer %>%
  ggplot(aes(y = CTDPRS, x = value, shape = SECT_ID)) + # 
  geom_point(alpha = 0.7, size = 1.5) +
  geom_line(orientation = 'y') +
  # geom_path() +
  scale_y_reverse() +  
  theme_bw(base_family = "GillSans", base_size = 18) +
  scale_x_continuous(position = "bottom") +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  labs(x = '', y = expression(~Profundidad~(dbar)))  +
  scale_shape_discrete(name = 'Coordenadas') +
  facet_grid(ocean~var, scales = 'free_x', labeller = labeller(var = label_parsed)) -> p

ggsave(p, path = '~/Documents/DOCTORADO/', filename = 'facet_vars_depth.png', width = 20, height = 10)

# 
# 
# y_longer %>%
#   ggplot(aes(y = value, x = SALNTY)) + # 
#   geom_point(alpha = 0.7, size = 1.5) +
#   facet_grid(var~ocean, scales = 'free', labeller = labeller(var = label_parsed))

# Graficar a una misma escala, H4SiO4 y PO4 de la coordenada 32 N del Pacifico.




x %>%
  filter(SECT_ID %in% 'P16N') %>%
  filter(round(LATITUDE) == 32) %>%
  mutate_at(vars(c('SILCAT', 'PHSPHT', 'CTDPRS', 'SALNTY')), as.numeric) %>% 
  filter(PHSPHT > 0 & CTDPRS <= 5000) %>%
  # mutate()
  ggplot(data = ., aes(y = CTDPRS)) +
  scale_y_reverse() + 
  geom_line(aes(x = PHSPHT), colour = 'red', orientation = 'y') +
  geom_line(aes(x = SILCAT/50), color = 'blue', orientation = 'y') +
  scale_x_continuous(expression(~PO[4]~(µmol~Kg^{-1})),
                     sec.axis = sec_axis(~.*50, 
                                         name = expression(~H[4]~SiO[4]~(µmol~Kg^{-1})))) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  labs(y = expression(~Profundidad~(dbar)), 
       caption = "Perfiles verticales de silicato y fosfato del Pacifico Norte (32°N)") +
  theme(axis.text.x.top = element_text(color = "blue"),
        axis.text.x.bottom = element_text(color = "red"),
        axis.ticks.x.top = element_line(color = "blue"),
        axis.ticks.x.bottom = element_line(color = 'red'),
        axis.title.x.top = element_text(color = "blue"),
        axis.title.x.bottom = element_text(color = 'red')) -> p3

ggsave(p3, path = '~/Documents/DOCTORADO/', filename = 'DSi_PO4_profiles.png', width = 5.5, height = 5.5)


 