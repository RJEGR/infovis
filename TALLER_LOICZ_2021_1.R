library(tidyverse)
library(geosphere)
library(patchwork)

rm(list = ls())


wd <- "~/Documents/DOCTORADO/"

setwd(wd)

t <- read.csv('TALLER_2021_1_DATOS_transecto_BSQ.csv') %>% as_tibble()



t %>% distinct(Lat) %>% pull() -> x 
distm(t(rbind(x,x)), fun = distCosine)[,1] -> dist

# meters to km 
dist/1000 -> dist

# 
# Caja A Caja C Oceano 
# 6      4      3

t  %>%
  mutate(cumdist = rep(dist, 2)) %>%
  mutate(reg = '') %>%
  mutate(reg = ifelse(Estacion %in% 1:5, 'Oceano', reg)) %>% 
  mutate(reg = ifelse(Estacion %in% 6:12, 'Caja C', reg)) %>%
  mutate(reg = ifelse(Estacion %in% 13:21, 'Caja A', reg)) %>%
  mutate(Escenario = ifelse(Escenario == 'A', 'Escenario 1', 'Escenario 2')) %>%
  mutate(Estacion = paste0('E', Estacion)) -> t

t %>% mutate(reg = factor(reg, levels = c('Oceano', 'Caja C', 'Caja A'))) -> t

# t %>% group_by(reg, Escenario) %>% 
#   summarise(mT = mean(TEMP), mS = mean(SALINIDAD)) %>% view()

# graficar los perfiles de sal y temperatura.. ----
xlab <- expression("Distancia relativa a la estación 1"~(Km))
color_values <- c("#F78482", "#02B8BB")

t  %>%
  ggplot(aes(y = TEMP, x = cumdist, color = Escenario)) +
  geom_vline(xintercept = c(8.5, 13), 
             colour = "grey75", alpha = 0.75, linetype = 'dashed') +
  # annotate("text", x = c(3.7, 10.7, 17), y = 5, 
  #          label = c("Ocean", "BoxC", "BoxA"), 
  #          colour = "grey70",
  #          size = 4) +
  geom_point(alpha = 0.7) +
  # geom_path() +
  # scale_y_reverse()+  
  scale_x_continuous(position = "top") +
  ggthemes::scale_color_calc(name ='') +
  theme_classic(base_family = "GillSans", base_size = 14) +
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1)) + 
  labs(x = xlab, y = expression(~Temperatura~('°C'))) +
  theme(legend.position = "none") -> p1

t %>% distinct(reg, .keep_all = T) %>% select(reg, cumdist)
#


t %>% select(TEMP, SALINIDAD, Escenario, reg, cumdist, Estacion) %>%
  pivot_longer(cols = c('TEMP', 'SALINIDAD')) %>%
  mutate(name = ifelse(name %in% 'SALINIDAD', 'Salinidad', 'Temperatura (°C)'))-> tLonger

tLonger %>% distinct(Estacion) %>% pull() -> sec.axis.labels
tLonger %>% distinct(cumdist) %>% pull() -> sec.axis.brakes

tLonger %>% 
  ggplot(aes(y = value, x = cumdist, color = reg)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_x_continuous(position = "bottom") +
  geom_point(alpha = 0, size = 0) +
  # scale_x_continuous('Estación',
  #                    sec.axis = sec_axis(~., name = xlab), 
  #                    labels = sec.axis.labels, 
  #                    breaks = sec.axis.brakes) +
  ggsci::scale_color_igv(name ='') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = xlab, y = '') +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 3)) +
  facet_grid(name ~ Escenario , scales = 'free_y') -> psave

ggsave(psave, path = wd, filename = 'transectos_T_S.png', width = 5, height = 5)


t  %>% 
  ggplot(aes(y = SALINIDAD, x = cumdist, color = Escenario)) +
  geom_vline(xintercept = c(8.5, 13), 
             colour = "grey75", alpha = 0.7, linetype = 'dashed') +
  annotate("text", x = c(3.7, 10.7, 17), y = 30, 
           label = c("Ocean", "BoxC", "BoxA"), 
           colour = "grey50",
           size = 4) +
  geom_point(alpha = 0.7) +
  # geom_path() +
  # scale_y_reverse()+  
  scale_x_continuous(position = "top") +
  ggthemes::scale_color_calc(name ='') +
  theme_classic(base_family = "GillSans", base_size = 14)+
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1)) + labs(x = expression(~Distancia~(Km)), y = expression(Salinidad)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top") +
  guides(color = guide_legend(ncol = 1)) -> p2

p1 / p2 -> psave



# ggsave(psave, path = wd, filename = 'transectos_T_S.png', width = 5, height = 5)

tLonger %>%
  ggplot(aes(x = reg, y = value, color = Escenario)) +
  facet_grid(name ~., scales = 'free_y') + 
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6)) +
  # geom_point(aes(color = Escenario)) +
  # ggrepel::geom_text_repel(aes(color = reg, label = Estacion)) +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  ggthemes::scale_color_calc(name ='') -> psave

psave + theme(legend.position = "top",
              strip.background = element_blank(),
              axis.line.y =  element_line(colour = "black"),
              panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 1)) -> psave

ggsave(psave, path = wd, filename = 'boxplot_T_S.png', width = 3, height = 4.5)


shape <- c(0, 3)


# graficar el diagrama de mezcla del aicdo silicico de cada escenario ----

t %>% group_by(Escenario) %>% filter(Estacion %in% 'E1') %>% pull(SILICATO) -> a

t %>% group_by(Escenario) %>% filter(Estacion %in% 'E1') %>% pull(SALINIDAD) -> minS
t %>% group_by(Escenario) %>% filter(Estacion %in% 'E21') %>% pull(SALINIDAD) -> maxS

(maxS/minS)-> b
Lteorica <- a*b

# Example. 0.5 unidades de una sustancia disuelta en el oceano y un 10% de incremento de una sustancia conservativa. (10*0.5/100)+0.5

data.frame(Salinidad_Oc = minS , Salinidad_Y = maxS,
           Factor = round(b, digits = 2),
           AcS_Oc = a, AcS_Y = round(Lteorica, digits = 2),
           Escenario = c("Escenario 1", "Escenario 2"))

data.segm <- data.frame(x = minS ,y = a, xend = maxS, yend = Lteorica,
                      Escenario = c("Escenario 1", "Escenario 2"))

t %>%
  ggplot(aes(y = SILICATO, 
             x = SALINIDAD, color = reg)) +
  theme_bw(base_family = "GillSans", base_size = 18) + 
  labs(y = expression(~H[4]~SiO[4]~(µmol~Kg^{-1})), x = 'Salinidad') +
  ggsci::scale_color_igv(name ='') +
  facet_grid(~ Escenario) +
  geom_segment(data=data.segm,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE,
               color = 'grey15', size = 0.9, linetype='dashed') + #
  # geom_smooth(color = 'black', alpha = 0.4) +
  geom_point(size = 2) +
  # ggpubr::stat_cor(color = 'black', label.y = 35, label.x = 34) +
  # ggrepel::geom_text_repel(aes(label = Estacion)) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 3)) -> psave

# width = 7, height = 7.3
ggsave(psave, path = wd, filename = 'diagMezcla_Si.png', width = 7.5, height = 4.5)

# Balances ---- 

t %>% group_by(Escenario, reg) %>%
  filter(!reg %in% 'Oceano') %>%
  summarise(mean(SALINIDAD))

# salinidad del sistema (caja A)
t %>% group_by(Escenario) %>% 
  filter(reg %in% 'Caja A') %>%
  summarise(mean(SALINIDAD))

# salinidad del flujo residual (Caja A : Caja C)
t %>% group_by(Escenario) %>% 
  filter(!reg %in% 'Oceano') %>%
  summarise(SR = mean(SALINIDAD))

# salinidad del flujo adyacente (Caja A : Caja C)
t %>% group_by(Escenario, reg) %>% 
  filter(reg %in% 'Caja C') %>%
  summarise(mean(SALINIDAD))

# determinar la estructura de las comunidades
nitr.label <- quote(NO[3]~(µmol~L^{-1})) 
fos.label <- quote(PO[4]~(µmol~L^{-1}))
tem.label <- quote(Temperatura~('°C'))


t %>% 
  select(SILICATO, NITRATO, FOSFATO, TEMP,
         Escenario, reg, cumdist, Estacion) %>%
  pivot_longer(cols = c('NITRATO', 'FOSFATO')) %>%
  mutate(name = factor(name, 
                       labels = c(nitr.label, fos.label))) %>%
  ggplot(aes(x=SILICATO, y = value, group = name, 
             color = reg)) +
  facet_grid(name~Escenario, scales = 'free',
             labeller = labeller(name = label_parsed)) +
  geom_smooth(method = "lm", se = F, color = 'black', formula = y ~ x) +
  geom_point() +
  ggpmisc::stat_poly_eq(aes(label = stat(eq.label)), formula = y ~ x,
                        parse = T) +
  # scale_color_viridis_c(option = 'magma') # if  color = TEMP
  ggsci::scale_color_igv(name ='') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = expression(~H[4]~SiO[4]~(µmol~Kg^{-1})), y = '') +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) -> psave

ggsave(psave, path = wd, filename = 'nutrientes_vs_DSi.png', width = 6.5, height = 6.5) 

t %>% 
  select(SILICATO, TEMP,
         Escenario, reg, cumdist, Estacion) %>%
  ggplot(aes(x=SILICATO, y = TEMP)) +
  facet_grid(Escenario~., scales = 'free',
             labeller = labeller(name = label_parsed)) +
  geom_point(aes(color = reg), size = 3) +
  ggpubr::stat_cor(label.x = 25, p.digits = 2) +
  ggsci::scale_color_igv(name ='') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = expression(~H[4]~SiO[4]~(µmol~Kg^{-1})), 
       y = expression(Temperatura~('°C'))) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) 
  # ggforce::geom_mark_ellipse(aes(group = reg))





# 10) transectos de nitrato y fosfato del océano hasta la estación 21 para cada una de los escenarios. ----

# 

nitr.label <- quote(DIN~(µmol~L^{-1}))  # quote(NO[3]~(µmol~L^{-1})) 
fos.label <- quote(DIP~(µmol~L^{-1})) # quote(PO[4]~(µmol~L^{-1}))
facet.labs <- c(`P` = fos.label, `N` = nitr.label)

# nitr.label <- expression(~NO[3](µmol L^{-1}))
# fos.label <- expression(~PO[4](µmolL^{-1}))

         
t %>% 
  select(NITRATO, FOSFATO, Escenario, reg, cumdist, Estacion) %>%
  pivot_longer(cols = c('NITRATO', 'FOSFATO')) %>%
  mutate(name = ifelse(name %in% 'NITRATO', 'N', 'P')) -> tLonger


tLonger %>%
  group_by(reg, Escenario, name) %>%
  summarise(m = mean(value)) 

tLonger %>%
  # mutate(name = ifelse(name %in% 'NITRATO', nitr.label, fos.label)) %>%
  mutate(name = factor(name, labels = c(nitr.label, fos.label))) -> tLonger

tLonger %>%
  ggplot(aes(y = value, x = cumdist, color = reg)) +
  # geom_path() +
  scale_x_continuous(position = "bottom") +
  # geom_smooth(color = 'black', method = "lm", se = F, linetype='dashed') +
  geom_point(alpha = 0.7, size = 2) +
  # ggpubr::stat_cor(color = 'black') +
  ggsci::scale_color_igv(name ='') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = xlab, y = '') +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 3)) +
  facet_grid(name~Escenario, scales = 'free', 
             # switch = 'y',
             labeller = labeller(name = label_parsed)) -> psave

ggsave(psave, path = wd, filename = 'transectos_N_P.png', width = 5, height = 5)

tLonger %>%
  ggplot(aes(x = reg, y = value, color = Escenario)) +
  facet_grid(name ~., scales = 'free_y', 
             labeller = labeller(name = label_parsed)) + 
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6)) +
  # geom_point(aes(color = Escenario)) +
  # ggrepel::geom_text_repel(aes(color = reg, label = Estacion)) +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  ggthemes::scale_color_calc(name ='') -> psave

psave + theme(legend.position = "top",
              strip.background = element_blank(),
              axis.line.y =  element_line(colour = "black"),
              panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 1)) -> psave

ggsave(psave, path = wd, filename = 'boxplot_N_P.png', width = 3, height = 4.5)


# (11) Dibujar los diagramas de mezcla para N y P para cada escenario. ----
# Fosfato 
# t %>% group_by(Escenario) %>% filter(Estacion %in% 'E1') %>% pull(FOSFATO) -> a
# t %>% group_by(Escenario) %>% filter(Estacion %in% 'E1') %>% pull(SALINIDAD) -> minS
# t %>% group_by(Escenario) %>% filter(Estacion %in% 'E21') %>% pull(SALINIDAD) -> maxS

# use mean istead of 1st station

t %>% group_by(Escenario) %>% filter(reg %in% 'Oceano') %>% summarise(x=mean(FOSFATO)) %>% pull(x) -> a
t %>% group_by(Escenario) %>%  filter(reg %in% 'Oceano') %>% summarise(x=mean(SALINIDAD)) %>% pull(x) -> minS
t %>% group_by(Escenario) %>% filter(Estacion %in% 'E21') %>% pull(SALINIDAD) -> maxS

(maxS/minS)-> b
Lteorica <- a*b


data.frame(Salinidad_Oc = minS , Salinidad_Y = maxS,
           Factor = round(b, digits = 2),
           x_Oc = a, x_Y = round(Lteorica, digits = 2),
           Escenario = c("Escenario 1", "Escenario 2")) 

data.segm <- data.frame(x = minS ,y = a, xend = maxS, yend = Lteorica,
                        Escenario = c("Escenario 1", "Escenario 2"))

t %>%
  ggplot(aes(y = FOSFATO, 
             x = SALINIDAD, color = reg)) +
  theme_bw(base_family = "GillSans", base_size = 18) + 
  labs(y = fos.label, x = 'Salinidad') +
  ggsci::scale_color_igv(name ='') +
  facet_grid(~ Escenario) +
  geom_segment(data=data.segm,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE,
               color = 'grey15', size = 0.9, linetype='dashed') + # 
  # geom_smooth(color = 'black', alpha = 0.4) +
  geom_point(size = 2) +
  # ggpubr::stat_cor(color = 'black') +
  # ggrepel::geom_text_repel(aes(label = Estacion)) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 3)) -> psave

# width = 7, height = 7.3
ggsave(psave, path = wd, filename = 'diagMezcla_P.png', width = 7.5, height = 4.5)

# Nitrato

# t %>% group_by(Escenario) %>% filter(Estacion %in% 'E1') %>% pull(NITRATO) -> a
# t %>% group_by(Escenario) %>% filter(Estacion %in% 'E1') %>% pull(SALINIDAD) -> minS
# t %>% group_by(Escenario) %>% filter(Estacion %in% 'E21') %>% pull(SALINIDAD) -> maxS

# use mean istead of 1st station

t %>% group_by(Escenario) %>% filter(reg %in% 'Oceano') %>% summarise(x=mean(NITRATO)) %>% pull(x) -> a


(maxS/minS)-> b
Lteorica <- a*b


data.frame(Salinidad_Oc = minS , Salinidad_Y = maxS,
           Factor = round(b, digits = 2),
           x_Oc = a, x_Y = round(Lteorica, digits = 2),
           Escenario = c("Escenario 1", "Escenario 2")) 

data.segm <- data.frame(x = minS ,y = a, xend = maxS, yend = Lteorica,
                        Escenario = c("Escenario 1", "Escenario 2"))

t %>%
  ggplot(aes(y = NITRATO, 
             x = SALINIDAD, color = reg)) +
  theme_bw(base_family = "GillSans", base_size = 18) + 
  labs(y = nitr.label, x = 'Salinidad') +
  ggsci::scale_color_igv(name ='') +
  facet_grid(~ Escenario) +
  geom_segment(data=data.segm,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE,
               color = 'grey15', size = 0.9, linetype='dashed') + # 
  # geom_smooth(color = 'black', alpha = 0.4) +
  geom_point(size = 2) +
  # ggpubr::stat_cor(color = 'black') +
  # ggrepel::geom_text_repel(aes(label = Estacion)) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 3)) -> psave

# width = 7, height = 7.3
ggsave(psave, path = wd, filename = 'diagMezcla_N.png', width = 7.5, height = 4.5)

# balance de nutrientes: ----
# masa molec
no3 <- 62.0049 # g/mol
po4 <- 94.9714 # g/mol * 1000 -> mg/mol

# Concentración promedio de elementos en agua subterránea (YG)

DIPg <- c(0.2, 0.2) # mg/L
DINg <- c(15, 15) # mg/L

# de L a m3
DIPg/0.001 -> DIPg # mg/m3
DINg/0.001 -> DINg

# conversion de unidades
# Para convertir de gramos a moles se divide la concentración de masa entre la concentración molar del compuesto

DIPg/(po4*1000) # mol/m3
DINg/(no3*1000)

# ratios
# 
t %>% 
  pivot_longer(cols = c('NITRATO', 'FOSFATO')) %>%
  mutate(ratio = 1 - value/SILICATO) %>%
  ggplot(aes(y = ratio, x = reg, color = Escenario)) +
  facet_wrap(~name, scales = 'free_y') +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6)) +
  # geom_point(aes(color = Escenario)) +
  # ggrepel::geom_text_repel(aes(color = reg, label = Estacion)) +
  labs(x = '', y = 'Y:Si Ratios') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  ggthemes::scale_color_calc(name ='')

# 
# 'SILICATO', 'SALINIDAD', 'TEMP')

t %>%
  pivot_longer(cols = c('NITRATO', 'FOSFATO')) %>%
  ggplot(aes(y = value, 
             x = SILICATO, color = reg)) +
  facet_grid(name ~ Escenario, scales = 'free') +
  theme_bw(base_family = "GillSans", base_size = 18) + 
  # labs(y = nitr.label, x = 'Salinidad') +
  ggsci::scale_color_igv(name ='') +
  # geom_smooth(color = 'black', alpha = 0.4) +
  geom_smooth(method = "lm", se = F, color = 'black', formula = y ~ x) +
  geom_point(size = 2) +
  # ggpubr::stat_cor(color = 'black') +
  ggpubr::stat_cor(color = 'black') +
  # ggrepel::geom_text_repel(aes(label = Estacion)) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  guides(color = guide_legend(ncol = 3)) 

cols <- c('NITRATO', 'FOSFATO', 'SILICATO', 'SALINIDAD',
          'TEMP')

t %>%
  pivot_longer(cols = cols) %>%
  group_by(Escenario, name) %>%
  summarise(cor(value, cumdist))

M <- t %>% filter(Escenario %in% 'Escenario 1') %>% select(c(cols, 'cumdist')) %>% cor()
library(corrplot)
corrplot(M, method = c("number"), type = 'upper')
M<- t %>% filter(Escenario %in% 'Escenario 2') %>% select(c(cols, 'cumdist')) %>% cor()
corrplot(M, method="number", type = 'upper')
