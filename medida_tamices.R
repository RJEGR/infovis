x <- read.csv('~/Documents/DOCTORADO/tamices_completos.csv')

library(tidyverse)

findc <- function(a, b) { sqrt(a^2 + b^2) }


x %>% mutate(n = 1:nrow(.)) %>%
  pivot_longer(cols = names(x)) %>%
  mutate(id = substr(name, nchar(name), nchar(name))) %>%
  mutate(name = substr(name, 1, nchar(name)-1)) %>%
  pivot_wider(names_from = id, values_from = value) %>%
  mutate(Theorical = findc(A,B)) %>%
  ggplot(aes(Theorical, C)) +
  facet_wrap(~name, scales = 'free') +
  geom_point(aes(color = name)) +
  ggpubr::stat_cor(aes(group = name)) +
  geom_smooth(method = "lm", 
    linetype="dashed", size = 0.5, alpha=0.5, 
    se = F, na.rm = TRUE, aes(group = name)) +
  labs(x = 'c = sqrt(a^2 + b^2)')

Lname <- c('M35', 'M46', 'M55', 
  'M75', 'M105', 'M150', 
  'M200', 'M300')

x %>% mutate(n = 1:nrow(.)) %>%
  pivot_longer(cols = names(x)) %>%
  mutate(id = substr(name, nchar(name), nchar(name))) %>%
  mutate(name = substr(name, 1, nchar(name)-1)) %>% 
  pivot_wider(names_from = id, values_from = value) %>%
  mutate(Theorical = findc(A,B)) %>%
  select(name, C, Theorical) %>%
  pivot_longer(-name, names_repair = 'minimal', names_to = 'g') %>%
  mutate(g = ifelse(g %in% 'C', 'Hipotenusa','Teorico')) %>%
  # mutate(name = substr(name, 5, nchar(name))) %>%
  mutate(name = gsub(c('X4'), '', name)) %>%
  mutate(name = gsub(c('X10X'), '', name)) %>%
  mutate(name = gsub(c('X10'), '', name)) -> df_clean 

df_clean %>%
  filter(g %in% 'Teorico') %>%
  mutate(name = factor(name, levels = Lname)) %>%
  ggplot(aes(name, value)) +
  # geom_boxplot() +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6), 
    outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
    size=1, position = position_dodge(0.6), color = 'red') +
  scale_color_hue('') +
  facet_wrap(~name, scales = 'free') +
  labs(caption = 'Los valores teoricos corresponden a los calculados con la ecuacion c = sqrt(a^2 + b^2)', x = 'Tipo de malla', y = 'Tamaño (um)') +
  theme_bw()

df_clean %>%
  drop_na() %>%
  filter(g %in% 'Teorico') %>%
  group_by(name) %>%
  summarise(Promedio = mean(value), n = length(g)) %>%
  arrange(match(name, Lname))
 