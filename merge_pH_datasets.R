rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)

path <- '~/Documents/DOCTORADO/pH_measures/2022_/'

namesL <- c('Canal-1', 'Canal-2', 'Canal-3', 'Canal-4')

# RColorBrewer::display.brewer.pal(12, 'Paired')

getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

getPalette <- structure(getPalette, names = namesL) 


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


# 
read.files <- function(path, pattern_f = '.dat$') {
  
  # Final Output generated is data.frame with cols:
  # [1] "date"    "Canal-1" "Canal-2" "Canal-3" "Canal-4"
  
  
  file <- list.files(path, pattern = pattern_f,  full.names = TRUE)
  
  df <- lapply(file, read.table)
  
  head(df <- do.call(rbind, df))
  
  n1 <- ncol(df) - 3
  
  head(out <- cbind(df[1], df[n1:ncol(df)]))
  
  names(out) <- c('date', namesL)
  
  return(out)
}

dff1 <- read.files(path)

head(dff1)

ocount <- function(df, zscore) {
  
  df %>%
    pivot_longer(cols = all_of(namesL), values_to = 'Obs') %>%
    drop_na(Obs) %>%
    group_by(name) %>%
    summarise(qqfun(Obs)) %>%
    mutate(z = zfun(y)) -> odf
  
  out <- odf %>%
    mutate(outlier = ifelse(abs(z) >= zscore, 'Y', 'N')) %>%
    count(name, outlier)
  
  out %>%
    pivot_wider(values_from = n, 
      names_from = 'outlier') %>%
    mutate(g = paste0(zscore, 'sc'))

}

ocount(dff1, 2) 

l <- list()

for(i in 5:1) {
  l[[i]] <- ocount(dff1, i)
}


subtitle <- expression('Outlier detection by the z Standard Deviations from the Mean')
caption <- expression('z = (x - mean(x)) / sd(x)')


do.call(rbind,l) %>%
  ungroup() %>%
  mutate(total = N+Y, pct = 100*(Y/total)) %>%
  # filter(name %in% c('Canal-3', 'Canal-4')) %>%
  ggplot(aes(x = g, y = pct, color = name)) +
  geom_line(aes(group = name), orientation = "x", linetype = 'dashed') +
  geom_point() +
  labs(x ='Z score', y = '% Outliers', caption = caption) +
  scale_color_manual(values = getPalette) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.grid = element_line(size = rel(0.5)),
    panel.grid.minor = element_line(size = rel(0)),
    panel.border = element_blank()) +
  guides(colour = guide_legend("")) -> psave

# psave

dff1 %>%
  pivot_longer(cols = namesL, values_to = 'Obs') %>%
  drop_na(Obs) %>%
  group_by(name, date) %>%
  summarise(qqfun(Obs)) %>%
  mutate(z = zfun(y)) %>%
  mutate(outlier = ifelse(abs(z) >= 3, 'Y', 'N')) -> df_summ

# to hard to plot !!!

# df_summ %>%
#   ggplot(aes(x, z)) +
#   geom_rug() +
#   geom_point(size = 0.2, alpha = 0.5) + 
#   theme_bw(base_family = "GillSans", base_size = 14) 


# do.call(rbind,l) %>% view()

dff1 %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = all_of(namesL), values_to = 'Obs') %>%
  drop_na(Obs) -> df_longer

df_longer %>%
  drop_na(Obs) %>%
  group_by(date, name) %>%
  mutate(value = Obs) %>%
  summarise(
    a = mean(value), sd = sd(value), IC = IC(value),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) -> out_stats

out_stats %>%
  # filter(name %in% 'Canal-2') %>%
  ggplot(aes(y = a, x = date, color = name)) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  # facet_grid(.~Buffer, scales = 'free_y') +
  geom_point() +
  geom_line(aes(group = name), orientation = "x", 
    linetype = 'dashed') +
  # geom_bar(stat="identity", fill = 'white', color = 'grey67') + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_color_manual('', values = getPalette)
  
# bind (ie. rbind data from larvae and poslarvae) 

path2 <- '~/Documents/DOCTORADO/pH_measures/2022_/asentamiento/'

dff2 <- read.files(path2)

df <- rbind(dff1, dff2)

df %>% separate(col = date, 
  into = c('D', 'M', 'Y'), sep = "/") %>%
  mutate(date = paste(M,D,Y, sep = '/')) %>%
  select(all_of(names(df))) -> df

df %>% arrange(date) %>% distinct(date) %>% pull(date) -> dateL

df %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = all_of(namesL), values_to = 'Obs') %>%
  drop_na(Obs) %>%
  mutate(date = factor(date, levels = dateL)) %>% 
  filter(name %in% c('Canal-3', 'Canal-4')) -> df_longer

df_longer

# l <- list()
# 
# for(i in 5:1) {
#   l[[i]] <- ocount(df, i)
# }
# 
# do.call(rbind,l) %>%
#   ungroup() %>%
#   mutate(total = N+Y, pct = 100*(Y/total)) %>%
#   # filter(name %in% c('Canal-3', 'Canal-4')) %>%
#   ggplot(aes(x = g, y = pct, color = name)) +
#   geom_line(aes(group = name), orientation = "x", linetype = 'dashed') +
#   geom_point() +
#   labs(x ='Z score', y = '% Outliers', caption = caption) +
#   scale_color_manual(values = getPalette) +
#   theme_bw(base_size = 10, base_family = "GillSans") +
#   theme(
#     legend.position = 'top',
#     legend.text = element_text(size = 7),
#     strip.background = element_blank(),
#     panel.grid = element_line(size = rel(0.5)),
#     panel.grid.minor = element_line(size = rel(0)),
#     panel.border = element_blank()) +
#   guides(colour = guide_legend("")) -> psave

df_longer %>%
  group_by(name, date) %>%
  summarise(qqfun(Obs)) %>%
  mutate(z = zfun(y)) %>%
  mutate(z = ifelse(abs(z) >= 3, NA, z)) %>%
  drop_na(z) -> df_summ

# No vale la pena hacer un scatterplot de todos los datos, pues son casi 2Kmillones. Asi que hacemos un summary

df_summ %>%
  group_by(date, name) %>%
  mutate(value = y) %>%
  summarise(
    a = mean(value), sd = sd(value), IC = IC(value),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) -> out_stats



# Include control data


df_longer_c <- readRDS(paste0(path2, '/df_longer_control.rds'))

df_longer_c %>%
  filter(name %in% 'Canal-2') %>%
  mutate(date = V1) %>%
  group_by(date, name) %>%
  mutate(value = Obs) %>%
  summarise(
    a = mean(value), sd = sd(value), IC = IC(value),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) %>%
  separate(col = date, 
    into = c('D', 'M', 'Y'), sep = "/") %>%
  mutate(date = paste(M,D,Y, sep = '/')) %>%
  select(all_of(names(out_stats))) -> out_stats2


dviz <- rbind(out_stats, out_stats2)

distinct(dviz, date) %>% pull(date) -> recode_g

struc_group <- c('Embryo', rep('Larvae', 5), rep('Settlement', 26))

level_key <- structure(struc_group, names = recode_g)

dviz %>% mutate(g = recode_factor(date, !!!level_key)) -> dviz

Labels <- c('Control', 'Experimental I', 'Experimental II')
  
dviz %>%
  # filter(name %in% 'Canal-2') %>%
  mutate(date = factor(date, levels = dateL)) %>%
  ggplot(aes(y = a, x = date, color = name)) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  # geom_point(aes(size = n), alpha = 0.5) +
  geom_line(aes(group = name), orientation = "x", linetype = 'dashed') +
  scale_color_manual('', values = getPalette,
    limits = namesL[-1], labels = Labels) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  scale_y_continuous(n.breaks = 10) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  # facet_grid(.~g, scales = 'free_x', space = 'free', switch = 'x') +
  # theme(legend.position = 'top',
  #   panel.border = element_blank(), 
  #   panel.spacing.x = unit(0,"line"),
  #   strip.text = element_text(size = 7),
  #   strip.background = element_rect(linetype ='dashed', colour="black",fill="white"),
  #   axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = 'top',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) +
  labs(y = 'Average pH', x = '') -> psave

# Anadimos un bottom plot del conteo de datos que se usaron por dia

dviz %>% 
  group_by(g, date) %>%
  summarise(l = length(name), N = sum(n), n = N/l) %>%
  # mutate(pct = n / sum(n)) %>%
  mutate(date = factor(date, levels = dateL)) %>%
  ggplot() +
  geom_point(aes(size = n, x = date, y = as.factor(1))) +
  # geom_col(aes(y = n, x = date), width = 0.5) +
  scale_fill_manual('', values = getPalette,
    limits = namesL[-1], labels = Labels) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  facet_grid(.~g, scales = 'free_x', space = 'free', switch = 'x') +
  theme(legend.position = 'bottom',
    axis.line.x = element_blank(),
    panel.border = element_blank(), 
    panel.spacing.x = unit(0,"line"),
    panel.spacing.y = unit(0,"line"),
    strip.text = element_text(size = 7),
    strip.background = element_rect(size = 0.7, linetype = 'dashed', 
      colour="black",fill="white"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  # theme(legend.position = 'none',
  #   axis.text.x = element_blank(),
  #   axis.ticks.x = element_blank()) +
  labs(y = 'Size', x = '') -> psave2
  # scale_y_reverse() 

library(patchwork) 

psave/psave2 + patchwork::plot_layout(heights = c(4,1)) -> pout

pout


ggsave(pout, 
  filename = 'pH_times_series_average.png', 
  path = path, width = 7, height = 5)

dviz %>%
  mutate(date = factor(date, levels = dateL)) %>%
  mutate(y = upper - lower) %>%
  ggplot(aes(y = sd, x = n)) + 
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = TRUE, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  geom_point() +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  scale_color_manual('', values = getPalette,
    limits = namesL[-1], labels = Labels) +
  scale_y_continuous(n.breaks = 10) 
  # facet_grid(~ g, scales = 'free_x')


# Posteriori test ----
# Buscamos diferencias entre las medias de pH en los sistemas. La intencion es evaluar si hubo o no diferencias entre los dias y/o en que grado. Una vez hecha la prueba, hacemos un promedio por sistema siempre y cuando no haya dif significativas entre las medias.

library(rstatix)

# Identify outliers
dviz %>% group_by(name) %>% identify_outliers(a)

# The normality assumption can be checked by computing the Shapiro-Wilk test. If the data is normally distributed, the p-value should be greater than 0.05.

dviz %>% group_by(name) %>% shapiro_test(a) %>% adjust_pvalue()

# A partir de la salida, si el valor p es mayor que el nivel de significación 0,05, indica que la distribución de los datos no es significativamente diferente de la distribución normal. En otras palabras, podemos asumir la normalidad.

# En vista de que encontramos outliers y no normalidad, evaluamos a traves de un test de wilcoxon

dviz %>%
  ungroup() %>%
  # t_test(a ~ name)
  pairwise_wilcox_test(a ~ name,  conf.level = 0.95) %>%
  adjust_pvalue() %>%
  add_significance("p") -> stats.test

stats.test %>% add_xy_position(x = "name") -> stats

# stats$y.position <- max(out_stats$a)+out_stats$sd

subtitle <- get_pwc_label(stats, "expression")

# and plot

dviz %>%
  ungroup() %>%
  mutate(name = factor(name, levels = c('Canal-3', 'Canal-4', 'Canal-2'))) %>%
  ggplot(aes(y = a, x = as.factor(name))) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  stat_boxplot(geom ='errorbar', width = 0.07) +
  geom_boxplot(width = 0.3, outlier.alpha = 0) +
  geom_point(aes(size = n), alpha = 0.5) +
  stat_summary(fun=mean, geom="point", shape=23, size=1, color = 'red') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(legend.position = 'bottom') + 
  scale_y_continuous(n.breaks = 10) +
  scale_x_discrete(labels = c("Experimental I", "Experimental II", 'Control')) +
  scale_size('Set Size', range = c(0,5)) +
  labs(y = 'Average pH', x = '') -> psave

# psave + coord_flip()  -> psave

psave + theme(panel.border = element_blank()) -> psave

psave + 
  ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
    remove.bracket = F, tip.length = 0.01, linetype = 'dashed') +
  labs(subtitle = subtitle) -> psave

psave

dviz %>%
  group_by(name) %>%
  summarise(a = mean(a))

ggsave(psave, filename = 'average_pH_boxplot.png', path = path, width = 5, height = 6)

# analyze chemistry ----

path <- '~/Documents/DOCTORADO/pH_measures/'

head(df <- read.csv(paste0(path, 'quimicas.csv')))

ylabs = expression(CO[3]^{-2}~(µmol~Kg^{-1}))
xlabs = expression(Omega["ara"])

# Lnames <- c(expression(Omega["ara"]),
#   expression(CO[3]^{-2}~(µmol~Kg^{-1})),
#   expression(pCO[2]~(uatm)),
#   expression(O[2]~(µmol~Kg^{-1})), "pH")


df %>%
  select(ID, c('CO3','DIC', 'TA', 'Ara')) %>%
  separate(col = ID, 
    into = c('hpf', 'pH', 'g'), sep = "-") %>%
  pivot_longer(cols = c('DIC', 'TA', 'Ara')) %>%
  ggplot(aes(value, `CO3`)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = F, na.rm = TRUE) +
  geom_point() +
  geom_rug() +
  facet_grid(~name, scales = 'free_x') +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'none') +
  labs(y = ylabs, x = '')

df %>%
  ggplot(aes(DIC, TA, group = name)) +
  geom_smooth(method = "lm", linetype="dashed", color = 'black', se = T) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", 
    p.accuracy = 0.001, label.y = 2300) +
  # ggrepel::geom_label_repel(aes(label = ID, color = pH)) +
  geom_point(size = 2.5, alpha = 0.5, aes(color = name)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  scale_color_manual('', values = structure(getPalette[-1], names =  Labels)) +
  facet_grid(~name)

# in comparison to 

dviz %>% 
  ungroup() %>%
  select(date, name, a) %>%
  pivot_wider(names_from = name, values_from = a) 

level_key <- structure(Labels, names = c('Canal-2', 'Canal-3', 'Canal-4'))

# como no hay pH de control al inicio, no vas a tener join con los datos de df
# sum(dviz$date %in% df$date)

dviz %>% ungroup() %>%
  mutate(name = recode_factor(name, !!!level_key)) %>%
  select(date, name, a) %>%
  inner_join(df) %>%
  ggplot(aes(a, DIC)) +
  geom_point(aes(color = name)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001)

# normalizar datos a la salinidad

df %>%
  mutate(nDIC = mean(df$Sal)*(DIC/Sal), nTA = mean(df$Sal)*(TA/Sal)) %>%
  select(name, nTA, nDIC, TA, DIC) %>%
  ggplot(aes(color = name)) +
  # geom_point(aes(nTA, TA))
  geom_point(aes(nDIC, DIC))

boxplot(df$Sal)
mean(df$Sal)



