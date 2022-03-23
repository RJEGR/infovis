# https://csdaw.github.io/ggprism/articles/web-only/palettes.html
# https://csdaw.github.io/ggprism/articles/web-only/palettes.html

rm(list = ls()) # Limpiar la memoria de la sesion de R

options(stringsAsFactors = FALSE) # 

#

# Define Color Palette: ----

# colNames <- c("Aspecto.1",	"Aspecto.2")

# RColorBrewer::display.brewer.all() # Revisamos todas las paletas de colores que tenemos

getPalette <- RColorBrewer::brewer.pal(3, 'Set1')

# Cvalues <- structure(getPalette, names = colNames)

namesL <- c('8.0', '8.0', '7.6', '7.8')


getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

getPalette <- structure(getPalette, names = namesL) 


#

library(tidyverse)


# 

path <- '~/Documents/DOCTORADO/Respirometry/'

file <- '.csv$'

file <- list.files(path, pattern = file, full.names = TRUE)

df <- lapply(file, read.csv) 

df <- do.call(rbind, df)

df %>% as_tibble()

df %>%
  mutate(N = `Aspecto.1`+ `Aspecto.2`) %>%
  mutate(across(c(5:6),
    .fns = ~./N)) %>% 
  pivot_longer(cols = all_of(c('Aspecto.1','Aspecto.2'))) %>%
  drop_na(value) -> df_longer

df_longer %>%
  filter(hpf == 24) %>%
  ggplot(aes(x = ID, y = value*100, fill = name)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual('', values = getPalette) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave

psave + 
  ggh4x::facet_nested(hpf + pH ~., 
    scales = "free", space = "free", switch = "y")

# hacer un barplot con barras de error para 24 hpf unicamente comparando el % de aspectos 

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

# Is gaussian?  ----
# test caldo de la bruja de todos los datos

df_longer %>% 
  mutate(x = value) %>%
  drop_na(x) %>%
  group_by(pH, hpf) %>%
  summarise(n= sum(x > 0),
    Zmax = max(Ztransform(x)),
    Zmin = min(Ztransform(x)),
    rCal = Rcalculate_coeff(x),
    rCri = Rcritical_coeff(n),
    gaussian = ifelse(rCal > rCri, TRUE, FALSE),
    outliers = sum(abs(Ztransform(x) > 3)))

# not! ----

  # facet_grid(hpf~pH)

# There are outliers? ----

subtitle <- expression('Outlier detection by the i Standard Deviations from the Mean')
caption <- expression('Zscore = (x - mean(x)) / sd(x)')

df %>%
  group_by(hpf, pH) %>%
  summarise(qqfun(N)) %>%
  mutate(z = zfun(y)) %>%
  mutate(outlier = ifelse(abs(z)>=3, TRUE, FALSE)) -> outliersdf

outliersdf %>%
  ggplot(aes(x, z, color = outlier)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  geom_point(size = 2.5, alpha = 0.8) + 
  facet_grid(hpf~pH) +
  ggpubr::stat_cor(aes(group = hpf), method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = "Expected", y = expression(Z[score]), caption = caption, subtitle = subtitle) +
  scale_color_manual(name = expression("Outlier-"~sigma), 
    values = c('black', 'red')) +
  theme(legend.position = "none")

outliersdf %>%
  ggplot(aes(x, z, color = outlier)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  geom_point(size = 2.5, alpha = 0.8) + 
  geom_rug(aes(color = outlier), length = unit(0.01, "npc")) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = "Expected", y = expression(Z[score]), caption = caption, subtitle = subtitle) +
  scale_color_manual(name = expression("Outlier-"~sigma), 
    values = c('black', 'red')) +
  theme(legend.position = "none")

# Normality ----
# After remove outliers (if there are) What is the degree of normally are our Observed data (multicolinearity)

# include corplot of Expected vs abundance relationship 
# library(corpcor)
# 
caption = c('Observed vs Expected Pearson Coefficients.\nIt highlight the degree of association between variables.\nTherefore we argue homogeneity through our sampling.')

outliersdf %>%
  filter(outlier == FALSE) %>%
  mutate(pH = as.factor(pH), hpf = as.factor(hpf)) %>%
  group_by(hpf, pH) %>%
  rstatix::cor_test(x,y) -> cor_df

# heatmap or
m <- cor_df %>% 
  select(hpf, pH, cor) %>%
  pivot_wider(values_from = cor, names_from = pH) %>%
  as.data.frame(row.names = hpf)

hc <- hclust(dist(m[-1]), method = "ward.D")

hc$labels <- m$hpf

cor_df$p

cor_df %>%
  # summarise(corr = cor(x,y), cov = cov(x,y)) %>%
  ggplot(aes(pH, hpf, fill = cor)) + 
  geom_tile(color = 'white', size = 0.5) + # aes(alpha = 1-p)
  geom_text(aes(label = cor), color = 'white') +
  theme_classic(base_size = 12, base_family = "GillSans") +
  ggsci::scale_fill_gsea(name="Corr", reverse = T) +
  scale_x_discrete(position = 'top') +
  labs(x = '', caption = caption) +
  ggh4x::scale_y_dendrogram(hclust = hc) 

# barp

cor_df %>%
  ggplot(aes(y = cor, hpf, fill = pH)) +
  geom_bar(stat = "identity", width = 0.7, 
    position = position_dodge(width = 0.75)) +
  theme_classic(base_size = 12, base_family = "GillSans") 

# not ----

# (omit) homocelasticidad or homogeneity of variance across groups

# df_longer %>%
#   # filter(hpf == 24) %>%
#   filter(name %in% 'Aspecto.1') %>%
#   drop_na(value) %>% 
#   levene_test(value ~ as.factor(pH))

# (omit) The p-value of the Levene’s test is significant, suggesting that there is a significant difference between the variances of the two groups. Therefore, we’ll use the Weltch t-test, which doesn’t assume the equality of the two variances.

# Lets test non parametric comparsion between pH across the times

# The Wilcoxon rank sum test is a non-parametric alternative to the independent two samples t-test for comparing two independent groups of samples, in the situation where the data are not normally distributed


# Make a barplot data based on the mean values and include significance differences as previous codes

# By N values

df %>%
  filter(hpf == 24) %>%
  group_by(hpf) %>%
  ggplot(aes(y = N,  x = as.factor(pH))) +
  geom_boxplot() -> psave


stats.test %>% add_xy_position(x = "pH") -> stats

psave + 
  ggpubr::stat_pvalue_manual(stats, label = "p.adj", 
    remove.bracket = F, tip.length = 0) +
  theme_bw(base_family = "GillSans", base_size = 14)

# By pct

# geom_bar(stat = "identity", width = 0.7, 
#   position = position_dodge(width = 0.4)) +
# ggh4x::facet_nested(Phylum  ~wrap, scales = "free", space = "free", switch = "y") +
# coord_flip() +
# geom_errorbar(aes(ymin = logFC - SE, 
#   ymax = logFC + SE), width = 0.2,
#   position = position_dodge(0.05), color = "black") + 
# geom_text(aes(y = logFC+2.5*sign(logFC), label=star), 
#   vjust=.7, color="black", position=position_dodge(width = .5)) +
# ggsci::scale_fill_uchicago()

# nos interesa comparar el  % de larvas que, por tratamiento, se observo un mayor/menor % competencia entre una fase y otra. Para ello, consideramos el Aspecto 1, a las 24 hpf, donde se evalue el numero de individuos presentes, en cada replica, que presentaban caracteristicas de larvas trocoforas. Los valores de Aspecto 1 son normalizados respecto al total de los individuos (N) y aplico prueba a priori y posteriori para considerar diferencias entre los tratamientos. 


# Asi que,

df_longer %>% filter(hpf == 24) %>% filter(name %in% 'Aspecto.1') -> df4viz

df4viz <- df_longer %>% filter(name %in% 'Aspecto.1')

df4viz %>%
  drop_na(value) %>%
  group_by(hpf, name, pH) %>%
  summarise(
    a = mean(value), sd = sd(value), IC = IC(value),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) -> out_stats


out_stats %>%
  mutate(name = gsub("[.]", " ", name)) %>%
  mutate(hpf = as.factor(hpf)) %>%
  ggplot(aes(y = a, x = as.factor(pH), fill = hpf)) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  facet_grid(.~ hpf, scales = 'free_y') +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  ylim(0,1) +
  # labs(y = expression("RE (" ~Log[10]~")"), x = "") +
  scale_color_manual(values = c('black', 'Blue')) +
  scale_fill_brewer() +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(legend.position = 'none') -> psave

psave

# En vez del valor N debe ser pct

df4viz %>%
  group_by(hpf) %>%
  mutate(pH = as.factor(pH)) %>%
  pairwise_wilcox_test(value ~ pH, 
    # ref.group = "8", 
    conf.level = 0.95) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p") -> stats.test

# stats.test %>% view()

stats.test %>% add_xy_position(x = "pH", group = '8') -> stats

stats$y.position <- max(out_stats$a)+out_stats$sd

subtitle <- get_pwc_label(stats)

psave + 
  ggpubr::stat_pvalue_manual(stats, hide.ns = T, label = "p.adj") +
  labs(subtitle = subtitle, y = 'Average', x = 'pH')

stats %>% view()
