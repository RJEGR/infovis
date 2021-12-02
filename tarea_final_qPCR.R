### Ref
# CAT	Catalasa		
# GGT1	Glutatión hidrolasa		
# NATD1	N-acetiltransferasa Dominio 1		

### Targets
# TUBA	Tubulina α		
# RPE	Ribulosa Fosfato Epimerasa		


rm(list = ls())

options(stringsAsFactors = FALSE, digits = 4, numerals = "allow.loss")

IC <- function(x, conf = 0.99) {
  a <- mean(x)
  sd <- sd(x)
  n <- length(x)
  err <- qnorm(p = conf, mean = a, sd = sd, lower.tail = F) * sd / sqrt(n)
  
  return(err)
}

path <- '~/Documents/DOCTORADO/expression_analysis'

file <- '2-CAT-GGT1.csv'

x <- list.files(path, pattern = file, full.names = TRUE)

df <- read.table(x, sep = '\t', header = T)

library(tidyverse)

df %>% as_tibble()

df %>% pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'Cq') -> df_longer


# Seleccion de los genes de refeencia mas estbles ----

## reffinder

file <- 'ranking.csv'

x <- list.files(path, pattern = file, full.names = TRUE)

rdf <- read.table(x, sep = '\t', header = T)

rdf %>% pivot_longer(-HKG, names_to = 'Método', values_to = 'Ranking') %>%
  filter(!Método %in% 'GM') %>%
  ggplot(aes_string(x = 'HKG', y = 'Ranking')) +
  geom_col() +
  geom_text(aes(label = Ranking), nudge_y = 0.1,
    family = "GillSans") +
  coord_flip() +
  facet_grid( Método ~., scales = 'free_y') +
  theme_bw(base_family = "GillSans", base_size = 14)

# rdf %>% pivot_longer(-HKG, names_to = 'Método', values_to = 'Ranking') %>%
#   filter(!Método %in% 'GM') %>%
#   ggplot(aes_string(x = 'HKG', y = 'Ranking', fill = 'Método')) +
#   geom_col(position = "dodge") +
#   theme_bw(base_family = "GillSans", base_size = 14)

# For best keeper ----

# df_longer %>%
#   group_by(Lagunas, Genes) %>%
#   summarise(min = min(Cq), max = max(Cq), 
#     a = mean(Cq), sd = sd(Cq), IC = IC(Cq),
#     upper = a+IC, lower = a-IC,
#     zupper = a+(3*sd), zlower = a-(3*sd), n = n())

# test if parametric
# no pasan la prueba:
# Laguna Cármen, RPE +
# Laguna Mecoacan, CAT
# Cero valores de alto ruido para todos los factores (ie. lagunas y genes)

source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

df %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'x') %>% 
  # group_by(Lagunas) %>%
  group_by(Genes) %>%
  summarise(n= sum(x > 0),
    Zmax = max(Ztransform(x)),
    Zmin = min(Ztransform(x)),
    rCal = Rcalculate_coeff(x),
    rCri = Rcritical_coeff(n),
    gaussian = ifelse(rCal > rCri, TRUE, FALSE),
    outliers = sum(abs(Ztransform(x) > 3)))

library(ggpubr)

# df_longer %>% 
#   # mutate(group = forcats::fct_reorder(Lagunas, Cq)) %>%
#   ggqqplot(x = "Cq", color = 'Genes', add.params = list(linetype = "dashed"),
#     conf.int.level = 0.95) +
#   theme_classic(base_family = "GillSans", base_size = 14) +
#   facet_wrap(~Lagunas, scales = "free_y")  +
#   theme(legend.position = "top") +
#   labs(x = "Expected", y = expression("Observed-Cq"), # ~Log[2]~(x~+1)
#     caption = "Correlation Coefficient (Pearson) between theorical normal distribution (expected) against (observed) distribution in the count data.\nThe observed distribution fit a linear relatioship to a normal distribution (alfa 5%)") 



# en base al rankeo, sabemos que TUBA (1), RPE (2) y NATD1 (3)
# calculate deltas (caldo de la bruja) 
# usando los primeros dos genes

# deltaq <- function(x) {mean(x) - x} # no sirve

deltadf <- df[-1]
deltadf <- colMeans(deltadf) - deltadf


# valor RQ 
# Elevamos el valor de eficiencia teorico al valor de los delta

Eff <- 2.00 # 100%

rqdf <- Eff^deltadf

gmean <- function(x) {exp(mean(log(x[x>0])))} # funcion de media geom.

refdf <- rqdf[c('TUBA', 'RPE')] # Genes de referencia estables

GM <- apply(refdf, 1, gmean) # Media geometrica de los genes de ref estables

redf <-  rqdf[c('CAT','GGT1')] # Genes blancos

out <- cbind(df[1], redf/GM)

# plot(apply(refdf, 1, gmean), apply(refdf, 1, mean))
# cor(apply(refdf, 1, gmean), apply(refdf, 1, mean))

out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  group_by(Lagunas, Genes) %>%
  mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
  mutate(RE = log10(RE+1)) %>%
  summarise(min = min(RE), max = max(RE), 
    a = mean(RE), sd = sd(RE), IC = IC(RE),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) -> out_stats

# RE Is normal???
# Not!!

out %>%
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'x') %>% 
  group_by(Genes, Lagunas) %>% # # Hacemos caldo de la bruja para los genes
  mutate(x = log10(x+1)) %>%
  ggplot(aes(x)) +
  geom_density(alpha = 0.5) +
  facet_grid(Lagunas~ Genes, scales = 'free')
 
out %>%
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'x') %>% 
  group_by(Genes, Lagunas) %>% # # Hacemos caldo de la bruja para los genes
  mutate(x = log10(x+1)) %>%
  summarise(n= sum(x > 0),
    Zmax = max(Ztransform(x)),
    Zmin = min(Ztransform(x)),
    rCal = Rcalculate_coeff(x),
    rCri = Rcritical_coeff(n),
    gaussian = ifelse(rCal > rCri, TRUE, FALSE),
    outliers = sum(abs(Ztransform(x) > 3)))


# out_stats %>%
#   ggplot(aes(y = a, x = Lagunas, fill = Genes)) + 
#   geom_bar(stat="identity") + 
#   geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.02) + 
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
#   facet_grid(~Genes)

# out %>% 
#   pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
#   group_by(Lagunas, Genes) %>%
#   summarise(a = mean(RE))

# Plot qqplot
library(rstatix)

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

out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'x') %>%
  mutate(x = log10(x+1)) %>%
  group_by(Genes) %>%
  summarise(qqfun(x)) %>%
  mutate(z = zfun(y)) %>%
  mutate(outlier = ifelse(abs(z)>3, TRUE, FALSE)) -> outliersdf

outliersdf %>%
  ggplot(aes(x, y, group = Genes, color = outlier)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  geom_point(size = 2.5, alpha = 0.8) + # , label.y = 2.5
  facet_wrap(~Genes) +
  ggpubr::stat_cor(aes(group = Genes), method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = "Expected", 
    y = expression("Observed "~Log[10]~(x~+1)),
    caption = "") +
  scale_color_manual(name = expression("Outlier-"~sigma), 
    values = c('black', 'red')) +
  theme(legend.position = "none")
# Homocelasticidad no es tan necesario probar pues ya se comprobo que no son normaleslas distribuciones, pero existen prueba basadas para datos no normales

# https://rpubs.com/Joaquin_AR/218466
# Fligner-Killeen’s test: a non-parametric test which is very robust against departures from normality.


out %>%
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>% 
  group_by(Lagunas) %>% # # Hacemos caldo de la bruja para los genes
  mutate(RE = log10(RE+1)) %>%
  fligner.test(RE ~ Genes, data = .)
  
# The test reveals a p-value lower than 0.05, indicating that there is  significant difference between the group variances in location.

# out %>% 
# pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
# mutate(LogRE = log10(RE+1)) %>%
# ggqqplot(x = "LogRE", color = 'Genes', add.params = list(linetype = "dashed"),
#   conf.int.level = 0.95) +
# theme_classic(base_family = "GillSans", base_size = 14) +
# facet_wrap(~Lagunas, scales = "free_y")  +
# theme(legend.position = "top")


# non parametric test for the RE values ----
# The Wilcoxon rank sum test is a non-parametric alternative to the independent two samples t-test for comparing two independent groups of samples, in the situation where the data are not normally distributed

# Hay diferencias entre los genes blanco a lo largo de las lagunas?
# 
out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  group_by(Genes, Lagunas) %>%
  # summarise(RE = mean(RE)) %>%
  mutate(RE = log10(RE+1)) %>%
  group_by(Lagunas) %>%
  kruskal_test(RE ~ Genes) %>%
  add_significance()

# Hay diferencias entre las lagunas?
# 
out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  group_by(Genes, Lagunas) %>%
  # summarise(RE = mean(RE)) %>%
  mutate(RE = log10(RE+1)) %>%
  group_by(Genes) %>%
  kruskal_test(RE ~ Lagunas) %>%
  add_significance()

# No se rechaza la Ho (A = B = C) y se propone la prueba a posteriori de contraste multiple w-w p/ muetras independientes (Ho -> Ta = Tb).

# Comparando las lagunas que presentan diferencias segun la prueba a priori":
comparisons <- list(c("Laguna Mecoacan", "contrast"), 
  c("Exposición Hidrocarburos", "contrast"),
  c("Laguna Mecoacan", "Exposición Hidrocarburos"))

# Comparando todas las lagunas vs el experimental:

# g <- c('Laguna Mecoacan', 'Exposición Hidrocarburos')
g <- c('Exposición Hidrocarburos')
comparisons <- list(c("contrast", "Exposición Hidrocarburos"))

out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  mutate(contrast = ifelse(!Lagunas %in% g , 'contrast', Lagunas)) %>%
  mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
  mutate(RE = log10(RE+1)) %>%
  group_by(Genes) %>%
  pairwise_wilcox_test(RE ~ contrast, paired = F, conf.level = 0.95, 
    comparisons = comparisons) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() -> stat.test

stat.test %>% view()

# comparando todas las lagunas y experimental, ninguna es significativa

out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  # mutate(contrast = ifelse(!Lagunas %in% g , 'contrast', Lagunas)) %>%
  mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
  mutate(RE = log10(RE+1)) %>%
  group_by(Genes) %>%
  pairwise_wilcox_test(RE ~ Lagunas, paired = F, conf.level = 0.95) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() -> stat.test

# out %>% 
#   pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
#   # mutate(contrast = ifelse(!Lagunas %in% g , 'contrast', Lagunas)) %>%
#   mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
#   mutate(RE = log10(RE+1)) %>%
#   group_by(Genes) %>%
#   pairwise_wilcox_test(RE ~ Lagunas, paired = F, conf.level = 0.95) %>%
#   add_significance() %>% view()
# # stat.test %>% view()
# 
# 
out %>%
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  ggplot(aes(log10(RE+1), Lagunas, color = Genes)) +
  geom_point(alpha = 0.5, position = position_dodge(0.6)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(width = 0.3, position = position_dodge(0.6),
    outlier.alpha = 0) +
  stat_summary(fun=mean, geom="point", shape=23,
    size=1, position = position_dodge(0.6)) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  theme(legend.position = "top") +
  labs(y = '', x = expression("RE (" ~Log[10]~")"))

library("ggpubr")
# 
# out %>% 
#   pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
#   mutate(RE = log10(RE+1)) %>%
#   # filter(RE < 3) %>%
#   mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
#   # mutate(Lagunas = gsub("\n$", "", Lagunas)) %>%
#   # mutate(Lagunas = factor(Lagunas, levels = LLagunas)) %>%
#   ggline(x = "Lagunas", y = "RE", 
#     color = 'Genes', 
#   add = c("mean_se", "jitter"), 
#   # order = c("ctrl", "trt1", "trt2"),
#   ylab = expression("RE (" ~Log[10]~")"), xlab = "") +
#   # facet_grid( Genes ~., scales = 'free_y') +
#   theme_bw(base_family = "GillSans", base_size = 14) +
#   theme(legend.position = 'top',
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12))

# Bar plots showing mean +/- SD

# out %>% 
#   pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
#   mutate(RE = log10(RE+1)) %>%
#   filter(RE < 3) %>%
#   mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
#   ggplot(aes(x = Lagunas, y = RE)) +
#   geom_boxplot(outlier.size = 0, alpha=0.2) + # outlier.shape=NA
#   # ylim(0, 20) +
#   theme_bw(base_family = "GillSans", base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
#   facet_grid( Genes ~., scales = 'free_y') 


out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  group_by(Lagunas, Genes) %>%
  mutate(RE = log10(RE+1)) %>%
  # mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
  filter(RE < 3) %>%
  summarise(min = min(RE), max = max(RE), 
    a = mean(RE), sd = sd(RE), IC = IC(RE),
    upper = a+IC, lower = a-IC,
    zupper = a+(3*sd), zlower = a-(3*sd), n = n()) %>%
  mutate(RE = a ) -> out_stats

# out_stats %>% distinct(Lagunas) -> LLagunas

out_stats %>%
  mutate(contrast = ifelse(!Lagunas %in% g , 'contrast', Lagunas)) %>%
  mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
  ggplot(aes(y = RE, x = Lagunas, fill = contrast)) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  facet_grid( Genes ~., scales = 'free_y') +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(y = expression("RE (" ~Log[10]~")"), x = "") +
  scale_fill_manual(values = c('black', 'Blue'))-> p

# 
# out_stats %>%
#   mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
#   ggplot(aes(y = RE, x = Lagunas, fill = Genes)) + 
#   theme_bw(base_family = "GillSans", base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, 
#     hjust = 1, vjust = 1, size = 10)) +
#   geom_bar(position="dodge", stat = "identity") + 
#   geom_errorbar(aes(ymin = lower, ymax = upper),
#     position = position_dodge(), colour="black")



# Add the p-value manually

out %>% 
  pivot_longer(-Lagunas, names_to = 'Genes', values_to = 'RE') %>%
  mutate(Lagunas = gsub(" ", "\n", Lagunas)) %>%
  mutate(RE = log10(RE+1)) %>%
  group_by(Genes) %>%
  pairwise_wilcox_test(RE ~ Lagunas, paired = F, conf.level = 0.95) %>%
  add_significance() -> stat.test

# stat.test %>% filter(grepl("Exposición", group1)) -> stats

stat.test %>% filter(grepl("^[*]", p.adj.signif)) -> stats

stats <- stats %>% add_xy_position(x = "Lagunas")

p + stat_pvalue_manual(stats, label = "p.adj.signif", 
  remove.bracket = F)

p
