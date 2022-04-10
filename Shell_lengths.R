# Morfoemtry

rm(list = ls())

options(stringsAsFactors = FALSE)

# Prepare functions and libraries

library(tidyverse)
library(rstatix)

path <- '~/Documents/DOCTORADO/Shell_length/'

pattern_f <- 'csv'


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

IC <- function(x, conf = 0.99) {
  a <- mean(x)
  sd <- sd(x)
  n <- length(x)
  err <- qnorm(p = conf, mean = a, sd = sd, lower.tail = F) * sd / sqrt(n)
  
  return(err)
}

# Load data ----

col_names <- c("Days","hpf","pH","Ind","group","L1","L2")

file <- list.files(path, pattern = pattern_f,  full.names = TRUE)

df <- lapply(file, read.csv)

head(df <- do.call(rbind, df))

df <- df %>% as_tibble() %>% select(col_names)

# Pal color

namesL <- c('8.0', '8.0', '7.6', '7.8')

getPalette <- RColorBrewer::brewer.pal(4, 'Paired')

getPalette <- structure(getPalette, names = namesL) 

# There are outliers? ----

subtitle <- expression('Outlier detection by the i Standard Deviations from the Mean')
caption <- expression('Zscore = (x - mean(x)) / sd(x)')


df %>%
  group_by(hpf, pH) %>%
  summarise(qqfun(L1)) %>%
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
  labs(x = "Theorical", y = expression(Z[score]), caption = caption, subtitle = subtitle) +
  scale_color_manual(name = expression("Outlier-"~sigma), 
    values = c('black', 'red')) +
  theme(legend.position = "none")

# Normally test :::: ----
# Visual inspection, described in the previous section, is usually unreliable. It’s possible to use a significance test comparing the sample distribution to a normal one in order to ascertain whether data show or not a serious deviation from normality. There are several methods for normality test such as Kolmogorov-Smirnov (K-S) normality test and Shapiro-Wilk’s test. Note that, normality test is sensitive to sample size. Small samples most often pass normality tests. Therefore, it’s important to combine visual inspection and significance test in order to take the right decision.

# From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from a normal distribution. In other words, we can assume the normality.

outliersdf %>% rstatix::shapiro_test(x) # Test it into the normal set
outliersdf %>% rstatix::shapiro_test(y) # Test it into the real sample set

# we processed to continue w/ a non parametric analysis.

# Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends the two-samples Wilcoxon test in the situation where there are more than two groups. It’s recommended when the assumptions of one-way ANOVA test are not met. T

# We want to know if there is any significant difference between the average lengths of larvae in the 3 pH conditions.

df %>%
  group_by(hpf) %>%
  rstatix::kruskal_test(L1 ~ pH) %>%
  adjust_pvalue(method = "none") %>%
  add_significance("p") -> kruskal.stats


# Interpretation: As the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the treatment groups.

df %>% ggplot(aes(y = L1, x = as.factor(pH))) + 
  # geom_boxplot(outlier.alpha = 0) +
  geom_point(alpha = 0.3) +
  stat_boxplot(geom ='errorbar', 
    width = 0.2, position = position_dodge(0.6)) -> psave

psave + labs(subtitle = get_test_label(kruskal.stats, detailed = TRUE))

# From the output of the Kruskal-Wallis test, we know that there is a significant difference between groups, but we don’t know which pairs of groups are different, so:

# Additionally, The Wilcoxon rank sum test is a non-parametric alternative to the independent two samples t-test for comparing two independent groups of samples, in the situation where the data are not normally distributed

df %>%
  group_by(hpf) %>%
  pairwise_wilcox_test(L1 ~ pH, 
    # ref.group = "8", 
    conf.level = 0.95) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p") -> stats.test


stats.test %>% add_xy_position() -> stats

stats$y.position <- stats$y.position+out_stats$sd

title <- get_pwc_label(stats)
subtitle <- get_test_label(kruskal.stats, detailed = TRUE)

psave + 
  ggpubr::stat_pvalue_manual(stats, hide.ns = T, tip.length = 0.001) +
  labs(title = title, subtitle = subtitle, 
    y = 'Shell length', x = 'pH') +
  theme_classic(base_size = 10, base_family = "GillSans")
  



# Despues de tener todos los datos de hpf (ie. 24, 38, 60, 108 y 26 dpa) hacer el siguiente grafico, de lineas

df %>%
  group_by(hpf, pH) %>%
  summarise(
    a = mean(L1), sd = sd(L1), IC = IC(L1),
    upper = a+sd, lower = a-sd, n = n()) -> out_stats

L <- df$L1
m <- round(min(L))
M <- round(max(L))

out_stats %>%
  ggplot(aes(y = a, x = as.factor(pH), fill = as.factor(pH))) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  # geom_line(aes(group = pH), orientation = "x", linetype = 'dashed') +
  # geom_point() +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) 
# scale_fill_manual('', values = getPalette[-3]) +
# scale_y_continuous(breaks = seq(m,M, by = 100)) 