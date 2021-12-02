rm(list = ls())

options(stringsAsFactors = FALSE, digits = 4, numerals = "allow.loss")


path <- '~/Documents/DOCTORADO/expression_analysis/'

file <- 'Cq'

x <- list.files(path, pattern = file, full.names = TRUE)

### Functions ----

# calculate the % of confidence intervals

IC <- function(x, conf = 0.99) {
  a <- mean(x)
  sd <- sd(x)
  n <- length(x)
  err <- qnorm(p = conf, mean = a, sd = sd, lower.tail = F) * sd / sqrt(n)

  return(err)
}

# IC(c(13.01,13.06,13.21)) # 0.075


# df <- lapply(file, function(x) read.table(x, sep = ',', header = T))
df <- read.table(x, sep = ',', header = T)

library(tidyverse)

# PROMEDIO	DESVESTA	INT DE CONF	MENOR	MAYOR
# = CONFIDENCE(0.01,E3,C3:C5)

df %>%
  group_by(Sample) %>%
  summarise(
    sd = sd(Cq), 
    a = mean(Cq),
    IC = IC(Cq),
    min = a - (IC* sd), max = a + (IC * sd)) -> df_stats

df_stats %>% arrange(desc(a)) %>% pull(Sample) -> Lsample

df %>%
  mutate(Sample = factor(Sample, levels = Lsample)) %>%
  ggplot(aes(Cq, Sample)) +
  geom_point(alpha = 0.5) +
  # stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(0.6)) +
  # geom_boxplot(width = 0.3, 
    # position = position_dodge(0.6),
    # outlier.alpha = 0) +
  stat_summary(fun = function(x) {c(mean(x) - (IC(x)*sd(x)), mean(x) + (IC(x)*sd(x)))}, 
    geom = "line", alpha = 0.5, #shape=23, 
    size=1, position = position_dodge(0.6))
  
df_stats %>%
  mutate(Sample = factor(Sample, levels = Lsample)) %>%
  ggplot(aes(y = a, x = Sample)) + 
  geom_bar(stat="identity") + geom_crossbar(aes(ymin = min, ymax = max), width = 0.02) + 
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.2)



# =IF(AND(C3>=G3,C3<=H3),C3,"fuera_rango")
