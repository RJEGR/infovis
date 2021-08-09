
path <- '~/Documents/DOCTORADO/'

library(ggplot2)
library(tidyverse)


x <- read.csv('~/Documents/DOCTORADO/pone.0104371.s007.csv', header = T)

mtd <- read.csv('~/Documents/DOCTORADO/pone.0104371.s007_mtd.csv', header = T)

x <- x[, grepl('RPM', names(x))]

c('Egg', 'Embryo', 'Trochophore', 'D-shape larvae', 
  'Umbo larvae', 'Pediveliger', 'Spat', 'Juvenile',
  "Gills","Digestive gland", "Mantle",
  "Hemocyte", "Labial palp","Female gonad", "Adductor muscle") -> stagesLev

sum(unique(mtd$Stage) %in% stagesLev)

library(ggsci)

# gcolor <- structure(c('#D33F3A', '#45B8DA'), names = c('Tissues', 'Developmental stages'))
# 
# mtd$Stage[mtd$Group %in% "Developmental stages"] -> c1
# 
# colorRampPalette(c(gcolor[2], '#2166ac','#f4a582'), alpha = F)(length(c1)) -> dPalette
# dPalette = structure(dPalette, names = c1) 
# 
# # factor colors for tissues
# mtd$Stage[mtd$Group %in% "Tissues"] -> c2
# 
# colorRampPalette(c(gcolor[1], '#D84F4B','#F1BEBD'), alpha = F)(length(c2)) -> tPalette
# tPalette = structure(tPalette, names = c2) 

x %>% pivot_longer(cols = names(x), names_to = 'ID') %>%
  group_by(ID) %>%
  summarise(n = sum(value)) %>%
  left_join(mtd) %>%
  mutate(Description = factor(Description, levels = mtd$Description)) %>%
  mutate(Stage = factor(Stage, levels = stagesLev)) %>%
  ggplot(aes(y = n, x = Description)) +
  geom_col(aes(fill = Group)) +
  labs(x = '', y = 'Reads Per Million (Precursor reads)') +
  scale_fill_manual('',values = gcolor) +
  theme_classic(base_family = "GillSans", base_size = 16) +
  theme(
    legend.position = 'top',
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 14),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave
  # theme(
  #   legend.position = 'top',
  #   legend.text = element_text(size = 7),
  #   axis.text.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   strip.background = element_blank(),
  #   panel.border = element_blank()) 


ggsave(psave, filename = 'precursor_reads_plot.png', path = path, width = 12, height = 7)
