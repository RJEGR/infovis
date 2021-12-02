
path <- '~/Documents/DOCTORADO/'

library(ggplot2)
library(tidyverse)


x <- read.csv('~/Documents/DOCTORADO/pone.0104371.s007.csv', header = T)

miR_df <- x[, c(1,2)]

mtd <- read.csv('~/Documents/DOCTORADO/pone.0104371.s007_mtd.csv', header = T)

x <- x[, grepl('RPM', names(x))]


pcadf <- cbind(miR_df, x)

prevelancedf = apply(X = x,
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

df <- data.frame(Prevalence = prevelancedf, 
                TotalAbundance = rowSums(x),
                miR_df) %>% 
  as_tibble(rownames = 'id') %>%
  mutate(arm = ifelse(grepl('-5p$', mature), '5p', '3p'))

# con la siguiente figura analizamos la prevalencia de mirnas en todas las muestras, tanto de desarrollo como de tejido, y revisamos que mirna puede ser el usado como valor de normalizacion para todas las muestras: 

df %>%
  filter(TotalAbundance > 0) %>%
  mutate(mature = ifelse(TotalAbundance > quantile(TotalAbundance, probs = 0.95 ), mature, ''),
         Prevalence = Prevalence/16 * 100) %>%
  ggplot(aes(Prevalence, log10(TotalAbundance + 1))) +
  geom_point(alpha = 0.7)  +# aes(color = arm)
  theme_classic(base_family = "GillSans", base_size = 16) +
  labs(x = "Prevalence (Frac. Samples)", 
       y = "Total Reads (log10 + 1)") + # color = "arm"
  ggrepel::geom_text_repel(aes(label = mature), 
                           max.overlaps = 30, size = 2.5) +
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_blank())

c('Egg', 'Embryo', 'Trochophore', 'D-shape larvae', 
  'Umbo larvae', 'Pediveliger', 'Spat', 'Juvenile',
  "Gills","Labial palp", "Female gonad", "Hemocyte",
  "Digestive gland", "Mantle", "Adductor muscle") -> stagesLev

sum(unique(mtd$Stage) %in% stagesLev)

whone <- df$TotalAbundance > quantile(df$TotalAbundance, probs = 0.99 )

# hkeepid <- df$mature[whone]
# which(df$mature %in% hkeepid)


# apply(x, 2, function(x) (x / x[43])*100) -> count

# rownames(x) <- df$mature

# membership.cols <- mtd$Group[mtd$ID %in% colnames(x)]

# superheat::superheat(x, pretty.order.rows = T, membership.cols = membership.cols)



# heatmap con patchwork para barplot
rownames(x) <- df$mature

hclust <- hclust(dist(x), "complete")

matureLev <- rownames(x)[hclust$order]

x %>%
  as_tibble(rownames = 'mature') %>%
  pivot_longer(cols = names(x), names_to = 'ID') %>%
  left_join(mtd) %>%
  mutate(tag = ifelse(mature %in% 'cgi-miR-33-5p', '*', '')) %>%
  mutate(value = ifelse(value > 0, value, NA)) %>%
  mutate(Stage = factor(Stage, levels = stagesLev)) %>% 
  arrange(match(Stage, stagesLev)) %>%
  mutate(Description = factor(Description, levels = unique(Description))) %>%
  mutate(mature = factor(mature, levels = matureLev)) %>%
  ggplot(aes(x = mature, y = Description, fill = log10(value+1))) +
  geom_tile(color = 'black', size = 0.25) + 
  geom_text(aes(label = tag),  vjust = 0.75, hjust = 0.5, size= 6) +
  scale_fill_viridis_c(name = "", na.value = "white") +
  labs(x = '', y = '') +
  guides(fill = guide_colorbar(barwidth = unit(3.5, "in"),
                               barheight = unit(0.15, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12)) -> pheat

# pheat + scale_y_discrete(position = 'right') -> pheat
# element_text(angle = 45, hjust = 1, vjust = 1, size = 7)


library(ggsci)

gcolor <- structure(c('#D33F3A', '#45B8DA'), 
                    names = c('Tissues', 'Developmental stages'))

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
  mutate(Stage = factor(Stage, levels = stagesLev)) %>% 
  arrange(match(Stage, stagesLev)) %>%
  mutate(Description = factor(Description, levels = unique(Description))) %>%
  mutate(Group = factor(Group, levels = c('Tissues', 'Developmental stages'))) %>%
  ggplot(aes(y = n, x = Description)) +
  geom_col(aes(fill = Group)) +
  labs(x = '', y = '') + #  Precursor miRNA reads (RPM)
  scale_fill_manual('',values = gcolor) +
  theme_classic(base_family = "GillSans", base_size = 12) +
  coord_flip() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 9),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) -> p

p + ggh4x::facet_nested(Group~., 
                            scales = 'free_y', space = 'free', 
                            switch = 'y') -> p


library(patchwork)

pheat + p + plot_layout(widths = c(2, 0.5))-> psave

ggsave(psave, filename = 'precursor_reads_plot.png', path = path, 
       width = 10, height = 4.5)


# PCA

df <- pcadf[, -c(1,2)]

rownames(df) <- pcadf$mature

# df <- t(df)
# df %>% 
#   as_tibble(rownames = 'ID') %>% 
#   left_join(mtd %>%
#       select(ID, Stage, Group)) %>%
#   select(-ID) %>%
#   as.data.frame() -> pcadf

# library(FactoMineR)

PCA <- prcomp(t(log2(df+1)), scale. = FALSE) # log2(count+1))
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], 
  PC2 = PCA$x[,2])

dtvis %>%
  mutate(ID = rownames(.)) %>%
  left_join(mtd %>% select(ID, Stage, Group)) %>%
  mutate(Stage = factor(Stage, levels = stagesLev)) %>% 
  arrange(match(Stage, stagesLev)) %>%
  mutate(Group = factor(Group, levels = c('Tissues', 'Developmental stages'))) %>%
  ggplot(., aes(PC1, PC2, color = Group)) +
  # geom_point(size = 0.5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggrepel::geom_label_repel(aes(label = Stage), alpha = 0.9) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') -> psave

ggsave(psave, filename = 'miRNAs_RPM_plots_PCA.png', path = path, width = 7, height = 5)
