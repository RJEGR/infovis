rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)
library(tidyverse)
library(ggh4x)
library(pheatmap)
# install_github("raivokolde/pheatmap")

obj <- read_xlsx(files[1])
mtd_file <- list.files(path = dir, pattern = "CIAD_MappingFile.txt$", full.names = TRUE)

mtd <- read.csv(mtd_file, sep = "\t") %>% 
  select(Sample_IDs, Tissue, Time) %>%
  rename(Index = Sample_IDs)


obj %>%
  select_if(is.double) %>%
  names() -> colNames

obj %>%
  pivot_longer(cols = all_of(colNames), 
               values_to = "RA", 
               names_to = "Index") %>% 
  inner_join(mtd) -> obj_longer

# arrange by factor


obj_longer %>%
  arrange(desc(Category)) %>%
  mutate(SuperPathway = factor(SuperPathway, 
                               levels = unique(SuperPathway))) -> dataViz

# view(dataViz)

# Which noised data we've? (look at x-axis)

dataViz %>%
  group_by(Index) %>%
  mutate(Z = RA - mean(RA) / sd(RA)) %>%
  ggplot(aes(RA, fill = cut(RA, 100))) +
  geom_histogram(show.legend = F, bins = 100) +
  facet_wrap(~ Tissue)


obj_longer %>%
  # filter(RA > 12000) %>%
  group_by(SuperPathway, Tissue, Category) %>%
  summarise(RA = sum(RA)) %>%
  group_by(Tissue) %>%
  mutate(Rank = rank(RA)) -> dataViz

LCategory <- dataViz %>% group_by(Category) %>% summarise(t = sum(RA)) %>% arrange(desc(t)) %>% pull(Category)

library(ggsci)

labels <- dataViz %>% pull(Category) %>% unique() %>% sort()
fillCol <- length(labels)

if(fillCol > 26) {
  getPalette <- colorRampPalette(pal_ucscgb()(26))(fillCol)
} else
  getPalette <- pal_ucscgb()(fillCol)

dataViz %>%
  mutate(Category = factor(Category, 
                               levels = LCategory)) %>%
  mutate(SuperPathway = forcats::fct_reorder(SuperPathway, Rank)) %>%
  ggplot() +
  geom_col(aes(x = SuperPathway ,
               y = RA, fill = Category)) +
  labs(y = "Relative Abundance of predicted gene (%)") +
  scale_fill_manual(values = getPalette) +
  facet_grid(Category ~ Tissue , scales = "free", space = "free_y") +
  coord_flip() +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank()) -> barPlot

# barPlot

ggsave(barPlot, filename = "barKegg.png", path = dir, 
       width = 18, height = 10)

# pheatmap ----
# en base a las 33 categorias de L2 ie. obj %>% group_by(SuperPathway)

obj %>% select_at(vars(!colNames)) %>% distinct(SuperPathway, .keep_all = T) -> left_j

obj %>% 
  group_by(SuperPathway) %>%
  summarise_at(vars(colNames), sum) %>%
  left_join(left_j) -> L3_agg

annotation_colors <- unique(L3_agg$Category)
names(annotation_colors) <- unique(L3_agg$Category)
annotation_colors <- list(Category = annotation_colors)

my_gene_col <- data.frame(row.names = L3_agg$SuperPathway, 
                          Category = factor(L3_agg$Category, levels = unique(L3_agg$Category)))

my_sample_col <- mtd[-1] %>% as.data.frame()
row.names(my_sample_col) <- colNames


# 
# n <- nrow(obj)
# pal <- colorRampPalette(ggsci::pal_material(palette = c("purple"),n = n, alpha = 1, reverse = F)(n))

L3_agg %>%
  select_at(vars(colNames)) %>%
  data.frame(row.names = rownames(my_gene_col), .) %>% 
  pheatmap(., 
           annotation_row = my_gene_col, 
           annotation_col = my_sample_col,
           cluster_rows = T,
           # color = pal,
           cluster_cols = T,
           fontsize_row = 7,
           show_rownames = T) -> my_pheatmap
           # annotation_colors = annotation_colors) 

save_pheatmap_png <- function(x, filename, width= 1500, height=1400, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_pheatmap,
                  paste0(dir, "/pheatmap_kegg.png"))
