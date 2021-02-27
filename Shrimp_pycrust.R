rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)
library(tidyverse)
library(ggh4x)
library(pheatmap)
# install_github("raivokolde/pheatmap")

obj <- read_xlsx(files[1])

mtd_file <- list.files(path = dir, pattern = "CIAD_MappingFile.csv$", full.names = TRUE)

mtd <- read.csv(mtd_file, sep = "\t")
  # select(Sample_IDs, Tissue, Time) %>%
  # rename(Index = Sample_IDs)

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
  filter(RA > 0) %>%
  group_by(Index) %>%
  mutate(Z = RA - mean(RA) / sd(RA)) %>%
  ggplot(aes(Z, fill = cut(RA, 100))) +
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
# # estadistico por tejido para tomar las significativos funcionales de cada tejido


obj %>% select_at(vars(!colNames)) %>% distinct(SuperPathway, .keep_all = T) -> left_j


# filtering by: ? ----

obj %>% select_if(is.double) %>% pull() %>% na.omit() %>% data.frame(x = .) %>% 
  ggplot() + geom_histogram(aes(x), bins = 70)

obj %>% 
  pivot_longer(cols = all_of(colNames), 
               values_to = "RA", 
               names_to = "Index") %>%
  left_join(mtd) %>%
  group_by(SuperPathway, Tissue) %>%
  summarise(mean_ab = mean(RA)) %>%
  filter(mean_ab > 1) %>%
  pivot_wider(names_from = Tissue, values_from = mean_ab) %>%
  left_join(left_j) %>%
  ungroup() %>%
  arrange(Category) %>%
  mutate(SuperPathway = forcats::fct_reorder(SuperPathway, Category)) -> L3_agg



# annotation_colors <- unique(L3_agg$Category)
labels <- L3_agg %>% pull(Category) %>% unique() %>% sort()

fillCol <- length(labels)

if(fillCol > 26) {
  annotation_colors <- colorRampPalette(pal_ucscgb()(26))(fillCol)
} else
  annotation_colors <- pal_ucscgb()(fillCol)



names(annotation_colors) <- labels

annotation_colors <- list(Category = annotation_colors)

my_gene_col <- data.frame(row.names = as.character(L3_agg$SuperPathway), 
                          Category = factor(L3_agg$Category, levels = unique(L3_agg$Category)))

my_sample_col <- mtd[-1] %>% as.data.frame()

# row.names(my_sample_col) <- colNames

# 
n <- nrow(L3_agg)

# pal <-  colorRampPalette(pal_material(palette = c("purple"))(10))(n)

pal <- viridis::viridis(n)

L3_agg %>%
  select_if(is.double) %>%
  data.frame(row.names = rownames(my_gene_col), .) %>% 
  pheatmap(., 
           annotation_row = my_gene_col, 
           annotation_colors = annotation_colors,
           cluster_rows = T,
           cluster_cols = T,
           cutree_rows = 3,
           cutree_cols = 3,
           na_col = "white",
           cellwidth = 30,
           cex = 1,
           border_color = T,
           angle_col = c("45"),
           # legend_breaks = -1:4
           # annotation_name_row = T,
           color = pal,
           fontsize = 14, 
           show_rownames = T) -> my_pheatmap

save_pheatmap_png <- function(x, filename, width= 2000, height=1400, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_pheatmap,
                  paste0(dir, "/pheatmap_kegg_cluster_true.png"))

library(grid)

grid.ls(grid.force())
grid.gedit("GRID.rect.17286", gp = gpar(col="white"))
grid.draw(my_pheatmap)
