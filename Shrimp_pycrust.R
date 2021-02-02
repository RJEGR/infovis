rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)

obj <- read_xlsx(files[1])


library(tidyverse)

mtd <- read.csv(list.files(path = dir, pattern = "txt$", full.names = TRUE), sep = "\t") %>% 
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

dataViz %>%
  mutate(Category = factor(Category, 
                               levels = LCategory)) %>%
  mutate(SuperPathway = forcats::fct_reorder(SuperPathway, Rank)) %>%
  ggplot() +
  geom_col(aes(x = SuperPathway ,
               y = RA, fill = Category)) +
  labs(y = "Relative Abundance of predicted gene (%)") +
  ggsci::scale_fill_ucscgb() +
  facet_grid(Category ~ Tissue , scales = "free", space = "free_y") +
  coord_flip() +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank()) -> barPlot

# barPlot

ggsave(barPlot, filename = "barKegg.png", path = dir, 
       width = 18, height = 10)

# TaxAb ----
  
obj <- read_xlsx(files[2])

obj %>%
  select_if(is.double) %>%
  names() -> colNames

# obj_longer %>% pull(Family) %>% unique() 

obj %>% 
  pivot_longer(cols = colNames, names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  inner_join(mtd) %>%
  group_by(Phylum, Tissue) %>%
  summarise(ab = sum(ab)) %>%
  group_by(Tissue) %>%
  mutate(RA = (ab / sum(ab)) * 100) %>%
  mutate(Phylum = ifelse(RA < 1, "ZLow", Phylum)) -> dataViz

labels <- sort(unique(dataViz$Phylum))
colourCount = length(labels)
getPalette <- pal_locuszoom()(colourCount)

getPalette[length(getPalette)] <- "Black"
labels[length(labels)] <- "Low abundance"

dataViz %>%
  ggplot(aes(x = Tissue, y = RA, fill = Phylum)) +
  geom_col() +
  coord_flip() +
  labs(x = "Tissue", y = "Relative Abundance (%)", 
       caption = "Low Abundance (Relative Ab < 1 %)") +
  scale_fill_manual(labels = labels, values = getPalette) +
  theme_classic(base_size = 16, base_family = "GillSans") -> splot

ggsave(splot, filename = "Bar_Phylum.png", path = dir, 
       width = 7, height = 5)

# by family

obj %>% 
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  group_by(Family, Index) %>%
  summarise(ab = sum(ab)) %>%
  inner_join(mtd) %>%
  filter(grepl('ceae', Family)) %>%
  group_by(Tissue) %>%
  mutate(RA = (ab / sum(ab)) * 100)  -> Fam_agg

Fam_agg %>% group_by(Tissue) %>% summarise(sum(RA))

Level <- "Family"

obj %>%
  select_at(vars(c(Level, colNames))) %>%
  filter(grepl('ceae', Family)) %>%
  group_by_at(vars(one_of(Level))) %>%
  summarise_at(vars(colNames), sum) %>%
  ungroup() -> agg_wide


pick_top <- function(x, y, top = 10) {
  
  # x <- vector of abundance
  # y <- vector of taxonomic name
  
  ordered <- order(x, decreasing = T)
  topPos <- head(ordered, top)
  taxPos <- y[topPos]
  
  return(taxPos)
}

apply(agg_wide[-1], 2, pick_top, top = 10, y = agg_wide$Family) %>%
  as_tibble() %>%
  pivot_longer(all_of(colNames), names_to = 'Index', 
               values_to = Level) %>%
  distinct_at(Level) %>%
  inner_join(agg_wide) -> fam_top

# make tax clustering 

library(superheat)
m <- data.frame(fam_top[-1])
rownames(m) <- fam_top$Family

raf <- function(x) x/sum(x) * 100

superheat(apply(m, 2, raf),
                     scale = F,
                     row.dendrogram = T,
                     col.dendrogram = T,
                     clustering.method = 'hierarchical',
                     dist.method = 'euclidean',
                     print.plot = F) -> sh

tax_hclust <- sh$order.rows
tax_hclust <- rownames(m)[tax_hclust]

obj %>%
  select_at(vars(ranks)) %>%
  distinct(Family, .keep_all = T) %>%
  inner_join(fam_top, by = Level) %>%
  mutate_at(colNames, raf) %>%
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ra') %>%
  filter(ra > 0) %>% inner_join(mtd) -> dataHeat

dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ra)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

dataHeat %>%
  mutate(Family = factor(Family, levels = tax_hclust),
         Phylum = factor(Phylum, levels = PhylumLevel),
                         ra = ifelse(ra == 0, NA, ra)) %>%
  ggplot(aes(y = Family, x = Index, fill = ra)) +
  geom_tile() +
  facet_grid( ~ Tissue, scales = "free", space = "free") +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white',
                             limits = c(0,100),
                             labels = scales::percent_format(scale = 1)) +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_family = "GillSans", base_size = 14) +
  guides(fill = guide_colorbar(barheight = unit(9, "in"), 
                               barwidth = unit(0.5, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 14))
  ) +
  
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1),
    strip.text.y = element_text(
      angle = 0, 
      size = 14),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.4),
    panel.spacing = unit(0.007, "lines")) -> heatPlot

heatPlot + facet_grid(Phylum + Class ~ Tissue, 
                      scales = "free", space = "free" ) -> heatPlot

ggsave(heatPlot, filename = "heatmap.png", path = dir, 
       width = 16, height = 14)

# Check redundant abreviation nomenclature 
# raw datavis

obj_longer %>% filter(!grepl('ceae', Family)) %>% arrange(desc(ab_agg))-> rawDv

library(stringr)

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

obj_longer %>% mutate_at(ranks, funs(str_replace_all(., c("^NA_[1-1000]" = NA)))) 

mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__" = "")))) -> obj_longer

obj_longer %>%
  separate(col = Replicate, into = c("Sample", "Rep"),
           sep = "_") -> obj_longer
  



# metacoder

# install.packages("devtools")
# devtools::install_github("grunwaldlab/metacoder")

library(metacoder)

obj %>%
  select_at(vars(ranks)) %>%
  distinct(Family, .keep_all = T) %>%
  inner_join(fam_top, by = Level) %>%
  parse_tax_data(class_cols = ranks, named_by_rank = TRUE) -> x

# getting per-taxon information
x$data$tax_abund <- calc_taxon_abund(x, "tax_data",
                                       cols = mtd$Index)

x$data$diff_table <- compare_groups(x, 
                                      dataset = "tax_abund",
   # What columns of sample data to use
                                      cols = mtd$Index,
   # What category each sample is assigned to
                                      groups = mtd$Tissue) 

set.seed(1)

heat_tree_matrix(x,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree_Shrimp.pdf")

# or


x$data$tax_occ <- calc_n_samples(x, "tax_abund", groups = mtd$Tissue, cols = mtd$Index)

set.seed(1) 
heat_tree(x, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Intestine, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          output_file = "Intestine_heat_tree_Shrimp.pdf") 

heat_tree(x, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Hepatopancreas, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          output_file = "Hepatopancreas_heat_tree_Shrimp.pdf") 

heat_tree(x, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Stomach, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          output_file = "Stomach_heat_tree_Shrimp.pdf") 
