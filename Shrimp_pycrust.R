rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)

obj <- read_xlsx(files[1])


library(tidyverse)
library(ggh4x)

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
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

obj <- read_xlsx(files[2])

# clean taxonomy
obj %>%
  mutate_at(ranks, funs(str_replace_all(., c("_[1-9]" = "")))) -> obj

obj %>%
  select_if(is.double) %>%
  names() -> colNames

barTax <- function(obj, colNames, agglom_lev = "Phylum", low_ab = 1) {
  obj %>% 
    pivot_longer(cols = colNames, names_to = "Index", values_to = 'ab') %>%
    filter(ab > 0) %>%
    inner_join(mtd) %>%
    rename( "Level" = agglom_lev) %>%
    group_by(Level, Tissue) %>%
    summarise(ab = sum(ab), Freq = length(Level > 0)) %>%
    group_by(Tissue) %>%
    mutate(RA = (ab / sum(ab)) * 100) %>%
    mutate(Level = ifelse(RA < low_ab, "ZLow", Level)) -> dataViz
  
  labels <- dataViz %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  if(colourCount > 7) {
    getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
  } else
    getPalette <- pal_locuszoom()(colourCount)

  
  getPalette[length(getPalette)] <- "Black"
  labels[length(labels)] <- "Low abundance"
  
  dataViz %>%
    ggplot(aes(x = Tissue, y = RA, fill = Level)) +
    geom_col() +
    coord_flip() +
    labs(x = "Tissue", y = "Relative Abundance (%)", 
         caption = paste0("Low Abundance (Relative Ab <", low_ab, " %)")) +
    scale_fill_manual(agglom_lev, labels = labels, values = getPalette) +
    theme_classic(base_size = 18, base_family = "GillSans")
  
}

PPhylum <- barTax(obj, colNames, "Phylum")
PClass <- barTax(obj, colNames, "Class", low_ab = 1)
POrder <- barTax(obj, colNames, "Order", low_ab = 1)

ggsave(PPhylum, filename = "Bar_Phylum.png", path = dir, 
       width = 10, height = 8)
ggsave(PClass, filename = "Bar_Class.png", path = dir, 
       width = 10, height = 8)
ggsave(POrder, filename = "Bar_Order.png", path = dir, 
       width = 10, height = 8)
# by family

agglom_lev <- "Family"

obj %>% 
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  rename("Level" = agglom_lev) %>%
  # mutate(Level = str_replace_all(Level, c("_[1-9]" = ""))) %>%
  group_by(Level, Index) %>%
  summarise(ab = sum(ab), Freq = length(Level > 0)) %>%
  inner_join(mtd) %>%
  filter(grepl('ceae', Level)) %>% 
  group_by(Tissue) %>%
  mutate(RA = (ab / sum(ab)) * 100)  -> Fam_agg

Fam_agg %>% group_by(Tissue) %>% summarise(sum(RA))

# prepare data to select top

obj %>%
  select_at(vars(c(agglom_lev, colNames))) %>%
  rename("Level" = agglom_lev) %>%
  # mutate(Level = str_replace_all(Level, c("_[1-9]" = ""))) %>%
  filter(grepl('ceae', Level)) %>%
  group_by(Level) %>%
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

apply(agg_wide[-1], 2, pick_top, top = 10, 
      y = agg_wide$Level) %>%
  as_tibble() %>%
  pivot_longer(all_of(colNames), names_to = 'Index', 
               values_to = "Level") %>%
  distinct(Level) %>%
  inner_join(agg_wide) -> fam_top

# make tax clustering 

library(superheat)

m <- data.frame(fam_top[-1])
rownames(m) <- fam_top$Level

raf <- function(x) x/sum(x) * 100

superheat(apply(m, 2, raf),
                     scale = F,
                     row.dendrogram = T,
                     col.dendrogram = T,
                     clustering.method = 'hierarchical',
                     dist.method = 'euclidean',
                     print.plot = T) -> sh

tax_hclust <- sh$order.rows
tax_hclust <- rownames(m)[tax_hclust]

obj %>%
  select_at(vars(ranks)) %>%
  distinct_at(agglom_lev, .keep_all = T) %>%
  rename("Level" = agglom_lev) %>%
  inner_join(fam_top, by = "Level") %>%
  mutate_at(colNames, raf) %>%
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ra') %>%
  filter(ra > 0) %>% inner_join(mtd) -> dataHeat

# sanity check

dataHeat %>% group_by(Tissue, Index) %>% summarise(sum(ra))

# set left-panel (Phylum) ordering 
dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ra)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

dataHeat %>%
  mutate(Level = factor(Level, levels = tax_hclust),
         Phylum = factor(Phylum, levels = PhylumLevel),
                         ra = ifelse(ra == 0, NA, ra)) %>%
  ggplot(aes(y = Level, x = Index, fill = ra)) +
  geom_tile() +
  # facet_grid( ~ as.factor(Tissue) , scales = "free", space = "free", switch = "x") +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white',
                             limits = c(0,100),
                             labels = scales::percent_format(scale = 1)) +
  # ggh4x::facet_nested(~ Tissue + Time, scales = "free", space = "free") +
  ggh4x::facet_nested(Phylum + Class ~ Tissue + Time,
                      scales = "free", space = "free" ) +
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
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text.y = element_text(
      angle = 0, 
      size = 14),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.4),
    panel.spacing = unit(0.007, "lines")) -> heatPlot

ggsave(heatPlot, filename = "heatmap.png", path = dir, 
       width = 22, height = 16)

# Check redundant abreviation nomenclature ----
# raw datavis

obj %>% filter(!grepl('ceae', Family)) %>% arrange(desc(ab_agg))-> rawDv

library(stringr)

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

obj_longer %>% mutate_at(ranks, funs(str_replace_all(., c("^NA_[1-1000]" = NA)))) 

mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__" = "")))) -> obj_longer

obj_longer %>%
  separate(col = Replicate, into = c("Sample", "Rep"),
           sep = "_") -> obj_longer
  
# test features diversity ----
source("~/Documents/GitHub/metagenomics/estimate_richness.R")

obj %>%
  select_at(vars(c(colNames))) %>%
  estimate_richness(., measures = c("Observed", "Shannon", "Chao1", "Fisher")) %>%
  as_tibble(rownames = "Index") %>%
  inner_join(mtd) -> diversityDat

diversityDat %>%
  ggplot(aes(x = Chao1, y = Observed)) +
  geom_point(data = diversityDat %>% select(-Time), 
             colour = "grey70", alpha = 0.3, size = 2.5) +
  geom_point(aes(colour = Tissue), size = 2.5) +
  xlab("Predicted (Chao1 Estimator)") +
  ylab("Observed") +
  facet_wrap(~ Time ) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  ggsci::scale_color_rickandmorty() -> div1 

ggsave(div1, filename = "obs_vs_chao_tissue_time.png", path = dir, 
       width = 12, height = 6)

diversityDat %>%
  ggplot(aes(x = Chao1, y = Observed)) +
  geom_point(data = diversityDat %>% select(-Tissue), colour = "grey70",
             alpha = 0.3) +
  geom_point(aes(colour = Shannon)) +
  facet_wrap(~ Tissue ) +
  ggrepel::geom_text_repel(max.overlaps = 30, size = 2.5,
                            aes(Chao1, Observed,
                                label = Index, color = Shannon)) +
  labs(xlab = "Predicted (Chao1 Estimator)",
         ylab = "Observed") +
  scale_colour_viridis_c(option = "magma", begin = 0.2, end = 0.8) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top") +
  guides(color = guide_colorbar(barheight = unit(0.2, "in"),
                                barwidth = unit(8, "in"),ticks.colour = "black", frame.colour = "black",label.theme = element_text(size = 10))) -> div2

ggsave(div2, filename = "obs_vs_chao_tissue_shannon.png", path = dir, 
       width = 12, height = 6)

diversityDat %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  # mutate(Time = ifelse(Time %in% "Farm", 0, Time)) %>%
  # mutate(Time = as.double(Time)) %>%
  # mutate(Time = factor(Time, levels = c(0,20,40,60,80))) %>%
  ggplot(aes(x = Time, y = Shannon)) +
  geom_boxplot(lwd = 0.5) +
  # stat_boxplot(geom ='errorbar', linetype = "dotted" ) +
  geom_point(aes(x = Time, y = Shannon, color = Tissue)) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  ggsci::scale_color_rickandmorty() -> div3
div3

ggsave(div3, filename = "Time_shannon_tissue.png", path = dir, 
       width = 10, height = 4)

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

# pos-clustering w/ Lulu ----

# devtools::install_github("tobiasgf/lulu") 
library(lulu)

min_r <- 1
min_match <- 98
min_cooccurence <- 0.95 

# Analysis composition of microbiomes w/ ANCOM-BC ----

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")

# http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#
library(ANCOMBC)
library(phyloseq)

sam <- data.frame(mtd, row.names = mtd$Index) %>% 
  arrange(match(Index, colNames)) %>% mutate_if(is.character, as.factor)

dat <- obj %>% select_at(colNames) %>% data.frame(row.names = obj$`Feature ID`)
tax <- obj %>% select_at(ranks)  %>% data.frame(row.names = obj$`Feature ID`)
identical(names(dat), rownames(sam))
identical(rownames(dat), rownames(tax))

# and parse
phyloseq = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(sam)) %>%
  subset_taxa(Family %in% tax_hclust) %>%
  prune_taxa(taxa_sums(.) > 0, .)

ancombc_data <- aggregate_taxa(phyloseq, "Family")

# tax_mat = as(tax_table(phyFam), "matrix")

# Run ancombc function

out <- ANCOMBC::ancombc(phyloseq = ancombc_data, formula = "Tissue",
              p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
              group = "Tissue", struc_zero = TRUE, neg_lb = TRUE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = TRUE)

res = out$res

# res_global = out$res_global
# res = cbind(taxon = rownames(out$feature.table), out$res)

df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id") %>%
  pivot_longer(-taxon_id, names_to = "group", values_to = "logFC")
df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id") %>% 
  pivot_longer(-taxon_id, names_to = "group", values_to = "SD")

# set the contrast (look at the contrasts)
mtd %>%
  with(., table(Tissue, Time))
# colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")

df_fig1 %>% left_join(df_fig2) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>%
  mutate(group = str_replace_all(group, c("Tissue" = ""))) %>%
  mutate(wrap = group) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig

df_fig$taxon_id = factor(df_fig$taxon_id, levels = unique(df_fig$taxon_id))

p = ggplot(data = df_fig, 
           aes(x = taxon_id, y = logFC, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logFC - SD, 
                    ymax = logFC + SD), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  scale_fill_discrete("Intercept") +
  guides(color = "none") +
  labs(x = NULL, y = "Log fold change", 
       title = "Analysis composition of microbiomes w/ bias correction (ANCOM-BC)\n",
       caption = "g1 (ie. Tissue in the facet) and\ng2 (Hepatopancreas as intercept)") + 
  theme_classic(base_size = 16, base_family = "GillSans") + 
  facet_wrap(~wrap) +
  coord_flip() +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(p, filename = "ANCOMBC.png", path = dir, 
       width = 10, height = 8)


