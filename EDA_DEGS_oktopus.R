
rm(list = ls())

# ==============
## Checking and Load packages ----
# ==============

.cran_packages <- c("wTO", "CoDiNA") # "tidyverse"
# .bioc_packages <- c("edgeR","DESeq2", "biomaRt", "topGO", "Rgraphviz", "wTO")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
# .inst <- .bioc_packages %in% installed.packages()
# if(any(!.inst)) {
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install(.bioc_packages[!.inst], ask = F)
# }

source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

# load count data

path <- "~/transcriptomics/oktopus_full_assembly/"
pattern <- "counts_table_length_ajus_gen_level-aproach2-filtered.txt"
countf <- list.files(path = path, pattern = pattern, full.names = T)

dim(count0 <- read.delim(countf, sep = "\t"))

mtd <- read.delim(paste0(path, "metadata.tsv"), sep = "\t") 


# load up/down regulated overlaps  (lists)

path <- "~/transcriptomics/oktopus_full_assembly/DGE_Pavel_2/"
pattern <- "just_upgenes.txt_add_annotation.txt.sort"
overlapsf <- list.files(path = path, pattern = pattern, full.names = T)

readGenes <- function(file) {
  
  require(tidyverse)
  
  df <- read.delim(file, sep = "\t") 
  group <- basename(file)
  group <- sapply(strsplit(group, "[.]"), `[`, 1)
  
  df %>% mutate(group = group) %>% separate(col = group, into = c("Tissue", "Intercept"), sep = "vs") %>%
    mutate(Develope = substr(Tissue, nchar(Tissue)-2, nchar(Tissue))) %>% 
    mutate(Tissue = substr(Tissue, 1, nchar(Tissue)-3)) %>%
    mutate(Tissue = recode(Tissue, Glandula_optica = "Optic\nGland",
                           Glandula_oviducal = "Oviducal\nGland", 
                           Lobulo_Optico = "Optic\nLobe"))

}

# readGenes(overlapsf[1])
# EDA 
# 1)
lapply(overlapsf, readGenes) %>% do.call(rbind, .) %>% as_tibble() -> dff

# saveRDS(dff, file = paste0(path, 'upgenes.rds'))

dff %>% group_by(Tissue, Intercept, Develope) %>% summarise(n = length(ID))

library(UpSetR)

colNames <- names(count)
  
# nota: la lista de genes up Glandula_oviducalDESvsLO_GO, corresponde a la comparacion de la glandula oviductal(en desove) vs el promedio de la expresión de Glandula optica (GO) con lobulo optico (LO).
# count %>%
#   filter(rownames(.) %in% unique(dff$ID)) %>%
#   # mutate_at(vars(colNames), function(x) {ifelse(x > 0, 1, 0) }) %>%
#   as_tibble(rownames = "ID") %>%
#   pivot_longer(cols = colNames, names_to = "id") %>%
#   left_join(mtd) %>%
#   group_by(ID, Tissue) %>%
#   summarise(value = sum(value)) %>%
#   mutate(value = ifelse(value > 0, 1, 0)) %>%
#   pivot_wider(names_from = Tissue, values_from = value) %>%
#   select(-ID) %>% as.data.frame() %>%
#   upset(number.angles = 0, point.size = 3.5, 
#         line.size = 2,  #nsets = 20,
#         nintersects = 40, 
#         mainbar.y.label = "Gene Intersections", sets.x.label = "Genes per Intercept", 
#         text.scale = c(1.3, 2, 1, 2, 2, 2),
#         order.by = "freq")


length(unique(dff$ID)) # 3744 DEGs 
length(dff$ID) # separated by interects of 6896

dff %>% with(., table(Tissue, Intercept, Develope))

dff %>% 
  mutate(Develope = factor(Develope, levels = c("PRE", "DES", "POS"))) %>%
  ggplot() +
  geom_histogram(aes(FDR), bins = 70) +
  geom_vline(color = "black", xintercept = 0.01, linetype="dashed", alpha=0.5) +
  labs(x = expression(~P[adjustable]), y = "Number of genes") +
  ggh4x::facet_nested(Tissue+Intercept ~., scales = "free") +
  # ggh4x::facet_nested(Tissue+Intercept ~ Develope, scales = "free") +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) -> p1

# ggsave(p1, filename = "significance_hist_tissueDev_vs_intercept.png", path = path, width = 8, height = 7)

ggsave(p1, filename = "significance_hist_tissue_vs_intercept.png", path = path, width = 7, height = 7)
  
# test outliers
colNames <- names(count0) 

count0 %>% 
  filter(rownames(.) %in% unique(dff$ID)) %>%
  pivot_longer(cols = colNames, values_to = "x", names_to = "id") %>% 
  group_by(id) %>% mutate(ecdf = ecdf(x)(x)) %>% arrange(desc(x)) %>%
  left_join(mtd) %>%
  filter(x > 0 ) %>%
  group_by(Tissue) %>%
  mutate(x = log2(x+1), cor = cor(x, ecdf, method = "pearson")) %>% 
  # ggplot() + geom_point(aes(x, ecdf , color = cor), alpha = 0.8) +facet_grid(~Tissue)
  ggplot() + stat_ecdf(aes(x , color = Tissue), size = 1) + # group = group
  # geom_smooth(aes(x, ecdf, color = Tissue), se = F, method = "lm")  +
  scale_y_reverse() +
  scale_color_brewer(palette = "Set1") +
  labs(y = expression(~f*italic("(x)")),  x= expression(~Log[10]~("x"~+1)), 
       caption = "Empirical Cumulative Distribution of degs") +
  theme_classic(base_family = "GillSans", base_size = 16) -> p2

ggsave(p2, filename = "ecdf_degs.png", path = path, width = 7, height = 7)

#

# heatmap of degs

degs <- unique(dff$ID)
# 
sum(rownames(count0) %in% degs) / length(degs)
# 
colNames <- names(count0)

count0 %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(all_of(colNames), names_to = "group", values_to = "value") %>%
  left_join(mtd, by = "group") -> count_longer

count_longer %>%
  filter(gene %in% degs) %>%
  select(gene, id, value) %>%
  mutate(value =  ifelse(value < 1, 0, value)) %>%
  pivot_wider(names_from = id, values_from = value, values_fill = NA) %>%
  mutate_if(is.double, function(x) log2(x+1)) %>%
  data.frame(row.names = 'gene') -> df

hclust <- hclust(dist(df), "complete")


count_longer %>%
  filter(gene %in% degs) %>%
  mutate(value =  ifelse(value < 1, NA, value)) %>%
  mutate(stage = factor(stage, levels = c("PRE", "Spawing", "POST"))) %>%
  mutate(group = forcats::fct_reorder(group, as.numeric(stage))) %>%
  ggplot(aes(y = gene, x = group, fill = log2(value+1))) +
  geom_tile(width = 1.5) +
  ggh4x::facet_nested(~ Tissue+stage, scales = "free", space = "free", nest_line = T) +
  scale_fill_viridis_c(name = expression(~Log[2]~("x"~+1)), na.value = "white") +
  labs(x = "", y = "Neuropeptide (blastp)") +
  theme_classic(base_family = "GillSans", base_size = 16) +
  theme(panel.spacing = unit(0, "lines"),
        axis.title.y = element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        strip.background=element_rect(fill=NA, color=NA)) +
  guides(fill = guide_colorbar(barheight = unit(4, "in"), 
                               barwidth = unit(0.3, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 12))) -> pheat

pheat + ggh4x::scale_y_dendrogram(hclust = hclust) -> pheat

ggsave(pheat, filename = paste0(path, "heatmap_degs.png"), width = 10, height = 10)

# 2) run net by file
# 1) choosing neuropeptide as overlap list of genes (peptideGenes)
# 2) subseting data count by degs
# 3) log2(x+1) transformation of count data

wd <-'~/transcriptomics/Diana_Lara/neuropeptides'
neuropep <- dir(path = wd, pattern = 'neuropeptide.xls', full.names = T)

require('trinotateR')

y <- read_trinotate(neuropep)

length(peptideGenes <- as.character(unique(y$gene_id)))

write(peptideGenes, file = paste0(path, "peptideGenes.list"))

# dim(df <- readGenes(overlapsf[1]))
# dim(df <- df %>%  filter(FDR < 0.005))

dim(df <- dff %>%  filter(FDR < 0.01 & Tissue %in% "Glandula_optica")) # table(dff$Tissue)

length(degs <- unique(df$ID))

# writeLines(degs, paste0(path, "/degs.lists"))


sum(degs %in% peptideGenes)
Overlap <- degs[degs %in% peptideGenes]

# length(Overlap <- unique(dff$ID))
# subset count-matrix based on significance and prevalence

sum(rownames(count) %in% degs)
dim(Data <- count[rownames(count) %in% degs,])
Data <- Data %>% mutate_at(vars(colNames), function(x) {log2(x + 1)})

# 
# geneSum <- 10
# prevalence <- 1
# dim(count <- count[rowSums(cpm(count) > geneSum) >= prevalence,])

# Note that, the Pearson’s correlation coefficient is sensitive to extreme values, and therefore it can exaggerate or under-report the strength of a relationship. The Spearman Rank Correlation is recommended when data is monotonically correlated, skewed or ordinal, and it is less sensitive to extreme outliers than the Pearson coefficient (Daysi et al 2018). I found extreme values skewed at the right side as well as an non-normal distribution, therefore we perform the net using spearman

Data %>% 
  pivot_longer(cols = colNames, values_to = "x", names_to = "id") %>% 
  left_join(mtd) %>%
  rename("g" = "Tissue") %>%
  is_parametric()

network = wTO.Complete(n = 100, k = 5,  Data = Data, 
                          method_resampling = 'Bootstrap', 
                          Overlap = Overlap, method = 's', plot = F, savecor = T) 

wTO = network$wTO
summary(wTO)

wTO$wTO = ifelse(wTO$Padj_sig<0.05, wTO$wTO_sign, 0 )

# after run wTO_complete.R in all the degs, lets merge it by:

require(CoDiNA)

