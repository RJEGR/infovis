
# red de interacion microbioma
# test https://www.nature.com/articles/s41522-018-0077-y
# (file:///Users/cigom/Downloads/v48i04.pdf qgraph method used in previous ref)
# network microbiome review at https://www.sciencedirect.com/science/article/abs/pii/S0966842X16301858
# The complexity of microbiomes motivates a movement from reductionist approaches that focus on individual pathogens in isolation to more holistic approaches that focus on interactions among members of the community and their hosts. 

# continue w/ Correlation-Based Methods inthe review

# https://academic.oup.com/bioinformatics/article/32/17/2611/2450750?login=true

rm(list = ls())

# library(wTO)
# require(CoDiNA)
library(microbiome)
library(phyloseq)
# library(Biostrings)

library(tidyverse)

library(ggraph)
library(igraph)
library(ggforce)
library(scatterpie)
library(RColorBrewer)

# fasta_file <- paste0(dir, "fastq/dada2_asv/", "features.fasta")
# dna <- Biostrings::readDNAStringSet(fasta_file)
# dna[names(dna) %in% df$Feature_IDs]

# Comensales/basales (C), asignacion no consistente (NC) y Patogena (P)

ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species', 'pplacer_sp')
sampleName <- c('Cantiles', 'Coloradito', 'Granito', 'Machos', 'Partido', 'Rasito')
dir <- '~/metagenomics/Loberas_MG/'

fileNames <- list.files(dir, pattern = "xls", full.names = T)

readXLS <- function(x) {
  group <- paste(sapply(strsplit(basename(x), "_"), `[`, 1),collapse = '_')
  readxl::read_xlsx(x, na = 'NA') %>% mutate(group = group)
}

featureDF <- lapply(fileNames, function(x) {
  y <- readXLS(x)
  # y %>% select(ranks, sampleName)
  return(y)})

featureDF <- do.call(rbind, featureDF)

names(featureDF)[9] <- 'pplacer_sp'

featureDF %>% 
  mutate(Relationship = recode_factor(Relationship, P = 'Pathogenic', 
                                      NC = 'Inconsistent', C = 'non-pathogenic')) -> featureDF


featureDF %>%  mutate_all(., funs(str_replace_all(., c("Bacteroides tectu$" = "Bacteroides tectus")))) -> featureDF
featureDF %>% mutate_at(sampleName, as.double) -> featureDF
# grep tectu = tectus

# clean some non pathogenic taxa
cleanTax <- c('Campylobacter insulaenigrae', 'Helicobacter muridarum',
              'Campylobacter coli', 'Sutterella stercoricanis', 'Campylobacter concisus')

featureDF %>% filter(!pplacer_sp %in% cleanTax) -> featureDF

# featureDF %>% filter(pplacer_sp %in% cleanTax) %>%
# select(Feature_IDs, group, pplacer_sp, sampleName)


saveRDS(featureDF, file = paste0(dir, '/featureDF.rds'))




# seleccionar asvs abundantes ----
# 

featureDF %>% select_at(ranks) -> tax

tax %>% select(-pplacer_sp) -> tax

max.rank <-  ncol(tax) -1 # Kingdom is not taking in account

tax$Resolution <- max.rank

for(t in max.rank:1){
  rr <- t-1  #rank real 
  rc <- t+1  #column-index (rank) para corregir
  if(t == max.rank){
    tax$Resolution[is.na(tax[,rc])] <- rr
  } else {
    tax$Resolution[is.na(tax[,rc]) & tax$Resolution <= t ] <- rr
  }
}

# table(tax$Resolution) - sum( table(tax$Resolution))

featureDF %>% select_at(sampleName) %>% mutate_at(sampleName, as.double)-> ab

prevelancedf = apply(X = ab,
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

# identical(rownames(tax), rownames(ab)) # sanity check

r <- c(0:6)
names(r) <- ranks[-8]

data.frame(Resolution = r) %>% as_tibble(rownames = "rank") -> r


df = data.frame(Prevalence = prevelancedf, 
                TotalAbundance = rowSums(ab),
                tax, group = featureDF$group, Relationship = featureDF$Relationship) %>% 
  as_tibble(rownames = 'id')

df %>% left_join(r) %>% select(-Resolution) %>% mutate(rank = factor(rank, levels = ranks)) -> df

# df %>% filter(group %in% 'V3') %>% drop_na(Order)

df %>% 
  group_by(rank, group) %>%
  tally() %>% 
  group_by(group) %>%
  mutate(nn = cumsum(n) - n, pct = 1 - nn/sum(n), cs = sum(n)-nn) %>% # filter(rank %in% 'Genus')
  ggplot(aes(x = rank, y = pct, color = cs, group = group)) +
  geom_path(size = 2, alpha=0.6) +
  geom_point(size=2, alpha=0.6) +
  geom_text(aes(label = cs), size=4, vjust = -1, color = 'black') +
  scale_color_viridis_c(option = "A", direction = -1) +
  labs(y = "Changes in the assignation (%)", x = 'rdp classification', color = "Number of\nFeatures") +
  facet_grid(group ~.) +
  theme_bw(base_size = 17) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave


df %>%
  mutate(Prevalence = Prevalence/ncol(ab)) %>%
  mutate(Relationship = ifelse(Relationship %in% 'Pathogenic', 'Pathogenic', "np")) %>%
  ggplot(aes(TotalAbundance, Prevalence)) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "R", 
                   p.accuracy = 0.001, color = 'black') +
  geom_smooth(color = 'blue',
              method = "lm",
              linetype="dashed", size = 0.5, alpha=0.5, 
              se = T, na.rm = TRUE) +
  # geom_point(aes(color = Relationship), size = 2, alpha = 0.7) +
  scale_x_log10() +
  labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Sp. resolution") +
  facet_wrap(~group) +
  theme_bw(base_size = 17) +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) 
  
#

df %>%
  # mutate(Relationship = ifelse(Relationship %in% 'Pathogenic', 'Pathogenic', "np")) %>%
  mutate(group= ifelse(Phylum %in% keepPhyla, group, "Low Taxa")) %>%
  # filter(!Phylum %in% keepPhyla)
  ggplot(aes(group, TotalAbundance, fill = Relationship)) +
  geom_col(position = position_dodge(0.6))
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), outlier.shape=NA) +
  stat_summary(fun = mean, geom="point", shape=20, 
               size = 3, color="red", fill="red", aes(group = Relationship),
               position = position_dodge(0.6)) +
  scale_y_log10() +
  labs(x = "", y = "Total Abundance (log10)")

df %>% 
  select(Prevalence, TotalAbundance, group, Relationship) %>%
  cbind(ab, .) %>% as_tibble() %>% 
  mutate(st)

# test uniformidad

kruskal.test(ab, subset = featureDF$group)
# expected
m <- rowSums(ab) / ncol(ab)
m <- matrix(m, nrow = nrow(ab), ncol = ncol(ab))
colnames(m) <- names(ab)
cor <- cor(m, ab, method = 'pearson')
corrplot::corrplot(cor, type = 'upper')

library(factoextra)
res.pca <- prcomp(ab,  scale = TRUE)
factoextra::fviz_pca(res.pca, label="var", col.ind = "contrib")
res.pca <- prcomp(m,  scale = TRUE)
factoextra::fviz_pca(res.pca, label="var", col.ind = "contrib")


# species detected by rdp -----
featureDF %>% 
  # filter(Relationship %in% 'Pathogenic') %>%
  group_by(Species, group) %>% 
  summarise_at(vars(sampleName), sum) %>% 
  pivot_longer(cols = sampleName, names_to = 'Sample', values_to = 'ab') %>%
  group_by(Species, group) %>% summarise(ab = sum(ab)) %>%
  pivot_wider(names_from = group, values_from = ab, values_fill = 0) %>%
  arrange(Species)

# species detected by pplacer ----

featureDF %>% 
  filter(Relationship %in% 'Pathogenic') %>%
  group_by(pplacer_sp, group) %>% 
  summarise_at(vars(sampleName), sum) %>% 
  pivot_longer(cols = sampleName, names_to = 'Sample', values_to = 'ab') %>%
  group_by(pplacer_sp, group) %>% summarise(ab = sum(ab)) %>%
  pivot_wider(names_from = group, values_from = ab, values_fill = 0) %>%
  arrange(pplacer_sp) # %>% view()


# Gracias Papa
# Gracias, de nuevo 03/06/21
ggsave(psave, filename = "resolution.png", path = dir, width = 16, height = 7)

# Prepare data count, solving otus sporious ----
# https://www.nature.com/articles/ismej2015235
# remove vary rare otus-asvs (sparse less 50%)


source("~/Documents/GitHub/metagenomics/estimate_richness.R")


# specnumber(ab, MARGIN = 2)
# vegan::diversity(ab, index = "invsimpson")


alphadf <- function(featureDF, g) {
  
  measures <- c("Observed", "Shannon","InvSimpson")
  
  featureDF %>% 
    filter(group == g) %>%
    select_at(sampleName) %>% 
    mutate_at(sampleName, as.double) -> ab
  
  Tab <- colSums(ab)
  
  ab %>%
    estimate_richness(., measures = measures) %>%
    as_tibble(rownames = "Index") %>%
    mutate(ab = Tab, group = g)
}

alphaDF <- lapply(unique(featureDF$group), function(x) {
  y <- alphadf(featureDF, g = x)
  return(y)})

alphaDF <- do.call(rbind, alphaDF)

alphaDF %>%
  group_by(group) %>%
  summarise(cor = cor(Observed, InvSimpson))

# multiple linear regression in r ggplot

alphaDF %>%
  ggplot(aes(Observed, InvSimpson, color = group)) +
  geom_point() +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", 
                   p.accuracy = 0.001) +
  geom_smooth(aes(),
              method = "lm", 
              linetype="dashed", size = 0.5, alpha=0.5, 
              se = F, na.rm = TRUE) +
  theme_bw(base_size = 16, base_family = "GillSans")

  # geom_mark_hull(aes(fill = group, label = group))

alphaDF %>%
  ggplot(aes(x = group, y = Observed)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), outlier.shape=NA) +
  stat_summary(fun = mean, geom="point", shape=20, 
               size = 3, color="red", fill="red") +
  coord_cartesian(ylim=c(0,NA)) +
  theme_bw(base_family = "GillSans", base_size = 14)

#
  
featureDF %>% 
  select_at(c(sampleName, 'group')) %>%
  pivot_longer(cols = sampleName) %>%
  group_by(group) %>%
  # filter(value == 0) %>%
  # count() %>% mutate(n = n/6)
  ggplot(aes(group, value)) +
  # geom_col()
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), outlier.shape=NA) +
  stat_summary(fun = mean, geom="point", shape=20, 
               size = 3, color="red", fill="red") +
  coord_cartesian(ylim=c(0,NA))
  

# hist(diversity(ab, index = "invsimpson"))

# -> diversityDat

nfeat <- function(x) {length(unique(x))}

df %>%
  group_by(Phylum, group) %>%
  summarise(n = nfeat(id), nR = nfeat(Family), divergence = 1-(1/nfeat(Family)), nTotal = sum(TotalAbundance)) %>%
  mutate(pct = nTotal/ sum(nTotal) * 100) %>%
  arrange(n) %>%
  mutate(cs = cumsum(n), csR = cumsum(nR)) -> dataV 

dataV %>%
  filter(divergence > 0.5) %>% # pull(Phylum)
  summarise(cor(nR, n, method = "pearson"),
            cor(divergence, nR, method = "pearson"),
            cor(divergence, nTotal))

dataV %>%
  mutate(facet = ifelse(csR > 10, "A", "B")) %>%
  ggplot(aes(divergence, csR)) +
  geom_point(aes(size = pct)) +
  scale_size(name = "Abundance (%)") +
  # ggsci::scale_color_gsea() +
  ggrepel::geom_text_repel(aes(label = Phylum), size = 4) +
  facet_grid(facet~group, scales = "free_y",space = "fixed") +
  # scale_color_manual(values = c('black', 'red')) +
  labs(x = "x", y = 'y') +
  theme_bw(base_size = 16) +
  theme(legend.position = "top",
        panel.border = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) -> p1

ggsave(p1, filename = "divergence.png", path = dir, 
       width = 12, height = 9)

# Now lets investigate low prevelance/abundance phylum and subset them out.

df %>%
  group_by(Phylum, group) %>%
  summarise(cor = cor(Prevalence, TotalAbundance, method = "pearson"),
            mean_prevalence = mean(Prevalence), 
            total_abundance = sum(TotalAbundance)) %>%
  drop_na(Phylum) %>%
  arrange(desc(total_abundance)) -> summary_prevalence

# OrderPhylumLevel <- summary_prevalence$Phylum 
# Using the table above, determine the phyla to filter based on the 0.001 threshold
threshold <- 0.001

sum(summary_prevalence$total_abundance)*threshold
table(summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= threshold)
keepPhyla <- summary_prevalence$Phylum[summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= threshold]

keepPhyla <- unique(keepPhyla)
summary_prevalence <- summary_prevalence[summary_prevalence$Phylum %in% keepPhyla,]

summary_prevalence 

# 

# Individual taxa filtering
# Subset to the remaining phyla by prevelance.

df %>%
  filter(TotalAbundance > 1) %>%
  filter(Phylum %in% keepPhyla) %>%
  group_by(Phylum) %>%
  summarise(cor = cor(Prevalence, TotalAbundance))


df %>%
  mutate(is_sp = ifelse(rank == 'Species', TRUE, FALSE)) %>%
  mutate(Phylum= ifelse(Phylum %in% keepPhyla, Phylum, "Low Taxa")) %>%
  mutate(Phylum = factor(Phylum, levels = c(keepPhyla, "Low Taxa"))) %>%
  mutate(Prevalence = Prevalence/ncol(ab)) %>%
  ggplot(aes(TotalAbundance, Prevalence, color = is_sp)) +
  geom_point(size = 2, alpha = 0.7) + 
  # geom_hline(yintercept = 0.10, alpha = 0.5, linetype = 2) +
  # geom_vline(xintercept = 0.1, alpha = 0.5, linetype = 2) +
  scale_x_log10() +
  # geom_smooth(se = F, color = "red") +
  labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Sp. resolution") +
  facet_wrap(~Phylum) +
  theme_bw(base_size = 17) +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> p2

ggsave(p2, filename = "Prevalence.png", path = dir, width = 8, height = 6)

# network!! ----

cocurrance_filtering_sp <- paste0(dir, 'cocurrance_filtering_sp.list') 

# filter by cocurrance, define wether or not use loberas nor marker type to make the network
# agglomerate rdp-genus rdp-sp and pplacer-sp

from_pyseq_to_wTO <- function(phyloseq, overlap = NULL) {
  
  phyloseq %>% prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
    otu_table() %>%
    as.data.frame() -> metagenomics_abundance 
  
  if (is.null(overlap)) {
    wTO <- wTO.fast(Data = metagenomics_abundance,
                    Overlap = row.names(metagenomics_abundance),
                    method = 's', sign = 'sign', n = 250, 
                    method_resampling = 'Bootstrap')
    
  } else
  
  wTO <- wTO.fast(Data = metagenomics_abundance,
                  Overlap = overlap,
                  method = 's', sign = 'sign', n = 250, 
                  method_resampling = 'Bootstrap')
  # method_resampling = 'BlockBootstrap', lag = 2
  
  return(wTO)
  
}

wTO_V67 <- phyloseq  %>% from_pyseq_to_wTO(.) # overlap = Pathogens, if only pathogens must be evaluated
wTO <- wTO_V67

saveRDS(wTO_V67, file = paste0(dir, '/wTO_V67.rds'))

# tax %>% as_tibble(rownames = 'id')


wTOcutoff <- function(wTO_out, cutoff = 0.01, tau = 0.5) {
  
  n <- length(unique(c(wTO_out$Node.1, wTO_out$Node.2)))
  
  wTO_out %>% mutate_if(is.factor, as.character) %>%
    mutate(wTO = ifelse(pval-wTO < cutoff, wTO, 0 )) %>%
    filter(wTO != 0 ) %>%
    filter(abs(wTO) > tau) %>%
    # as.data.frame() -> out
    as_tibble() -> out
  cat("Dimension of df:", nrow(out), "\n")
  cat("Number of significant nodes interacting: ", length(unique(c(out$Node.1, out$Node.2))))
  cat("\nProportion of : ", length(unique(c(out$Node.1, out$Node.2)))/ n, "\n")
  
  return(out)
  
}
exportNet <- function(WTO, cutoff = 0.01, tau = 0.5) {
  
  WTO %>% wTOcutoff(cutoff, tau) -> WTO
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  nodes <- data.frame(id = sort(unique(c(Node.1,Node.2))))
  
  # nodes <- nodes %>% left_join(df, by = c('id' = 'Feature_IDs')) %>% as_tibble()
  
  # phyloseq %>% phyloseq::tax_table() %>% as.data.frame() %>% as_tibble(rownames = 'id') -> tax
  # phyloseq %>% phyloseq::otu_table() %>% as.data.frame() %>% as_tibble(rownames = 'id') -> ab
  
  # nodes <- nodes %>% left_join(tax, by = c('id')) %>% as_tibble()
  
  # nodeInfo %>% select(id, sampleName) -> scatterNodes
  # pivot_longer(cols = sampleName, names_to = 'Sample', values_to = 'Abundance') %>%
  # group_by(Sample, id) %>% summarise(Abundance = sum(Abundance)) %>%
  # # group_by(OTU) %>% mutate(Abundance = Abundance/sum(Abundance)) %>%
  # pivot_wider(names_from = Sample, values_from = Abundance) -> scatterNodes
  
  graph <- tbl_graph(nodes = nodes, edges = WTO, directed = FALSE)
  
  graph %>% activate(nodes) %>%
    mutate(degree = centrality_degree()) %>%
    mutate(value = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(weight = value * 2 + 1) -> graph
  
  gr <- function(x) {ifelse(x > 0, '+', '-')}
  
  graph %>% activate("edges") %>%  mutate(color = gr(wTO)) -> graph
  
  
  return(graph)
  
}


wTO_V67 %>%
  mutate(delta = pval-wTO) %>%
  # mutate(wTO = ifelse(pval-wTO < 0.01, wTO, 0 )) %>%
  # wTOcutoff() %>%
  ggplot(aes(color = delta)) +
  geom_point(aes(wTO, pval)) +
  ggsci::scale_color_gsea()

summary(wTO_V67$wTO)

library(tidygraph)

exportNet(wTO_V67, cutoff = 1, tau = 0) -> graph

# Hub detection as Layeghifard, M.,et al 2019
# top-ranked network hub taxa:
# We applied the PageRank algorithm to the microbiome networks to identify key members of CF lung microbiome in each patient. PageRank is a link analysis algorithm with the underlying assumption that hubs are likely to be more connected to other nodes when compared to non-hub nodes. This algorithm assigns a numerical weight to each node of a network as a measure of its relative importance within the network. The numerical weights of the nodes (also called PageRanks of the nodes) were then sorted and the top five or ten taxa with highest PageRank were selected as microbiome network hubs (or key taxa)


graph %>% activate("nodes") %>% mutate(betweenness = betweenness(.),
                                       membership = igraph::cluster_louvain(.)$membership,
                                       pageRank = page_rank(.)$vector) -> graph


plot(V(graph)$betweenness, V(graph)$pageRank)
quantile(V(graph)$pageRank,  probs = 0.75) -> cutOff
# quantile(V(graph)$betweenness, probs = 0.75) -> cutOff

graph %>% activate(nodes) %>%
  filter(pageRank > cutOff) -> g
  # filter(betweenness > cutOff) %>%
  # mutate(membership = igraph::cluster_louvain(.)$membership) -> g

table(V(g)$membership)
quantile(abs(E(g)$wTO))
quantile(abs(E(g)$wTO), probs = 0.75) -> EdgeCut

g %>% activate("edges") %>%
  mutate(wTO = ifelse(abs(wTO) > EdgeCut, wTO, 0)) %>%
  filter(wTO != 0) -> g

# g %>% activate('nodes') %>% filter(Relationship %in% 'P') -> g

# scatterDat <- as_data_frame(g, "vertices") 
# data.frame(row.names = scatterDat$id, scatterDat[, names(scatterDat) %in% sampleName]) -> scatterDat

rshipCol = c('P'='#6a3d9a', 'NC'='#ffff99', 'C'='#4daf4a')

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

ggraph(layout) +
  ggforce::geom_mark_hull(
    aes(x, y, group = as.factor(membership), fill = as.factor(membership)), #  fill = as.factor(membership)
    color = NA,
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25)  +
  guides(fill = FALSE) -> psaveNet

psaveNet + # psaveNet
  geom_edge_link(aes(edge_alpha = abs(wTO), color = color), edge_width = 1)+
  geom_node_point(aes(size = weight, color = Relationship)) + # , alpha = GS, color = abs(GS), 
  geom_node_text(aes(label = Pplacer.classification), repel = TRUE, size = 2) +
  scale_edge_width(range = c(0.2, 2)) +
  scale_color_manual(values = rshipCol) +
  scale_edge_color_manual(values = c('#e41a1c', '#377eb8')) +
  # scale_color_viridis_c(name = 'weight', direction = -1) +
  # scale_color_continuous(name = 'Gene Correlation') + # values = nodeCol
  theme_graph(base_family = "GillSans") +
  # guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top")+
  coord_fixed() + 
  facet_nodes(~ membership)

# test later
