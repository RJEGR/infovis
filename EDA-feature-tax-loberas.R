
# red de interacion microbioma
# test https://www.nature.com/articles/s41522-018-0077-y
rm(list = ls())

library(wTO)
require(CoDiNA)
library(microbiome)
library(Biostrings)

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

ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'pplacer_sp') # 'Species'
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

saveRDS(featureDF, file = paste0(dir, '/featureDF.rds'))

names(featureDF)[9] <- 'pplacer_sp'

df <- readCSV(fileNames[2])
df %>% select_if(is.integer) -> ab

df %>% filter(Relationship %in% 'P') %>% distinct() %>% pull(Feature_IDs) -> Pathogens

# df %>% select_at(ranks) %>% data.frame(row.names = df$Feature_IDs) -> tax

df %>% select_at(all_of(ranks)) %>% 
  mutate(id = 1:nrow(df)) %>%
  pivot_longer(cols = ranks) %>% fill(value) %>%
  pivot_wider(names_from = name) %>%
  select(-id) %>%
  data.frame(row.names = df$Feature_IDs) -> tax

dat <- df %>% select_if(is.integer) %>% data.frame(row.names = df$Feature_IDs)

# identical(names(dat), rownames(sam))
identical(rownames(dat), rownames(tax))

# and parse
phyloseq = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix'))) 

# transform_sample_counts(function(x) sqrt(x / sum(x)))
saveRDS(phyloseq, file = paste0(dir, 'phyloseq_loberas.rds'))


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

featureDF %>% select_at(sampleName) -> ab

prevelancedf = apply(X = ab,
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

# identical(rownames(tax), rownames(ab)) # sanity check

r <- c(0:7)
names(r) <- ranks
data.frame(Resolution = r) %>% as_tibble(rownames = "rank") -> r


df = data.frame(Prevalence = prevelancedf, 
                TotalAbundance = rowSums(ab),
                tax, group = featureDF$group, Relationship = featureDF$Relationship) %>% as_tibble(rownames = 'id')

df %>% left_join(r) %>% select(-Resolution) %>% mutate(rank = factor(rank, levels = ranks)) -> df

df %>% filter(group %in% 'V3') %>% drop_na(Order)

df %>% 
  group_by(rank, group) %>%
  tally() %>% 
  group_by(group) %>%
  mutate(nn = cumsum(n) - n, pct = 1 - nn/sum(n), cs = sum(n)-nn) %>%
  ggplot(aes(x = rank, y = pct, color = cs, group = group)) +
  geom_path(size = 2, alpha=0.6) +
  geom_point(size=2, alpha=0.6) +
  geom_text(aes(label = cs), size=4, vjust = -1, color = 'black') +
  scale_color_viridis_c(option = "A", direction = -1) +
  labs(y = "Changes in the assignation (%)", x = 'rdp classification', color = "Number of\nFeatures") +
  facet_grid(~group) +
  theme_bw(base_size = 17) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave

# Gracias Papa
ggsave(psave, filename = "resolution.png", path = dir, width = 14, height = 7)

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
       width = 10, height = 7)

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

featureDF %>% 
  pivot_longer(sampleName, values_to = 'Abundance', names_to = 'Sample') %>%
  select(Feature_IDs, Abundance, Sample, ranks) %>%
  filter(Abundance > 0) %>%
  group_by()
  pivot_longer(ranks, names_to = 'Rank', values_to = 'name') -> out

out %>% group_by(Rank, Sample) %>% summarise(nreads = sum(Abundance)) -> nr 

nasv <- function(x) {length(unique(x))}

out %>% group_by(Rank, Sample) %>% count() -> nf
# nf <- aggregate(out[,'Feature_IDs'], by=list(out$Rank, out$Sample), FUN = nasv)
ng <- aggregate(out[,'name'], by=list(out$Rank, out$Sample), FUN = nasv)

n <- left_join(nr, nf)

ggplot() +
  theme(panel.background = element_rect(fill = "transparent")) +
  geom_point(data = n, 
             aes(x=Rank, y = n, shape = Sample, group = Sample),
             size = 3, alpha = 1) +
  geom_path(data = n, aes(x = Rank, y = n, group = Sample))

ggplot(data = n) +
  theme(panel.background = element_rect(fill = "transparent")) +
  geom_point(aes(x = Rank, y = n, color = nreads, group = Sample, shape = Sample),
             size = 3, alpha = 1) + 
  geom_path(data = n, aes(x = Rank, y = ngroup, group = Ship))

reads_n_features <- function(phyloseq, ...) {
  
  # psample <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)
  
  # out <- psmelt(psample)
  out <- filter(out, Abundance > 0) %>% select(OTU, Abundance, Sample, ranks)
  out <-reshape2::melt(out, 
                       measure.vars = ranks, 
                       variable.name = 'Rank', value.name='name')
  
  nr <- aggregate(out[,'Abundance'], 
                  by=list(out$Rank, out$Sample), FUN = sum)
  
  nasv <- function(x) {length(unique(x))}
  
  nf <- aggregate(out[,'Feature_IDs'], by=list(out$Rank, out$Sample), FUN = nasv)
  
  ng <- aggregate(out[,'name'], by=list(out$Rank, out$Sample), FUN = nasv)
  
  total <- sum(sample_sums(psample)) 
  
  pct <- round(nr[,3] / total * 100, 3)
  
  # sanity check 
  if(identical(nr[,1], nf[,1])) {
    
    n <- data.frame(Rank = nr[,1], Sample = nr[,2], 
                    nreads = nr[,3], pct = pct, nasvs = nf[,3], 
                    ngroup = ng[,3])
  } else
    n <- data.frame(Rank = nr[,1], Sample = nr[,2], nreads = nr[,3], 
                    pct = pct, Rank = nr[,1], nasvs = nf[,3], 
                    ngroup = ng[,3])
  
  return(n)
  
}

n <-reads_n_features(phyloseq)

nasv <- function(x) {length(unique(x))}


df %>%
  mutate(Resolution = as.factor(Resolution)) %>%
  select(id, TotalAbundance, Resolution, group) %>% 
  group_by(Resolution, group) %>%
  summarise(ab = sum(TotalAbundance), n = nasv(id)) %>%
  arrange(group) %>%
  ggplot() +
  theme(panel.background = element_rect(fill = "transparent")) +
  geom_point(aes(x = as.factor(Resolution), y = n, color = ab, group = group, shape = group),
             size = 3, alpha = 1) +  
  geom_path(aes(x = as.factor(Resolution), y = n, group = group, color = ab))

# Individual taxa filtering
# Subset to the remaining phyla by prevelance.

df %>%
  filter(TotalAbundance > 1) %>%
  filter(Phylum %in% keepPhyla) %>%
  group_by(Phylum) %>%
  summarise(cor = cor(Prevalence, TotalAbundance))


df %>%
  mutate(is_sp = ifelse(Resolution == 6, TRUE, FALSE)) %>%
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
  
  nodes <- nodes %>% left_join(df, by = c('id' = 'Feature_IDs')) %>% as_tibble()
  
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

