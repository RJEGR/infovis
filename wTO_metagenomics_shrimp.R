rm(list = ls())

library(wTO)
require(CoDiNA)
library(microbiome)

library(tidyverse)

dir <- '~/Documents/Shrimp_Estefany/'

ancom_res <- read_tsv(paste0(dir, "ANCOMBC_res.tsv"))

# load(paste0(dir, "euler_outputs.RData"))
load(paste0(dir, "objects.RData"))


tax %>% select_at(ranks[2:5]) %>% distinct(Family, .keep_all = T)

phyloseq <- readRDS(paste0(dir, 'phyloseq_ancom.rds'))
phyloseq <- phyloseq %>% subset_taxa(Phylum %in% keepPhyla) %>% 
  prune_taxa(taxa_sums(.) > 0, .) %>% aggregate_taxa(., "Order")
phyloseq

phyloseq %>% tax_table %>% as.data.frame() %>% distinct(unique, .keep_all = T) -> tax
# trabjar a nivel familia, es 'particionar' los taxones en nodos que , separados, no son significativamente estadisticos, la prueba wTO no es tan potente para clasificar esta red. Por lo que conviene trabajar a un nivel mas alto como lo es el filo o quiza el Order (aunque mas complicado de describir)

datax <- unique(ancom_res$taxon_id)
phyloseq %>% subset_taxa(Family %in% taxon_g$Hindgut) %>% subset_samples(Tissue == 'Hindgut') %>% prune_taxa(taxa_sums(.) > 0, .) %>% aggregate_taxa(., "Phylum")


  # subset_taxa(Family %in% keepPhyla) %>%
  # prune_taxa(taxa_sums(.) > 0, .)

from_pyseq_to_wTO <- function(phyloseq) {
  
  phyloseq %>% prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
    otu_table() %>%
    as.data.frame() -> metagenomics_abundance 
  
  wTO <- wTO.fast(Data = metagenomics_abundance,
                  Overlap = row.names(metagenomics_abundance),
                  method = 's', sign = 'sign', n = 250, 
                  method_resampling = 'BlockBootstrap', lag = 2)
  
  return(wTO)
  
}

# usando overlaps, no se consigen valores p significativos,
# usamos toda la matriz como overlap
# ----
ancom_res %>% 
  separate(., col = wrap, into = c('-group', '+group'), sep = '-') %>%
  mutate(group = ifelse(logFC > 0, `+group`, `-group`)) %>% 
  distinct(group, taxon_id) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_wider(names_from = group, values_from = taxon_id) -> taxon_g

# ancom_res %>% distinct(taxon_id) %>% pull() -> overlaps
overlaps <- taxon_g$Hindgut
overlaps <- overlaps[!is.na(overlaps)]

# ----
# If our package is used with metagenomics data, for instance for analyzing co-occurrence networks, we recommend the abundance data to be normalized per day/ sample.

# par(mfrow = c(3,3))
# for ( i in 1:nrow(metagenomics_abundance)){
#   acf(t(metagenomics_abundance[i,]))
# }


table(sample_data(phyloseq)$Tissue)

wTO_Hindgut <- phyloseq %>% subset_samples(Tissue == 'Hindgut') %>% from_pyseq_to_wTO(.)
wTO_Midgut <- phyloseq %>% subset_samples(Tissue == 'Midgut') %>% from_pyseq_to_wTO(.)
wTO_Foregut <- phyloseq %>% subset_samples(Tissue == 'Foregut') %>%from_pyseq_to_wTO(.)

summary(wTO_Hindgut$pval.adj);summary(wTO_Midgut$pval.adj);summary(wTO_Foregut$pval.adj)

sum(wTO_Hindgut$pval.adj < 0.01);sum(wTO_Midgut$pval.adj < 0.01);sum(wTO_Foregut$pval.adj < 0.01)

# We found that x out of y taxa had at least one significant interaction (padj-value <0.01)

rbind(mutate(wTO_Hindgut, group = 'Hindgut'),
      mutate(wTO_Midgut, group = 'Midgut'),
      mutate(wTO_Foregut, group = 'Foregut')) %>%
  # ggplot(aes(pval.adj, fill = group)) + geom_histogram(bins = 60) + geom_vline(xintercept = 0.05)
  mutate(wTO = ifelse(pval.adj-abs(wTO) < 0.01, wTO, 0 )) %>%
  # filter(pval < 0.05) %>%
  ggplot(aes(color = pval.adj)) +
  geom_point(aes(pval.adj, wTO)) +
  # geom_vline(xintercept = 0.05, color = '#084594', linetype = 'dashed') +
  # geom_hline(yintercept = 0, color = '#084594', linetype = 'dashed') +
  ggsci::scale_color_gsea() +
  facet_grid(~group)

# how to interpretate: The nodes within x cluster contains only negative interactions (green links), suggesting that the bacterial species in this cluster do not co-exist. We also notice, that many of the bacteria belonging to the same level_taxa are well connected by purple links, indicating that they co-exist and share interactions. However, the number of interactions among non-related bacteria demonstrate that interactions are not intra-level_taxa specific. Positive correlations in co-occurence networks may represent symbiotic or commensal relationships, while negative correlations may represent predator-prey interactions, allelopathy or competition for limited resources

wTOcutoff <- function(wTO_out, cutoff = 0.01) {
  
  
  wTO_out %>% mutate_if(is.factor, as.character) %>%
    mutate(wTO = ifelse(pval.adj-abs(wTO) < cutoff, wTO, 0 )) %>%
    filter(wTO != 0 ) %>%
    as.data.frame() -> out
  
 cat("Number of significant nodes interacting: ", length(unique(c(out$Node.1, out$Node.2))))
  
  return(out)
  
}
wTO_Hindgut %>% wTOcutoff() -> HindgutN
wTO_Midgut %>% wTOcutoff()  -> MidgutN
wTO_Foregut %>% wTOcutoff() -> ForegutN

sum(HindgutN$pval.adj < 0.01);sum(MidgutN$pval.adj < 0.01);sum(ForegutN$pval.adj < 0.01)

length(unique(c(HindgutN$Node.1, HindgutN$Node.2)))
length(unique(c(HindgutN$Node.1, HindgutN$Node.2)))
length(unique(c(HindgutN$Node.1, HindgutN$Node.2)))

rbind(mutate(HindgutN, group = 'Hindgut'),
      mutate(MidgutN, group = 'Midgut'),
      mutate(ForegutN, group = 'Foregut')) %>%
  ggplot(aes(color = pval.adj)) +
  geom_point(aes(pval.adj-wTO, wTO)) +
  ggsci::scale_color_gsea() +
  facet_grid(~group)

plotNet <- function(WTO, tau = 0.5) {
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  wTO <- WTO$wTO
  padj <- WTO$pval.adj
  cutoff <- list(kind = "Threshold", value = tau)
  # cutoff <- list(kind = "pval", value = 0.05)
  
  NetVis(Node.1, Node.2, wTO, padj = padj, cutoff = cutoff, MakeGroups = 'louvain')
  
  # MakeGroups should be FALSE or one of the following options: 'walktrap', 'optimal', 'spinglass', 'edge.betweenness', 'fast_greedy', 'infomap', 'louvain', 'label_prop', 'leading_eigen'. 
}

# plotNet(HindgutN, tau = 0.3)
# plotNet(ForegutN, tau = 0.3) # wTO_Foregut
# plotNet(MidgutN, tau = 0)
# 

aesthNet <- function(wTO, tau = 0.3) {
  
  library(ggraph)
  library(ggforce)
  
  
  g <- plotNet(wTO, tau = tau)
  
  g$Nodes %>% left_join(tax, by = c('id' = 'unique')) -> Nodes
  
  Nodes <- Nodes %>% mutate(size = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(size = size * 2 + 1)
  
  Edges <- wTO %>% filter(abs(wTO) > tau) %>% 
    mutate(wTOc = ifelse(wTO > 0, '+', '-'), 
           width = 0.5 + 5 * abs((wTO - min(wTO))/(max(wTO) -  min(wTO))))
  graph = graph_from_data_frame(Edges, directed = FALSE, Nodes)
  layout = create_layout(graph, layout = 'igraph', algorithm = 'kk')
  
  # test aesthetics https://www.r-bloggers.com/2020/03/ggraph-tricks-for-common-problems/
  
  ggraph(layout) + 
    geom_edge_link(aes(edge_colour = wTOc, edge_alpha = width), width = 1.2,
                   arrow = arrow(
                     angle = 10,
                     length = unit(0.1, "inches"),
                     ends = "last",
                     type = "closed"
                   )) +  # aes(alpha = Score)
    geom_node_point(aes(size = degree * 2 + 1)) + # color = as.character(group)), color = Phylum
    geom_node_text(aes(label = name), repel = TRUE) +
    geom_mark_hull(
      aes(x, y, group = group, label=group),
      fill = "grey", color = NA,
      concavity = 4,
      con.size = 0.3,
      con.linetype = 2,
      expand = unit(2, "mm"),
      alpha = 0.25) +
    theme_void() + 
    theme(legend.position="top")
# 
#   pg + geom_edge_link(
#     # aes(end_cap = circle(node2.degree + 0.5, "pt"),
#         # edge_colour = wTOc, width = abs(wTO)),
#     arrow = arrow(
#       angle = 10,
#       length = unit(0.1, "inches"),
#       ends = "last",
#       type = "closed"
#     )) 
  
}

phyloseq %>% tax_table() %>% as.data.frame() %>% pull(Phylum) %>% sort() %>% unique() -> labels 

colourCount = length(labels)

library(ggsci)

if(colourCount > 7) {
  getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
} else
  getPalette <- pal_locuszoom(alpha = 1)(colourCount)  

names(getPalette) <- labels

p1 <- aesthNet(HindgutN, tau = 0.3) + geom_node_point(aes(color = Phylum)) + 
  guides(col = guide_legend(nrow = 3)) + scale_color_manual(values = getPalette)
p2 <- aesthNet(ForegutN, tau = 0.3) + geom_node_point(aes(color = Phylum)) + 
  guides(col = guide_legend(nrow = 3)) + scale_color_manual(values = getPalette)
p3 <- aesthNet(MidgutN, tau = 0) + geom_node_point(aes(color = Phylum)) + 
  scale_color_manual(values = getPalette)

ggsave(p1, filename = "HindgutWTO.png", path = dir, width = 12, height = 8)
ggsave(p2, filename = "ForegutWTO.png", path = dir, width = 12, height = 8)
ggsave(p3, filename = "MidgutWTO.png", path = dir, width = 12, height = 8)

# library(patchwork)
# p1+p2
# if differential net ----


DiffNet <- MakeDiffNet(Data = list(HindgutN, MidgutN,ForegutN),
                       Code = c('Hindgut','Midgut', 'Foregut'))

CoDiNA::plot.CoDiNA(DiffNet, Cluster = TRUE, sort.by.Phi = TRUE)

# continue w/ https://deisygysi.github.io/rpackages/Pack-2

# Clustering the nodes into Φ and Φ̃ categories
# based on median

int_C = quantile(DiffNet$Score_internal, 0.3)
ext_C = quantile(DiffNet$Score_Phi, 0.3)

Nodes_Groups = ClusterNodes(DiffNet = DiffNet, 
                            cutoff.external = ext_C, 
                            cutoff.internal = int_C)
table(Nodes_Groups$Phi_tilde)

Graph = plot(DiffNet, cutoff.external = ext_C, cutoff.internal = int_C, layout = 'layout_components')

g = CoDiNA::as.igraph(Graph) 

library(igraph)
library(RColorBrewer)
vgroup <- as.factor(Graph$Nodes$Phi_tilde)
color_pal <- brewer.pal(4, "Set1")
vertex.color <- color_pal[as.numeric(vgroup)]
names(vertex.color) <- Graph$Nodes$Phi_tilde

plot(g, layout = layout.fruchterman.reingold(g), vertex.color=vertex.color)
legend("bottomleft", legend= levels(vgroup)  , 
       col = color_pal , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col= color_pal , horiz = FALSE, inset = c(0.1, 0.1))

library(ggraph)

g = graph_from_data_frame(Graph$Edges, directed = FALSE, Graph$Nodes)
layout = create_layout(g, layout = 'igraph', algorithm = 'nicely')

ggraph(layout) + 
  # geom_edge_density(aes(fill = Score)) +
  geom_edge_link(alpha = 0.25, edge_colour = "grey") +  # aes(alpha = Score)
  geom_node_point(aes(size = Degree_Total, alpha = Degree_Total, color = Phi_tilde)) +
  geom_node_text(aes(label = name)) +
  scale_color_manual(values = color_pal) -> graphSave
  # facet_edges(~Group)
  # facet_nodes(~ Phi)
ggsave(graphSave, filename = "codina.png", path = dir, 
       width = 8, height = 8)

# or chord

rbind(HindgutN, MidgutN, ForegutN)


df <- wTO_Hindgut %>% filter(abs(wTO) > 0.5)

library(circlize)

set.seed(260220)

samGroup <- phyloseq  %>% sample_data() %>% pull(Tissue) %>% unique()

n <- length(samGroup)
grid.col <- ggsci::pal_rickandmorty(alpha = 0.8)(n)
names(grid.col) <- c("Hindgut", "Midgut", "Foregut")


df %>%
  mutate_if(is.factor, as.character) %>%
  # mutate(from = ifelse(wTO < 0, Node.2, Node.1)) %>%
  # mutate(to = ifelse(wTO > 0, Node.2, Node.1))
  select(Node.1, Node.2) %>%
  mutate(color = 'black') %>%
  as.data.frame() -> arr.col

circos.clear()
circos.par(start.degree = 0, gap.degree = 4, 
           track.margin = c(-0.01, 0.01), 
           points.overflow.warning = FALSE)
df %>% 
  select(Node.1, Node.2, wTO) %>%
  # with(., table(Node.1, Node.2)) %>%
  chordDiagram(
    #grid.col = c(grid.col, getPalette),
               directional = -1,
               diffHeight = mm_h(5), target.prop.height = mm_h(4),
               direction.type = "arrows",
               link.arr.col = arr.col, 
               link.arr.length = 0.2,
               preAllocateTracks = 1,
               small.gap = 10, big.gap = 15)
