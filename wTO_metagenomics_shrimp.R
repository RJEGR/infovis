rm(list = ls())

library(wTO)
require(CoDiNA)

dir <- '~/Documents/Shrimp_Estefany/'

ancom_res <- read_tsv(paste0(dir, "ANCOMBC_res.tsv"))

load(paste0(dir, "euler_outputs.RData"))

phyloseq <- readRDS(paste0(dir, 'phyloseq_ancom.rds'))

phyloseq <- phyloseq %>% aggregate_taxa(., "Family") 
  # subset_taxa(Family %in% intersected_taxa) %>%
  # prune_taxa(taxa_sums(.) > 0, .)



from_pyseq_to_wTO <- function(phyloseq) {
  
  phyloseq %>% prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
    otu_table() %>%
    as.data.frame() -> metagenomics_abundance 
  
  wTO <- wTO.fast(Data = metagenomics_abundance,
                  Overlap = row.names(metagenomics_abundance),
                  method = 'p', sign = 'sign', n = 250, 
                  method_resampling = 'BlockBootstrap', lag = 2)
  
  return(wTO)
  
}


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

sum(wTO_Hindgut$pval.adj < 0.05);sum(wTO_Midgut$pval.adj < 0.05);sum(wTO_Foregut$pval.adj < 0.05)

wTO_Foregut %>%
  mutate(wTO = ifelse(pval.adj-wTO < pval.adj, wTO, 0 )) %>%
  # filter(pval < 0.05) %>%
  ggplot(aes(color = pval.adj)) +
  geom_point(aes(pval, abs(wTO)))


Node.1 = as.character(wTO_Foregut$Node.1)
Node.2 = as.character(wTO_Foregut$Node.2)
wTO <- wTO_Foregut$wTO

NetVis(Node.1 = Node.1, Node.2 = Node.2, wTO)

wTO_Hindgut %>% mutate(wTO = ifelse(pval.adj-wTO < pval.adj, wTO, 0 )) -> HindgutN
wTO_Midgut %>% mutate(wTO = ifelse(pval.adj-wTO < pval.adj, wTO, 0 )) -> MidgutN
wTO_Foregut %>% mutate(wTO = ifelse(pval.adj-wTO < pval.adj, wTO, 0 )) -> ForegutN

DiffNet = MakeDiffNet(Data = list(HindgutN, MidgutN, ForegutN),
                       Code = c('Hindgut','Midgut', 'Foregut'))

plot(DiffNet)

# continue w/ https://deisygysi.github.io/rpackages/Pack-2

# Clustering the nodes into Φ and Φ̃ categories
# based on median

int_C = quantile(DiffNet$Score_internal, 0.5)
ext_C = quantile(DiffNet$Score_Phi, 0.5)

Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = ext_C, cutoff.internal = int_C)
table(Nodes_Groups$Phi_tilde)

Graph = plot(DiffNet, cutoff.external = ext_C, cutoff.internal = int_C, layout = 'layout_components')

g = CoDiNA::as.igraph(Graph) 

library(igraph)
library(RColorBrewer)
vgroup <- as.factor(Graph$Nodes$Phi_tilde)
color_pal <- brewer.pal(4, "Set1")
vertex.color <- color_pal[as.numeric(vgroup)]
names(vertex.color) <- Graph$Nodes$Phi_tilde

plot(g, layout = layout.fruchterman.reingold(g), vertex.color=my_color)
legend("bottomleft", legend= levels(vgroup)  , 
       col = color_pal , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col= color_pal , horiz = FALSE, inset = c(0.1, 0.1))

library(ggraph)



g = graph_from_data_frame(Graph$Edges, directed = FALSE, Graph$Nodes)
layout = create_layout(g, layout = 'igraph', algorithm = 'nicely')

ggraph(layout) + 
  geom_edge_density(aes(fill = Score)) +
  geom_edge_link(alpha = 0.25, edge_colour = "grey") +  # aes(alpha = Score)
  geom_node_point(aes(size = Degree_Total, alpha = Degree_Total, color = Phi_tilde)) +
  geom_node_text(aes(label = name)) +
  scale_color_manual(values = color_pal) -> graphSave
  # facet_edges(~Group)
  # facet_nodes(~ Phi)
ggsave(graphSave, filename = "codina.png", path = dir, 
       width = 8, height = 8)
