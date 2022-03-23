# From file run018_test01_ASVs_animalia.taxonomy start practicing networks

path <- "~/metagenomics/MG_COI/run018/"

# countf <- '^run018_test01_ASVs_animalia_count.table'

# taxf <- 'run018_test01_ASVs_animalia.taxonomy'

psf <- 'run018_test01_ASVs_animalia_physeq.rds'

psf <- list.files(path, pattern = psf,  full.names = TRUE)


library(tidyverse)
library(phyloseq)

ps <- readRDS(psf)

taxName <- colnames(tax_table(ps))[3]
# Which are a major groups?

ps.ord <- ordinate(ps, "NMDS","bray")

p1 = plot_ordination(ps, ps.ord, type="samples", color = 'Tipo', title="taxa")

print(p1)

# PCA
count_feat <- ps %>% 
  transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
  otu_table()

PCA <- prcomp(t(count_feat), scale. = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAviz <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

ps %>% sample_data() %>% as_tibble(rownames = 'id') -> mtd

g <- mtd %>% filter(grepl('Developmetal', Sample.type )) %>% pull(stage)

n <- length(unique(g))
grid.col <- c(ggsci::pal_d3()(10), ggsci::pal_aaas()(n-10))
grid.col <- structure(grid.col, names = stagesLev)

PCAviz %>%
  mutate(id = rownames(.)) %>%
  left_join(mtd) %>%
  ggplot(., aes(PC1, PC2, color = Tipo)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # geom_text(aes(label = Description), alpha = 0.9) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') 
# scale_color_manual(values = grid.col) +
# coord_fixed(ratio = sd_ratio) +
# ggforce::geom_mark_ellipse(aes(group = stage, label = stage))

p1$data %>%
  ggplot(., aes(NMDS1, NMDS2, fill = Tipo, color = Tipo)) +
  # ggrepel::geom_text_repel(aes(label = Estación), alpha = 0.9)
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(legend.position = 'top') +
  ggforce::geom_mark_ellipse(
    expand = unit(2, "mm"),
    alpha = 0.25) +
  geom_point(size = 3, alpha = 0.7) 


library(wTO)

require(CoDiNA)



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

colnames(tax_table(ps))[6] -> taxName

ps %>% microbiome::aggregate_taxa(., taxName) -> ps_agg

wTO <- ps_agg %>% from_pyseq_to_wTO(.)


# wTOL <- ps_agg %>% subset_samples(Tipo == 'Larvas') %>% from_pyseq_to_wTO(.)
# wTOH <- ps_agg %>% subset_samples(Tipo == 'Huevos') %>% from_pyseq_to_wTO(.)

# rbind(mutate(wTOL, g = 'Larvas'), mutate(wTOH, g = 'Huevos')) %>%

wTO %>% ggplot(aes(pval.adj)) + geom_histogram(bins = 100)

wTO %>%
  mutate(wTO = ifelse(pval.adj-abs(wTO) < 0.01, wTO, 0 )) %>%
  ggplot(aes(color = pval.adj)) +
  geom_point(aes(wTO, pval.adj)) +
  geom_vline(xintercept = 0, color = '#084594', linetype = 'dashed') +
  geom_hline(yintercept = 0.01, color = '#084594', linetype = 'dashed') +
  ggsci::scale_color_gsea() 
# facet_grid(~ g)

wTOcutoff <- function(wTO_out, cutoff = 0.01) {
  
  
  wTO_out %>% mutate_if(is.factor, as.character) %>%
    mutate(wTO = ifelse(pval.adj-abs(wTO) < cutoff, wTO, 0 )) %>%
    filter(wTO != 0 ) %>%
    as.data.frame() -> out
  
  cat("Number of significant nodes interacting: ", length(unique(c(out$Node.1, out$Node.2))))
  
  return(out)
  
}

wTO %>% wTOcutoff(cutoff = 0.01) -> wTO_out

ps %>% tax_table() %>% as.data.frame() %>% pull(Filo) %>% sort() %>% unique() -> labels

library(ggsci)

c <- length(labels)

if(c > 7) {
  getPalette <- colorRampPalette(pal_locuszoom()(7))(c)
} else
  getPalette <- pal_locuszoom(alpha = 1)(c)  

getPalette <- structure(getPalette, names = labels)

plotNet <- function(WTO, tau = 0.05, type = 'Threshold') {
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  wTO <- WTO$wTO
  
  padj <- WTO$pval.adj
  
  cutoff <- list(kind = type, value = tau)
  
  # It's a list containing the kind of cutoff to be used (pval, Threshold or pval.adj)and it's value. Example: cutoff= list(kind = "Threshold", value = 0.5) || list(kind = "pval", value = 0.05)
  
  NetVis(Node.1, Node.2, wTO, padj = padj, cutoff = cutoff, MakeGroups = 'louvain')
  
  # MakeGroups should be FALSE or one of the following options: 'walktrap', 'optimal', 'spinglass', 'edge.betweenness', 'fast_greedy', 'infomap', 'louvain', 'label_prop', 'leading_eigen'. 
}

aesthNet <- function(wTO, ps_agg, tau = 0.05) {
  
  library(igraph)
  library(ggraph)
  library(ggforce)
  
  tax <- ps_agg %>% tax_table() %>% data.frame()
  
  g <- plotNet(wTO, tau = tau)
  
  g$Nodes %>% left_join(tax, by = c('id' = 'unique')) -> Nodes
  
  Nodes <- Nodes %>% mutate(size = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(size = size * 2 + 1)
  
  Edges <- wTO %>% filter(abs(wTO) > tau) %>% 
    mutate(wTOc = ifelse(wTO > 0, '+', '-'), 
      width = 0.5 + 5 * abs((wTO - min(wTO))/(max(wTO) -  min(wTO))))
  
  # graph = igraph::graph_from_data_frame(Edges, directed = FALSE, Nodes)
  
  # graph <- delete.edges(graph, E(graph)[wd[,3] < 0.2])
  
  library(tidygraph)
  
  graph <- tidygraph::tbl_graph(nodes = Nodes, edges = Edges, directed = FALSE)
  
  # hist(Edges$width)
  
  # Edges %>% filter(width >= quantile(width, 0.5))
  
  quant50 <- quantile(Edges$width, 0.5)
  
  graph <- graph %>% activate("edges") %>% mutate(width = ifelse(width > quant50, width, 0))
  
  graph <- graph %>% activate("edges") %>% filter(width > 0)
  
  layout = create_layout(graph, layout = 'igraph', algorithm = 'kk')
  
  # test aesthetics https://www.r-bloggers.com/2020/03/ggraph-tricks-for-common-problems/
  
  ggraph(layout) + 
    geom_edge_link(, width = 1.2, # aes(edge_colour = wTOc, edge_alpha = width)
      arrow = arrow(
        angle = 10,
        length = unit(0.1, "inches"),
        ends = "last",
        type = "closed"
      )) +  # aes(alpha = Score)
    geom_node_point(aes(size = degree * 2 + 1, color = Filo)) + 
    geom_node_text(aes(label = id), repel = TRUE) +
    geom_mark_hull(aes(x, y, group = group),
      fill = "grey", color = NA,
      expand = unit(2, "mm"),
      alpha = 0.25) +
    theme_void() + 
    theme(legend.position="top")
  
}

# plotNet(wTO_out, tau = 0.01, type = 'Threshold')

aesthNet(wTO_out, ps_agg, tau = 0.01) + 
  geom_node_point(aes(color = Filo)) + 
  guides(col = guide_legend(nrow = 3)) + scale_color_manual(values = getPalette)


tax <- ps_agg %>% tax_table() %>% data.frame()

g <- plotNet(wTO_out, tau = 0.01)

g$Nodes %>% left_join(tax, by = c('id' = 'unique')) -> Nodes

Nodes <- Nodes %>% mutate(size = (degree - min(degree))/(max(degree) - min(degree))) %>%
  mutate(size = size * 2 + 1)

Edges <- wTO %>% filter(abs(wTO) > 0.01) %>% 
  mutate(wTOc = ifelse(wTO > 0, '+', '-'), 
    width = 0.5 + 5 * abs((wTO - min(wTO))/(max(wTO) -  min(wTO))))

Nodes %>% as_tibble() %>% group_by(Filo) %>% 
  summarise(N = sum(degree), count = length(Familia)) %>% 
  arrange(N) %>% pull(Filo) -> FiloL

Nodes %>% mutate(Filo = factor(Filo, levels = FiloL)) %>% 
  ggplot(aes(x = Filo, y = degree)) + geom_col() + coord_flip()
