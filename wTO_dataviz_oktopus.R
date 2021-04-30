
# 1. se probo la net-aesthetic, estamos solo usando datos up/anotados, funciona bien, 
# 2. Estamos evaluando sin significancia ajustada, pero usando el valor wTO_abs-pval
# 3. Based on testing spearman/pearson (method) and vst/raw (inputs), we choose, raw/spearman as parameter for up/down degs network analysis

rm(list = ls())

library(tidyverse)

library(wTO)

library(ggraph)
library(igraph)
library(ggforce)
library(scatterpie)
library(RColorBrewer)

# path <- "~/transcriptomics/oktopus_full_assembly/"
# path <- "~/transcriptomics/oktopus_full_assembly/wTO/spearman_out/"
# path <- "~/transcriptomics/oktopus_full_assembly/wTO/wTO_subset/"

path <- "~/transcriptomics/oktopus_full_assembly/wTO/spearman_by_stages/"

pattern <- "DES"

rds <- list.files(path = path, pattern = pattern, full.names = T)

n1 <- readRDS(rds[1])$wTO
n2 <- readRDS(rds[2])$wTO
n3 <- readRDS(rds[3])$wTO

# qqnorm(n1$wTO_sign, pch = 1, frame = FALSE)
# qqline(n1$wTO_sign, col = "steelblue", lwd = 2)


rbind(mutate(n1, group = 'net1'),
      mutate(n2, group = 'net2'),
      mutate(n3, group = 'net3')) -> datViz


# testing method ----
# hist(n3$pval_sig);hist(n2$pval_sig);hist(n1$pval_sig)
# hist(n3$wTO_sign);hist(n2$wTO_sign);hist(n1$wTO_sign)

# rbind(mutate(n1, group = 'raw-pearson'),
#       mutate(n2, group = 'raw-spearman'),
#       mutate(n3, group = 'vst-pearson')) -> datViz

datViz %>%
  group_by(group) %>%
  sample_n(100) %>%
  ggplot(aes(pval_sig-wTO_abs, color = group)) +
  geom_vline(xintercept = 0.01, color = '#084594', linetype = 'dashed') +
  geom_density()

datViz %>%
  group_by(group) %>%
  sample_n(100) %>%
  ggplot(aes(wTO_sign, color = group)) +
  geom_density() + facet_grid(~group)

datViz %>% 
  mutate(wTO_sign = ifelse(pval_sig-wTO_abs < 0.01, wTO_sign, 0 )) %>%
  filter(wTO_sign != 0) %>%
  group_by(group) %>%
  sample_n(1000) %>%
  ggplot(aes(pval_sig, wTO_sign, color = pval_sig-wTO_abs)) + 
  ggsci::scale_color_gsea() +
  geom_point(alpha = 0.5) + 
  geom_hline(yintercept = 0.01, color = '#084594', linetype = 'dashed') +
  facet_grid(~group)

datViz %>% 
  mutate(wTO_sign = ifelse(pval_sig-wTO_abs < 0.01, wTO_sign, 0 )) %>%
  filter(wTO_sign != 0) %>%
  group_by(group) %>%
  sample_n(1000) %>%
  ggplot(aes(x = group, y = wTO_sign)) +
  geom_boxplot()

datViz %>%
  mutate(wTO_sign = ifelse(pval_sig-wTO_abs < 0.01, wTO_sign, 0 )) %>%
  filter(wTO_sign != 0 ) %>%
  group_by(group) %>%
  sample_n(1000) %>%
  ggplot(aes(color = pval_sig)) +
  geom_point(aes(wTO_sign, -log10(pval_sig))) +
  ggsci::scale_color_gsea() +
  facet_grid(~group)

wTOcutoff <- function(wTO_out, cutoff = 0.01, tau = 0.5) {
  
  n <- length(unique(c(wTO_out$Node.1, wTO_out$Node.2)))
              
  wTO_out %>% mutate_if(is.factor, as.character) %>%
    mutate(wTO_sign = ifelse(pval_sig-wTO_abs < cutoff, wTO_sign, 0 )) %>%
    filter(wTO_sign != 0 ) %>%
    filter(abs(wTO_sign) > tau) %>%
    # as.data.frame() -> out
    as_tibble() -> out
  cat("Dimension of df:", nrow(out), "\n")
  cat("Number of significant nodes interacting: ", length(unique(c(out$Node.1, out$Node.2))))
  cat("\nProportion of : ", length(unique(c(out$Node.1, out$Node.2)))/ n)
  
  return(out)
  
}

n1 %>% wTOcutoff(cutoff = 0.01, tau = 0.7) -> n1Cutoff
n2 %>% wTOcutoff(cutoff = 0.01, tau = 0.7) -> n2Cutoff
n3 %>% wTOcutoff(cutoff = 0.01, tau = 0.5) -> n3Cutoff

# hist(n1Cutoff$pval_sig)

gr <- function(x) {ifelse(x > 0, '+', '-')}

n1Cutoff %>% as_tibble() %>% group_by(Node.1) %>%
  # filter(Node.1 %in% 'TRINITY_DN56991_c0_g1') %>%
  summarise(pos = sum(gr(wTO_sign) != '+'), 
            neg = sum(gr(wTO_sign) != '-')) %>%
  mutate(type = ifelse(pos > neg, 'pos', 'neg')) %>%  view()
  # right_join(dff, by = c("Node.1"="ID")) %>% view()

plotNet <- function(WTO, tau = 0.5) {
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  wTO <- WTO$wTO_sign
  pval <- WTO$pval_sig
  cutoff <- list(kind = "Threshold", value = tau)
  # cutoff <- list(kind = "pval", value = 0.05)
  
  NetVis(Node.1, Node.2, wTO, pval = pval, cutoff = cutoff, MakeGroups = 'FALSE')
  
  # MakeGroups should be FALSE or one of the following options: 'walktrap', 'optimal', 'spinglass', 'edge.betweenness', 'fast_greedy', 'infomap', 'louvain', 'label_prop', 'leading_eigen'. 
}

g <- plotNet(n1Cutoff)

aesthNet <- function(wTO, tau = 0.5) {
  
  library(ggraph)
  library(ggforce)
  
  
  g <- plotNet(wTO, tau = tau)
  
  # g$Nodes %>% left_join(kegg, by = c('id' = 'gene')) -> Nodes # pseudo-code
  
  g$Nodes -> Nodes
  
  Nodes <- Nodes %>% mutate(size = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(size = size * 2 + 1)
  
  names(wTO)[3] <- 'wTO'
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
                     # ends = "last",
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
  
}


require(CoDiNA)

# Code <- c('Optic-Gland','Oviducal-Gland', 'Optic-Lobe')
Code <- c('Net1','Net2', 'Net3')
DiffNet = MakeDiffNet(Data = list(n1Cutoff, n2Cutoff, n3Cutoff),
                       Code = Code)

DiffNet

# test second


path <- "~/transcriptomics/oktopus_full_assembly/wTO/outputs/"
pattern <- "rds"

rds <- list.files(path = path, pattern = pattern, full.names = T)

n1 <- readRDS(rds[1])$wTO
n2 <- readRDS(rds[2])$wTO
n3 <- readRDS(rds[3])$wTO


rbind(mutate(n1, group = 'DES'),
      mutate(n2, group = 'POS'),
      mutate(n3, group = 'PRE')) -> datViz


n1 %>% wTOcutoff(cutoff = 0.01) -> n1Cutoff
n2 %>% wTOcutoff(cutoff = 0.01) -> n2Cutoff
n3 %>% wTOcutoff(cutoff = 0.01) -> n3Cutoff

g <- plotNet(n1Cutoff, tau = 0.5)

aesthNet <- function(wTO, tau = 0.5) {
  
  library(ggraph)
  library(ggforce)
  
  
  g <- plotNet(wTO, tau = tau)
  
  # g$Nodes %>% left_join(kegg, by = c('id' = 'gene')) -> Nodes # pseudo-code
  
  g$Nodes -> Nodes
  
  Nodes <- Nodes %>% mutate(size = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(size = size * 2 + 1)
  
  names(wTO)[3] <- 'wTO'
  Edges <- wTO %>% filter(abs(wTO) > tau) %>% 
    mutate(wTOc = ifelse(wTO > 0, '+', '-'), 
           width = 0.5 + 5 * abs((wTO - min(wTO))/(max(wTO) -  min(wTO))))
  
  graph = graph_from_data_frame(Edges, directed = FALSE, Nodes)
  
  create_layout(graph, layout = 'igraph', algorithm = 'kk') %>%
    ggraph(.) + 
    geom_edge_link(aes(edge_colour = wTOc, edge_alpha = width), width = 1.2) +
    geom_node_point(size = 4) + 
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
    theme(legend.position="top") -> gplot
  
  return(gplot)
  
  # return(graph)

}

aesthNet(n1Cutoff, tau = 0.5)
aesthNet(n2Cutoff, tau = 0.5)
aesthNet(n3Cutoff, tau = 0.5)


# diffNet ----
library(CoDiNA)


Graph = CoDiNA::plot.CoDiNA(DiffNet, layout = 'layout_components', Cluster = TRUE)


library(ggraph)

Edges <- Graph$Edges
Nodes <- Graph$Nodes

Edges %>% mutate_at('Group', funs(str_replace_all(., c("^[a-z][.]"="")))) -> Edges
# Nodes %>% left_join(kegg, by = c('id' = 'unique')) -> Nodes

graph <- graph_from_data_frame(Edges, directed = FALSE, Nodes)

group <- igraph::cluster_louvain(graph)$membership

nodes = plyr::join(Nodes, data.frame(id = igraph::V(graph)$name, 
                                     group = group))

graph = graph_from_data_frame(Edges, directed = FALSE, nodes)


# V(graph)$facet_node <- facet_nodes[match(V(graph)$name, facet_nodes$id), 'Group']

layout = create_layout(graph, layout = 'igraph', algorithm = 'kk')

ggraph(layout) + 
  geom_edge_link(aes(edge_colour = Group, edge_alpha = Score)) +
  # facet_edges(~Phi, scales = "free") +
  # facet_graph(Phi~Group)
  geom_node_point(aes(size = Degree_Total)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() -> graphSave

graphSave + geom_mark_hull(
  aes(x, y, group = group, label=group),
  fill = "grey", color = NA,
  concavity = 4,
  con.size = 0.3,
  con.linetype = 2,
  expand = unit(2, "mm"),
  alpha = 0.25)
