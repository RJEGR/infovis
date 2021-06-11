
# Biological systems have to regulate when, under which conditions, and how much of a particular protein or RNA is expressed. Interestingly, the molecules that perform this regulation (gene regulatory factors) by binding to DNA are themselves encoded in the DNA, thus creating gene regulatory networks (figure 2c). The effect of the binding can be deduced from analysing expression changes of target genes using RNA-Seq

library(tidygraph)
library(tidyverse)
library(igraph)
library(ggraph)

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
  cat("\nProportion of : ", length(unique(c(out$Node.1, out$Node.2)))/ n, "\n")
  
  return(out)
  
}

path <- "~/transcriptomics/oktopus_full_assembly/wTO/spearman_by_stages/"
pattern <- "DES"
rds <- list.files(path = path, pattern = pattern, full.names = T)


exportNet <- function(fileName, cutoff = 0.01, tau = 0.5) {
  
  group <- paste(sapply(strsplit(basename(rds[1]), "_"), `[`, 1:3),collapse = '_')
  
  WTO <- readRDS(fileName)$wTO
  WTO %>% wTOcutoff(cutoff, tau) -> WTO
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  wTO <- WTO$wTO_sign
  
  input_vis = data.frame(Node.1 = Node.1, Node.2 = Node.2, 
                         wTO = as.numeric(wTO))
  input_vis = input_vis[!is.na(input_vis$wTO), ]
  input_vis = plyr::arrange(input_vis, input_vis$Node.1, input_vis$Node.2)
  
  nodes <- data.frame(id = sort(unique(c(as.character(input_vis$Node.1), 
                                         as.character(input_vis$Node.2)))), group)
  
  graph <- tbl_graph(nodes = nodes, edges = input_vis, directed = FALSE)
  
  graph %>% activate(nodes) %>%
    mutate(degree = centrality_degree()) %>%
    mutate(value = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(weight = value * 2 + 1) -> graph
  
  gr <- function(x) {ifelse(x > 0, '+', '-')}
  
  graph %>% activate("edges") %>%  mutate(color = gr(wTO)) -> graph

  
  return(graph)
  
}

exportNet(rds[1]) -> graph

graph %>% activate("nodes") %>% mutate(betweenness = betweenness(.)) -> graph

graph %>% arrange(desc(degree))
# membership = igraph::cluster_louvain(.)$membership
quantile(V(graph)$betweenness)
quantile(V(graph)$betweenness, probs = 0.75) -> cutOff

graph %>% activate(nodes) %>% filter(betweenness >= cutOff) -> g


quantile(E(g)$wTO, probs = 0.75) -> EdgeCut

g %>% activate("edges") %>%
  mutate(wTO = ifelse(wTO > EdgeCut, wTO, 0)) %>%
  filter(wTO != 0) -> g

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

# ggraph(layout) +
#   ggforce::geom_mark_hull(
#     aes(x, y, group = as.factor(membership)), #  fill = as.factor(membership)
#     color = NA,
#     concavity = 4,
#     con.size = 0.3,
#     con.linetype = 2,
#     expand = unit(2, "mm"),
#     alpha = 0.25)  -> psaveNet

ggraph(layout) + # psaveNet
  geom_edge_link(aes(edge_alpha = abs(wTO), color = color), edge_width = 1)+
  geom_node_point(aes(size = weight)) + # , alpha = GS, color = abs(GS), 
  # geom_node_text(aes(label = id), repel = TRUE, size = 2) +
  scale_edge_width(range = c(0.2, 2)) +
  # ggsci::scale_color_gsea(name = 'Gene Correlation', reverse = T) +
  # scale_color_viridis_c(name = 'weight', direction = -1) +
  # scale_color_continuous(name = 'Gene Correlation') + # values = nodeCol
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top")+
  coord_fixed() 
  # facet_nodes(~ group)
