# Ricardo Gomez-Reyes, 05/07/2022
# Based on https://github.com/scwatts/fastspar toyset
# Prepare scripts for merge pvalues, median correlation and features into a S3 object for network

# Most networks show a power-law degree distribution, where a few nodes have a very large number of connections, while other nodes have no or few connections [53]. The highly connected nodes are called hubs, and networks following a power-law degree distribution are often called scale-free networks (Figure ID). Cellular networks, genetic regulatory networks, and protein–protein interaction networks are biological examples of scale-free networks 
# https://icon.colorado.edu/#!/

# https://www.nature.com/articles/s41467-019-08746-5

# A. Clean canvas and memory


#SampleID	BioflocType

# Control1	Control
# Control2	Control
# Control3	Control
# OatBran1	Oat bran
# OatBran2	Oat bran
# OatBran3	Oat bran
# WheatBran1	Wheat bran
# WheatBran2	Wheat bran
# WheatBran3	Wheat bran
# Amaranth1	Amaranth
# Amaranth2	Amaranth
# Amaranth3	Amaranth

rm(list = ls());

if(!is.null(dev.list())) dev.off()

# B. Prepare work directory and file paths

path <- '~/Downloads/datos_marce/'

options(readr.show_col_types = FALSE, stringsAsFactors = FALSE)


pvals_f <- list.files(path = path, pattern = 'pvalues.tsv', full.names = T)

cor_f <- list.files(path = path, pattern = 'median_correlation.tsv', full.names = T)

mtd_f <- list.files(path = path, pattern = 'tabla_final_para_R.tsv', full.names = T)

tax_f <- list.files(path = path, pattern = 'taxonomy.tsv', full.names = T)

# C. load key functions and r package for merge data and prepare input for data visualization

prep_net <- function(cor_f, pvals_f) {
  
  cor_data <- read_tsv(cor_f, comment = "#", col_names = F) %>% select(-X1)
  p_data <- read_tsv(pvals_f, comment = "#", col_names = F) %>% select(-X1)
  
  # Correlation to adjacency Conversion file ----
  # 0) replace upper triangle of the matrix with NA for filtering step
  
  to_adj <- function(data, values_to = 'values') {
    
    adj_mat <- data
    
    adj_mat[upper.tri(adj_mat)] <- NA
    
    
    # 1) parse colnames and rownames
    # 2) filter out the upper matrix values
    # 3) filter out the node loops (self correlations)
    # 4) create adjacency list
    
    colnames(adj_mat) <- paste0('sp-',1:nrow(adj_mat))
    
    adj_mat %>%
      mutate(fromNode = paste0('sp-',1:nrow(adj_mat))) %>%
      pivot_longer(-fromNode, names_to = 'toNode', values_to = all_of(values_to)) %>%
      drop_na(values_to) %>%
      filter(fromNode != toNode) -> adj_mat
    
    return(adj_mat)
    
  }
  
  to_adj(p_data, 'pvalue') -> adj_mat
  
  to_adj(cor_data, 'weight') %>%
    left_join(adj_mat, by = c("fromNode", "toNode")) %>%
    mutate_at(vars('weight', 'pvalue'), as.numeric)
  
}

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)

# D. Merge step

Edges <- prep_net(cor_f, pvals_f)

# D.1 Prepare metadata

Node <- read_tsv(cor_f, comment = "#", col_names = F) %>% pull(X1)

Node <- data.frame(Node = paste0('sp-',1:length(Node)), Family = Node)

# D.2 And node Taxonomy

tax <- read_tsv(tax_f, col_names = T)

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  
tax %>% 
  mutate(Taxon = gsub("D_[0-9]__", "", Taxon, perl = T)) %>%
  separate(col = Taxon, into = ranks, sep = ';') %>%
  # rename("Node" = "Family") %>%
  select(-`Feature ID`, -Confidence) %>%
  distinct(Family, .keep_all = T) -> tax

# D.3 Merge taxonomy labels to node

Node %>% left_join(tax) %>% as_tibble() -> Node

# (Omit) Warning messages: Expected 6 pieces. Additional pieces ...

# E. Separate correlation by group

hist(Edges$weight)
hist(Edges$pvalue)

Edges %>% 
  mutate(type = NA) %>% 
  mutate(type = ifelse(weight > 0, '+', type)) %>%
  mutate(type = ifelse(weight < 0, '-', type)) %>%
  mutate(weight = abs(weight)) -> Edges


# E.1 Edge viz

Edges %>% 
  mutate(weight = ifelse(type == "-", -weight, weight)) %>%
  mutate(pvalue = as.numeric(pvalue)) %>%
  mutate(g = ifelse(pvalue < 0.05, 'sig', 'no.sig')) %>%
  group_by(g, type) %>%
  tally()

Edges %>% 
  mutate(weight = ifelse(type == "-", -weight, weight)) %>%
  mutate(pvalue = as.numeric(pvalue)) %>%
  ggplot(aes(weight, -log10(pvalue), color = pvalue)) +
  ggsci::scale_color_gsea() +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red')


# F. add feature info into nodes (for example using freq and ab of feature count)

feature_count <- read_tsv(mtd_f, comment = "#", col_names = T) %>% 
  as.data.frame()

rownames(feature_count) <- feature_count$OTU_ID

feature_count$OTU_ID <- NULL

prevelancedf = apply(X = feature_count,
  MARGIN = 1,
  FUN = function(x){sum(x > 0)})

df = data.frame(Prevalence = prevelancedf, 
  Prev_pct = prevelancedf/max(prevelancedf),
  TotalAbundance = rowSums(feature_count)) %>% 
  as_tibble(rownames = "Node") %>%
  mutate(Node = paste0('sp-',1:nrow(.))) 

Nodes <- unique(c(Edges$fromNode, Edges$toNode))

Nodes <- df %>% filter(Node %in% Nodes) %>% 
  left_join(Node) 

# G. Create network object

# net <- graph.data.frame(adj_mat, directed = FALSE) # FOR igraph

g <- tidygraph::tbl_graph(edges = Edges, nodes = Nodes,  directed = F) # instead of directed = T

# plot(g, edge.curved=T)

# H. Filter edges by weight or pvals

# g %>% activate('edges') %>% filter(abs(weight) > 0.1) -> g
g %>% activate('edges') %>% filter(pvalue < 0.05) -> g

# plot(g, edge.curved=T)

# Calculate additional node features based on graph theory

degree <- degree(g)

# in most real networks, the degree distribution is highly asymmetric (or skewed)
# http://users.dimi.uniud.it/~massimo.franceschet/ns/plugandplay/scale-free/scale-free.html#4

skewness = function(x) mean( (x - mean(x))^3 ) / sd(x)^3

skewness(degree)


# I. Using ggplot grammar for visualization edition

w_filter <- quantile(Edges$weight, probs = c(0.5))

g %>% activate("edges") %>%  mutate(weight = ifelse(weight > w_filter, weight, NA)) -> g
g %>% activate("edges") %>% filter(!is.na(weight)) -> g


# La centralidad de grado («degree centrality»)
# La intermediación («betweenness»)
# PageRank algorithm to the microbiome networks to identify key members 

g %>% activate("nodes") %>%  
  mutate(
    betweenness = betweenness(.), 
    degree = centrality_degree(),
    pageRank = page_rank(.)$vector,
    membership = igraph::cluster_louvain(.)$membership # if directed = F
  ) -> g

g %>% activate("nodes") %>% filter(degree > 0) -> g


layout = ggraph::create_layout(g, layout = 'igraph', algorithm = 'kk')

# ggraph(layout) +
  # geom_edge_fan(aes( color = type), curvature = 0.2)

ggraph(layout) +
  geom_edge_arc(aes(edge_alpha = weight, color = type), strength = 0.1,) + # edge_width
  geom_node_point(aes(color = Prev_pct, size = degree)) + 
  geom_node_text(aes(label = Family), repel = TRUE, size = 2) +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top")+
  coord_fixed() -> psave

psave + facet_edges(~ type) -> psave2


# J. Save as png

ggsave(psave, filename = "network.png", path = path, width = 16, height = 7, dpi = 300)

ggsave(psave2, filename = "network_by_edges.png", path = path, width = 16, height = 7, dpi = 300)


# K. test layouts for better viz
# 
# ggraph(g, layout = 'linear', circular = TRUE, sort.by = type) + 
#   geom_edge_arc(aes(edge_alpha = weight, color = type)) +
#   geom_node_point(aes(color = Prev_pct, size = degree))

# c.layout <- ggraph::create_layout(g, layout = "linear", 
#   circular=TRUE, sort.by = Prev_pct)
# 
# ggraph(c.layout) + 
#   geom_node_point(aes(size = degree)) +
#   geom_edge_arc(aes(edge_alpha = weight, color = type,
#     circular=TRUE)) +
#   facet_edges(~ type)

# L. Analyze node interactions
g %>% activate("nodes") %>% as_tibble() -> Nodes

g %>% activate("edges") %>% 
  filter(type %in% '+') %>% 
  activate("nodes") %>%
  mutate(
  degree = centrality_degree(),
  pageRank = page_rank(.)$vector) %>% 
  as_tibble() -> positive_nodes

g %>% activate("edges") %>% 
  filter(type %in% '-') %>% 
  activate("nodes") %>%
  mutate(
    degree = centrality_degree(),
    pageRank = page_rank(.)$vector) %>% 
  as_tibble() -> negativa_nodes

write.csv(x = positive_nodes, file = paste0(path,'positive_interactions.csv'))
write.csv(x = negativa_nodes, file = paste0(path,'negative_interactions.csv'))

# ggsave(psave, filename = "fastspar_nodes.png", path = path, width = 5, height = 5)

# Calculate nodes features (omit) ----

# comp <- components(g, mode = 'weak')
# table(comp$membership)
# group <- igraph::cluster_louvain(g)$membership


g %>% activate("nodes") %>%  
  mutate(betweenness = betweenness(.), 
    degree = centrality_degree(),
    centrality = components(.)$membership,
    transitivity = transitivity(.),
    membership = igraph::cluster_louvain(g)$membership) -> g

hist(V(g)$degree)

# 
# La centralidad de grado («degree centrality»)
# La cercanía («closeness»)
# La intermediación («betweenness»)

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

ggraph(layout) +
  ggforce::geom_mark_hull(
    aes(x, y, fill = as.factor(membership)),
    color = NA,
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25) -> psaveNet

