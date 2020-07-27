# 
# Explore	interesing	modules	
# - hub	genes	
# - funcional	characterization	

library(WGCNA)
library(flashClust)
library(tidyverse)

path <- '~/transcriptomics/Diana_Lara/Diana_results/'

fileName <- 'Network_signed_nomerge_RLDfiltered.RData'
datExpr_file <- 'datExpr.rds'
  
datExpr <- readRDS(file = paste0(path, '/', datExpr_file))

load(paste0(path,'/', fileName))

# how modules where obtained:
nm <- table(moduleColors)
cat('Number of mudules obtained\n :', length(nm))
print(nm)

plotDendroAndColors(geneTree, moduleColors, c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

# Load neuropeptide

wd <-'~/transcriptomics/Diana_Lara/neuropeptides'
# options(stringsAsFactors = F)

neuropep <- dir(path = wd, pattern = 'neuropeptide.xls', full.names = T)

annot <- dir(path = wd, pattern = 'Trinotate.xls', full.names = T)

require('trinotateR')

y <- read_trinotate(annot)

pfam <- split_pfam(y)
go <- split_GO(y)
blastx <- split_blast(y, "sprot_Top_BLASTX_hit")

goBP <- go %>% 
  as_tibble() %>%
  filter(ontology == 'biological_process') %>%
  distinct(transcript, go, .keep_all = T) %>%
  select(transcript, go, name)
  

# Merge neuropeptide List

signifGenes <- unique(y$transcript_id)

which_annot <- which(colnames(datExpr) %in% signifGenes)

moduleColors[which_annot]
colnames(datExpr)[moduleColors == "brown"]

colnames(datExpr)[which_neuropep]

datExpr_net <- t(datExpr) %>%
  as_tibble(rownames = 'transcript') %>%
  pivot_longer(-transcript) %>%
  mutate(sample = paste0("LOF_", sapply(strsplit(name, "_"), `[`, 2))) %>%
  group_by(transcript, sample) %>%
  summarise(mean = mean(value)) %>%
  pivot_wider(names_from = sample, values_from = mean,  
              values_fill = list(mean = 0)) %>%
  ungroup() %>%
  mutate_if(is.numeric, round) %>%
  mutate(moduleColors = moduleColors)

sam <- names(datExpr_net)[grepl('LOF', names(datExpr_net))]

datExpr_net %>%
  pivot_longer(sam, values_to = 'count', names_to = 'sam') %>%
  group_by(moduleColors, sam) %>%
  summarise(count = sum(count), n = n()) %>%
  ggplot(aes(x = sam, y = moduleColors))

datExpr_net[sam] %>%
  superheat::superheat(., 
                     scale = T,
                     membership.rows = moduleColors,
                     membership.cols = 1:length(sam),
                     left.label.text.size = 2.0,
                     left.label.col = 'white',
                     bottom.label = 'variable',
                     left.label = 'cluster',
                     grid.hline = T,
                     grid.vline = T,
                     smooth.heat = T,
                     grid.vline.col = 'white',
                     grid.hline.col = 'white',
                     bottom.label.text.size = 3.0)

datExpr_net %>%
  inner_join(blastx) %>%
  group_by(transcript) %>%
  select(name, moduleColors)

# net ----

library(igraph)

adjacency <- readRDS(paste0(path,'/', 'adjacency.rds'))
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify

adj <- TOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
g <- graph.adjacency(adj)
g <- simplify(g)

V(g)$color <- moduleColors

g <- delete.vertices(g, degree(g)==0)

g <- as.undirected(g)

library(tidygraph)

g <- as_tbl_graph(g, directed = FALSE)

g %>% 
  activate(nodes) %>% 
  mutate(degree = centrality_degree()) -> g
  # activate(edges) %>% 
  # mutate(betweenness = centrality_edge_betweenness()) %>% 
  # arrange(betweenness) -> g

g %>% 
  activate(nodes) %>% 
  left_join(age)
