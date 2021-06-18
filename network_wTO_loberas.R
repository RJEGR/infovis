#

rm(list = ls())

# 


cols <- ggpubr::get_palette(palette = 'default', 3)
names <- c("P", "NP", "NC")
cols <- ggpubr::get_palette(palette = 'default', 3)

cols = structure(cols, names = names)

dir <- '~/metagenomics/Loberas_MG/'
# ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species', 'pplacer_sp')
sampleName <- c('Cantiles', 'Coloradito', 'Granito', 'Machos', 'Partido', 'Rasito')

fileNames <- list.files(dir, pattern = "xls", full.names = T)

readFeatures <- function(filename) {
  
  ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species', 'pplacer_sp')
  
  readXLS <- function(x) {
    group <- paste(sapply(strsplit(basename(x), "_"), `[`, 1),collapse = '_')
    readxl::read_xlsx(x, na = 'NA') %>% mutate(group = group)
  }
  
  features <- readXLS(filename)
  
  names(features)[9] <- 'pplacer_sp'
  
  # brief step, replace ab matriz w/ rarefacted
  
  abFile <- list.files('~/metagenomics/Loberas_MG/Rarefied_tables/', pattern = "tsv", full.names = T)
  
  features %>% distinct(group) %>% pull() ->  prefixF
  
  read.delim(abFile[grepl(prefixF, abFile)], header = F, comment.char = '#') -> ab
  
  colnames(ab) <- c('Feature_IDs', sampleName)
  
  left_join(features %>% select(!sampleName), ab) -> features
  
  features %>% 
    mutate(Relationship = recode_factor(Relationship, P = 'P', 
                                        NC = 'NC', C = 'NP')) %>%
    mutate_all(., funs(str_replace_all(., 
                                       c("Bacteroides tectu$" = "Bacteroides tectus",
                                            'Escherichia/Shigella' = 'Enterobacteriaceae',
                                         'Neisseriaceae'='Neisseria')))) %>%
    mutate_at(sampleName, as.double) -> features
  
  dat <- features %>% select_at(sampleName) %>% mutate(id = 1:nrow(features))
  
  # curate taxonomy

  features %>%
    mutate(id = 1:nrow(features)) %>%
    select(id, Relationship) -> pplacer
  
  features %>% 
    select_at(all_of(ranks)) %>%
    mutate(id = 1:nrow(features)) %>%
    pivot_longer(cols = ranks) %>% fill(value) %>% 
    pivot_wider(names_from = name) %>%
    left_join(pplacer) %>% 
    mutate(Species = ifelse(Relationship %in% 'P', pplacer_sp, Species)) %>%
    left_join(dat, by = 'id')

}

prepare_ps <- function(filename, agg = T, agg_level = 'Species') {
  
  ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species')
  
  features <- readFeatures(filename)
  
  # features %>% select(pplacer_sp, Relationship) 
  
  dat <- features %>% select_at(sampleName) %>% data.frame(row.names = 1:nrow(features))
  

  
  # hay problemas con la taxonomia que francesco curo, por tanto tener cuidado al usar datos aglomerados por taxonomia
  
  features %>% 
    mutate(Species = paste0(Species, '__' ,Relationship)) %>%
    select_at(all_of(ranks)) %>%
    data.frame(row.names = 1:nrow(features)) -> tax
  
  # view(tax)
  
  #  identical(rownames(dat), rownames(tax))
  
  # and parse
  ps = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                tax_table(as(tax, 'matrix'))) 
  
  if(agg) {
    microbiome::aggregate_taxa(ps, level = agg_level)
    # tax_glom(ps, taxrank = 'Species')
  } else
    return(ps)
  
  
}

from_pyseq_to_wTO <- function(phyloseq, overlap = NULL, n = 250) {
  
  phyloseq %>% prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
    otu_table() %>%
    as.data.frame() -> metagenomics_abundance 
  
  if (is.null(overlap)) {
    wTO <- wTO.fast(Data = metagenomics_abundance,
                    Overlap = row.names(metagenomics_abundance),
                    method = 's', sign = 'sign', n = n, 
                    method_resampling = 'Bootstrap')
    
  } else
    
    wTO <- wTO.fast(Data = metagenomics_abundance,
                    Overlap = overlap,
                    method = 's', sign = 'sign', n = n, 
                    method_resampling = 'Bootstrap')
  # method_resampling = 'BlockBootstrap', lag = 2
  
  return(wTO)
  
}

# dataList <- readRDS(paste0(dir, '/phyloseqList.rds'))

library(wTO)
library(tidyverse)
library(phyloseq)
library(tidygraph)


library(ggraph)
library(igraph)
library(ggforce)
library(scatterpie)
library(RColorBrewer)

ps <- prepare_ps(fileNames[1], agg = T)

# ps %>% tax_table() %>% as.data.frame() %>%
#   mutate_if(is.character, as.factor) %>% 
#   filter(grepl('__P', unique)) %>% pull(Species) -> overlaps

wTOout <- ps  %>% from_pyseq_to_wTO(., n = 20)

head(wTOout)

# a medida que aumenta el valor P, disminuye la relacion, 
# usaremos una medida de filtrado basado en pval-abs(wTO) < umbral (0.01/0.05)
# wTOout %>% ggplot(aes(abs(wTO), pval)) + geom_point()

wTOout %>%
  separate(col = 'Node.1', into = c('Node.1', 'r1'), sep = '__') %>%
  separate(col = 'Node.2', into = c('Node.2', 'r2'), sep = '__') %>%
  mutate(facet = paste0(r1,'-' ,r2)) %>%
  filter(grepl('^P-|-P$', facet)) -> wTOdv 


wTOdv %>%
  ggplot() + # aes(color = group)
  geom_point(aes(wTO, pval-abs(wTO)), alpha = 0.5) +
  # ggsci::scale_color_gsea(reverse = T) + 
  facet_grid(r1~r2) +
  labs(x = 'Correlation',   y = expression(~Log[10]~('Pvalue')),caption = '# Interactions (a)') -> p1

dat_text <- wTOdv %>% group_by(r1, r2) %>% tally() 
dat_text <- cbind(dat_text, pval = 1, wTO = 1)

p1 + geom_text(
  data    = dat_text, family = "GillSans", color = 'black',
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, "a")),
  hjust   = -0.5,
  vjust   = -2
) + theme_bw(base_size = 16, base_family = "GillSans")


# Most biological networks show a power-law degree distribution, where a few nodes have a very large number of connections, while other nodes have no or few connections

summary(wTOout$wTO)
summary(wTOout$pval)

# version para analizar una sola red -----
library(tidygraph)

wTOcutoff <- function(wTO_out, cutoff = 0.01, tau = 0.5) {
  
  n <- length(unique(c(wTO_out$Node.1, wTO_out$Node.2)))
  
  wTO_out %>% mutate_if(is.factor, as.character) %>%
    mutate(wTO = ifelse(pval-abs(wTO) < cutoff, wTO, 0 )) %>%
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
  
  nodes %>% 
    mutate_if(is.character, as.factor) -> nodes
  
  ps %>% tax_table() %>% as.data.frame() -> infoNodes
  
  ps %>%
    transform_sample_counts(., function(x) x / sum(x)) %>%
    otu_table() %>% 
    as.data.frame() %>% cbind(unique = rownames(.)) %>%
    left_join(infoNodes, by = 'unique') %>% 
    mutate_if(is.character, as.factor) %>% as_tibble() -> infoNodes
  
  # sum(infoNodes$unique %in% nodes$id)
  
  nodes %>%
    as_tibble() %>%
    left_join(infoNodes, by = c('id' = 'Species')) %>%
    separate(col = 'id', 
             into = c('Species', 'Relationship'), sep = '__') -> nodes
  
  
  WTO %>%
    mutate_if(is.character, as.factor) %>%
    separate(col = 'Node.1', into = c('Node.1', 'r1'), sep = '__') %>%
    separate(col = 'Node.2', into = c('Node.2', 'r2'), sep = '__') %>%
    mutate(typeEdge = paste0(r1, '-', r2)) %>%
    select(-r1, -r2) -> edges
    
  edges %>%
    filter(grepl('^P-|-P$', typeEdge)) -> edges
  
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = F)
  
  graph %>% activate(nodes) %>%
    mutate(degree = centrality_degree()) %>%
    mutate(value = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(weight = value * 2 + 1) -> graph
  
  gr <- function(x) {ifelse(x > 0, '+', '-')}
  
  graph %>% activate("edges") %>%  mutate(color = gr(wTO)) -> graph
  
  
  return(graph)
  
}

# en vez de usar cutoff, usar solo tao y basarse en el cuartil 75 para recortar los edges

# wTOout %>% wTOcutoff(cutoff, tau) -> WTO

exportNet(wTOout, cutoff = 1, tau = 0) -> graph

graph %>%  
  activate("edges") %>% 
  mutate(typeEdge = case_when(
  grepl('^P',typeEdge) ~ 'P',
  grepl('^NC',typeEdge) ~ 'NC',
  TRUE ~ 'NP')) -> graph

components(graph)


graph %>% 
  activate("nodes") %>% 
  mutate(betweenness = betweenness(.),
         pageRank = page_rank(.)$vector) %>%
  as_data_frame(., "vertices") %>% 
  select(Species, Relationship, degree, betweenness, pageRank ) %>%
  pivot_longer(cols = c('degree', 'betweenness')) %>%
  ggplot(aes(value, pageRank)) + 
  geom_point() +
  facet_grid(Relationship~name, scales = 'free')


# 1

# scales::show_col(cols)

graph %>%
  activate("edges") %>% # pval-abs(wTO) < cutoff, wTO, 0
  mutate(wTO = ifelse(pval-abs(wTO) < 0.05, wTO, 0)) %>% filter(wTO > 0) -> g

g %>% activate("nodes") %>% 
  mutate(degree = degree(.)) %>% filter(degree > 0) -> g

# Hub detection as Layeghifard, M.,et al 2019
# top-ranked network hub taxa:
# We applied the PageRank algorithm to the microbiome networks to identify key members of CF lung microbiome in each patient. PageRank is a link analysis algorithm with the underlying assumption that hubs are likely to be more connected to other nodes when compared to non-hub nodes. This algorithm assigns a numerical weight to each node of a network as a measure of its relative importance within the network. The numerical weights of the nodes (also called PageRanks of the nodes) were then sorted and the top five or ten taxa with highest PageRank were selected as microbiome network hubs (or key taxa). For the network hub taxa, we recorded both the top five/top ten taxa for further analysis. 

g %>% 
  activate("nodes") %>% 
  mutate(betweenness = betweenness(.),
         membership = igraph::cluster_louvain(.)$membership,
         pageRank = page_rank(.)$vector) -> g



g %>% activate("edges") %>%
  mutate(betweenness = edge.betweenness(.)) -> g


# hive works better when is directed network
g %>%
  ggraph(., 'hive', axis = Relationship, sort.by = degree) + 
  geom_edge_hive(aes(edge_alpha = abs(wTO), color = color)) +
          #        arrow = arrow(
          # angle = 10,
          # length = unit(0.1, "inches"),
          # ends = "last",
          # type = "closed")) + # factor(color)
  geom_axis_hive(aes(colour = Relationship), size = 2, label = FALSE) + 
  coord_fixed() +
  facet_edges(~color) +
  theme_graph(base_family = "GillSans") +
  theme(legend.position = "top") +
  scale_edge_color_manual(values = cols) +
  scale_color_manual(values = cols)

# 2


layout = create_layout(g, layout = 'igraph', algorithm = 'kk')


ggraph(layout) +
  ggforce::geom_mark_hull(
    aes(x, y, group = as.factor(membership), 
        fill = as.factor(membership)), #  fill = as.factor(membership)
    color = NA,
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25)  +
  guides(fill = FALSE) -> psaveNet

psaveNet +
  geom_edge_link(aes(color = typeEdge, alpha = abs(wTO)), # edge_width = betweenness
                 arrow = arrow(
                   angle = 10,
                   length = unit(0.1, "inches"),
                   ends = "last",
                   type = "closed")) +
  geom_node_point(aes(size = pageRank, color = Relationship)) +
  geom_node_text(aes(label = Species), repel = T) +
  coord_fixed() +
  # scale_size('Hub') +
  facet_edges(~color) +
  theme_graph(base_family = "GillSans") +
  theme(legend.position = "top") +
  scale_edge_colour_manual(values = cols, guide = "none") +
  scale_color_manual(values = cols, name = '') +
  facet_edges(~color)

# 3 subgroup taxa per memberships

# que % de la ab. relativa de los nodos corresponden a interacciones patogenicas?? :

as_data_frame(g, "vertices") %>% as_tibble() %>%
  select(all_of(sampleName)) %>% colSums()

# Ej. para v2, al parecer cantiles, granito y rosarito resultan con mayor interaccioenes de Patog

as_data_frame(g, "edges") %>% as_tibble() %>%
  ggplot(aes(wTO, betweenness, color = typeEdge)) +
  geom_point()

as_data_frame(g, "vertices") %>% as_tibble() %>%
  # select(Species, all_of(sampleName),Relationship, pageRank, membership, degree) %>%
  arrange(desc(pageRank)) %>%
  # pivot_longer(cols = sampleName) %>%
  ggplot(aes(y = pageRank, x = Relationship)) +
  geom_boxplot() + labs(caption = 'Hub taxa detection')
  # facet_grid(~)

as_data_frame(g, "vertices") %>% as_tibble() %>%
  filter(pageRank > quantile(V(g)$pageRank,  probs = 0.75)) %>%
  ggplot(aes(y = pageRank, x = degree, size = betweenness)) +
  geom_point() + ggrepel::geom_text_repel(aes(label = Species), size = 4) +
  labs(y = 'Hubs taxa detected')

# loop the network into lists to diffNet ----
# recuperar solo las interacciones P
# Gracias Padre
WTO <- lapply(fileNames[1:5], function(x) {
  
  g <- strsplit(basename(x), "_")
  g <- paste(sapply(g, `[`, 1),collapse = '_')
  
  prepare_ps(x, agg = T) %>%
    from_pyseq_to_wTO(., n = 250) -> WTO
  
  WTO %>% 
    separate(col = 'Node.1', into = c('Node.1', 'r1'), sep = '__') %>%
    separate(col = 'Node.2', into = c('Node.2', 'r2'), sep = '__') %>%
    mutate(pair = paste0(r1,'-' ,r2)) %>%
    select(-r1, -r2) %>%
    # filter(grepl('^P-|-P$', pair)) %>%
    cbind(., group = g) -> WTO
  
  return(WTO)
})

# revisar todas las redes como independientes y comparar vs las diffnets
# saveRDS(WTO, file = paste0(dir, 'MakeDiffNet_in.rds'))
# filter(grepl('^P-|-P$', pair)) para quedarse solo con las interacciones P-NC y P-NP

WTO <- readRDS(paste0(dir, 'MakeDiffNet_in.rds'))

exportNet_ <- function(WTO, cutoff = 0.01) {
  
  WTO %>%
    mutate(wTO = ifelse(pval-abs(wTO) < cutoff, wTO, 0)) %>%
    filter(abs(wTO) > 0) -> WTO
  
  WTO %>% separate(col = 'pair', into = c('r1','r2'), sep = '-') %>%
    mutate(Node.1 = paste0(Node.1, '__',r1),
           Node.2 = paste0(Node.2, '__',r2)) -> WTO
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  Node.1 = paste0(Node.1, '__', as.character(WTO$group))
  Node.2 = paste0(Node.2, '__', as.character(WTO$group))
  
  nodes <- data.frame(id = sort(unique(c(Node.1,Node.2))))
  
  nodes %>% mutate_if(is.character, as.factor) -> nodes
  
  rownames(nodes) <- as.character(nodes$id)
  
  into <- c('Species', 'Relationship', 'group')
  
  nodes %>%
    as_tibble() %>%
    separate(col = 'id', 
             into = into, sep = '__') -> nodes
  
  WTO %>%
    mutate_if(is.character, as.factor) %>%
    separate(col = 'Node.1', into = c('Node.1', 'r1'), sep = '__') %>%
    separate(col = 'Node.2', into = c('Node.2', 'r2'), sep = '__') %>%
    mutate(typeEdge = paste0(r1, '-', r2)) %>%
    select(-r1, -r2) -> edges
  
  # edges %>% filter(grepl('^P-|-P$', typeEdge)) -> edges
  
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = F)
  
  graph %>%
    activate("nodes") %>%
    mutate(degree = centrality_degree(),
           betweenness = betweenness(.),
           membership = igraph::cluster_louvain(.)$membership,
           pageRank = page_rank(.)$vector) %>%
    mutate(weight = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(weight = weight * 2 + 1)  -> graph
  
  gr <- function(x) {ifelse(x > 0, '+', '-')}
  
  graph %>% activate("edges") %>%  mutate(color = gr(wTO)) -> graph
  
  return(graph)
}

# exportNet_(WTO[[1]]) 
#   as_data_frame(., what = 'both') %>% str()

out <- list()

V <- list()
E <- list()

for(i in 1:length(WTO)) {
  j <- i
  g <- exportNet_(WTO[[j]])
  out[[j]] <- g
  
  V[[j]] <- as_data_frame(out[[j]], "vertices") 
  E[[j]] <- as_data_frame(out[[j]], "edges") 
  

}


# no es posible unir, pues los edges estan organizados de modo independiente
do.call(rbind, V) -> nodes
do.call(rbind, WTO) -> WTOout

# write files!

nodes %>% as_tibble() %>%
  group_by(group) %>%
  mutate(quantile75 = ifelse(pageRank > quantile(nodes$pageRank,  probs = 0.75), TRUE, FALSE)) -> nodedf


WTOout %>%
  group_by(group) %>%
  mutate(wTO = ifelse(pval-abs(wTO) < 0.01, wTO, 0)) %>%
  filter(abs(wTO) > 0) -> WTOout
  # filter(abs(wTO) > quantile(abs(wTO),  probs = 0.75)) -> WTOout
  # mutate(quantile75 = ifelse(abs(wTO) > quantile(abs(wTO),  probs = 0.75), TRUE, FALSE)) -> WTOout



length(unique(c(WTOout$Node.1, WTOout$Node.2)))
length(nodedf %>% pull(Species) %>% unique())

write_excel_csv(nodedf, file = paste0(dir, "nodesInfo.csv"))
write_excel_csv(WTOout, file = paste0(dir, "edgesInfo.csv"))

# do.call(rbind, E) -> edges

# g <- tbl_graph(nodes = nodes, edges = edges, directed = F)

nodes %>% as_tibble() %>%
  arrange(desc(pageRank)) %>%
  ggplot(aes(y = pageRank, x = group, color = Relationship)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), 
               outlier.shape=NA) +
  scale_color_manual("",values = cols) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank()) +
  labs(caption = 'Key members detection (Layeghifard, M.,et al 2019)', x = '') -> p1

ggsave(p1, filename = "pageRankBox.png", path = dir, width = 4, height = 4)

nodes %>% as_tibble() %>%
  group_by(group, Relationship) %>%
  tally(sort = T) %>%
  mutate(pct = n/sum(n), Relationship = as.factor(Relationship)) %>%
  ggplot(aes(x = group, y = pct, fill = Relationship)) +
  geom_col() +
  scale_fill_manual("",values = cols) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank()) +
  labs(caption = 'Low-correlated nodes filtered. Relationships\npresented in the nodes-nodes interactions', 
       y = '%', x = '') -> p2

ggsave(p2, filename = "ProportionNodes.png", path = dir, width = 4, height = 4)

nodes %>% as_tibble() %>%
  group_by(group) %>%
  filter(pageRank > quantile(nodes$pageRank,  probs = 0.75)) -> hub_nodes 

caption = 'Low-correlated nodes filtered. Relationships\npresented in the nodes-nodes interactions'

nodes %>% as_tibble() %>%
  group_by(group, Relationship) %>%
  tally(sort = T) %>% 
  ggplot(aes(x = group, y = n, fill = Relationship)) +
  geom_col(position = position_dodge2()) +
  scale_fill_manual("",values = cols) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank()) +
  labs(caption = caption,  y = 'Nodes', x = '') -> p3

ggsave(p3, filename = "NumberOfNodes.png", path = dir, width = 4, height = 4)


# upset -----
hub_nodes %>% 
  select(Species, Relationship, group) %>%
  mutate(val = 1) %>% 
  pivot_wider(names_from = group, values_from = val, values_fill = 0) -> ocurrance

ocurrance %>% count(Relationship)

cols <- names(ocurrance)[-c(1,2)]

euler <- as.data.frame(ocurrance[, -1])

rownames(euler) <- ocurrance$Species

library(UpSetR)

# no funciona tanto un diagrama de venn para esta serie de datos (nodos), ya que la resolucion taxonomica esta ajustada en cada marcador. Por tanto, el id de nodo puede variar entre marcador. Asi que nos saltamos 

upset(euler, number.angles = 0, 
      point.size = 3.5, line.size = 0, sets = cols, 
      keep.order = T, 
      mainbar.y.label = "Intersections Size", 
      sets.x.label = "Set Size", 
      text.scale = c(1.3, 1.5, 1, 1.5, 2, 2),
      order.by = "degree")

library(VennDiagram)


# alluvial ? ----
# nada interesante aqui
library(ggalluvial)

alluv_long <- to_lodes_form(ocurrance, key = "x", axes = 3:7)

ggplot(data = alluv_long,
       aes(x = x, stratum = stratum, alluvium = alluvium,
           color = Relationship)) +
  # stat_stratum(geom = "errorbar", na.rm = TRUE, alpha = 0.5) +
  geom_line(stat = "alluvium") +
  stat_alluvium(geom = "pointrange", cement.alluvia = T) +
  # scale_color_manual(values=scale) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            color = 'black') +
  # scale_x_discrete(limits = c("complete", "incomplete")) +
  theme_minimal()

# graficar redes 17/06/21 -----
# si quieres ver la red compelta, usa el objeto WTO[[i]], si deseas la red filtrada usar out[[i]]
# Gracias Papa Dios

cols <- ggpubr::get_palette(palette = 'default', 3)
names <- c("P", "NP", "NC")
cols <- ggpubr::get_palette(palette = 'default', 3)

cols = structure(cols, names = names)

plotNet <- function(WTO, cutoff = 0.01) {
  g <- exportNet_(WTO, cutoff) %>%
    activate(nodes) %>% filter(degree > 0)
  
  # layout = create_layout(g, layout = 'igraph', algorithm = 'kk')
  layout = create_layout(g, layout = 'stress')
  
  ggraph(layout) +
    ggforce::geom_mark_hull(
      aes(x, y, group = as.factor(membership)), #  fill = as.factor(membership)
      color = NA, fill = "grey76",
      concavity = 4,
      con.size = 0.3,
      con.linetype = 2,
      expand = unit(2, "mm"),
      alpha = 0.25)  +
    guides(fill = FALSE) -> psaveNet
  
  psaveNet +
    # facet_nodes(~group) +
    geom_edge_link() + # edge_colour = color, edge_alpha = abs(wTO), 
    geom_node_point(aes(size = degree, color = Relationship)) +
    coord_fixed() +
    scale_size('Degree', range = c(0, 5)) +
    theme_graph(base_family = "GillSans") +
    theme(legend.position = "top") +
    scale_edge_colour_manual(values = c('black', 'blue'), guide = "none") +
    scale_color_manual(values = cols, name = '') +
    facet_edges(~as.factor(color))
  
}

# exportNet_(WTO[[1]], cutoff = 0.01) %>% activate(edges)

# plotNet(WTO[[2]], cutoff = 0.01)  + facet_edges(~as.factor(color))

ggsave(plotNet(WTO[[1]]), filename = paste0(dir, "network_V2.png"), width = 8.5)
ggsave(plotNet(WTO[[2]]), filename = paste0(dir, "network_V3.png"), width = 8.5)
ggsave(plotNet(WTO[[3]]), filename = paste0(dir, "network_V4.png"), width = 8.5)
ggsave(plotNet(WTO[[4]]), filename = paste0(dir, "network_V67.png"), width = 8.5)
ggsave(plotNet(WTO[[5]]), filename = paste0(dir, "network_V8.png"), width = 8.5)

# important nodes interactions labeled

WTO[[1]] %>%
  exportNet_(cutoff = 0.01) %>%
  activate(nodes) %>% filter(degree > 0) %>%
  filter(pageRank > quantile(pageRank,  probs = 0.75)) %>%
  activate(edges)  -> g
  # mutate(wTO = ifelse(abs(wTO) > quantile(abs(wTO), probs = 0.75), wTO, NA)) %>%
  # mutate(Species = ifelse(pageRank > quantile(pageRank,  probs = 0.75), Species, '')) -> g

layout = create_layout(g, layout = 'stress')

ggraph(layout) +
  # facet_nodes(~group) +
  geom_edge_link(aes(edge_colour = color)) + # edge_colour = color, edge_alpha = abs(wTO), 
  geom_node_point(aes(size = degree, color = Relationship)) +
  geom_node_text(aes(label = Species), repel = T) +
  coord_fixed() +
  scale_size('Degree', range = c(0, 5)) +
  theme_graph(base_family = "GillSans") +
  theme(legend.position = "top") +
  scale_edge_colour_manual("", values = c('black', 'blue')) + # guide = "none"
  scale_color_manual(values = cols, name = '')

# from net to adj matrix ----
# as_adjacency_matrix(out[[1]])
as_adjacency_matrix(g, sparse = F, type = "upper", attr="wTO") -> adj
colnames(adj) <- g %>% as_data_frame("vertices") %>% pull(Species)
rownames(adj) <- g %>% as_data_frame("vertices") %>% pull(Species)
superheat::superheat(adj, bottom.label.text.angle = 90)
                     # row.dendrogram = T, col.dendrogram = T)

# bipartite net ----

# https://users.dimi.uniud.it/~massimo.franceschet/ns/syllabus/learn/graph/graph.html
types = c(rep(TRUE,7), rep(FALSE,4))
edges = c(8,1, 8,2, 8,3, 9,2, 9,3, 9,4, 9,5, 10,4, 10,6, 11,5, 11,6, 11,7)
bg = make_bipartite_graph(types, edges, directed=FALSE)

bg %>% as_tbl_graph()

lay = layout.bipartite(bg)
plot(bg, layout=lay[,2:1], vertex.color=c("skyblue","green")[V(bg)$type + 1], vertex.size = 20)

(B = as_incidence_matrix(bg))

# edges -----

# plot hives

do.call(rbind, E) -> edges


edges %>% as_tibble() %>%
  ggplot(aes(y = abs(wTO), x = group, color = color)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), 
               outlier.shape=NA) +
  scale_color_manual("",values = c("blue", "black")) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank()) +
  labs(x = '', "wTO") -> p1


edges %>% as_tibble() %>%
  filter(grepl('^P-|-P$', typeEdge)) %>%
  filter(abs(wTO) > quantile(abs(wTO),  probs = 0.75)) %>%
  group_by(typeEdge, color, group) %>%
  tally() %>% 
  # separate(col = typeEdge, into = c("f1", "f2")) %>%
  ggplot(aes(x = typeEdge, y = n, fill = color)) +
  geom_col(position = position_dodge2()) +
  # geom_text(aes(label = group), position = position_jitter()) +
  # scale_fill_manual("",values = cols) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank())
  # facet_grid(typeEdge ~.) 

# ggsave(p2, filename = "interactions.png", path = dir, width = 4, height = 3)

ggsave(p1, filename = "wtoBoxplot.png", path = dir, width = 4, height = 3)

plotHive <- function(g) {
  g %>%
    activate("edges") %>%
    filter(abs(wTO) > quantile(abs(wTO),  probs = 0.75)) %>%
    # mutate(type = case_when(
    #   grepl('^P',typeEdge) ~ 'P',
    #   grepl('^NC',typeEdge) ~ 'NC',
    #   grepl('^NP',typeEdge) ~ 'NP')) %>%
    ggraph(., 'hive', axis = Relationship, sort.by = degree) + 
    geom_edge_hive(aes(edge_alpha = abs(wTO))) + # color = typeEdge
    geom_axis_hive(aes(colour = Relationship), size = 2, label = FALSE) + 
    coord_fixed() +
    facet_edges(~color) +
    theme_graph(base_family = "GillSans") +
    theme(legend.position = "top") +
    scale_color_manual(values = cols, name = '')
}

ggsave(plotHive(out[[1]]), filename = paste0(dir, "hiveplot_V2.png"), width = 8.5)
ggsave(plotHive(out[[2]]), filename = paste0(dir, "hiveplot_V3.png"), width = 8.5)
ggsave(plotHive(out[[3]]), filename = paste0(dir, "hiveplot_V4.png"), width = 8.5)
ggsave(plotHive(out[[4]]), filename = paste0(dir, "hiveplot_V67.png"), width = 8.5)
ggsave(plotHive(out[[5]]), filename = paste0(dir, "hiveplot_V8.png"), width = 8.5)



# para cada marcador ----
# Cuantas veces las patogenas interaccionan con otras? 
# es negativa o positiva la interaccion?

wTOout <- do.call(rbind, WTO)

wTOout %>%
  mutate(int = ifelse(wTO > 0, 1, -1)) %>%
  select(Node.1, Node.2, pair, group, int) %>%
  separate(col = 'pair', into = c('from', 'to'), sep = '-') -> data_long

data_long %>%
  group_by(from, to, group, int) %>%
  tally() %>%
  filter(group %in% 'V3') %>%
  mutate(value = int*n) %>%
  # mutate(from = ifelse(int < 0,  paste0('-', from), from)) %>%
  # mutate(to = ifelse(int < 0,  paste0('-', to), to)) %>%
  ungroup() %>%
  select(from, to, value) -> adjacencyData

adjacencyData
  
library(circlize)

# int <- adjacencyData$int
# group = structure(int, names = adjacencyData$from)

names <- c("P", "NP", "NC")
grid.col <- ggpubr::get_palette(palette = 'default', 3)

grid.col = structure(grid.col, names = names)
# grid.col = structure(rep(grid.col, 2), names = c(names, paste0("-", names)))
                  
circos.clear()

circos.par(start.degree = 0, 
           gap.degree = 1.2, 
           track.margin = c(mm_h(1), 0),
           # track.margin = c(-0.01, 0.01),
           track.height = mm_h(4), # la distancia entre el highlight.sector y el chord
           points.overflow.warning = FALSE)

# Make the circular plot

adjacencyData %>% mutate(value = ifelse(value > 0, 'red', 'black')) -> arr.col
# link.arr.col = arr.col
arr.col %>% distinct()

chordDiagram(adjacencyData,
             # link.visible = adjacencyData$value < 0,
             directional = 1,
             grid.col = rep(grid.col, 2),
             diffHeight = -mm_h(2),
             target.prop.height = mm_h(2),
             annotationTrack = c("grid", "axis"),
             annotationTrackHeight = c(0.03, 0.01),
             direction.type = c("arrows"), # , "diffHeight"
             # link.arr.type = "big.arrow",
             link.decreasing = T,
             link.sort = TRUE, symmetric = T,
             scale = T, preAllocateTracks = list(
               track.height = mm_h(4),
               track.margin = c(mm_h(4), 0)), transparency = 0.4,
             link.arr.col = arr.col)

# , "-NC"
highlight.sector(c("NC"), track.index = 1, col = "#619CFFFF", 
                 text = "NC", cex = 0.8, text.col = "white", niceFacing = TRUE)
# c("P", "-P")
highlight.sector(c("P"), track.index = 1, col = "#F8766DFF", 
                 text = "P", cex = 0.8, text.col = "white", niceFacing = TRUE)
# c("NP", "-NP")
highlight.sector(c("NP"), track.index = 1, col = "#00BA38FF", 
                 text = "NP", cex = 0.8, text.col = "white", niceFacing = TRUE)
wTOdv %>% 
  mutate(wTO = ifelse(pval-abs(wTO) < 0.05, wTO, 0)) %>%
  filter(grepl('^P-|-P$', facet) & wTO != 0) %>% 
  tail() %>%
  pivot_wider(names_from = 'group', values_from = 'wTO')

wTOdv %>% 
  filter(grepl('^P-|-P$', facet)) %>%
  ggplot(aes(y = wTO, x = group)) + # fill = facet
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), 
               outlier.shape=NA, notch = T)
  # stat_summary(fun = mean, geom="point", shape=20,
  #              size = 3, color="red", fill="red") 



library(CoDiNA)

Code <- c('V2', 'V3', 'V4', 'V67') # , 'V4', 'V67', 'V8')

w1 <- WTO[[1]]
w2 <- WTO[[2]]
w3 <- WTO[[3]]
w4 <- WTO[[4]]
# w5 <- WTO[[5]]

Data <- list(w1, w2, w3, w4)

DiffNet <- MakeDiffNet(Data = Data,
                       Code = Code)


# continue w/ https://deisygysi.github.io/rpackages/Pack-2

# Clustering the nodes into Φ and Φ̃ categories
# based on median

int_C = quantile(DiffNet$Score_internal, 0.3) # the closer to zero, the better.
ext_C = quantile(DiffNet$Score_Phi, 0.75) # the closer to 1, the better.

Nodes_Groups = ClusterNodes(DiffNet = DiffNet, 
                            cutoff.external = ext_C, 
                            cutoff.internal = int_C)

table(Nodes_Groups$Phi_tilde)

Graph = CoDiNA::plot.CoDiNA(DiffNet, 
                            # cutoff.external = ext_C, cutoff.internal = int_C, 
                            layout = 'layout_components', Cluster = TRUE)


library(ggraph)
library(igraph)
library(ggforce)

Edges <- Graph$Edges
Nodes <- Graph$Nodes

Edges %>% mutate_at('Group', 
                    funs(str_replace_all(., c("^[a-z][.]"="")))) -> Edges

Nodes %>% left_join(tax, by = c('id' = 'unique')) -> Nodes

graph <- graph_from_data_frame(Edges, directed = FALSE, Nodes)

group <- igraph::cluster_louvain(graph)$membership

nodes = plyr::join(Nodes, data.frame(id = igraph::V(graph)$name, 
                                     group = group))

graph = graph_from_data_frame(Edges, directed = FALSE, nodes)

layout = create_layout(graph, layout = 'igraph', algorithm = 'kk')

# test aesthetics https://www.r-bloggers.com/2020/03/ggraph-tricks-for-common-problems/

ggraph(layout) + 
  geom_edge_link(aes(edge_colour = Phi, edge_alpha = Score), width = 1.2,
                 arrow = arrow(
                   angle = 10,
                   length = unit(0.1, "inches"),
                   ends = "last",
                   type = "closed"
                 )) +  # aes(alpha = Score)
  geom_node_point(aes(size = Degree_Total * 2 + 1)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  ggforce::geom_mark_hull(
    aes(x, y, group = group, label=group),
    fill = "grey", color = NA,
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25) +
  theme_void() + 
  theme(legend.position="top") +
  geom_node_point(aes(color = Phylum)) + 
  scale_color_manual(values = getPalette) -> graphSave

ggraph(layout) + 
  geom_edge_link(aes(edge_colour = Group, edge_alpha = Score)) +
  facet_edges(~Phi, scales = "free") +
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

ggsave(graphSave, filename = "codina_facet.png", path = dir, 
       width = 18, height = 8)
