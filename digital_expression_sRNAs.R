
# 16/07/21, Ricardo Gore
# Digital Expression values of the miRNA biogenesis genes measured by RNA-seq analysis for Mg and Cg (expressed as percentage of El1a).
#  Based on total mapped reads, we computed TPM values and we used elongation factor 1 α (El1α) as normalizer housekeeping gene to compare the expression level of the different genes in each sample

path <- '~/Documents/DOCTORADO/'


x <- read.csv('~/Documents/DOCTORADO/digital_expression_sRNA_biogenesis_umbertoEtAl2016.csv', header = T)
mtd <- read.csv('~/Documents/DOCTORADO/Supplementary_File_2_sRNA_biogenesis_umbertoEtAl2016.csv', header = T) 

mtd %>% mutate_all(., funs(str_replace_all(., c("SRR" = "")))) %>% as_tibble() -> mtd

battery <- c('DROSHA', 'DGCR8', 'XPO5', 'DICER', 'TARBP2')


library(ggplot2)
library(tidyverse)

x %>% select_if(is.double) %>% names() -> samNames

count <- x %>% select_at(samNames)
rownames(count) <- x$ID

count <- apply(count, 2, function(x) (x / max(x))*100)
colSums(count)

sapply(strsplit(colnames(count), "_"), `[`, 1) -> membership.cols

# superheat::superheat(count[-1,], row.dendrogram = T, membership.cols = membership.cols,
                     # bottom.label.text.angle = 90) # log2(count+1)

annotdf <- x %>% select(ID, Annotation)

cbind(annotdf, count) %>% 
  pivot_longer(cols = samNames, names_to = 'sample', values_to = 'de') %>%
  separate(col = sample, into = c('group', 'SRA.ID'), sep = '_') %>%
  mutate(sam = paste0(group, SRA.ID)) -> xlonger

xlonger %>% left_join(mtd) -> xlonger

# xlonger %>%
#   mutate(Relationship = recode_factor(Relationship, P = 'Pathogenic', 
#                                       NC = 'Inconsistent', C = 'non-pathogenic'))

# annotdf <- xlonger %>% distinct(ID, Annotation)

#  most of the miRNA biogenesis genes were expressed at remarkable levels during the early stages of the oyster development: mainly from two cells to the rotary movement and, for some genes, also in the next developmental stages until D-shaped larvae, with no detectable signals afterward in spat and juveniles. Hence, these genes are particularly active in the early development, in particular one AGO (CGI_10020511) and two PIWI transcripts from the egg to trocophora (Fig. 5). In the same developmental stages we also noticed a remarkable expression of the key miRNA genes, with the co-expression of DROSHA and DGCR8 evident in all the analyzed samples.


xlonger %>% group_by(Annotation) %>% summarise(n=sum(de)) %>% arrange(desc(n)) 

# encontramos que durante el desarrollo larval hay mayor actividad de la maquinaria de la biogenesis de sRNAs

xlonger %>%
  filter(Annotation != 'EL1a') %>%
  # mutate(Description = str_replace(Description, "developmental stage ", "")) %>%
  ggplot(aes(x = sam, y = de, fill = Sample.type)) +
  geom_col() +
  # facet_grid(~group, scales = 'free_x', space = 'free_x') +
  labs(y = '% El1a (TPM)', x = 'Sample ID') +
  theme_bw(base_family = "GillSans") +
  theme(
    legend.position = 'top',
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 7),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave


ggsave(psave, filename = 'digitalExp_sRNA_biogenesis.png', path = path, width = 13, height = 7)

# only developmental stage

xlonger %>% 
  filter(Annotation != 'EL1a' & group == 'DEV') %>%
  mutate(Description = str_replace(Description, "developmental stage ", "")) -> devdf

devdf %>% group_by(Description) %>% summarise(n=sum(de)) %>% arrange(desc(n)) %>% 
  pull(Description) -> levels

# basado en cgigas_embryo_larva_development_SEM.pdf

c('Egg', 'Two cells', 'Four Cells', 'Morula', 'Blastula','Rotary movement' , 
  'free swimming', 'Gastrula', 'Trocophore',
  'Dshaped', 'Umbo', 'Spat', 'Juvenile') -> stagesLev


devdf %>%
  arrange(match(stage, stagesLev)) %>%
  # mutate_if(is.character, as.factor) %>%
  # mutate(stage = factor(stage, levels = stagesLev)) %>%
  mutate(Description = factor(Description, levels = unique(Description))) %>%
  ggplot(aes(x = Description, y = de, fill = stage)) +
  geom_col() + 
  labs(y = '% El1a (TPM)') +
  theme_bw(base_family = "GillSans") +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.background = element_blank(),
    panel.border = element_blank()) -> psave

ggsave(psave, filename = 'developmental_sRNA_biogenesis.png', path = path, width = 12, height = 7)

# x %>% select_if(Negate(is.double))


dim(countD <- count[,grepl('DEV', samNames)])

PCA <- prcomp(t(log2(countD+1)), scale. = FALSE) # log2(count+1))
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], 
                    PC2 = PCA$x[,2])
                    # mtd %>% distinct(group, .keep_all = T))


dtvis %>%
  mutate(id = rownames(.)) %>%
  separate(col = id, into = c('group', 'SRA.ID'), sep = '_') %>%
  left_join(mtd) %>% 
  mutate(Description = str_replace(Description, "developmental stage ", "")) %>%
  mutate(subdesc = str_replace(Description, "[1-9]$", "")) %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # geom_text(aes(label = Description), alpha = 0.9) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top')+
  # coord_fixed(ratio = sd_ratio) +
  ggforce::geom_mark_ellipse(aes(label = subdesc, group = subdesc)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") 
  # guides(color = FALSE, fill = FALSE)



# test WGCNA network using log2 transformed data ----

library(WGCNA)
library(flashClust)

power_pct <- 0.9

count <- x %>% select_at(samNames)
rownames(count) <- annotdf$ID

# count <- count[-1,] # Considera meter EL1a como control negativo, ya que no interactua con alguna proteina de la biogenesis

datExpr <- t(log2(count+1)) # log2(count+1) # 

str(datExpr)

cat("\n:::::\n")

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

# detect max power

max_power <- 50

powers = c(c(1:10), seq(from = 10, to = max_power, by=1)) 

#powers = unique(powers)
allowWGCNAThreads()

sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 5, 
                        networkType = "signed")


# Construct a gene co-expression matrix and generate modules ----

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

power_pct <- 0.6

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit,signed R^2'
title2 = 'Mean Connectivity'
sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  geom_text(aes(label = Power)) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
       caption = "lowest power for which the scale free topology index reaches the 90%") +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  # scale_x_continuous(position = "top") +
  theme_light()

#


enableWGCNAThreads()

# specify network type

adjacency <- adjacency(datExpr, 
                       power = softPower, 
                       type = "signed")

# translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:


TOM <- TOMsimilarity(adjacency, TOMType = "signed") # specify network type

heatmap(TOM)

dissTOM = 1 - TOM

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", 
     main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module

minClusterSize <- 0

dynamicMods <- cutreeDynamic(dendro= geneTree, 
                             distM = dissTOM,
                             method = "hybrid",
                             deepSplit = 2, 
                             cutHeight = 0.8,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minClusterSize)

dynamicColors = labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)

MEList = moduleEigengenes(datExpr, 
                          colors= dynamicColors,
                          softPower = softPower)

MEs = MEList$eigengenes

MEDiss= 1 - cor(MEs)

METree = flashClust(as.dist(MEDiss), method= "average")

MEDissThres = 0.0

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, dynamicColors, c("Modules"), 
                    dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

#plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic\n(cutHeight=0.2)"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1


MEs <- mergedMEs


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))
print(nm)

# Plot TOM
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
TOMplot(plotTOM, dendro = geneTree, Colors = moduleColors)
# Relate gene expression modules to traits ----

# Define number of genes and samples

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

sam <- rownames(datExpr)
x <- sapply(strsplit(sam, "_"), `[`, 1)

g <- unique(x)

binSamples <- data.frame(G = as.numeric(g[1] == x),
                         DG =  as.numeric(g[2] == x),
                         HAE =  as.numeric(g[3] == x),
                         MUS =  as.numeric(g[4] == x),
                         MAN =  as.numeric(g[5] == x),
                         GON =  as.numeric(g[6] == x),
                         DEV =  as.numeric(g[7] == x))

datTraits <- data.frame(row.names = sam,  binSamples)

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, 
                        moduleColors)$eigengenes
MEs = orderMEs(MEs0)


moduleTraitCor = cor(MEs, datTraits, use= "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits

textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix) == dim(moduleTraitCor)

# sizeGrWindow(5,5)
# sets the bottom, left, top and right margins respectively
par(mar= c(3, 30, 3, 3))

png(paste0(path, "/eigengenes_digitalExpression_sRNA_biogenesis.png"), res = 120)

# Display the corelation values with a heatmap plot
pal <- ggsci::pal_gsea(palette = c("default"), n = 30, 
                       alpha = 1, reverse = T)(30)
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= str_replace_all(names(MEs), '^ME', ''),
               colorLabels= FALSE,
               colors= pal,
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.65,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))

#


dev.off()

# Intramodular connectivity
# We would find modules containing genes w/ high positive / negative correlations in spite of a variable (for example HC concetration) previosly correlated with plotEigengeneNetworks
# Calculate the correlations between modules

geneModuleMembership <- as.data.frame(WGCNA::cor(datExpr, MEs, use = "p"))
superheat::superheat(geneModuleMembership, pretty.order.rows = T)

#datKME = signedKME(datExpr, MET)

# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# What's the correlation for the trait?
geneTraitSignificance <- as.data.frame(cor(datExpr,datTraits, use = "p"))

# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))
names(geneTraitSignificance) <- paste("GS.", names(datTraits), sep = "")
names(GSPvalue) <- paste("p.GS.", names(datTraits), sep = "")


GSpval <- GSPvalue %>%
  tibble::rownames_to_column(var = "ids")

gMM_df <- geneModuleMembership %>%
  tibble::rownames_to_column(var = "ids") %>%
  gather(key = "moduleColor", value = "moduleMemberCorr", -ids) 

gMM_df %>% 
  mutate(moduleColor = str_replace_all(moduleColor, '^ME', '')) %>%
  # mutate(moduleColor = factor(moduleColor, levels = yLabels)) %>% 
  as_tibble() -> gMM_df

# Prepare gene significance df
GS_bacprod_df <- geneTraitSignificance %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ids")


# Put everything together 
allData_df <- gMM_df %>%
  left_join(GS_bacprod_df, by = "ids") %>%
  left_join(GSpval, by = "ids") 
  # mutate(degs = ifelse(ids %in% updegs, 'up', 
  #                      ifelse(ids %in% dwndegs, 'down', 'ns')))

allData_df %>% group_by(moduleColor) %>% count() -> dat_text

# considere Development results (GS.DEV)

allData_df %>%
  mutate(moduleMemberCorr = abs(moduleMemberCorr), GS = abs(GS.DEV)) %>%
  ggplot(aes(moduleMemberCorr, GS, color = p.GS.DEV)) +
  geom_point(size = 0.7, alpha = 0.7) +
  scale_color_viridis_c(name = "Significance",
                        breaks = c(0, 0.25 ,0.5, 0.75, 1), 
                        # labels=c("Minimum",0.5,"Maximum"),
                        limits= c(0, 1)) +
  facet_wrap(~ moduleColor) +
  labs(x = 'Module membership', y = 'Gene Significance') +
  guides(color = guide_colorbar(barheight = unit(3, "in"),
                                ticks.colour = "black", 
                                frame.colour = "black",
                                label.theme = element_text(size = 12)))

##

ggplot() + geom_point(data = allData_df,
                      aes(moduleMemberCorr, abs(GS.DEV), alpha = p.GS.DEV, 
                          color = p.GS.DEV), size = 0.7) + 
  facet_wrap(~ moduleColor) + theme_classic(base_size = 16, base_family = "GillSans") +
  labs(x = 'Module membership', y = 'Gene Significance') 

# continue w/ heatmap of module of clusters 
# prapare data for network w/ ggplot

moduleTraitCor %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'moduleTraitCor') -> df1
moduleTraitPvalue %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

df1 %>%
  # filter(name %in% c('HC', 'Ctrl')) %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
                       ifelse(corPvalueStudent <.01, "**",
                              ifelse(corPvalueStudent <.05, "*", "")))) -> df1

df1 %>%
  # filter(module %in% psME & name != 'Time') %>%
  # mutate(moduleTraitCor = ifelse(abs(moduleTraitCor) < 0.5, NA, moduleTraitCor)) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  geom_tile(color = 'black', size = 0.5) + 
  geom_text(aes(label = star),  vjust = 0.75, hjust = 0.5, size= 6) +
  ggsci::scale_fill_gsea(name = "Membership", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  labs(x = '', y = 'Module') +
  guides(fill = guide_colorbar(barwidth = unit(3, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top") -> p1

#



moduleColors %>% as_tibble(., rownames = 'transcript') %>%
  group_by(value) %>% count(sort = T) %>% mutate(module = paste0('ME',value))-> ModuleDF

ModuleDF %>% mutate(module = factor(module, levels = hclust$labels[hclust$order])) -> ModuleDF

ModuleDF %>% ungroup() %>% mutate(pct = n / sum(n)) -> ModuleDF

ModuleDF %>% 
  # mutate(degs = ifelse(degs == 'ns', '', degs)) %>%
  ggplot(aes(x = module, y = n)) + # , fill = degs
  labs(y = 'Number of transcripts') +
  geom_col() + coord_flip() +
  # geom_col(color = 'black', size = 0.25) + coord_flip() +
  scale_fill_manual(name = '', values = c("grey30", "#EE4141", "#2428D0")) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
        axis.title.y = element_blank(), axis.text.y= element_blank(),
        axis.ticks.y=element_blank(), axis.line.y = element_blank()) -> p2

library(patchwork)

caption = c("*** corPvalueStudent <.001\n** corPvalueStudent <.01\n* corPvalueStudent < 0.05")
p1 + p2 + plot_layout(widths = c(0.5, 1)) +
  labs(caption = caption) -> psave

ggsave(psave, filename = 'ModuleTraitRelationship_digitalExpression_sRNA_biogenesis.png', path = path, width = 12, height = 7)


hclust <- hclust(dist(geneTraitSignificance), "complete")

geneTraitSignificance %>% 
  as_tibble(rownames = 'nodeName') %>% 
  left_join(annotdf, by = c('nodeName' = 'ID')) %>%
  left_join(moduleColors %>% as_tibble(., rownames = 'nodeName')) %>%
  mutate(nodeName = factor(nodeName, levels = hclust$labels[hclust$order])) %>%
  mutate(Module = value) %>% select(-value) -> heatDF

heatDF %>%
  pivot_longer(cols = names(geneTraitSignificance)) %>%
  ggplot(aes(y = nodeName, x = name, fill = value)) +
  geom_tile(color = 'black', size = 0.5) + 
  ggsci::scale_fill_gsea(name = "Membership", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  # ggh4x::scale_y_dendrogram(hclust = hclust) +
  facet_grid(Module~., scales = 'free_y', space = 'free') +
  labs(x = '', y = 'Gene') +
  guides(fill = guide_colorbar(barwidth = unit(3, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top") -> psave

ggsave(psave, filename = 'digital_expression_sRNAs_biogenesis-geneTraitSignificance.png', 
       path = path, height = 9, width = 7) 


# export net

library(ggraph)
library(igraph)
library(tidygraph)
library(WGCNA)
library(tidyverse)

nodeNames <- names(moduleColors)
nodeAttr <- moduleColors

altNodeNames <- annotdf %>% filter(Annotation != 'EL1a')

identical(nodeNames, altNodeNames$ID)

Net = exportNetworkToCytoscape(TOM,
                               # edgeFile = edgeFile,
                               # nodeFile = nodeFile,
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = nodeNames, nodeAttr = nodeAttr, 
                               altNodeNames = altNodeNames$Annotation)

Nodes <- Net$nodeData
names(Nodes) <- c('nodeName', 'altName', 'Module')

g <- tbl_graph(nodes = Nodes, 
                   edges = Net$edgeData, directed = FALSE)

hist(E(g)$weight)
boxplot(E(g)$weight)

quantile(E(g)$weight, probs = c(0.5)) -> qq

g %>% activate(edges) %>% filter(weight > qq) -> g

g %>% activate(edges) %>% as.data.frame() %>% view()

g %>% activate("nodes") %>%  
  mutate(betweenness = betweenness(.), 
         degree = centrality_degree(),
         membership = components(.)$membership,
         transitivity = transitivity(.)) -> g

geneTraitSignificance %>% 
  as_tibble(rownames = 'nodeName') %>% 
  # left_join(annotdf, by = c('nodeName' = 'ID')) %>%
  arrange(match(nodeName, igraph::V(g)$nodeName)) -> nodeInfo

g %>% left_join(nodeInfo) -> g

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

ggraph(layout) +
  geom_edge_link(aes(edge_alpha = weight)) + # edge_width = weight
  geom_node_point(aes(color = GS.MAN, size = degree)) + # , alpha = GS
  geom_node_text(aes(label = altName), repel = TRUE, size = 2) +
  scale_edge_width(range = c(0.2, 2)) +
  # ggsci::scale_color_gsea(name = 'Gene Correlation', reverse = T) +
  scale_color_viridis_c(name = 'Gene Significance', option = 'magma', direction = -1) +
  # scale_color_continuous(name = 'Gene Correlation') + # values = nodeCol
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top")+
  coord_fixed() -> pnet
  # facet_nodes(~ Module)

# library(ggforce)

# pnet +
#   geom_mark_hull(
#     aes(x, y, group = membership),
#     fill = "grey", color = NA,
#     concavity = 4,
#     con.size = 0.3,
#     con.linetype = 2,
#     expand = unit(2, "mm"),
#     alpha = 0.25)

ggsave(pnet, filename = 'digital_expression_sRNAs_biogenesis-Network-dev.png', 
       path = path, height = 7, width = 12) 
