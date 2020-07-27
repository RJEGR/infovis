# Ricardo Gomez-Reyes, 2020, July
# Use:  Rscript --vanilla adjacency_WGCNA.R sft.rds datExpr.rds 0.8 <power_pct> nThreads <interger>

system('module load R-3.5.0')

library(WGCNA)
library(flashClust)

# library(edgeR)
# library(DESeq2)
# devtools::install_github("krlmlr/ulimit")
# library(ulimit)

options(stringsAsFactors = FALSE)

# Loadding data
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  sft_file <- 'sft.rds'
  datExpr_file <- 'datExpr.rds'
  power_pct <- 0.9
  nThreads <- 24
  
} else {
  sft_file = args[1]
  datExpr_file = args[2]
  power_pct = as.numeric(args[3])
  nThreads = as.numeric(args[4])
}

path <- getwd()

# path <- '~/transcriptomics/Diana_Lara/Diana_results/'
# 
sft <- readRDS(file = paste0(path, '/', sft_file))
datExpr <- readRDS(file = paste0(path, '/', datExpr_file))


str(datExpr)
cat("\n:::::\n")
str(sft)
cat("\n:::::\n")

# From this plot, we would choose the intersect power because it's the lowest power for which the scale free topology index reaches power_pct


# Construct a gene co-expression matrix and generate modules ----

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 3)

# power_pct <- 0.9
  
softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]


# if(sum(soft_values >= 0.9) < 1) {
#   softPower <- 3 # random value from 
# } else
#   softPower <- min(softPower)

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')

# Build a adjacency "correlation" matrix


# ulimit::memory_limit(2000)
# enableWGCNAThreads()

enableWGCNAThreads()
# specify network type
adjacency <- adjacency(datExpr, power = softPower, type = "signed")

saveRDS(adjacency, file = paste0(path,'/', 'adjacency.rds'))

# heatmap(adjacency, labRow=FALSE, labCol=FALSE)

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:

# allowWGCNAThreads() # nThreads = 24

adjacency <- readRDS(paste0(path,'/', 'adjacency.rds'))

TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1 - TOM

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module

minClusterSize <- 3

dynamicMods <- cutreeDynamic(dendro= geneTree, 
                             distM = dissTOM, 
                             deepSplit = 2, 
                             cutHeight = 0.8,
                             pamRespectsDendro= FALSE,
                             minClusterSize = minClusterSize)

dynamicColors = labels2colors(dynamicMods)

MEList = moduleEigengenes(datExpr, 
                          colors= dynamicColors,
                          softPower = softPower)

MEs = MEList$eigengenes

MEDiss= 1 - cor(MEs)

METree = flashClust(as.dist(MEDiss), method= "average")

# Set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0

MEDissThres = 0.0

merge = mergeCloseModules(datExpr, 
                          dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1

MEs <- mergedMEs

fileName <- 'Network_signed_nomerge_RLDfiltered.RData'


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))

print(nm)

save(MEs, moduleLabels, 
     moduleColors, geneTree, datExpr,
     file = paste0(path,'/', fileName))

# Because the position in dataExpr is the same in moduleColors try:
# 
# colnames(datExpr)[moduleColors == "brown"]

quit(save = 'no')

# Relate gene expression modules to traits ----
# Correlate traits -----

# Define number of genes and samples

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datTraits <- rownames(datExpr)
rep <- sapply(strsplit(datTraits, "_"), `[`, 2)
datTraits <- data.frame(datTraits = as.numeric(as.factor(datTraits)),rep = as.numeric(as.factor(rep)))

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits

textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix) == dim(moduleTraitCor)

par(mar= c(6, 8.5, 3, 3))


# Display the corelation values with a heatmap plot

labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))

# Network ----
library(igraph)

adj <- TOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
g <- graph.adjacency(adj)
g <- simplify(g)  # removes self-loops

# results <- blockwiseModules(datExpr, 
#                             power = softPower, TOMType="signed", 
#                             networkType="signed")

V(g)$color <- moduleColors

par(mar=c(0,0,0,0))

# remove unconnected nodes
g <- delete.vertices(g, degree(network)==0)


# Analicemos propiedades globales de la red
# Longitud promedio de caminos más cortos
mis_caminos <- average.path.length(g)
# clustering coefficient global
mi_clustering <- transitivity(g, type = "global")
# Diametro
mi_diametro <- diameter(g)
# propiedades de nodos ----

# signar nuevas propiedades a los nodos


V(g)$grado <- degree(g)
V(g)$agrupamiento <- transitivity(g, 
                                  type = "local", 
                                  isolates = "zero")

# Propiedades de los enlaces (Edges) ----

E(g)$intermediacion <- edge.betweenness(g)
g1 <- as.undirected(g)
comm.louvain <- cluster_louvain(g1)
V(g)$comm.louvain <- membership(comm.louvain)

mi_df_nodos <- get.data.frame(x = g, 
                              what = "vertices")
mi_df_nodos %>% head()

plot(g,
     # change size of labels to 75% of original size
     vertex.label.cex = .5, 
     vertex.size = degree(g),
     edge.curved=.25, # add a 25% curve to the edges
     edge.color="grey20",
     layout=layout_nicely)
