# Ricardo Gomez-Reyes, 2020, July
# Use:  Rscript --vanilla WGCNA.R iso.counts.matrix MIN_CPM <round> MIN_REPS <round>

# https://ramellose.github.io/networktutorials/wgcna.html

system('module load R-3.5.0')

library(WGCNA)
library(DESeq2)
library(flashClust)
library(edgeR)

options(stringsAsFactors = FALSE)

# Loadding data
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  count_arg <- "iso.counts.matrix$"
  MIN_CPM <- 1
  MIN_REPS <- 2 

} else {
  count_arg = args[1]
  MIN_CPM = as.numeric(args[2])
  MIN_REPS = as.numeric(args[3])

}

path <- getwd()

count_file <- list.files(path = path, pattern = count_arg, 
                         full.names = T)

count <- read.table(count_file, header=T, com='', 
                    row.names=1, check.names=F, sep='\t', 
                    stringsAsFactors = FALSE)

count <- count[rowSums(edgeR::cpm(count) > MIN_CPM) >= MIN_REPS, ]

conditions <- sapply(strsplit(names(count), "_"), `[`, 2)


conditions <- data.frame(conditions=factor(conditions))

rownames(conditions) <- colnames(count)

count <- round(count)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = conditions,
  design = ~ conditions )

dds <- estimateSizeFactors(ddsFullCountTable)


datExpr <- counts(dds, normalized = TRUE)
datExpr <- t(datExpr)

# datExpr <- t(head(datExpr, n = 100))

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

# filter if !TRUE ----
# #If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
      if (sum(!gsg$goodSamples)>0)
          printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
      datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
      }

saveRDS(datExpr, file = paste0(path, '/', 'datExpr.rds'))

# choosing a set of soft-thresholding powers

max_power <- 30

powers = c(c(1:10), 
           seq(from = 10, to = max_power, by=1)) 

#powers = unique(powers)
allowWGCNAThreads()

sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 5, networkType = "signed") # call network topology analysis function

saveRDS(sft, file = paste0(path,'/', 'sft.rds'))


sizeGrWindow(9, 5)
par(mfrow= c(1,2))

cex1 = 0.9
h_line = 0.9

plot(sft$fitIndices[,1], - sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex = cex1, col="red")
abline(h = h_line, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

#from this plot, we would choose the intersect power because it's the lowest power for which the scale free topology index reaches 0.90


# Construct a gene co-expression matrix and generate modules ----

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

softPower <- sft$fitIndices[,1][which(round(soft_values) >= cex1)]

# if(sum(soft_values >= 0.9) < 1) {
#   softPower <- 3 # random value from 
# } else
#   softPower <- min(softPower)

softPower <- min(softPower)

#build a adjacency "correlation" matrix
# enableWGCNAThreads()
allowWGCNAThreads()

# specify network type
adjacency = adjacency(datExpr, power = softPower, type = "signed") 

saveRDS(adjacency, file = paste0(path,'/', 'adjacency.rds'))

# heatmap(adjacency)

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:

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

# plots tree showing how the eigengenes cluster together

plot(METree, main = "Clustering of module eigengenes", 
     xlab= "", sub = "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0

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
     moduleColors, geneTree, 
     file = paste0(path,'/', fileName))

# Because the position in dataExpr is the same in moduleColors try:
# 
# colnames(datExpr)[moduleColors == "brown"]

quit(save = 'no')