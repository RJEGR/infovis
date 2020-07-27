# Ricardo Gomez-Reyes, 2020, July
# Use:  Rscript --vanilla WGCNA.R iso.counts.matrix MIN_CPM <round> MIN_REPS <round>

system('module load R-3.5.0')

library(WGCNA)

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

# path <- c('~/transcriptomics/Diana_Lara/Diana_results/')
count_file <- list.files(path = path, pattern = count_arg, 
                         full.names = T)

count <- read.table(count_file, header=T, com='', 
                    row.names=1, check.names=F, sep='\t', 
                    stringsAsFactors = FALSE)

count <- count[rowSums(edgeR::cpm(count) > MIN_CPM) >= MIN_REPS, ]

# if log2 transformation 

count <- readRDS(paste0(path, '/', 'count_prepared.rds'))

conditions <- sapply(strsplit(names(count), "_"), `[`, 2)


conditions <- data.frame(conditions=factor(conditions))

rownames(conditions) <- colnames(count)

count <- round(count)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = count,
  colData = conditions,
  design = ~ conditions )

dds <- estimateSizeFactors(ddsFullCountTable)


# if log2 transform ----

# ex <- count
# qntl <- c(0., 0.25, 0.5, 0.75, 0.99, 1.0)
# qx <- as.numeric(quantile(ex, qntl, na.rm=T))
# 
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# if (LogC) { 
#   ex[which(ex <= 0)] <- NaN
#   ex <- log2(ex) }
# 
# datExpr <- t(ex)

# else ----

datExpr <- counts(dds, normalized = TRUE)
datExpr <- t(log2(datExpr + 1))

# datExpr <- t(head(datExpr, n = 100))

datExpr <- readRDS(paste0(path, '/', 'datExpr.rds'))
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

# quit(save = 'no') 

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

rm(list = ls())

quit(save = 'no')
