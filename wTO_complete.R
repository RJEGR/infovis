#!/usr/bin/env Rscript

# wTO_complete.R 
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)

options(stringsAsFactors = FALSE)

# paralellize

require(parallel)
require(doParallel)

cores <- makeCluster(detectCores()-2, type='PSOCK')

print(cores)

registerDoParallel(cores)

Sys.setenv("MC_CORES" = length(cores))
options("mc.cores"=length(cores))

stopCluster(cores)


# ==============
## Checking and Load packages ----
# ==============

.cran_packages <- c("wTO", "CoDiNA") # "tidyverse"

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
} 
  
# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)


degsf <- args[1]
cuttof <-  0.05 # args[2]

fileOut <- basename(degsf)
fileOut <- sapply(strsplit(fileOut, "[.]"), `[`, 1)

# load count data ----

path <- getwd()
pattern <- "counts_table_length_ajus_gen_level-aproach2-filtered.txt"
countf <- list.files(path = path, pattern = pattern, full.names = T)

dim(count <- read.delim(countf, sep = "\t"))

# load overlaps ----

peptideGenes <- readLines(paste0(path, "/peptideGenes.list"))

# load up/down regulated overlaps  (lists) ----

readGenes <- function(file) {
  
  df <- read.delim(file, sep = "\t") 
  return(df)
}

df <- readGenes(degsf)

df <- subset(df, FDR < cuttof)

degs <- unique(df$ID)

cat("\nnumber of overlaps in the data: ", sum(degs %in% peptideGenes))

Overlap <- degs[degs %in% peptideGenes]

cat("\nNumber of genes to input into the network: ", sum(rownames(count) %in% degs))

# subset count-matrix based on significance and prevalence


dim(Data <- count[rownames(count) %in% degs,])

Data <- log2(Data+1)


network <- wTO::wTO.Complete(n = 100, k = 24,  Data = Data, method_resampling = 'Bootstrap', Overlap = Overlap, method = 's', pvalmethod = "bonferroni", plot = F, savecor = T) 


saveRDS(network, file = paste0(path, "/", fileOut, "_network.rds"))

quit(save = "no")
