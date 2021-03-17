# lets make wTO network per tissue / stage (n = 9 networks)
# using rownaes(count) as overlaps, 
# then, make a consensus network

rm(list = ls())


.cran_packages <- c("wTO", "CoDiNA") # "tidyverse"

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

library(tidyverse)

path <- "~/transcriptomics/oktopus_full_assembly/"
countf <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst.txt"
countf <- list.files(path = path, pattern = countf, full.names = T)

#
# Clean low abundance genes ----
#


dim(count <- read.delim(countf, sep = "\t"))
sum(rowSums(edgeR::cpm(count)) > 3) / nrow(count)
dim(count <- count[rowSums(edgeR::cpm(count)) > 3, ])

# mtd <- read.delim(paste0(path, "metadata.tsv"), sep = "\t")

#
# Select sample group ----
# 

samGroup <- 'GLO'

count %>% select_at(vars(contains(samGroup))) -> Data
dim(Data <- Data[rowSums(Data) > 3, ])
overlaps <- rownames(Data)
network <- wTO.Complete(n = 100, k = 1,  
                             Data = Data, method_resampling = 'Bootstrap', 
                             Overlap = overlaps, method = 'p', 
                             pvalmethod = "bonferroni", plot = F, 
                             savecor = T)




