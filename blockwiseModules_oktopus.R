
library(tidyverse)

path <- "~/transcriptomics/oktopus_full_assembly/"
countf <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst.txt"
countf <- list.files(path = path, pattern = countf, full.names = T)

dim(count <- read.delim(countf, sep = "\t"))

#
# Clean low abundance genes
#

sum(rowSums(edgeR::cpm(count)) > 3) / nrow(count)
dim(count <- count[rowSums(edgeR::cpm(count)) > 3, ])

outLierdf <- readRDS(paste0(path, 'outliersdf.rds'))
# outLierdf %>% ggplot() + geom_boxplot(aes(x = id, y = y, color = outlier))
# outLierdf %>% filter(abs(z) < 3) %>% ggplot(aes(y,group = id)) + geom_density()

library(WGCNA)


# Choose a set of soft-thresholding powers

max_power <- 30
powers = c(c(1:10),seq(from = 10, to = max_power, by=1))

# Call the network topology analysis function

allowWGCNAThreads()
sft = pickSoftThreshold(t(count), powerVector = powers, verbose = 5, networkType = "signed")

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])
soft_values <- round(soft_values, digits = 2)
power_pct <- 0.8
softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]
softPower <- min(softPower)

saveRDS(sft, paste0(path, 'sft.rds'))

# stop here

# cat("\nsoftPower value", softPower, '\n')

# run at nigth

# using the iterativeWGCNA
# /Users/cigom/transcriptomics/oktopus_full_assembly/WGCNA
# iterativeWGCNA -i counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst_orfs.txt --wgcnaParameters maxBlockSize=5000,corType=bicor,power=12