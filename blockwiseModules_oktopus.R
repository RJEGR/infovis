rm(list = ls())

library(tidyverse)

path <- "~/transcriptomics/oktopus_full_assembly/annotation"
filepath = paste0(path, "/Trinity_SAM_extract_properly_mapped_pairs.min1TPM.fasta.transdecoder.pep")



dnaset <- Biostrings::readAAStringSet(filepath)
AAnames <- sapply(strsplit(names(dnaset), "::"), `[`, 1)
AAfreq <- Biostrings::alphabetFrequency(dnaset) 
AAcons <- Biostrings::consensusMatrix(dnaset) 
tAAcons=t(AAcons)
matplot(tAAcons[,-1], type="l", lwd=2, xlab="Read Length", ylab= "Base frequency at each position")
legend(legend = colnames(tAAcons)[-1],"topright",col=1:20, lty=1:20, lwd=2)

barplot(AAfreq[1,])

tAAcons[,-1] %>%
  as_tibble(rownames = 'length') %>%
  pivot_longer(cols = colnames(tAAcons)[-1]) %>%
  # filter(value > 0) %>%
  # sample_n(1000) %>%
  ggplot(aes(x = length, y = value, group = name, color = name)) +
  geom_line() +
  scale_x_discrete(breaks = c(0, 1000, 2000, 3000, 4000))


path <- "~/transcriptomics/oktopus_full_assembly/"
# countf <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst.txt"
countf <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst_orfs.txt"
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

#allowWGCNAThreads()
#sft = pickSoftThreshold(t(count), powerVector = powers, verbose = 5, networkType = "signed")

#soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])
# soft_values <- round(soft_values, digits = 2)
# power_pct <- 0.8
# softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]
#softPower <- min(softPower)

# saveRDS(sft, paste0(path, 'sft.rds'))

# stop here

# cat("\nsoftPower value", softPower, '\n')

# run at nigth

# using the iterativeWGCNA
# /Users/cigom/transcriptomics/oktopus_full_assembly/WGCNA
# iterativeWGCNA -i counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst_orfs.txt --wgcnaParameters maxBlockSize=5000,corType=bicor,power=12

# Dataviz results
dir <- '/Users/cigom/transcriptomics/oktopus_full_assembly/WGCNA/'

dim(wgcna <- read.delim(paste0(dir, 'merged-0.05-membership.txt')))

hist(wgcna$kME)

ggplot(data = blastp) +
  geom_point(aes(identity, evalue), alpha = 0.6) +
  facet_grid(~domain)

blastp %>% distinct(gene, uniprot, .keep_all = T) -> blastp

# dim(y[match(y$gene_id, wgcna$Gene)])
# blastp[blastp$gene %in% wgcna$Gene]
wgcna %>%
  left_join(blastp, by = c('Gene'='gene')) %>% 
  as_tibble()


# load(paste0(dir, 'pass37/i1/wgcna-blocks.RData'))
# tree <- blocks$dendrograms[[1]]

