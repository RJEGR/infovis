# lets make wTO network per tissue / stage (n = 9 networks)
# using rownaes(count) as overlaps, 
# then, make a consensus network
# use the non-average dataset

rm(list = ls())


.cran_packages <- c("wTO", "CoDiNA") # "tidyverse"

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

library(tidyverse)

path <- "~/transcriptomics/oktopus_full_assembly/"
# countf <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst.txt"
countf <- "counts_table_length_ajus_gen_level-aproach2-filtered.txt"
countf <- list.files(path = path, pattern = countf, full.names = T)

dff <- readRDS(paste0(path, "upgenes.rds")) %>% mutate(Develope = recode(Develope, DES = "Spawing"))
dff %>% group_by(Tissue, Intercept, Develope) %>% summarise(n = length(ID))

dff %>% ggplot(aes(logFC, FDR)) + facet_grid(Develope~Tissue) + geom_point()

upsetDegs <- function(dff, names_from = 'Tissue') {
  
  library(UpSetR)
  # value_fun <- function(x) sum(x)/length(x)
  value_fun <- function(x) ifelse(sum(x) > 0, 1, 0)
  
  dff %>% 
    select_at(vars('ID', names_from, 'logFC')) %>% 
    pivot_wider(names_from = names_from, values_from = logFC, 
                values_fn = value_fun, values_fill = 0) %>%
    select(-ID) %>% as.data.frame() -> booleandf
  
  upset(booleandf, number.angles = 45, point.size = 2.3, line.size = 1.5, 
        order.by = "freq", scale.intersections = 'identity', mainbar.y.label = "Gene Intersections", sets.x.label = "Genes per sample", text.scale = c(1.3, 1, 1.3, 1, 1, 1))
}
upsetDegs(dff)
upsetDegs(dff, names_from = 'Develope')

#
# match DEGs into the count file ----
#

dim(count <- read.delim(countf, sep = "\t"))


# sum(rowSums(edgeR::cpm(count)) > 3) / nrow(count)
# dim(count <- count[rowSums(edgeR::cpm(count)) > 3, ])

mtd <- read.delim(paste0(path, "metadata.tsv"), sep = "\t")

library(DESeq2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = round(count),
  colData = mtd, 
  design = ~ 1 ) # if not rep use design = ~ 1

dds <- estimateSizeFactors(ddsFullCountTable)
vsd <- varianceStabilizingTransformation(dds)
vsd_df <- as.data.frame(assay(vsd))

file_out <- "counts_table_length_ajus_gen_level-aproach2-filtered_vsd.txt"
write.table(vsd_df, file = paste0(path, '/',file_out), sep = "\t", quote = F)

#
# Select sample group ----
# 

# table(mtd$Tissue);table(mtd$stage)
samGroup <- 'GL' # groups: GOV(3), LOP(9), and GL*(4)

mtd %>% filter(grepl(samGroup, group)) %>% pull(Tissue) %>% unique() -> sam
degs <- dff %>% filter(grepl(sam, Tissue)) %>% pull(ID) %>% unique()

length(degs)

dim(vsd_df %>% select_at(vars(contains(samGroup))) -> Data)
dim(Data %>% filter(rownames(.) %in% degs) -> Data)
overlaps <- rownames(Data)

# run within the cluster !!

network <- wTO.Complete(n = 1000, k = 2,  
                        Data = Data, method_resampling = 'Bootstrap', 
                        Overlap = overlaps, method = 'p', 
                        pvalmethod = "bonferroni", plot = F, 
                        savecor = T)




