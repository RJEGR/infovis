wd <-'~/transcriptomics/Diana_Lara/neuropeptides'
# options(stringsAsFactors = F)

# count <- dir(path = wd, pattern = 'iso.counts.matrix', full.names = T)
path <- "~/transcriptomics/oktopus_full_assembly/"
countf <- list.files(path = path, pattern = "counts_table_length_ajus_gen_level-aproach2-filtered.txt", full.names = T)
neuropep <- dir(path = wd, pattern = 'neuropeptide.xls', full.names = T)

# 


data <- read.delim(countf, sep = "\t")

# get colors ----
cond <- names(data)
rep <- paste0("LOF_", sapply(strsplit(cond, "_"), `[`, 2))
sam <- unique(rep)
n <- length(sam)

sample_color <- setNames(viridis::viridis(n), sam)

#

library(edgeR)
library(DESeq2)
library(tidyverse)

x <- round(data)
countData <- x
#countData <- x[rowSums(cpm(x) > 1) >= 2,]

colData <- names(countData)
colData <- data.frame(conditions=factor(colData))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ conditions)

dds <- estimateSizeFactors(dds)

# annotation ----
require('trinotateR')

y <- read_trinotate(neuropep)
summary_trinotate(y)

pfam <- split_pfam(y)
spfam <- summary_pfam(pfam)

go <- split_GO(y)
gos <- summary_GO(go)

blastx <- split_blast(y, "sprot_Top_BLASTX_hit")
sblastx <- summary_blast(blastx)

blastp <- split_blast(y, "sprot_Top_BLASTP_hit")
sblastp <- summary_blast(blastp)

data(cogs)
download.file("https://raw.githubusercontent.com/RJEGR/infovis/master/NOG.annotations.tsv", "NOG.annotations.tsv")
egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")
names(egg) <- c("db", "nog", "proteins", "species", "class", "description")



plot_NOGs(y, "transcript_id")

system("rm NOG.annotations.tsv")

htmlwidgets::saveWidget(widget, paste0(wd, "/pfam.html"))


# parse ----
sum(rownames(data) %in% y$gene_id)
# signifGenes <- unique(y$transcript_id)
signifGenes <- as.character(unique(y$gene_id))

counts(dds, normalized = TRUE)[signifGenes, ] %>%
  as_tibble(rownames = 'transcript') %>% 
  pivot_longer(-transcript) %>%
  mutate(sample = paste0("LOF_", sapply(strsplit(name, "_"), `[`, 2))) %>%
  group_by(transcript, sample) %>%
  summarise(mean = mean(value)) %>%
  # ungroup() %>%
  # filter(mean > 0) %>%
  pivot_wider(names_from = sample, values_from = mean,  
              values_fill = list(mean = 0)) %>%
  #filter_at(any_vars(sam), rowSums(. !=0) >2)
  inner_join(blastp) %>%
  ungroup() -> datavis
  #inner_join(pfam) -> datavis

write.table(blastp, file = paste0(wd,'/neuropeptide_blastp.txt'))
write.table(datavis, file = paste0(wd,'/neuropeptide_counts.txt'))

datavis %>% select(sam, name) %>%
  group_by(name) %>%
  summarise_at(vars(sam), sum) %>%
  ungroup() %>%
  as.data.frame() -> dataheat
  #mutate_at(vars(sam), function(x) {x / sum(x) * 100 })

dataheat %>% select(-name) %>% rowSums() -> rowsms
plot(density(rowsms / ncol(dataheat)-1))
abline(v=median(rowsms / 6), col = 'red')
abline(v=mean(rowsms / 6), col = 'blue')

x_order <- dataheat[order(dataheat$LOF_24DES), 'name']

barp <- reshape2::melt(dataheat) %>%
  mutate(Temp = substr(variable, 5, 6)) %>%
  mutate(state = substr(variable, 7, length(variable))) %>%
  mutate(state = factor(state, levels = c('PRE','DES','POST'))) %>%
  mutate(name = factor(name, levels = x_order)) %>%
  ggplot(aes(x = as.numeric(name), y = log2(value+0.5), fill = Temp))+
  geom_bar(stat="identity", color = 'black',
           position=position_dodge()) +
  scale_fill_manual(values = c("#2b83ba", "#d7191c"), name = '') +
  theme_bw() +
  scale_x_discrete(limits = x_order) +
  theme(axis.text.x = element_text(angle = 70,
                                   hjust = 1, size = 7
  )) + facet_grid(state~.) +
  labs(x = 'Neuropeptide associated terms (Uniprot-ORFS)',
       y = 'Log2(Normalized Counts + 0.5)')

# library(ggforce)
# barp +
#   facet_zoom(x = value > 1000) +
#   scale_x_continuous(
#     breaks = 1:length(x_order),
#     label = x_order)


rownames(dataheat) <- dataheat$name
dataheat$name <- NULL
xt.plot <- rowSums(dataheat)

png(paste0(wd,"/abundant_neuropeptides.png"), 
    height = 11, width = 14, res = 300,
    units = 'in')

superheat::superheat(dataheat, 
                     scale = T,
                     yr = xt.plot,
                     yr.plot.type = 'bar',
                     yr.axis.name = 'Normalized\nCount',
                     pretty.order.rows = T, 
                     col.dendrogram = T ,
                     left.label.text.size = 2.5,
                     left.label.col = 'white',
                     bottom.label = 'variable',
                     bottom.label.text.size = 3.0)

dev.off()
