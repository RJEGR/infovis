path <- "~/transcriptomics/oktopus_full_assembly/"
count0 <- read.delim(paste0(path, 'counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst.txt'), sep = "\t")
dim(count0)

dim(count <- count0[rowSums(count0) > 10, ])
# overlaps <- rownames(Data)

sum(rownames(count) %in% degs)/length(degs)

file_out <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst_rowSum_10.txt"
write.table(count, file = paste0(path, '/',file_out), sep = "\t", quote = F)


PCA <- prcomp(t(count), scale. = FALSE) # log2(count+1))
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], 
                    PC2 = PCA$x[,2],
                    mtd %>% distinct(group, .keep_all = T))


dtvis %>%
  mutate(Tissue = ifelse(Tissue %in% c("GLA","GLB"), "GLO", Tissue)) %>%
  mutate(stage = factor(stage, levels = c("PRE", "Spawing", "POST"))) %>%
  # mutate(group = ifelse(SampleType != "P", group, "Pool")) %>%
  ggplot(., aes(PC1, PC2)) +
  # geom_point(aes(color = Tissue, shape = stage), size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_label(aes(label = stage, fill = Tissue), alpha = 0.9) +
  labs(caption = 'Applying a variance stabilizing transformation (VST)\nto the count data from the fitted dispersion-mean;\nthen the data count is transformed (normalized by division by the size factors') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top')+
  coord_fixed(ratio = sd_ratio) +
  ggforce::geom_mark_ellipse(aes(label = Tissue, group = Tissue, fill = Tissue)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  guides(color = FALSE, fill = FALSE)


# 

library(Biostrings)
path <- "~/transcriptomics/oktopus_full_assembly/annotation/"
filepath = paste0(path, "Trinity_SAM_extract_properly_mapped_pairs.min1TPM.fasta.transdecoder.cds")
dnaset <- Biostrings::readDNAStringSet(filepath)
head(names(dnaset))

library(trinotateR)

y <- trinotateR::read_trinotate("~/transcriptomics/oktopus_full_assembly/Trinotate.xls")

# split_pfamdf <- trinotateR::split_pfam(y)
blastp <- split_blast(y, "sprot_Top_BLASTP_hit")


orfs <- sapply(strsplit(names(dnaset), "::"), `[`, 1)
orfs <- unique(orfs)

length(orfs) # 20486

sum(rownames(count0) %in% orfs) # cuantos genes (curados por Oscar y Pavel) presentaron orf
sum(orfs %in% unique(blastp$gene)) # cuantos orf presentaron anotacion funcional
sum(orfs %in% degs) # cuantos de estos orfs son degs??
dim(count0[rownames(count0) %in% orfs, ])

# guardar todos los orfs 

file_out <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps_vst_orfs.txt"
write.table(count0[rownames(count0) %in% orfs, ], file = paste0(path, '/',file_out), sep = "\t", quote = F)



blastp <- blastp[blastp$gene %in% rownames(count),]

table(blastp$domain)

blastp %>% 
  dplyr::count(domain, genus, sort = TRUE) %>%
  mutate(pct = n / sum(n) * 100) %>%
  mutate(genus =  ifelse(pct < 0.1, domain, genus)) %>%
  # mutate(genus = factor(genus, levels = c("PRE", "Spawing", "POST"))) %>%
  mutate(genus = forcats::fct_reorder(genus, n)) %>%
  ggplot(aes(genus, pct)) + geom_col() + coord_flip()

count %>%
  edgeR::cpm(.) %>% 
  as_tibble(rownames = 'gene') %>%
  pivot_longer(cols = colnames(count), names_to = 'sam', values_to = 'exp') -> count_longer


count_longer %>% group_by(sam) %>% dplyr::summarise(exp = sum(exp))

count_longer %>%
  right_join(blastp %>% dplyr::count(gene, domain, sort = TRUE))%>%
  group_by(sam, domain) %>%
  dplyr::summarise(exp = sum(exp)) %>% mutate(pct = exp / 1E6 * 100) %>%
  ggplot(aes(x = sam, y = pct)) + geom_col(aes(fill = domain)) + coord_flip()
  

count_longer %>% filter(gene %in% 'TRINITY_DN70503_c2_g1')
blastp %>% filter(gene %in% 'TRINITY_DN70503_c2_g1')

dnaset[order(width(dnaset), decreasing = T),]

dnaset[grepl('TRINITY_DN70503_c2_g1', names(dnaset))]
