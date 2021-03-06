
library(tidyverse)
library(RColorBrewer)
# Test count data metrics
# Transform count data 
# Dimension reduction
# 

path <- "~/transcriptomics/oktopus_full_assembly/"
countf <- list.files(path = path, pattern = "counts_table_length_ajus_gen_level-aproach2-filtered.txt", full.names = T)
count <- read.delim(countf, sep = "\t")
mtd <- read.delim(paste0(path, "metadata.tsv"), sep = "\t") 

# barplot(colSums(edgeR::cpm(count)))
# barplot(colSums(count))

colNames <- names(count)

# count <- edgeR::cpm(count)
count %>% 
  as_tibble(rownames = "gene") %>%
  pivot_longer(all_of(colNames), names_to = "id") %>%
  left_join(mtd, by = "id") -> count_longer

topGenes <- sort(rowSums(edgeR::cpm(count)), decreasing = T)
topGenes <- names(head(topGenes, n = 100))

count_longer %>%
  filter(gene %in% topGenes) %>%
  select(gene, id, value) %>%
  pivot_wider(names_from = id, values_from = value, values_fill = NA) %>%
  mutate_if(is.double, function(x) log2(x+1)) %>%
  data.frame(row.names = 'gene') %>%
  superheat::superheat(., row.dendrogram = T, membership.cols = mtd$Tissue, heat.na.col = 'white')
  # mutate(id = name) %>%
  # mutate(group = substr(name, 1, 14)) %>%
  # separate(col = name, into = c("Tissue", "Sex", "x", "y", "Temp", "SampleType"), sep = "_") -> count_longer


# count_longer %>%
#   select(-value) %>%
#   distinct() %>%
#   select(Tissue, SampleType, id, group, x) %>%
#   mutate(group = gsub("GL[A-B]","GLO" ,group)) %>%
#   mutate(Tissue = ifelse(Tissue %in% c("GLA","GLB"), "GLO", Tissue)) %>%
#   mutate(Tissue = recode(Tissue, GLO = "Optic\nGland",GOV = "Oviducal\nGland", 
#                          LOP = "Optic\nLobe")) %>%
#   mutate(x = recode(x, PR = "PRE",DE = "Spawing", PO = "POST")) %>%
#   rename("stage" = x) -> mtd


# write_delim(mtd, file = paste0(path, "metadata.tsv"), delim = "\t", col_names = T)

table(mtd$group)

library(rstatix)
library(ggpubr)
groups <- mtd$group
table(groups)
groups[which(table(groups) >= 2)]
mtd %>% filter(grepl("LOP", id) | grepl("GL[A-B]", id)) %>% pull(id) -> unpairedReps

# is normal?
source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

count_longer %>% 
  rename("g" = "id", x = value) %>%
  # mutate(x = log2(x +1)) %>%
  is_parametric() %>% arrange(desc(outliers))


qqfun <- function(x) {
  x <- x[x > 0]
  qq <- qqnorm(x, plot.it = F) 
  qq %>% as_tibble()
}

zfun <- function(x) {
  x <- x[x > 0]
  z <- c((x - mean(x)) / sd(x))
  # z %>% as_tibble()
  return(z)
}

xx <- log2(count$GLA_H_PR_AD_24_P+1)
qqfun(xx) %>% rstatix::cor_test(x, y)
plot(zfun(xx), qqfun(xx)$y)

count_longer %>%
  filter(value > 0) %>%
  group_by(id) %>%
  mutate(x = log2(value +1)) %>% # 
  summarise(qqfun(x)) %>%
  mutate(z = zfun(y)) %>%
  mutate(outlier = ifelse(abs(z)>3, TRUE, FALSE)) %>%
  left_join(mtd) %>%
  # rstatix::cor_test(x, y) %>%
  # filter(id %in% unpairedReps[1:2]) %>%
  ggplot(aes(x, y, group = id, shape = outlier)) +
  geom_point(aes(color = z)) +
  scale_shape_manual(name = expression("Outlier-"~sigma), values=c(3,1)) +
  facet_wrap(~id) +
  stat_cor(aes(group = id), method = "pearson", cor.coef.name = "R", p.accuracy = 0.001,
           label.y = 18) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  labs(x = "Expected", y = expression("Observed-"~Log[2]~(x~+1)),
       caption = "Correlation Coefficient (Pearson) between theorical normal distribution (expected) against (observed) distribution in the count data.\nThe observed distribution doesn't fit a linear relatioship to a normal distribution (alfa 5%)") +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
              se = TRUE, na.rm = TRUE) +
  # stat_regline_equation(label.y = 15) +
  scale_color_viridis_c(name = expression(sigma)) +
  guides(color = guide_colorbar(barheight = unit(5, "in"), 
                               barwidth = unit(0.3, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 14))) -> sigmaP



ggsave(sigmaP, filename = "outliers.png", path = path, 
       width = 10, height = 10)

count_longer %>% 
  # filter(id %in% unpairedReps[1:2]) %>%
  filter(value > 0) %>%
  mutate(value = log2(value+1)) %>%
  mutate(stage = factor(stage, levels = c("PRE", "Spawing", "POST"))) %>%
  mutate(group = forcats::fct_reorder(group, as.numeric(stage))) %>%
  ggqqplot(., x = "value", color = 'stage', add.params = list(linetype = "dashed"),
           conf.int = TRUE) +
  labs(x = "Expected", y = expression("Observed-"~Log[2]~(x~+1))) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  facet_wrap(~group, scales = "free_y")  +
  theme(legend.position = "top") -> qqP

ggsave(qqP, filename = "ggqqplot.png", path = path, 
       width = 10, height = 10)

# The Wilcoxon rank sum test is a non-parametric alternative to the independent two samples t-test for comparing two independent groups of samples, in the situation where the data are not normally distributed

count_longer %>% filter(id %in% unpairedReps) %>% 
  filter(value > 0) %>%
  mutate(SampleType = ifelse(grepl("GLA", id),"1", SampleType)) %>%
  mutate(SampleType = ifelse(grepl("GLB", id),"2", SampleType)) %>%
  # filter(Tissue %in% c("Optic\nLobe")) %>%
  # group_by(id) %>% get_summary_stats(value)
  group_by(group) %>%
  # wilcox_effsize(value ~ SampleType, paired = FALSE, conf.level = 0.95) %>%
  wilcox_test(value ~ SampleType, paired = FALSE, conf.level = 0.95) %>%
  add_significance() %>%
  mutate(group1 = paste0(group,"_",group1),
         group2 = paste0(group,"_",group2)) %>%
  mutate(group1 = ifelse(grepl("GLO_H_PR_AD_24_1", group1),"GLA_H_PR_AD_24_P", group1)) %>%
  mutate(group2 = ifelse(grepl("GLO_H_PR_AD_24_2", group2),"GLB_H_PR_AD_24_P", group2)) -> stat.test

# count_longer %>%
#   mutate(SampleType = ifelse(grepl("GLA", id),"1", SampleType)) %>%
#   mutate(SampleType = ifelse(grepl("GLB", id),"2", SampleType)) %>%
#   filter(id %in% unpairedReps)
#   # filter(Tissue %in% c("Optic\nLobe")) %>%
#   filter(value > 0) %>%
#   group_by(group) %>%
#   t_test(value ~ SampleType, paired = FALSE, conf.level = 0.95) %>%
#   mutate(group1 = paste0(group,"_",group1),
#          group2 = paste0(group,"_",group2)) %>%
#   mutate(group1 = ifelse(grepl("GLO_H_PR_AD_24_1", group1),"GLA_H_PR_AD_24_P", group1)) %>%
#   mutate(group2 = ifelse(grepl("GLO_H_PR_AD_24_2", group2),"GLB_H_PR_AD_24_P", group2)) -> stat.test

# Create a box plot
count_longer %>%
  filter(Tissue %in% c("Optic\nLobe")) %>%
  # filter(id %in% unpairedReps[1:2]) %>%
  filter(value > 0) %>%
  mutate(stage = factor(stage, levels = c("PRE", "Spawing", "POST"))) %>%
  # mutate(id = forcats::fct_reorder(id, as.numeric(SampleType))) %>%
  mutate(value = log2(value+1)) %>%
  ggplot(aes(id, value)) + 
  geom_boxplot(outlier.size = 0, alpha=0.2) + # outlier.shape=NA
  ylim(0, 20) +
  theme_classic(base_family = "GillSans", base_size = 14) +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  labs(x = "", y = expression(~Log[2]~(x~+1))) +
  ggh4x::facet_nested( ~ group, scales = "free", space = "free", nest_line = TRUE, switch = "x") -> p



# Add the p-value manually

stat.test %>% filter(grepl("LOP", group)) -> stats
p + stat_pvalue_manual(stats, label = "p.adj.signif", remove.bracket = F, y.position = c(15, 18, 20)) -> p


# stat.test %>% filter(grepl("GLO", group)) -> stats

#p1 + stat_pvalue_manual(stats, label = "p.adj.signif", remove.bracket = F, y.position = c(15)) -> p1
# p1  + theme(axis.title.y=element_blank(),
#              axis.text.y=element_blank(),
#              axis.ticks.y=element_blank(),
#             axis.line.y = element_blank()) -> p1

# p1 + scale_x_discrete(labels= c("GLA", "GLB")) -> p1
# 
p + scale_x_discrete(labels=rep(paste0("Replicate ", 1:3), 3)) -> p

library(patchwork)
p + p1 + plot_layout(widths = c(4, 1.2), heights = c(4, 1)) -> p2
p2 + labs(caption = "For replicates a Mann-Whitney U test was evaluated (confidence = 0.95)") -> p2

ggsave(p2, filename = "boxplot_replicates_samples_Mann-Whitney.png", path = path, 
       width = 8, height = 6)
  # limma::normalizeQuantiles()

# mean replicates ----

count0 %>%
  as_tibble(rownames = "ID") %>%
  pivot_longer(cols = colNames, names_to = "id") %>%
  left_join(mtd) %>%
  group_by(ID, group) %>%
  summarise(value = mean(value)) %>%
  pivot_wider(names_from = "group", values_from = "value") %>%
  data.frame(row.names = .$ID) %>% select(-ID) -> count

count <- round(count, digits = 3)
file_out <- "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps.txt"
write.table(count, file = paste0(path, '/',file_out), sep = "\t", quote = F)

# PCA ----
PCA <- prcomp(t(log2(count+1)), scale. = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], 
                    PC2 = PCA$x[,2],
                    mtd %>% distinct(group, .keep_all = T))

dtvis %>%
  mutate(Tissue = ifelse(Tissue %in% c("GLA","GLB"), "GLO", Tissue)) %>%
  # mutate(group = ifelse(SampleType != "P", group, "Pool")) %>%
  ggplot(., aes(PC1, PC2)) +
  # geom_point(aes(color = Tissue), size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_label(aes(label = x, fill = Tissue), alpha = 0.9) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  ggforce::geom_mark_ellipse(aes(label = Tissue, group = Tissue, fill = Tissue)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  guides(color = FALSE, fill = FALSE) +
  theme_bw(base_family = "GillSans", base_size = 16) -> p1

factoextra::fviz_eig(PCA, choice = "variance", ncp = 5, geom = "line", addlabels = T, main = "") + theme_bw(base_family = "GillSans", base_size = 18) -> p2


ggsave(p1, filename = "PCA_mean_rep.png", path = path, 
       width = 8, height = 8)
ggsave(p2, filename = "PCA_mean_rep_eig.png", path = path, 
       width = 6, height = 6)
# library(patchwork)
# p1 + inset_element(p2, left = 0, bottom = 0.6, right = 0.4, top = 1)

# compare vs previous method

countf <- list.files(path = path, pattern = "RSEM.gene.counts.matrix", full.names = T)
count0 <- read.delim(countf, row.names = 1)

# names(count0)[!names(count0) %in% colNames] %>% substr(., 1,3) %>% table(.)

dim(count0 <- count0[rownames(count0) %in% rownames(count),]) # select genes filtered
dim(count0 <- count0[, names(count0) %in% colNames]) # select samples filtered

identical(dim(count0), dim(count))
identical(names(count0), names(count))
# PCA0 <- prcomp(t(log2(count0+1)), scale. = FALSE)
# 
# PCA0$x %>% as_tibble(rownames = "id") %>% pivot_longer(-id) %>% mutate(group = "RSEM_gene") -> df1
# PCA$x %>% as_tibble(rownames = "id") %>% pivot_longer(-id) %>% mutate(group = "length_ajus_gen") -> df2
# 
# rbind(df1, df2) %>%
#   pivot_wider(names_from = group, values_from = value) %>%
#   group_by(id) %>%
#   summarise(cor(RSEM_gene, length_ajus_gen))


count %>% as_tibble(rownames = "gene") %>% pivot_longer(colNames) %>% mutate(method = "length_ajus_gen") -> df2 
count0 %>% as_tibble(rownames = "gene") %>% pivot_longer(colNames) %>% mutate(method = "RSEM_gene") -> df1


rbind(df1, df2) %>%
  mutate(group = substr(name, 1, 14)) %>%
  mutate(Tissue = substr(name, 1, 3)) -> df

df %>% 
  mutate(value = log10(value + 1)) %>%
  pivot_wider(names_from = method, values_from = value) %>%
  # sample_n(1000) %>%
  group_by(name) %>%
  mutate(cor = cor(RSEM_gene, length_ajus_gen)) %>%
  arrange(desc(cor)) %>%
  ungroup() %>%
  mutate(name = forcats::fct_reorder(name, cor, .desc = T)) %>%
  ggplot() +
  geom_point(aes(RSEM_gene, length_ajus_gen , color = cor), alpha = 0.8) +
  geom_abline(color = "black", slope = 1, linetype="dashed", alpha=0.5) +
  labs(caption = expression(~Log[10]~("x"~+1))) +
  scale_color_viridis_c(name = "Correlation\n(Pearson)") +
  facet_wrap(~name) +
  theme_classic(base_family = "GillSans", base_size = 16) +
  guides(color = guide_colorbar(barheight = unit(4, "in"), 
                               barwidth = unit(0.3, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 12))) -> pscatter

ggsave(pscatter, filename = "scatter_cor.png", path = path, 
       width = 12, height = 12)

# then, neuropep?

wd <-'~/transcriptomics/Diana_Lara/neuropeptides'
neuropep <- dir(path = wd, pattern = 'neuropeptide.xls', full.names = T)

require('trinotateR')

y <- read_trinotate(neuropep)

length(peptideGenes <- as.character(unique(y$gene_id)))

summary_trinotate(y)

blastp <- split_blast(y, "sprot_Top_BLASTP_hit")
blastx <- split_blast(y, "sprot_Top_BLASTX_hit")
go <- split_GO(y)
pfam <- split_pfam(y)

superheat::superheat(count[rownames(count) %in% peptideGenes, ], row.dendrogram = T, col.dendrogram = T)

sum(peptideGenes %in% rownames(count))  # son 109 genes con neurpep
# dim(peptide <- count[rownames(count) %in% peptideGenes, ])

count %>%
  as_tibble(rownames = "gene") %>% 
  left_join(blastp, by = "gene") %>%
  # left_join(pfam, by = "gene") %>% # filter(ontology %in% "biological_process") %>%
  drop_na() %>%
  group_by(name) %>%
  summarise_at(vars(colNames), sum) -> df

dfm <- data.frame(df[-1], row.names = df$name)
dim(dfm <- dfm[rowSums(cpm(dfm) > 1) >= 2,])
superheat::superheat(dfm, row.dendrogram = T, col.dendrogram = T) -> sheat


# sheat$order.cols


df %>% 
  pivot_longer(colNames, names_to = "id", values_to = "exp") %>%
  left_join(mtd) %>%
  group_by(name) %>%
  # mutate(Normalized = edgeR::cpm(exp)) %>%
  filter(exp > 1) -> df_longer

df_longer %>%
  mutate(name = factor(name, levels = df$name[sheat$order.rows])) %>%
  ggplot(aes(y = name, x = id, fill = log2(exp+1))) +
  geom_tile() +
  geom_tile(color = "black") +
  facet_grid(~ Tissue, scales = "free", space = "free" ) +
  scale_fill_viridis_c(name = expression(~Log[2]~("x"~+1))) +
  labs(x = "", y = "Neuropeptide (blastp)") +
  theme_classic(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
  guides(fill = guide_colorbar(barheight = unit(4, "in"), 
                                barwidth = unit(0.3, "in"),
                                ticks.colour = "black", 
                                frame.colour = "black",
                                label.theme = element_text(size = 12))) -> pheat

ggsave(pheat, filename = "peptide_heatmap.png", path = path, 
       width = 12, height = 12)


# library(DESeq2)
# x <- round(count)
# 
# countData <- x
# 
# colData <- names(countData)
# colData <- data.frame(conditions=factor(colData))
# 
# dds <- DESeqDataSetFromMatrix(
#   countData = countData,
#   colData = colData,
#   design = ~ conditions)
# 
# dds <- estimateSizeFactors(dds)
# colSums(count)
# colSums(counts(dds, normalized = TRUE))


