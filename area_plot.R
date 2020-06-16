# Ricardo Gomez Reyes, 2020, May
#

library(edgeR)
library(DESeq2)
library(tidyverse)

options(stringsAsFactors = F)

wd <- '~/transcriptomics/Diana_Lara/Diana_results/'

setwd(wd)

# Load count.matrix and filter ----
data = read.table("iso.counts.matrix", header=T, row.names=1, com='')

# get colors ----
cond <- names(data)
rep <- paste0("LOF_", sapply(strsplit(cond, "_"), `[`, 2))
sam <- unique(rep)
n <- length(sam)

sample_color <- setNames(viridis::viridis(n), sam)

# work separate to visualze area:

annot <-  read.csv("../relacion_ids.csv", header=T) %>% as_tibble()

# aggregate(annot$transcript_id, by = list(annot$Temp, annot$state, annot$Exp), length)
# parse colors to genes

annot[annot$state == 'Desove', 'state'] <- 'DES'
annot[annot$state == 'Post', 'state'] <- 'POST'

annot %>%
  as_tibble() %>%
  mutate(sample = paste0('LOF_', Temp, state)) %>%
  mutate(genes = ifelse(sprot_Top_BLASTX_hit == '.', 
                        transcript_id, sprot_Top_BLASTX_hit)) -> annot

  #filter(sprot_Top_BLASTX_hit != '.') -> annot
  # distinct(sprot_Top_BLASTX_hit, .keep_all = T) -> annot2

as_tibble(sample_color, rownames = 'sample') %>%
  inner_join(annot) %>%
  select(-sprot_Top_BLASTX_hit, -sample) %>%
  dplyr::rename('color' = 'value') -> annot
  # distinct(transcript_id, .keep_all = T) # -> annot

# annotation ----
# 1) transform read-count data
# 2) 
# filtrar por Exp y Temp
f <- annot %>% filter(Temp == '24' & Exp == 'Up')
signifGenes <- as.character(unique(f$transcript_id))

x <- round(data)

countData <- x[rowSums(cpm(x) > 1) >= 2,]

colData <- names(countData)
colData <- data.frame(conditions=factor(colData))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ conditions)

dds <- estimateSizeFactors(dds)

counts(dds, normalized = TRUE)[signifGenes, ] %>%
  as_tibble(rownames = 'transcript_id') %>% 
  pivot_longer(-transcript_id) %>%
  mutate(sample = paste0("LOF_", sapply(strsplit(name, "_"), `[`, 2))) %>%
  group_by(transcript_id, sample) %>%
  summarise(mean = mean(value), sum = sum(value), n = sum(value > 0)) %>%
  ungroup() %>%
  filter(mean > 0) -> normalized_sig_counts

# normalized_sig_counts %>%
#   group_by(transcript_id, sample) %>%
#   mutate_at(vars(mean), function(x) {x / sum(x) * 100 })
 


normalized_sig_counts %>%
  select(-sum, -n) %>%
  pivot_wider(names_from = sample, values_from = mean,  
              values_fill = list(mean = 0)) %>%
  mutate_at(vars(sam), function(x) {x / sum(x) * 100 }) %>% 
  pivot_longer(-transcript_id, values_to = 'RA', names_to = 'sample') %>%
  filter(RA > 0) -> normalized_sig_counts

f %>%
  #select(transcript_id, genes, value, Exp) %>%
  inner_join(normalized_sig_counts) %>%
  mutate(state = sapply(strsplit(sample, "_"), `[`, 2)) %>%
  mutate(Temp = substr(state, 1, 2)) %>%
  mutate(state = substr(state, 3, length(state))) %>%
  mutate(state = factor(state, levels = c('PRE','DES','POST'))) %>%
  filter(Temp == '24') -> datavis

#
x.color <- setNames(f$color, f$transcript_id)

datavis %>%
  #mutate(mean = log2(mean+0.5))
  ggplot(aes(x = genes, 
           y = RA,
           fill = sample,
           color = sample )) +
  geom_bar(stat="identity") +
  xlab("Genes (Up)") +
  ylab("Normalized Counts (mean)") +
  scale_fill_manual(values = sample_color, name = '') +
  scale_color_manual(values = sample_color, name = '') +
  theme_bw() +
  # scale_x_discrete(levels = x_order) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 5,
                                   color = x.color
  )) + facet_grid(state~., scales = "free")

datavis %>%
  #mutate(mean = log2(mean+0.5))
  ggplot(aes(y=sample, x=RA,  fill=sample)) +
  geom_density_ridges(alpha=0.6, stat="binline", bins=20) +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8))


#  length(unique(datavis$transcript_id)) # 340

# 3) plot


# counts(dds, normalized = TRUE)[signifGenes, ] %>%
#   as_tibble(rownames = 'transcript_id') %>%
#   filter(rownames(.) %in% names(yr_plot_max))
#   # mutate_at(vars(sam), function(x) {x / sum(x) * 100 }) 



# sample_color 

# normalized_sig_counts %>%
#   mutate(s = ifelse(pa != 0, id_sample, NA))

datavis %>%
  filter(Exp == 'Up') -> Up

x.color <- Up$value

Up %>%
  filter(Exp == 'Up') %>%
  # mutate(percentage = mean / sum(mean))
  ggplot(aes(x = genes, 
               y = log2(mean+0.5),
               fill = Temp,
               color = Temp,
               group = sample
             )) +
  geom_area(alpha=0.6 , size=.5, colour="white") +
  #geom_bar(stat="identity") +
  xlab("Genes (Up)") +
  ylab("Normalized Counts (mean)") +
  #scale_fill_manual(values = sample_color, name = '') +
  theme_bw() +
  # scale_x_discrete(levels = x_order) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 5
                                   #color = x.color
  )) + facet_grid(~state, scales = "free")
# ridges ----
# library
library(ggridges)
library(ggplot2)

# Diamonds dataset is provided by R natively
#head(diamonds)

x.color <- datavis$value

Up %>%
  ggplot(aes(x = mean, y = genes)) + # fill  = 
  geom_density_ridges() +
  theme_ridges(font_size = 7) +
  scale_fill_manual(values = sample_color, name = '') +
  theme(legend.position = "top", 
        axis.text.y = element_text(hjust = 1, size = 5, color = x.color)) +
  facet_grid(~Temp, scales = "free")



# Heatmap ---
 
library(superheat)
library(ggheatmap)

# make categorical groups
ggpubr::ggboxplot(normalized_sig_counts, y = 'mean', x = 'sample', yscale = "log2")

normalized_sig_counts %>%
  # convert state to factor and reverse order of levels
  mutate(mean=log2(mean+1)) %>%
  #mutate(mean=factor(mean,levels=rev(sort(unique(mean))))) %>%
  # create a new variable from count
  mutate(meanfactor = cut(mean,breaks=c(-1,0,1,10,100,500,1000, max(mean, na.rm=T)),
                         labels=c("0","0-1","1-10","10-100","100-500","500-1000",">1000"))) %>%
  # change level order
  mutate(countfactor=factor(as.character(meanfactor),levels=rev(levels(meanfactor)))) %>%
  ggplot(aes(x=sample,y=transcript_id,fill=meanfactor))+
  geom_tile(colour="white",size=0.2) +
  guides(fill=guide_legend(title=" "))+
  labs(x="",y="",title=" ")+
  scale_y_discrete(expand=c(0,0))+
  #scale_x_discrete(expand=c(0,0),breaks=c("1930","1940","1950","1960","1970","1980","1990","2000"))+
  scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"),na.value = "grey90") +
  theme_grey(base_size=7)

# or

library(dendextend)
library(gplots)

# normalized by column median (direction: 2) and log2.
# F_m2 <- (apply(F_m, 2, function(x){log2(x/median(x, na.rm = T))}))

normalized_sig_counts %>%
  select(-sum, -n) %>%
  # log2-transformed. 0.5 was added to prevent infinite values resulting from log2(0).
  mutate(mean=log2(mean+0.5)) %>%
  #mutate(mean=log2(mean/median(mean, na.rm = T))) %>%
  pivot_wider(names_from = sample, values_from = mean,  
              values_fill = list(mean = 0)) -> wider_nsc

wider_nsc <- data.frame(wider_nsc[-1], row.names = wider_nsc$transcript_id)

sam <- c( "LOF_24PRE", "LOF_24DES","LOF_24POST",
          "LOF_30PRE", "LOF_30DES","LOF_30POST")

wider_nsc %>% select(sam[4:6]) %>%
  superheat(col.dendrogram = T,
            row.dendrogram = T,
            left.label.text.size = 4.0,
            left.label.col = 'white',
            bottom.label = 'variable',
            bottom.label.text.size = 4.0)

  #mutate(s = ifelse(pa != 0, id_sample, NA)) %>%
  # pivot_wider(-pa, names_from = id_sample,   values_from = s) -> m_wider



normalized_sig_counts %>%
  ggplot(aes(x = sample, y = as.factor(transcript_id), 
             fill = log2(mean+0.5), color = log2(mean+0.5)))+
  geom_tile(color = "grey")

# https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/

