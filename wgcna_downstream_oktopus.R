# test zscore 
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001057
library(tidyverse)
library(WGCNA)

rm(list = ls())

path <- "~/transcriptomics/oktopus_full_assembly/all_deg_ricardo/"

load(paste0(path, 'wgcna_out.RData'))

mtd <- read.delim(paste0(path, "metadata.tsv"), sep = "\t")

# sft <- readRDS(file = paste0(path, '/', 'sft.rds'))
# adjacency <- readRDS(paste0(path,'/', 'adjacency_signed.rds'))

datExpr <- readRDS(paste0(path,'/', 'datExpr.rds'))
degsdf <- readRDS(paste0(path,'/', 'updown_degs.rds'))

degsdf %>% mutate(Exp = ifelse(logFC > 0, 'UP', 'DOWN')) -> degsdf
# updegs <- degs %>% filter(Exp == 'UP') %>% pull(ID) %>% unique()
# dwndegs <- degs %>% filter(Exp == 'DOWN') %>% pull(ID) %>% unique()
# sum(updegs %in% dwndegs)

degsdf %>% distinct(ID) %>% pull() -> degs

rownames(moduleTraitCor) <- str_replace_all(rownames(moduleTraitCor), '^ME', '')


# ---- moduleTraitRelationship

mtd %>% distinct(stage) %>% pull() -> stages
mtd %>% distinct(Tissue) %>% pull() -> tissues

data.frame(A = as.numeric(mtd %>% pull(stage) %in% stages[1]),
           B = as.numeric(mtd %>% pull(stage) %in% stages[2]),
           C = as.numeric(mtd %>% pull(stage) %in% stages[3]),
           D = as.numeric(mtd %>% pull(Tissue) %in% tissues[1]),
           E = as.numeric(mtd %>% pull(Tissue) %in% tissues[2]),
           G = as.numeric(mtd %>% pull(Tissue) %in% tissues[3])) -> datTraits

names(datTraits) <- c(stages, tissues)
rownames(datTraits) <- mtd$id

# 

moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, ncol(datTraits))

moduleTraitCor %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'moduleTraitCor') -> df1
moduleTraitPvalue %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

df1 %>%
  # filter(name %in% c('HC', 'Ctrl')) %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
                       ifelse(corPvalueStudent <.01, "**",
                              ifelse(corPvalueStudent <.05, "*", "")))) -> df1

df1 %>%
  # filter(module %in% psME & name != 'Time') %>%
  # mutate(moduleTraitCor = ifelse(abs(moduleTraitCor) < 0.5, NA, moduleTraitCor)) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  geom_tile(color = 'black', size = 0.5) + 
  geom_text(aes(label = star),  vjust = 0.75, hjust = 0.5, size= 6) +
  ggsci::scale_fill_gsea(name = "Membership", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  labs(x = '', y = 'Module') +
  guides(fill = guide_colorbar(barwidth = unit(3, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top") -> p1

#

moduleColors %>% 
  as_tibble(., rownames = 'transcript') %>%
  # mutate(degs = ifelse(transcript %in% updegs, 'up', 
                       # ifelse(transcript %in% dwndegs, 'down', ''))) %>%
  # dplyr::rename('module' = 'value') -> moduleColors

moduleColors %>% 
  as_tibble(., rownames = 'transcript') %>%
  # filter(degs != 'ns') %>% 
  group_by(module) %>% count(sort = T) -> ModuleDF

ModuleDF %>% mutate(module = factor(module, levels = hclust$labels[hclust$order])) -> ModuleDF

ModuleDF %>% ungroup() %>% mutate(pct = n / sum(n)) -> ModuleDF

ModuleDF %>% 
  # mutate(degs = ifelse(degs == 'ns', '', degs)) %>%
  ggplot(aes(x = module, y = n)) + # , fill = degs
  labs(y = 'Number of transcripts') +
  geom_col() + coord_flip() +
  # geom_col(color = 'black', size = 0.25) + coord_flip() +
  scale_fill_manual(name = '', values = c("grey30", "#EE4141", "#2428D0")) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
        axis.title.y = element_blank(), axis.text.y= element_blank(),
        axis.ticks.y=element_blank(), axis.line.y = element_blank()) -> p2

library(patchwork)

p1 + p2 + plot_layout(widths = c(0.5, 1)) +
  labs(caption = '* corPvalueStudent < 0.05 ') -> psave

ggsave(psave, filename = 'ModuleTraitRelationship.png', path = path, width = 8, height = 7)


