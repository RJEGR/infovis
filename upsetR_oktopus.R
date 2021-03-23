
library(tidyverse)

path <- "~/transcriptomics/oktopus_full_assembly/"
countf <- list.files(path = path, pattern = "counts_table_length_ajus_gen_level-aproach2-filtered_mean_reps.txt", full.names = T)

degs <- readLines(paste0(path, "/degs.lists"))
dff <- readRDS(file = paste0(path, 'upgenes.rds'))

dim(count0 <- read.delim(countf, sep = "\t"))

sum(rownames(count0) %in% degs)

dim(count <- count0[rownames(count0) %in% degs, ])
dim(count <- count0[rowSums(edgeR::cpm(count0) > 3 ) >= 2, ] )
colNames <- names(count)

apply(count, 2, function(x) ifelse(x > 0, 1, 0)) %>% as.data.frame() -> booleandf

table(rowSums(booleandf)) -> prevalence

apply(count0, 2, function(x) ifelse(x > 0, 1, 0)) %>% 
  as.data.frame() -> booleandf

table(rowSums(booleandf))


dim(booleandf <- booleandf[rowSums(booleandf) != ncol(booleandf),])

read.delim(paste0(path, "metadata.tsv"), sep = "\t")  %>% 
  distinct(group, .keep_all = T) %>% group_by(Tissue) %>%
  arrange(Tissue, match(stage, c("PRE", "Spawing", "POST"))) -> meta

# identical(names(booleandf), meta$group)
# 
# sets[match(meta$group, names(booleandf))]
# booleandf <- booleandf[, match(meta$group, names(booleandf))]
# 
# identical(names(booleandf), meta$group)
# 
# names(booleandf) <- paste0(substr(meta$group, 1, 3), '-',meta$stage)

sets <- names(booleandf)

labels <- meta$Tissue

n <- length(unique(labels))

grid.col <- ggsci::pal_rickandmorty()(n)

names(grid.col) <- unique(labels)

sets.bar.color <- grid.col[match(labels, names(grid.col))]

names(sets.bar.color) <- sets


# Boolean matrix for all intersections of degs genes, sorted by size. Dark circles in the matrix indicate sets that are part of the intersection


library(UpSetR)

png(filename = paste0(path, "upset_degs.png"), width = 780, height = 780, res = 150)

upset(booleandf, number.angles = 45, point.size = 2.3, line.size = 1.5, 
      order.by = "freq", nsets = 9, nintersects = 30, 
      sets = sets, sets.bar.color = sets.bar.color, scale.intersections = 'identity',
      mainbar.y.label = "Gene Intersections", sets.x.label = "Genes per sample",
      # attribute.plots=list(gridrows=60, plots=list(list(plot=histogram, x="GOV_H_PR_AD_24"))),
      text.scale = c(1.3, 1, 1.3, 1, 1, 0.7))

dev.off()

# en base a que criterio podriamos definir que los genes se estan co-expresando?
# si queremos comparar relevancia por tejido, quiza convenga conservar aquellos genes que prevalecen en todas las muestras, 
# prevalence of genes::

Tcount <- edgeR::cpm(count0)
# Tcount <- log2(count0+1)
prevelancedf = apply(Tcount, 1, function(x) sum(x > 0))
mean_se = apply(Tcount, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
                TotalAbundance = rowSums(Tcount),
                mean_se) %>% 
  as_tibble(rownames = "id") %>%
  arrange(desc(TotalAbundance)) -> df

df %>% 
  # filter(id %in% degs) %>%
  group_by(Prevalence) %>% summarise(mean_se(TotalAbundance)) %>%
  ggplot(aes(as.character(Prevalence), y)) + 
  geom_boxplot() + labs(y = 'mean(TotalAbundance)') +
  geom_errorbar(aes(ymin = ymin, 
                    ymax = ymax), width = 0.2,
                position = position_dodge(0.05), color = "black")

df %>% 
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  # filter(id %in% degs) %>%
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>%
  # ggpubr::ggqqplot(., x = "TotalAbundance" ,facet.by = 'Prevalence')
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- df %>%  filter(id %in% degs) %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n), Prevalence = paste0(1:9, ' Samples'))

p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " genes")),
  hjust   = -1,
  vjust   = -2
) + theme_classic(base_size = 16, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1))) -> p1

ggsave(p1, filename = "Prevalence.png", path = path, 
       width = 8, height = 6)

df %>% 
  filter(id %in% degs) %>% 
  # filter(Prevalence == 1) %>%
  left_join(dff, by = c('id'='ID')) %>%
  ggplot(aes(logFC, PValue, 
             color = Tissue, size = log2(TotalAbundance+1))) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~as.character(Prevalence))

# lo anterior es un claro ejemplo de una sobre-estimacion (inflacion) de valores P. el teste de degs tiene por objetivo comparar si la expresion del gen1 es diferencial entre la muestra A y B, el supuesto es que, a pesar de una abundancia en ambas muestras, es necesrio identificar si una abundancia u otra es significativamente diferencial en x mustra. Pero si la abundancia es zero, no es necesario aplicar ninguna prueba estadistica, tal como la de degs.
# verifiquemos la prevalencia
# podemos establecer una matriz con los genes que prevalencen a lo largo de todos las muestras (n = 9), estos genes (n = 20826) podrian catalogarse como basales para el funcionamiento de los tejidos del pulpo. Verifiquemos esto:

df %>% 
  filter(Prevalence < 4) %>% 
  arrange(desc(TotalAbundance)) %>%
  pull(id) -> prevGenes

dim(prev <- Tcount[rownames(Tcount) %in% prevGenes, ])
colSums(prev)
names(prev) <- substr(names(prev), 1,8)

hclust <- hclust(dist(prev), "complete")
sclust <- hclust(dist(t(prev)), "complete")

apply(prev, 2, function(x) ifelse(x == 0, NA, x)) %>% as.data.frame() -> prev

superheat::superheat(t(prev), 
                     order.cols = hclust$order, 
                     order.rows = sclust$order,
                     # membership.rows = meta$stage,
                     heat.na.col = "white",
                     bottom.label.text.angle = 45,
                     bottom.label.text.size = 4)

colNames <- colnames(Tcount)

Tcount %>%
  as_tibble(rownames = 'gene') %>%
  pivot_longer(cols = colNames, names_to = 'group') %>%
  filter(value > 0) %>% 
  left_join(df %>% select(Prevalence, id), by = c('gene'='id')) %>%
  group_by(group,Prevalence) %>%
  summarise(value = sum(value)) %>% arrange(desc(value)) %>%
  left_join(mtd %>% distinct(group, .keep_all = T), by = 'group') %>%
  mutate(stage = factor(stage, levels = c("PRE", "Spawing", "POST"))) %>%
  mutate(group = forcats::fct_reorder(group, as.numeric(stage))) %>%
  # filter(Prevalence > 8) %>% group_by(group) %>% summarise(sum(value)) %>% arrange() # min(value) to f
  # group_by(group) %>% summarise(sum(value)) # sanity check
  # mutate(value = ifelse(value >= 1, value, NA)) %>% 
  mutate(wrap = ifelse(Prevalence > 8, 'B', 'A')) %>%
  ggplot(aes(x = group, y = value)) +
  geom_col(aes(fill = as.character(Prevalence))) + coord_flip() +
  geom_abline(slope = 0, intercept = 829365, linetype="dashed", alpha=0.5) +
  ggh4x::facet_nested(Tissue+stage ~ wrap, scales = 'free')
  # facet_grid(Prevalence ~., scales = 'free_y')

# so, select the prevalence of 9 to network?

levs <- unique(unlist(dt, use.names = FALSE))
table(lapply(dt, factor, levs))

library(circlize)

getPalette <- RColorBrewer::brewer.pal(3, "Set1")
getPalette <- rep(getPalette, 1, each = 3)
names(getPalette) <- colNames

stages <- mtd %>% distinct(group, .keep_all = T)
group <- factor(stages$stage, levels = c("PRE", "Spawing", "POST"))
names(group) <- colNames

# sets.bar.color <- grid.col[match(labels, names(grid.col))]
circos.clear()
booleandf %>%  
  as_tibble(rownames = 'from') %>%
  pivot_longer(cols = colNames, names_to = 'to', values_to = 'pa') %>%
  filter(pa > 0) %>% 
# to make sample-to-sample chord ----
  group_by(from) %>%
  mutate(from = replace(from, duplicated(from), NA)) %>%
  mutate(from = ifelse(is.na(from), to, from)) %>%
  mutate(to = ifelse(grepl('^TR', from), to, NA)) %>%
  ungroup() %>% fill(to) %>% filter(!grepl('^TR', from)) %>%
  with(., table(from, to)) %>% 
  chordDiagram(grid.col = getPalette,
               directional = -1, 
               # group = group,
               diffHeight = mm_h(5), target.prop.height = mm_h(4))
# -----
  # left_join(df %>% select(Prevalence, id, TotalAbundance), by = c('from'='id')) %>%
  # left_join(mtd %>% distinct(group, .keep_all = T), by = c('to'='group')) %>%
  # filter(Prevalence < 2) %>% select(Prevalence, to, TotalAbundance) %>% chordDiagram()
  # mutate(from = 'transcripts') %>% filter(Prevalence < 5)  %>% 
  # with(., table(Prevalence, to)) %>% chordDiagram()
  # select(Prevalence, to, pa) %>% chordDiagram()
circos.clear()
circos.par(start.degree = 0, gap.degree = 4, 
           track.margin = c(-0.01, 0.01), 
           points.overflow.warning = FALSE)
