
library(tidyverse)

malacos <- read.csv('/Users/cigom/metagenomics/db/bold/BOLD_public_trim/malacostraca.csv', sep = ' ', header = F, stringsAsFactors = F)

tax.split <- strsplit(malacos[, ncol(malacos)], ";")

# Using the max rank assignation to names the taxonomy object
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
taxonomy <- as.data.frame(t(taxonomy))
tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)

tax <-data.frame(row.names = rownames(tax),
                 mutate_all(data.frame(tax), 
                            funs(str_replace_all(., c("_"=" ", "NA"=NA)))))


tax %>%
  as_tibble() %>%
  drop_na() %>%
  group_by(V6) %>%
  summarise(n = length(unique(V9)), divergence = 1-(1/length(unique(V9)))) %>%
  arrange(desc(n)) %>%
  mutate(pct = n / sum(n), 
         group = 'Malacostraca') %>%
  mutate(facet = ifelse(n > 499, "A", "B")) %>%
  ggplot(aes(divergence, n)) +
  geom_point(size = 1) +
  ggrepel::geom_text_repel(aes(label = V6, 
                               color = ifelse(V6 %in% "Euphausiacea", 'red','black')), size = 4) +
  facet_grid(facet~., scales = "free_y",space = "fixed") +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = "Divergence", y = '# Species') +
  
  # mutate(Order = ifelse(pct > 0.001, V6, "Low")) %>%
  # mutate(V6 = factor(V6, levels = unique(V6))) %>%
  # ggplot(aes(x = group, y = pct, fill = V6)) +
  # geom_col() + coord_flip() +
  # ggsci::scale_fill_d3() +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

# filter, then,
table(tax[,6])

eu <- tax %>%
  filter(V7 %in% 'Euphausiidae') %>%
  filter_at(.vars = vars(V9), all_vars(!grepl('sp',.))) %>%
  drop_na()

eu <- eu[, 7:9]

library(ggalluvial)

eu %>%
  as_tibble() %>%
  drop_na() %>%
  group_by(V8) %>%
  summarise(n = length(unique(V9))) %>%
  mutate(nn = n/21, group = 'Euphausiidae') %>%
  arrange(n) %>%
  mutate(cumsum = cumsum(n)) %>%
  # mutate(V8 = factor(V8, levels = unique(V8))) %>%
  ggplot(aes(x = nn, y = n, size = cumsum)) +
  geom_point() +
  ggrepel::geom_text_repel(arrow = NULL, aes(label = V8, 
                               color = ifelse(V8 %in% "Euphausia", 
                                              'red','black')), size = 4) +
  scale_color_manual(values = c('black', 'red')) +
  labs(x = "Fraction", y = '# Species') +
  theme_bw(base_size = 18) +
  theme(legend.position = "none")

rest <- tax[,6:7] %>%
  filter_at(.vars = vars(V6), all_vars(!grepl('sp',.))) %>% 
  drop_na()

library(tidyr)

data0 <- data.frame(unite(rest, sep = ";", remove = TRUE, col = 'Taxonomy'))


data <- data.frame(unite(eu, sep = ";", remove = TRUE, col = 'Taxonomy'))

library(taxa)


obj <- parse_tax_data(data, #rbind(data,data0),
                      class_cols = "Taxonomy",
                      class_sep = ";")
library(metacoder)

set.seed(1) # This makes the plot appear the same each time it is run

obj %>%
  # filter_taxa(supertaxa = TRUE) %>%
  heat_tree(., 
            node_label = taxon_names,
            node_label_size = log10(n_obs+1), 
            # node_size_range = c(0.01, 0.04),
            edge_size_range = c(0.005, 0.005),
            # node_size = n_obs / nrow(data),
            tree_label = taxon_names,
            node_color = n_obs,
            node_size_axis_label = "Index",
            node_color_axis_label = "Counts",
            layout = "da",
            initial_layout = "re")

# map of distribution
library("rnaturalearthdata")
library('rnaturalearth')
library('ggspatial')
library('sf')
library(tidyverse)

# lme <- sf::read_sf("~/Downloads/MarineRealmsShapeFile/MarineRealms.shp")
lme <- sf::read_sf("~/Downloads/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp")


lme <- lme %>% 
  # st_transform(., 54032) %>%
  sf::st_simplify(dTolerance = 0.01)

world <- ne_countries(scale = "medium", returnclass = "sf")

gg <- ggplot(data = world) +
  # geom_sf(fill = 'antiquewhite') +
  geom_sf(size = .2, fill = "gray80", col = "gray90")
# annotation_scale(location = "bl", width_hint = 0.5)
# annotation_north_arrow(location = "bl", which_north = "true", 
#                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
#                        style = north_arrow_fancy_orienteering)
# coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97))


unique(lme$TYPE)


coords <- read.csv('/Users/cigom/Downloads/bold_data.txt', 
                   sep = '\t', stringsAsFactors = F)

coords %>% 
  as_tibble() %>%
  select(genus_name,species_name, lon, lat) %>%
  mutate(genus_name = ifelse(genus_name == '', NA, genus_name)) %>%
  # group_by(species_name) %>%
  drop_na() -> coods
coods %>%
  filter(species_name %in% "Euphausia superba") -> coods
  # group_by(species_name)

library(RColorBrewer)

colourCount <- length(unique(coods$genus_name))
getPalette = colorRampPalette(brewer.pal(11, "BrBG"))
fvalues = getPalette(colourCount)

pal <- setNames(fvalues, unique(coods$genus_name))

coli.fun = colorRampPalette(brewer.pal(11, "Paired"))



gg +  
  geom_sf(data = lme, aes(fill = REALM), size = .2, col = 0, alpha=0.7)+
  ggtitle(paste("LME - Large Marine Ecosystems of the World - ", length(unique(lme$REALM)),"ecosystems")) +
  geom_point(data = coods,
             aes(x = lon, y = lat), color = "red",
             alpha = 0.9, shape = 20, size = 5, fill = 'NA') +
  scale_color_manual(values = pal) +
  scale_fill_manual(values=coli.fun(length(unique(lme$REALM)))) +
  # geom_sf_text(data = lme, aes(label = REALM), colour = "grey40", check_overlap =TRUE)+
  coord_sf(expand = FALSE) +
  theme(legend.position="top")

ggsave("LME_EUFASIDES2.png", dpi = 300,
         units = 'in', width = 17, height = 11)

gg +  
  geom_point(data = coods,
             aes(x = lon, y = lat, color = genus_name), 
             alpha = 0.9, shape = 20, size = 5, fill = 'NA') +
  coord_sf(expand = FALSE) +
  theme(legend.position="top") +
  scale_color_manual(values = pal)

ggsave("LME_EUFASIDES.png", dpi = 300,
       units = 'in', width = 17, height = 11)
# test https://semba-blog.netlify.app/12/01/2019/plotting-the-spatial-distribution-of-chlorophyll-from-modis/
# https://rstudio-pubs-static.s3.amazonaws.com/485396_7d5f60e87225469fb0c0c04684a0cf31.html

#  -----
.cran_packages <- c("dplyr", "tidyverse")
.bioc_packages <- c("Biostrings", "DECIPHER", "phangorn", 'dada2')
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

coords %>% 
  as_tibble() %>%
  select(genus_name,species_name, nucleotides) %>%
  mutate(nucleotides = str_replace_all(nucleotides, '-', '')) %>%
  distinct(nucleotides, .keep_all = T) -> tbl
tbl %>% 
  pull(nucleotides) -> seqs

# seqs <- DNAStringSet(seqs)
# seqs <- DECIPHER::RemoveGaps(seqs)
# perform the alignment
names(seqs) <- seqs

alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)

treeNJ <- NJ(dm) # Note, tip order != sequence order

fit = pml(treeNJ, data=phangAlign)

fitGTR <- update(fit, k=4, inv=0.2)

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

tree <- fitGTR$tree

library(dplyr)
library(ggtree)
library(RColorBrewer)

x <- as_data_frame(fitGTR$tree)
y <- full_join(x, tbl, by = 'label')
tree <- NULL
tree <- y %>% tidytree::as.treedata()



colourCount = length(unique(tbl$Orden))
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))

ggtree(tree, branch.length='none', layout='slanted', aes(color=Orden), na.rm = TRUE, size=1) +
  theme(legend.position="top") +
  scale_color_manual(values = c(getPalette(colourCount)), na.value = "grey", guide = guide_legend(ncol=3)) +
  #geom_tiplab(size=4, aes(label=paste0('italic(', label, ')~bolditalic(', Query.ID, ')~', Best.ID)), parse=FALSE)
  xlim(0, 10) +
  geom_tiplab(size=4, aes(label=Query.ID), offset=0.5) +
  geom_nodepoint(alpha=1/4, size=4)