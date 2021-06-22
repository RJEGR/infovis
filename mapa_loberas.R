# mapa con coordenadas y pieplot de 
# barplot sin incluir low_abundance
# heatmap con generos a la derecha y especies (patogenas) a la izq. incluir todos los marcadores
# los porcentajes de abundancia de los 4 filos mas abundantes (para incluir en el texto del paper)

library("rnaturalearthdata")
library('rnaturalearth')
library('ggspatial')
library('sf')
library(tidyverse)
library(geosphere)
library(maps)
library(ggplot2)
library(ggsci)
library(phyloseq)

# Barplot ----

# prepare_ps(fileNames[1], agg = T, agg_level = 'Phylum') %>% 
#   otu_table(.) %>% 
#   cbind() %>% as_tibble(rownames = 'Level') %>%
#   pivot_longer(cols = sampleName, names_to = 'site') %>%
#   left_join(df) -> dfLong



# Mapa ----
# Rookery name	Geographic coordinates
# Coloradito	30°02’ N - 114°29’ W
# Granito	29°33’N - 113°32’ W
# Cantiles	29°30’N - 113°27’ W
# Machos	29°17’N - 113°29’W
# Partido	28°54’N - 113°02’ W
# Rasito	28°50’N - 112°59’W

site <- c('Coloradito', 'Granito', 'Cantiles', 'Machos', 'Partido', 'Rasito')
lat <- c(30.02, 29.33, 29.30, 29.17, 28.54, 28.50)
long <- -c(114.29, 113.32, 113.27, 113.29, 113.02, 112.59)

df <- data.frame(long, lat, site)

prepare_ps(fileNames[1], agg = T, agg_level = 'Phylum') %>%
  otu_table(.) %>%
  cbind() %>% as_tibble(rownames = 'Level') %>%
  select_at(vars(ranks)) %>%
  pivot_longer(cols = sampleName, names_to = 'site') %>%
  pivot_wider(names_from = 'Level', values_from = value) -> dfLong

apply(dfLong[-1], 1, function(x) sum(x > 0)) -> ratio

ratio <- ratio/length(taxNames)

dfLong %>% select_if(is.double) %>% names(.) -> taxNames

colourCount <- length(taxNames)

library(ggsci)
library(RColorBrewer)
if(colourCount > 7) {
  # getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
  getPalette <- colorRampPalette(brewer.pal(12, 'Paired'))(colourCount)
} else
  # getPalette <- pal_locuszoom()(colourCount)
  getPalette <- brewer.pal(colourCount, 'Set1')

getPalette = structure(getPalette, names = taxNames)

# test extended color palette

agg_level = 'Order'
color_by = 'Phylum'

# borrar facet)Bar
facet_bar <- function(featureDF, agglom_lev = ranks[2], low_ab = 1, color_by = ranks[2]) {
  
  ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species')
  
  featureDF %>%
    pivot_longer(all_of(sampleName), values_to = 'ab', names_to = 'Sample') %>%
    filter(ab > 0) %>%
    group_by(Sample, group) %>%
    drop_na(agglom_lev) %>%
    mutate(RA = (ab / sum(ab)) * 100) %>% # summarise(sum(RA))
    filter(RA > low_ab) %>%
    mutate(RA = (ab / sum(ab)) * 100) -> dataViz
  
  labels <- dataViz %>% pull(agglom_lev) %>% unique() %>% sort()
  colourCount = length(labels)
  
  library(ggsci)
  
  if(identical(agg_level, color_by)) {
    if(colourCount > 7) {
      getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
    } else {
      getPalette <- pal_locuszoom()(colourCount)
      getPalette = structure(getPalette, names = labels)
    }
  } else {
    dataViz %>% 
      ungroup() %>% 
      distinct_at(vars(c(agg_level, color_by))) %>%
      arrange_at(color_by) -> labelsdf 
    
    labelsdf %>%
      group_by_at(color_by) %>% tally() %>% pull(n) -> n
      if(length(n) > 7) {
        getPalette <- colorRampPalette(pal_locuszoom()(7))(length(n))
        } else
          getPalette <- pal_locuszoom()(length(n))
        
        coul <- list()
        
        for(i in 1:length(getPalette)) {
          j <- i
          c <- colorRampPalette(c(getPalette[j], 'white'), alpha = F)(n[j]+1)
          c <- c[1:n[j]]
          coul[[j]] <- c
          }
    
    coul <- unlist(coul)
    labels <- labelsdf %>% pull(agglom_lev) %>% unique()
    coul = structure(coul, names = labels)
    
  }
  
  
  
  dataViz %>% group_by(group) %>%
    summarise(ab = sum(ab)) %>% 
    pivot_wider(names_from = 'group', values_from = ab) -> pout
  
  dataViz %>%
    arrange_at(color_by) %>%
    mutate_at(vars(agglom_lev), as.factor) %>%
    ggplot(aes_string(y = 'RA',x = 'Sample', fill = agglom_lev)) +
    geom_col() +
    labs(x = "", y = "Relative Abundance (%)", 
         caption = paste0("Low Abundance (Relative Ab <", low_ab, " %)")) +
    scale_fill_manual(agglom_lev, labels = labels, values = coul) +
    theme_classic(base_size = 18, base_family = "GillSans") +
    facet_grid(~group) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.background = element_blank(),
          panel.border = element_blank()) -> p
  
  
  print(pout)
  return(p)
  
}



# now, replace the color fill scale from the facet_bar function! 

facet_bar(featureDF, agglom_lev = agglom_lev, low_ab = 1, color_by = color_by)
facet_bar(featureDF, agglom_lev = color_by, low_ab = 1, color_by = color_by)

# tax %>% cbind(.,  col = unlist(coul))

pie(rep(1, length(coul)), col = coul , main="", labels = names(coul)) 

dfLong %>% right_join(df) %>% 
  select(names(df), taxNames) %>%
  cbind(., ratio) -> dfLong

# ylim <- c(min(lat)-2, max(lat)+2)
# xlim <- -c(min(lat)-2, max(lat)+2)

xlim=c(-115, -110);ylim = c(27, 32) # <- estas coords, funcionan bien para golfo de cali

world <- ne_countries(scale = "medium", returnclass = "sf")

gg <- ggplot(data = world) +
  geom_sf(fill = 'antiquewhite') +
  coord_sf(xlim = xlim, ylim = ylim) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)

# gg

library(scatterpie)


gg +
  geom_scatterpie(aes(x = long, y = lat, group = site), #  r = ratio
  cols = taxNames, data = dfLong, colour = NA, pie_scale = 2) + 
  scale_fill_manual(values = getPalette)

  # geom_scatterpie_legend(dfLong$ratio, x=-112, y=29)
  # geom_label(data = dfLong, aes(x = long, y = lat, 
  #                               label = site),  size = 2.5, fontface = "bold")
