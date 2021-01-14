# start beepr::beep(3)

# # # # # #
# Clean memory ----
# # # # # #

rm(list = ls())

# # # # # #
# Set path and files ----
# # # # # #

dir <- '~/Downloads/archivosarticulosloberas/'
files <- list.files(path = dir, pattern = ".csv", full.names = T)

# files es un ejemplo de un vector de tipo caracteres typeof(files)

obj <- read.csv(files[1])

typeof(obj)
class(obj)

# Trabajando con objetos

# # # # # # # # # # #
# Explorar dimensiones ----
# # # # # # # # # # #

dim(obj)
nrow(obj)
ncol(obj)

names(obj) # tambien se puede usar colnames

# # # # # # # # # # #
# Explorar estructura ----
# # # # # # # # # # #

head(obj)
tail(obj)
str(obj)
# View(obj)

# # # # # #
# ordenar tabla de datos ----
# # # # # #

# start beepr::beep(8)

# # # # # # # # # #
# Cargamos programas ----
# # # # # # # # # #

library(tidyr)
library(dplyr)
library(ggplot2)

# # # # # # #
# uso de tuberias
# # # # # # #

obj %>% head()
obj %>% tibble::glimpse()

# re-organizando datos:

obj <- obj %>% as_tibble()

# Metodos de seleccion de columnas

obj %>% 
  select(index, Samples_Name,  Sampling_Area, Replicate) %>% 
  View()

obj %>% 
  select_if(is.double)  %>% 
  View()

obj %>% 
  select_at(vars(contains("D_0"))) %>%
  View()

# 

library(stringr)

colNames <- obj %>% select_at(vars(contains("D_0"))) %>% names()

obj %>% 
  pivot_longer(cols = colNames, 
               names_to = "linaje", 
               values_to = "ab") -> obj_longer 

obj_longer %>% View()

# datavis
# 1.
obj_longer %>%
  ggplot(aes(x = linaje, y = ab)) +
  geom_col() +
  coord_flip()

# 2.
obj_longer %>%
  ggplot(aes(x = linaje, y = ab, fill = Sampling_Area)) +
  geom_col() +
  coord_flip()

# Expresiones regulares
ranks <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family')

obj_longer %>% 
  separate(col = linaje, into = ranks, sep = '[.]') %>% 
  mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__"=""))))  -> obj_longer

obj_longer %>%
  separate(col = Replicate, into = c('Sample','Rep'), sep = '_') -> obj_longer
  
obj_longer %>% tibble::glimpse()

# dataVis

obj_longer %>%
  ggplot(aes(x = Sampling_Are, y = ab, fill = Phylum)) +
  coord_flip() +
  geom_bar(stat = "identity", 
           position = "stack", color = "black")

raT <- function(x) {(x / sum (x) ) * 100}

obj_longer %>%
  group_by(index) %>%
  mutate(RA = raT(ab)) %>%
  filter(Rep %in% "1") %>%
  ggplot(aes(x = Sampling_Area, y = RA, fill = Phylum)) +
  coord_flip() +
  geom_bar(stat = "identity", 
           position = "stack", color = "black") -> plot

# Aesthetic

plot + ggsci::scale_fill_aaas()
plot + ggsci::scale_fill_cosmic()
plot + ggsci::scale_fill_rickandmorty()

plot + ggsci::scale_fill_rickandmorty() +
  theme_classic(base_size = 14)

# facet

obj_longer %>%
  group_by(index) %>%
  mutate(RA = raT(ab)) %>%
  ggplot(aes(x = index, y = RA, fill = Phylum)) +
  coord_flip() + labs(x = 'Replicate', y = 'Rel. Ab (%)') +
  geom_bar(stat = "identity", 
           position = "stack", color = "black") +
  facet_grid( Sampling_Area ~., scales = "free_y") +
  ggsci::scale_fill_rickandmorty() +
  theme_classic(base_size = 14)

# create my own method

my_tidy_function <- function(file) {
  
  obj <- read.csv(file)
  
  colNames <- obj %>% select_at(vars(contains("D_0"))) %>% names()
  
  obj %>% pivot_longer(cols = colNames, 
                 names_to = "linaje", 
                 values_to = "ab") %>%
    mutate(fileName = basename(file)) %>%
    select(fileName, Replicate, index, linaje, ab) %>%
    mutate(fileName = str_replace_all(fileName, "_level-5.csv", ""))
  
  
  
}

my_tidy_function(files[1])

# Concatenate tidy 

obj_longer <- lapply(files, my_tidy_function) %>% do.call(rbind, .)

obj_longer %>% 
  separate(col = linaje, into = ranks, sep = '[.]') %>% 
  mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__"="")))) %>%
  separate(col = Replicate, into = c('Sample','Rep'), sep = '_') %>%
  mutate(Amplicon = str_replace_all(fileName, "_level-5.csv", "")) -> obj_longer

obj_longer %>% tibble::glimpse()


obj_longer %>% View()

obj_longer %>%
  ggplot(aes(x = Sample, y = RA, fill = Phylum)) +
  coord_flip() + labs(x = 'Amplicon', y = 'Rel. Ab (%)') +
  geom_bar(stat = "identity", 
           position = "stack", color = "black") +
  facet_grid(Amplicon ~ Rep, scales = "free_y") +
  ggsci::scale_fill_rickandmorty() +
  theme_classic(base_size = 14)

#

# elaboracion de heatmap con facet de amplicones


obj_longer %>%
  filter(Rep %in% "1") %>%
  ggplot(aes(y = Family, x = Sample, fill = RA)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative\nAbundance") +
  guides(fill = guide_colorbar(barheight = unit(7, "cm"),  
                               ticks.colour = "black", 
                               frame.colour = "black")) +
  facet_grid(Phylum+Class ~ Amplicon, scales="free", space = "free") +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1),
    strip.text.y = element_text(
      angle = 0, 
      size = 5),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.1),
    panel.spacing = unit(0.001, "lines")
  )

# grouped barplot

obj_longer %>%
  ggplot(aes(y = ab, x = Sample)) +
  geom_bar(aes(fill = Amplicon), 
           position = "dodge", stat = "identity") +
  labs(y = "N reads") +
  facet_wrap(~ Phylum, scales = "free_y") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1))

# Probar bajo tu propio riesgo!!!

# Charge the circlize library
library(circlize)
# https://jokergoo.github.io/circlize_book/book/advanced-usage-of-chorddiagram.html


set.seed(100121)

grid.col <- rep("grey", 6)
names(grid.col) <- unique(obj_longer$Amplicon)

# establecer paleta de colores de filos

g <- ggplot_build(plot)
col_palette <- c(unique(g$data[[1]]["fill"]))$fill
names(col_palette) <- unique(g$plot$data$Phylum)

circos.clear()

obj_longer %>% 
  filter(ab > 0) %>%
  select(Amplicon, Phylum) %>%
  with(., table(Phylum, Amplicon)) %>%
  chordDiagram(grid.col = c(grid.col, col_palette),
               annotationTrack = "grid", 
               preAllocateTracks = 1,
               small.gap = 5)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
              CELL_META$sector.index, 
              facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

abline(h = 0.05, lty = 2, col = "#00000080")

#treemap

library(treemap)


# Make input for an alluvial plot of amplicons

obj_longer %>%
  group_by(Family, Sampling_Area) %>%
  filter(ab > 0) %>%
  summarise(n = length(Family), N = sum(ab)) %>%
  # inner_join(obj_longer %>% distinct(Family, .keep_all = T)) %>%
  # select_at(vars(ranks, n, N)) %>%
  treemap(index = c("Sampling_Area","Family"),
          vSize = "n", type="index",
          border.col=c("black","white"),
          palette = "Set1", 
          title="",
          fontsize.labels = 10)
# input for alluvial plot

library(ggalluvial)

obj_longer %>%
  filter(Rep %in% "1") %>%
  filter(Sample %in% "Cantiles") %>%
  # group_by(Amplicon, Phylum) %>% # Sanity check w/o Phylum grouping
  filter(ab > 0) %>%
  count(Amplicon, Phylum) %>%
  mutate_if(is.character, as.factor) %>%
  tibble::rowid_to_column() -> alluvInput

table(alluvInput$Amplicon)

alluvInput %>%
  ggplot(aes(x = Amplicon, stratum = Phylum)) +
  geom_stratum()





# 

dir <- '~/Downloads/archivosarticulosloberas/'
files <- list.files(path = dir, pattern = ".xlsx", full.names = T)

obj <- readxl::read_xlsx(files[1])

obj %>% 
  select_if(is.double) %>% names() -> colNames

obj %>% 
  pivot_longer(cols = colNames, values_to = "RA", names_to = "Index") %>%
  separate(Index, into = c("Sample", "index", "Rep"), sep = "-") %>%
  select(-pathway, -SubPathway, -index) -> obj_longer

obj_longer %>%
  filter(Rep %in% "1") %>%
  filter(RA > 1) %>%
  ggplot(aes(x = SuperPathway, y = RA, fill = Category)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ Sample) +
  theme_bw()

# arrange by factors and facets

# https://trinkerrstuff.wordpress.com/2016/12/23/ordering-categories-within-ggplot2-facets/


obj_longer %>%
  filter(Rep %in% "1") %>%
  group_by(Sample) %>%
  ggplot(aes(y = reorder_within(SuperPathway, RA, Sample), 
             x = RA, fill = Category)) +
  geom_col() +
  facet_wrap(~ Sample, scales = "free_y")


library(babynames)
top_names <- babynames %>% 
  filter(year >= 1950,
         year < 1990) %>% 
  mutate(decade = (year %/% 10) * 10) %>% 
  group_by(decade) %>%
  count(name, wt = n) %>% 
  group_by(decade) %>% 
  top_n(10)

top_names %>% 
  ggplot(aes(y = reorder_within(name, n, decade), 
             x = n)) + 
  geom_col(show.legend = FALSE) + 
  facet_wrap(~ decade, scales = "free_y") +
  labs(x = NULL, y = NULL)
