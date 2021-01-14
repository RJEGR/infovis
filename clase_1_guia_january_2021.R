# # # # # #
# Clean memory ----
# # # # # #

# para conocer mas informacion de las funciones se usa el signo ?
# ej:

?ls
?rm

rm(list = ls())

# # # # # #
# Set path and files ----
# # # # # #

dir <- "~/Downloads/archivosarticulosloberas/"
files <- list.files(path = dir, pattern = ".csv", full.names = T)

# 

obj <- read.csv(files[1])

# # # # # # # # # # #
# Explorar dimensiones ----
# # # # # # # # # # #

dim(obj)
nrow(obj)
ncol(obj)

names(obj)

# # # # # # # # # # #
# Explorar estructura ----
# # # # # # # # # # #

head(obj)
tail(obj)
str(obj)
View(obj)

# # # # # #
# ordenar tabla de datos ----
# # # # # #

library(tidyr)
library(dplyr)


# # # # # # #
# uso de tuberias
# # # # # # #

obj %>% head()

# re-organizando datos:

obj <- obj %>% as_tibble()

# Metodos de seleccion de columnas
# consultas de tablas

obj %>% 
  select(index, Samples_Name,  Sampling_Area, Replicate) %>% 
  View()

# is.double(10)

obj %>% select_if(is.double) %>% select(-Sampling_Date)

colNames <- obj %>% select_at(vars(contains("D_0"))) %>% names()

str(colNames)

obj_longer <- obj %>% pivot_longer(cols = colNames, names_to = "linaje", values_to = 'ab') 

library(ggplot2)

obj_longer %>%
  filter(ab > 100) %>% 
  ggplot(aes(x = linaje, y = ab)) + geom_col() + 
  coord_flip()

# Expresiones regulares ----

library(stringr)

string <- "D_0__Bacteria.D_1__Firmicutes.D_2__Clostridia.D_3__Clostridiales.D_4__Family.XI"

str_replace_all(string, c("D_[0-4]__" = ""))




ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family")

obj_longer %>% separate(col = linaje, into = c("Kingdom", "Phylum", "Class", "Order", "Family"), sep = "[.]") %>% mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__" = "")))) -> obj_longer

obj_longer %>%
  separate(col = Replicate, into = c("Sample", "Rep"),
           sep = "_") -> obj_longer


# Explore datavis

obj_longer %>%
  ggplot(aes(x = Sampling_Area, y = ab, fill = Phylum)) +
  geom_col() +
  coord_flip()

# tranform abundance to RA
library(tidyverse)

obj_longer %>%
  group_by(index) %>%
  mutate(RA = (ab / sum(ab)) * 100) %>%
  ggplot(aes(x = Rep, y = RA, fill = Phylum)) +
  geom_col() +
  coord_flip() +
  labs(x = "Sample_name", y = "Relative Abundance (%)",
       title = "V2 amplicon", subtitle = "Este es mi subtitulo", 
       caption = "Figura 1. Caption") +
  facet_grid(Sampling_Area ~., scales = "free_y") -> plot

plot + theme_classic(base_size = 16, base_family = "GillSans") -> plot

# Change color

plot + ggsci::scale_fill_rickandmorty() -> plot
# Facets

# 1)
plot + facet_grid(~ Sampling_Area ) # Hacen paneles como columnas (ie. Verticalmente)
# 2)
plot + facet_grid(Sampling_Area ~ . ) # Hacen paneles como filas (ie. Horizontalmente)

plot + facet_grid(Sampling_Area ~., scales = "free_y")




ggsave(plot, filename = "Figure_1.png", path = dir, height = 10, width = 7)

# reciclando nuestro trabajo

# 1. Cargue los archivos csv
# 2. arregle el orden de la tabla "cruda"
# 3. Edite la informacion taxomomica
# 4. deje una tabla lista para grafica

my_tidy_function <- function(file) {
  
  # 1. Cargue los archivos csv
  
  obj <- read.csv(file, stringsAsFactors = FALSE)
  
  # 2. arregle el orden de la tabla "cruda"
  
  colNames <- obj %>% select_at(vars(contains("D_0"))) %>% names()
  
  obj %>% pivot_longer(cols = colNames, names_to = "linaje", 
                       values_to = "ab") %>%
    mutate(fileName = basename(file)) %>%
    group_by(index) %>%
    mutate(RA = (ab / sum(ab)) * 100)
  
}


my_tidy_function(files[2]) %>% view()

lapply(files, my_tidy_function) -> obj_longer

# Necesitamos apilar nuestras tablas que se encuentran en la lista, ya que tenemos el mismo numero de columnas y nombres de columnas, lo podemos hacer con la funcion rbind.

obj_longer %>% do.call(rbind, .) -> obj_longer

obj_longer %>% view()

# to continue tomorrow

obj_longer %>%
  separate(col = linaje, into = ranks, sep = "[.]") %>% 
  mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__" = "")))) %>%
  separate(col = Replicate, into = c("Sample", "Rep"),
           sep = "_") %>% 
  mutate(fileName = str_replace_all(fileName, "_level-5.csv", "")) -> obj_longer

# grafica de barras por replica y marcador

obj_longer %>%
  ggplot(aes(x = Sample, y = RA, fill = Phylum)) +
  coord_flip() + 
  labs(x = "Sample Area", y = "Rel. Ab (%)") +
  geom_col() + 
  facet_grid(fileName ~ Rep, scales = "free_y") +
  theme_classic(base_size = 14) +
  ggsci::scale_fill_rickandmorty() -> barPlot
  # scale_fill_brewer(palette = "Paired")
  # ggsci::scale_fill_rickandmorty()

ggsave(barPlot, filename = "barPlot.png", 
       path = dir, width = 12, height = 14)


# barplot agrupado

obj_longer %>%
  ggplot(aes(y = ab, x = Sample, fill = fileName)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(y = "N reads", x = "Sample Area") +
  facet_wrap(~ Phylum, scales = "free_y") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   vjust = 1)) -> barPlot2

ggsave(barPlot2, filename = "barPlot2.png", 
       path = dir, width = 16, height = 12)

# heatmap

obj_longer %>%
  ggplot(aes(y = Order, x = Sample, fill = log10(ab + 1))) +
  geom_tile() +
  facet_grid(Phylum ~ fileName, scales = "free_y", space = "free") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   vjust = 1),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 12)) + 
  guides(fill = guide_colorbar(barheight = unit(14, "cm"),
                               barwidth = unit(1, "cm"))) +
  scale_fill_viridis_c(name = "Relative\nAbundance (log2)") -> heatPlot

ggsave(heatPlot, filename = "heatmap.png", path = dir, 
       width = 12, height = 14)

# Visualizar parte funcional

rm(list = ls())

dir <- '~/Downloads/archivosarticulosloberas/'
files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)

library(readxl)

obj <- read_xlsx(files[1])

library(tidyverse)

obj %>%
  select_if(is.double) %>%
  names() -> colNames

obj %>%
  pivot_longer(cols = colNames, values_to = "RA", names_to = "Index") %>%
  separate(Index, into = c("Rookery", "Year", "Sample"), sep = "-") -> obj_longer


obj_longer %>%
  group_by(Rookery, Category, SuperPathway) %>%
  summarise(mean = mean(RA)) -> dataViz

dataViz %>%
  ggplot(aes(x = SuperPathway, y = mean, fill = Category)) +
  geom_col() +
  labs(y = "Relative Abudance of predicted gene (mean %)") +
  coord_flip() +
  facet_grid(~ Rookery) +
  theme_bw() 
  
# arrange by factor

dataViz %>%
  arrange(desc(Category)) %>%
  mutate(SuperPathway = factor(SuperPathway, 
                               levels = unique(SuperPathway))) -> dataViz


dataViz %>%
  filter(mean > 10) %>%
  ggplot(aes(x = SuperPathway, y = mean, fill = Category)) +
  geom_col() +
  labs(y = "Relative Abundance of predicted gene (mean %)") +
  coord_flip() +
  facet_grid(~ Rookery) +
  theme_bw(base_size = 16) -> barPlot

ggsave(barPlot, filename = "barKegg.png", path = dir, 
       width = 18, height = 10)

# crear una funcion que:
# 1) Leer archivo
# 2) Reordenar columnas de loberas
# 3) sacar la media 
# 4) Hacer figura

my_2nd_function <- function(file) {
  
  obj <- read_xlsx(file)
  
  obj %>% select_if(is.double) %>% names() -> colNames
  
  obj %>%
    pivot_longer(cols = colNames, values_to = "RA", names_to = "Index") %>%
    separate(Index, into = c("Rookery", "Year", "Sample"), sep = "-") %>%
    group_by(Rookery, Category, SuperPathway) %>%
    summarise(mean = mean(RA)) %>%
    mutate(Amplicon = basename(file)) %>%
    mutate(Amplicon = str_replace_all(Amplicon, 
                                      "_path_abun_unstrat_descrip.xlsx", "")) %>%
    arrange(desc(Category)) %>%
    mutate(SuperPathway = factor(SuperPathway, 
                                 levels = unique(SuperPathway)))
}

# comprobamos 

my_2nd_function(files[1]) %>% 
  # z transform to clean outliers (datos de alto ruido) x - mean(x) / sd(x)
  mutate(zT = (mean - mean(mean)) / sd(mean)) %>%
  filter(mean > 10) %>%
  ggplot(aes(x = SuperPathway, y = mean, fill = Category)) +
  geom_col() +
  labs(y = "Relative Abundance of predicted gene (mean %)") +
  coord_flip() +
  facet_grid(~ Rookery) +
  theme_bw(base_size = 16)

# metemos todas los files en una sola tabla

lapply(files, my_2nd_function) %>% do.call(rbind, .) -> dataViz

# dataViz %>%
#   group_by(Amplicon) %>%
#   mutate(zT = (mean - mean(mean)) / sd(mean)) %>%
#   summarise(min  = min(zT), max = max(zT))

dataViz %>% 
  filter(mean > 100) %>%
  ggplot(aes(x = SuperPathway, y = mean, fill = Category)) +
  geom_col() +
  labs(y = "Relative Abundance of predicted gene (mean %)") +
  coord_flip() +
  facet_grid(~ Amplicon, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   vjust = 1,
                                   size = 10)) -> ampliconBarKegg

ggsave(ampliconBarKegg, filename = "ampliconBarKegg.png", path = dir,width = 18, height = 10)

# function to plot

my_plot_function <- function(file, min_abund = 10) {
  
  my_2nd_function(file) %>%
    filter(mean > min_abund) %>%
    ggplot(aes(x = SuperPathway, y = mean, fill = Category)) +
    geom_col() +
    labs(y = "Relative Abundance of predicted gene (mean %)") +
    coord_flip() +
    facet_grid(~ Rookery) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1,
                                     vjust = 1,
                                     size = 12))

  ggsave(filename = paste0(basename(file), ".png"), 
           path = dir, width = 18, height = 10)
}



lapply(files, my_plot_function)



