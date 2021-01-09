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
  ggplot(aes(x = Sampling_Area, y = ab, fill = Phylum)) +
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
obj_longer %>% tibble::glimpse()

obj_longer %>% 
  separate(col = linaje, into = ranks, sep = '[.]') %>% 
  mutate_at(ranks, funs(str_replace_all(., c("^D_[0-9]__"="")))) %>%
  separate(col = Replicate, into = c('Sample','Rep'), sep = '_') -> obj_longer

obj_longer %>% View()

obj_longer %>%
  filter(Rep %in% "1") %>%
  group_by(index, fileName) %>%
  mutate(RA = raT(ab)) %>%
  ggplot(aes(x = fileName, y = RA, fill = Phylum)) +
  coord_flip() + labs(x = 'Amplicon', y = 'Rel. Ab (%)') +
  geom_bar(stat = "identity", 
           position = "stack", color = "black") +
  facet_grid( Sample ~., scales = "free_y") +
  ggsci::scale_fill_rickandmorty() +
  theme_classic(base_size = 14)
