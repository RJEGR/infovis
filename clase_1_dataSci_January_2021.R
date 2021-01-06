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
View(obj)

# # # # # #
# ordenar tabla de datos ----
# # # # # #

# start beepr::beep(3)

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



# previs




obj %>% 
  as_tibble() %>%
  # select_if(Negate(is.double)) %>%
  select(Samples_Name,  Sampling_Area, Replicate )

ranks <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

obj %>% 
  as_tibble() %>%
  select_if(is.double)  %>% 
  select(-Sampling_Date) %>%
  t() %>% as_tibble(rownames = 'linaje') %>% 
  separate(col = linaje, into = ranks, sep = '[.]') %>%
  mutate_all(., funs(str_replace_all(., c("^D_[0-9]__"="")))) 

# Expresiones regulares

