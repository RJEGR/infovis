# Ricardo Gomez Reyes, 2021
# ir directo a # zoom map 
# Visualising temp, pCO2 and pH profiles
# Transecto Pacifico Norte
# Coastal Ocean Data Analysis Product in North America (CODAP-NA, Version 2021) (NCEI Accession 0219960)
# NOAA Ocean Acidification Program (OAP);
# https://www.ncei.noaa.gov/data/oceans/ncei/ocads/metadata/0219960.html

# SPATIAL COVERAGE:
# NORTH: 62.78
# WEST: -179.44
# EAST: -65.45
# SOUTH: 18.83

rm(list=ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

path <- '~/Documents/DOCTORADO/'
files <- list.files(pattern = 'CODAP_NA_v2021', path = path, full.names = T)

datesc <- c('Year_UTC',	'Month_UTC',	'Day_UTC')
coordc <- c('Latitude','Longitude', 'CTDPRES', 'Depth')
var1c <- c('recommended_Salinity_PSS78', 'recommended_Oxygen', 'DIC', 'TALK', 'CTDTEMP_ITS90')

# vars that are calculated from CO2SYS with DIC and TALK
vars2c <- c('pH_TS_insitu_calculated', 'fCO2_insitu_calculated', 'Carbonate_insitu_calculated', 'Aragonite', 'Calcite')
# [umol/kg] -> carbonate, recommended_Oxygen,DIC, TALK 
# fCO2 -> [uatm]

cols <- c(datesc, coordc, var1c, vars2c)

x <- lapply(files, function(x) {
  y <- read.csv(x, skip = 0, comment.char = "#", stringsAsFactors = F)
  y <- data.frame(y, GROUP = basename(x))
  y <- y %>% select_at(vars(all_of(cols)))
  # y <- y %>% rename_with(~ toupper(gsub("_SWS", "_TOT", .x, fixed = TRUE)))
  return(y)})

x <- do.call(rbind, x) %>% 
  as_tibble() %>% drop_na() 

# names(x)

# x %>% mutate(Year_UTC = as.factor(Year_UTC)) -> x

x %>% mutate(reg = ifelse(between(abs(Longitude), 100, 180), 'Pacific', 'Atlantic')) -> x

  # pivot_longer(cols = all_of(c(var1c, vars2c)), names_to = 'var') %>%
  # filter(value > 0) %>%
  # mutate(var = factor(var, levels = cols, labels = labels)) -> y_longer

x %>% pivot_longer(cols = all_of(c(var1c, vars2c)), names_to = 'var') %>% filter(value > 0) -> x_longer

# Atlantic
# U.S. East Coast; Gulf of Mexico; North Atlantic Ocean;


# Pacific

x %>% filter(reg %in% 'Pacific') %>% pull(Latitude) %>% summary()

# Bering Sea; Latitud = x,58
# Golfo de Alaska = 57
# U.S. West Coast; Latitud = 39
# North Pacific Ocean; Latitud = 26

x %>% filter(reg %in% 'Pacific') %>%
  mutate(oce = ifelse(between(Latitude, 24, 26), 'Oc. Pacifico Norte', 
                      ifelse(between(Latitude, 26.01, 39), 'U.S. Costa Oeste', 
                             ifelse(between(Latitude, 39.01, 57), 'Golfo de Alaska', 
                                    ifelse(between(Latitude, 57.01, 58), 'Mar de Bering', NA))))) -> pacific_df

names(pacific_df)[which(grepl('calculated', names(pacific_df)))] <- c('pH', 'fCO2', 'Carbonate')
names(pacific_df)[which(grepl('recommended', names(pacific_df)))] <- c('Salinity', 'Oxygen')
names(pacific_df)[which(grepl('CTD', names(pacific_df)))] <- c('PRES', 'Temp')


vars <- c('Temp','Oxygen','fCO2', 'pH', 'Carbonate','Aragonite', 'Calcite', 'Salinity','TALK', 'DIC')
oceL <- c('Oc. Pacifico Norte', 'U.S. Costa Oeste', 'Golfo de Alaska', 'Mar de Bering')

pacific_df %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>% 
  mutate(var = factor(var, levels = vars)) %>%
  mutate(oce = factor(oce, levels = oceL)) %>%
  drop_na(oce) %>%
  filter(value > 0) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite', 'Salinity', 'Oxygen')) %>%
  # group_by(Depth, oce, var) %>%
  # summarise(value = mean(value)) %>%
  filter(Depth > 10) %>%
  filter(oce %in% c('Oc. Pacifico Norte', 'U.S. Costa Oeste')) %>%
  ggplot(aes(y = Depth, x = value)) + # 
  geom_line(orientation = 'y', aes(color = oce)) +
  # geom_point(alpha = 0.5, size = 1.5, aes(color = oce)) +
  scale_y_reverse() +  
  theme_bw(base_family = "GillSans", base_size = 18) +
  scale_x_continuous(position = "bottom") +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) +
  labs(x = '', y = expression(~Profundidad~(dbar)))  +
  scale_shape_discrete(name = 'Coordenadas') +
  facet_grid(as.factor(Year_UTC)~var, scales = 'free_x', labeller = labeller(var = label_parsed))

# 


pacific_df %>% 
  pivot_longer(cols = all_of(vars), names_to = 'var') %>%
  group_by(Year_UTC,oce, var) %>% 
  filter(value > 0) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite', 'Salinity', 'Oxygen')) %>%
  summarise(y = mean(value)) 

library(ggpubr)

pacific_df %>%
  pivot_longer(cols = all_of(vars), names_to = 'var', values_to = 'y') %>%
  mutate(var = factor(var, levels = vars)) %>%
  mutate(Year_UTC = factor(Year_UTC, levels = unique(Year_UTC))) %>%
  filter(Depth <= 100) %>%
  drop_na(oce) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  filter(y > 0) %>%
  ggstripchart(., 
             x = "Year_UTC", y = "y",facet.by = 'var',
             color = "oce",
             # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "mean_sd")


# boxplot

pacific_df %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>% 
  mutate(var = factor(var, levels = vars)) %>%
  drop_na(oce) %>%
  filter(value > 0) %>%
  filter(Depth <= 100) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  ggplot(aes(x = as.factor(Month_UTC), y = value)) +
  facet_grid(var~oce, scales = 'free_y', labeller = labeller(var = label_parsed)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  # geom_boxplot(width = 0.4, position = position_dodge(0.6),
  #              outlier.alpha = 0.4) +
  stat_summary(fun=mean, geom="point", shape=23, 
               size=1, color = 'red', position = position_dodge(0.6)) +
  labs(y = '') +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) 

# group_by(oce, Year_UTC, var) %>% summarise(m = mean(value))

# or complete pacific

x %>%  filter(oce %in% 'Pacific') -> pacific_df
  # filter(Latitude <= 40) -> pacific_df

names(pacific_df)[which(grepl('calculated', names(pacific_df)))] <- c('pH', 'fCO2', 'Carbonate')
names(pacific_df)[which(grepl('recommended', names(pacific_df)))] <- c('Salinity', 'Oxygen')
names(pacific_df)[which(grepl('CTD', names(pacific_df)))] <- c('PRES', 'Temp')

pacific_df %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>% 
  mutate(var = factor(var, levels = vars)) %>%
  drop_na(oce) %>%
  filter(value > 0) %>%
  filter(Depth <= 100) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  ggplot(aes(x = as.factor(Year_UTC), y = value)) +
  facet_grid(var~., scales = 'free_y', labeller = labeller(var = label_parsed)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_violin(aes(fill = oce)) +
  # geom_boxplot(width = 0.4, position = position_dodge(0.6)) +
  stat_summary(fun=mean, geom="point", shape=23, size=1, color = 'red', position = position_dodge(0.6)) +
  labs(y = '') +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) 

# anualmente ----

pacific_df %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>% 
  mutate(var = factor(var, levels = vars)) %>%
  drop_na(oce) %>%
  filter(value > 0) %>%
  filter(Depth <= 100) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  # group_by(Year_UTC, var) %>%
  # summarise(value = mean(value)) %>%
  ggplot(aes(y = value, x = as.factor(Year_UTC), group = var)) + 
  facet_grid(var~., scales = 'free_y', labeller = labeller(var = label_parsed)) +
  geom_point(alpha = 0.7, size = 1.5) + # aes(color = as.factor(Month_UTC)
  geom_line(orientation = 'x') +
  labs(y = '') +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) 

pacific_df %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>% 
  mutate(var = factor(var, levels = vars)) %>%
  drop_na(oce) %>%
  filter(value > 0) %>%
  filter(Depth <= 100) %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  group_by(Month_UTC, var) %>%
  summarise(value = mean(value)) %>%
  ggplot(aes(y = value, x = as.factor(Month_UTC), group = var)) + 
  facet_grid(var~., scales = 'free_y', labeller = labeller(var = label_parsed)) +
  geom_point(alpha = 0.7, size = 1.5) + 
  geom_line(orientation = 'x') +
  labs(y = '') +
  theme_minimal() +
  theme(strip.background = element_blank(),
        axis.line.y =  element_line(colour = "black"),
        panel.border = element_blank()) 

# maps ----

x %>% filter(oce %in% 'Pacific') %>% pull(Latitude) %>% summary()

pacific_df

library("rnaturalearthdata")
library('rnaturalearth')
library('ggspatial')
library('sf')

world <- ne_countries(scale = "medium", returnclass = "sf")

gg1 <- ggplot(data = world) +
  geom_sf(fill = 'antiquewhite') +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-180, -100),  ylim = c(18, 76))

gg <- ggplot(data = world) +
  geom_sf(fill = 'antiquewhite') +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  coord_sf(xlim = c(-180, -100),  ylim = c(18, 76))

pacific_df %>%
  filter(Depth <= 100) %>%
  filter_at(vars(vars), ~ (. > 0)) %>% 
  group_by(Latitude, Longitude) %>%  # Year_UTC,
  summarise_at(vars(vars), mean) -> map_in

map_in %>% pull(Temp) %>% max() -> limit

gg1 +
  geom_point(data = map_in, 
             aes(Longitude, Latitude,
                 color = Temp), 
             alpha = 1, size = 0.5, fill = 'NA') + # shape = 15
  ggsci::scale_color_gsea(name = '') + # Temp (°C)
  # scale_size_continuous(name = '', guide = 'none') + 
  theme_bw() +
  theme(text = element_text(size = 5),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "top") +
  guides(color = guide_colorbar(barheight = unit(0.1, "in"), barwidth = unit(2, "in"), ticks.colour = "black", frame.colour = "black",label.theme = element_text(size = 10))) -> ptemp

# pH

gg +
  geom_point(data = map_in, 
             aes(Longitude, Latitude,
                 color = pH), 
             alpha = 1, size = 0.5, fill = 'NA') + # shape = 15
  scale_color_viridis_c(name = '') + # pH
  # scale_size_continuous(name = '', guide = 'none') + 
  theme_bw() +
  theme(text = element_text(size = 5),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "top") +
  guides(color = guide_colorbar(barheight = unit(0.1, "in"), barwidth = unit(2, "in"), ticks.colour = "black", frame.colour = "black",label.theme = element_text(size = 10))) -> ppH

# aragonite

gg +
  geom_point(data = map_in, 
             aes(Longitude, Latitude,
                 color = Aragonite), 
             alpha = 1, size = 0.5, fill = 'NA') + # shape = 15
  scale_colour_gradient(low = "yellow", high = "red",  name = '') + #Ara
  # scale_size_continuous(name = '', guide = 'none') + 
  theme_bw() +
  theme(text = element_text(size = 5),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "top") +
  guides(color = guide_colorbar(barheight = unit(0.1, "in"), barwidth = unit(2, "in"), ticks.colour = "black", frame.colour = "black",label.theme = element_text(size = 10))) -> pAr

library(patchwork)

ptemp + ppH + pAr -> psave

ggsave(psave, path = path, filename = 'ara_pH_temp_map.png', width = 10, height = 3.5)

# zoom map ----

library(gganimate)

gg1 <- ggplot(data = world) +
  geom_sf(fill = 'antiquewhite') +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), 
                         pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-130, -100),  ylim = c(18, 35))

gg <- ggplot(data = world) +
  geom_sf(fill = 'antiquewhite') +
  coord_sf(xlim = c(-130, -100),  ylim = c(20, 35))

pacific_df %>%
  filter(Depth <= 200) %>%
  filter(Depth > 10) %>%
  filter_at(vars(vars), ~ (. > 0)) -> map_in
  # group_by(Latitude, Longitude, Year_UTC) %>%  # Year_UTC,
  # summarise_at(vars(vars), mean) -> map_in

# COUNTOUR PLOT ----
# library(isoband)
library(metR)
# # https://eliocamp.github.io/metR/articles/Visualization-tools.html

map_in %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  ggplot(aes(value, Latitude)) +
  geom_density2d_filled(aes(color = Depth), contour_var = 'ndensity' ) +
  facet_grid(~ var, scales = 'free_x') +
  theme_classic() +
  labs(x = '') +
  theme(legend.position = "none",
        strip.background = element_blank(),
        # axis.line.y =  element_blank(),
        panel.border = element_blank()) -> psave

ggsave(psave, filename = 'density2d_seassons.png', path = path,width = 10,height = 2)

# hay algo raro aqui
map_in %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  filter(abs(Longitude) < 130 ) %>%
  ggplot(aes(value, Longitude)) +
  geom_density2d_filled(aes(color = value), contour_var = 'ndensity' ) +
  facet_wrap(~ var, scales = 'free_y',  strip.position = 'left') +
  theme_classic() +
  labs(y = '') +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.border = element_blank()) +
  coord_flip()


map_in %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>%
  filter(var %in% c('Temp', 'pH', 'Aragonite')) %>%
  ggplot(aes(x = as.factor(Month_UTC), y = value)) +
  facet_wrap(var~., scales = 'free_y') +
  geom_boxplot(width = 0.4, position = position_dodge(0.6), 
               outlier.shape = NA, color = 'grey') +
  stat_boxplot(geom ='errorbar', width = 0.6, color = 'grey') +
  stat_summary(fun = mean, geom = "line", aes(group = 1), col = "red") +
  labs(y = '', x = 'Mes') +
  theme_classic() +
  theme(strip.background = element_blank(),
        # axis.line.y =  element_blank(),
        panel.border = element_blank()) -> psave
  
ggsave(psave, filename = 'boxplot_seassons.png', path = path,width = 10,height = 2)

gg1 +
  geom_point(data = distinct(map_in, Latitude, Longitude), 
             aes(Longitude, Latitude),
                 # color = Year_UTC), 
             alpha = 1, size = 1, fill = 'NA', shape = 3) + # shape = 15
  # ggsci::scale_color_gsea(name = '') + # Temp (°C)
  theme_bw() +
  theme(text = element_text(size = 5),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "top") +
  guides(color = guide_colorbar(barheight = unit(0.1, "in"), barwidth = unit(2, "in"), ticks.colour = "black", frame.colour = "black",label.theme = element_text(size = 10))) -> psave

ggsave(psave, path = path, filename = 'pacific_map.png', width = 3.5, height = 3.5)

# COUNTOUR PLOT



# ptemp + transition_time(Year_UTC) -> ptemp
# anim_save(animation = ptemp, filename = "myanimation.gif", path = path)

# magick::image_write(anim, path= paste0(path, "/myanimation.gif"))


# library(sf)

crs = st_crs("EPSG:4326")
no2.sf = st_as_sf(map_in, coords = c("Latitude", "Longitude"), crs = "EPSG:4326") %>%
  st_transform(crs)

# build a grid over Germany:
library(stars)

st_bbox(gg1) %>%
  starst_as_stars(dx = 10000) %>%
  st_crop(gg1) -> grd

grd

# interpolate
# https://keen-swartz-3146c4.netlify.app/interpolation.html

library(gstat)
i = idw(Temp~1, no2.sf, grd)
