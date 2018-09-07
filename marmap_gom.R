# good ref at https://rstudio-pubs-static.s3.amazonaws.com/301056_e188ebc4c8644410b4abbc4ae98b6c98.html

library(marmap)
# Import bathymetry
# bat <- getNOAA.bathy(-99, -84, 18, 26, res = 1, keep = TRUE)

# Load location of individuals 
loc <- read.csv("~/Downloads/metadata_xiximi04.csv", header=TRUE, stringsAsFactors = FALSE)
loc <- loc[-c(17,19,21),] # remove duplicates in Y
loc[c(16:19),2] <- substr(loc[c(16:19),2], 1,2 ) # rename Ys
rownames(loc) <- loc$Station
loc$Presence <- "False"
presence <- c("A1", "A6", "B13", "C25", "E35", "E33", "E34", "F37", "F38", "H46", "H47","Y1", "Y3", "Y4")

# define position in presence 
which(loc$Station %in% presence)
loc[which(loc$Station %in% presence),8] <- "True"


# tr <- trans.mat(bat, min.depth = -5, max.depth = -300)
# cost <- lc.dist(tr, loc[, c(3:4)], res="path")
# blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
# greys <- c(grey(0.6), grey(0.93), grey(0.99))

# Plot map with isobaths every 1000m

# plot(bat, image = TRUE, 
#      land = TRUE, 
#      deep=-4000, 
#      shallow=-1000, 
#      step=1000, 
#      drawlabels = FALSE, 
#      bpal = list(c(min(bat,na.rm=TRUE), 0, blues), 
#                  c(0, max(bat, na.rm=TRUE), greys)), lwd = 0.1)
# 
# 
# plot(bat, image = FALSE)
# 
# scaleBathy(bat, deg = 2, x = "bottomleft", inset = 5)
# points(loc[, c(3:4)], bg = "orange", cex = 0.8, pch = 21)
######

# El segundo approach es mas rapido y sencillo de elaborar 
# Estrategia 1
# 
library(ggmap)

gom <- c(left = -99, bottom = 18, right = -84, top = 26)

map2 <- get_map(gom, maptype = "satellite")

ggmap(map2) +
  geom_point(aes(x = Longitude, y = Latitude, color = Presence, 
                 group = Presence), data = loc, alpha = .5) + 
  geom_text(data = loc, aes(x = Longitude, y = Latitude, label = ifelse(Presence=="True", Station,'')), 
            size = 3, vjust = 0, hjust = -0.5, check_overlap = FALSE) +
  scale_color_manual(values =c('#1F78B4', '#33A02C'))

#'#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C']
#'

# Estrategia 2
library(ggrepel) # https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
library(ggplot2)

world <- map_data("world", region = "Mexico") 

gg <- ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(xlim = c(-99,-84), ylim = c(18,26)) # #xlim and ylim can be manipulated to zoom in or out of the map

#lon1,lon2,lat1,lat2
#-99, -84, 18, 26

gg + 
  geom_point(data = loc, aes(Longitude, Latitude), 
             color = ifelse(loc$Presence == "True",'red','black'), size=1.5, alpha = .8) + 
  ggtitle("") +
  geom_text_repel(data = subset(loc, Presence=="True"), 
                  aes(Longitude, Latitude, label=Station),
                  segment.color = "grey50"
                  ) + 
  theme(text = element_text(size=12),
        panel.background = element_blank(), legend.position = "none") 
