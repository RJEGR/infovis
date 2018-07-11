# Import bathymetry
#bat <- getNOAA.bathy(-100, -80, 18, 31, res = 1, keep = TRUE)
bat <- getNOAA.bathy(-99, -84, 18, 26, res = 1, keep = TRUE)
# Load location of individuals (these are NOT from Viricel 2012)
# loc <- data.frame( x = c(-96.92707, -96.60861, -96.86875, -96.14351, -92.82518, -90.86053, -90.14208, -84.64081, -83.81274, -81.13277, -80.33498, -88.52732, -94.46049), y = c(25.38657, 25.90644, 26.57339, 27.63348, 29.03572, 28.16380, 28.21235, 26.71302, 25.12554, 24.50031, 24.89052, 30.16034, 29.34550) )

loc <- read.csv("~/Downloads/metadata_xiximi04.csv", header=TRUE, stringsAsFactors = FALSE)
loc <- loc[-c(17,19,21),] # remove duplicates in Y
loc[c(16:19),2] <- substr(loc[c(16:19),2], 1,2 ) # rename Ys
rownames(loc) <- loc$Station
#loc <- loc[order(loc$Station),]

loc$presence <- 0
presence <- c("A1", "A6", "B13", "C25", "E35", "E33", "E34", "F37", "F38", "H46", "H47","Y1", "Y3", "Y4")
which(loc$Station %in% presence)# define position in presence 

loc[which(loc$Station %in% presence),8] <- 1

# coor loc[, c(3:4)]
tr <- trans.mat(bat, min.depth = -5, max.depth = -300)
cost <- lc.dist(tr, loc[, c(3:4)], res="path")
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

# Plot map with isobaths every 1000m
plot(bat, image = TRUE, 
        land = TRUE, 
        deep=-4000, 
        shallow=-1000, 
        step=1000, 
        drawlabels = FALSE, 
        bpal = list(c(min(bat,na.rm=TRUE), 0, blues), 
        c(0, max(bat, na.rm=TRUE), greys)), lwd = 0.1)


plot(bat, image = FALSE)

scaleBathy(bat, deg = 2, x = "bottomleft", inset = 5)
points(loc[, c(3:4)], bg = "orange", cex = 0.8, pch = 21)

#text(sites[,1], sites[,2], lab = rownames(sites),
text(loc[which(loc$Station %in% presence), c(3:4)], lab = rownames(loc))


ggplot2::autoplot(atl, geom=c("raster", "contour"), colour="white", size=0.1) + scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen")



pos = c(3, 4, 1, 2), col = "blue")
lapply(out1, lines, col = "orange", lwd = 5, lty = 1) -> dummy
lapply(out2, lines, col = "black", lwd = 1, lty = 1) -> dummy



