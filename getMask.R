# Create a maks of x metres distance around a set of non-NA cells created from a bathymetry map
# The non-NA cells are createed by setting all values above and below 5 m to NA

library(raster)
library(rasterVis)
library(xts)
library(maptools)
library(gridExtra)
library(grid)
library(rgdal)


# setup bathymetry map
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")
locbb <- matrix(c(9,-40,41,-20),nrow=2,ncol=2)
x <- extent(locbb)
bathy.cont <- crop(bathy,x)
bathy.loc <- crop(bathy,x)      # crop bathy map
bathy.land <- crop(bathy,x)
rm(bathy,locbb,x)

borders <- readShapePoly("SAfrica-coast.shp")
borders <- SpatialPolygons(borders@polygons)
minlat <- -40;minlon <- 9;maxlat <- -20;maxlon <- 41

buf.dist <- 370400                               # coastal boudary distance metres
bathy.land[bathy.land >= 0] = NA                    # set land values to NA
bathy.loc[bathy.loc > 5] <- NA                      # select 5 m height as coastal limits
bathy.loc[bathy.loc < -5] <- NA
bathy.buf <- buffer(bathy.loc, width=buf.dist, doEdge=TRUE)
bathy.mask.temp <- mask(bathy.buf,bathy.land)       # remove the inland buffer distance

p <- rasterToPolygons(bathy.mask.temp, dissolve=T)

setwd("~/Desktop/")
save(p, file = "~/Desktop/200nm_mask.RData")

