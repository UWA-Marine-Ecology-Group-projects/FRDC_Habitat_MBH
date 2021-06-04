

###     ###     ###     Test comparison of GRTS and MBH sampling design     ###     ###     ###

###     ###     ###     Based on fine unsup class of Sentinel2   ###     ###     ###


library( rgdal)
library( sp)
library( raster)
library( rgdal)
library( rgeos)
#install.packages("MBHdesign")
#library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)
library( stringr)

# clear environment ----
rm(list = ls())

# Set working directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
p.dir <- paste(w.dir,"plots",sep="/")
o.dir <- paste(w.dir,"outputs",sep="/")
d.dir <- paste(w.dir,"data",sep="/")
s.dir <- paste(w.dir,"shapefiles",sep="/")
r.dir <- paste(w.dir,"rasters",sep="/")
x.dir <- paste(w.dir,"GRTS_Sent2_Fine",sep="/")

# list rasters ----
list1 <- dir(r.dir)
listMBH <-  dir(o.dir, pattern = ".shp")
listGRTS <- dir(x.dir, pattern = ".shp")


# Read raster of classes ----
r <- raster(paste(r.dir, list1[[1]], sep ='/'))
plot(r)

# Read sample design outputs ----
mbh <- readOGR(paste(o.dir, listMBH[[1]], sep='/'))
mbh
grts <- readOGR(paste(x.dir, listGRTS[[1]], sep='/'))
grts


# extract classes from sample designs ----

