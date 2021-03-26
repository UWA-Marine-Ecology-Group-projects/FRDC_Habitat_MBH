# Get inclusion probabilities ---

## Create inclusion probabilities ####

library( rgdal)
library( sp)
library( raster)
library( rgdal)
library( rgeos)
library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)

# clear environment ----
rm(list = ls())

# Set working directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
p.dir <- paste(w.dir,"plots",sep="/")
o.dir <- paste(w.dir,"outputs",sep="/")
d.dir <- paste(w.dir,"data",sep="/")
s.dir <- paste(w.dir,"shapefiles",sep="/")
r.dir <- paste(w.dir,"rasters",sep="/")


# set site ----
site <- "c1_AbEast"

classby <- "SentFine"


# Read raster of classes ----
# Classes based on sentinel 2 unsupervised classification of sentinel --

#c1 <- readOGR(paste(s.dir, "Sen2_c1_unsup.shp", sep='/'))
rc <- raster(paste(r.dir, "Sen2_c9_unsup.tif", sep='/'))
rc
plot(rc)
rcrs <- proj4string(rc)

# Read FRDC polygons ----
ps <- readOGR(paste(s.dir, "icoast_frdc_polygons.shp", sep='/'))

# transform to lat long --
psu <- spTransform(ps, CRSobj = rcrs)

c9 <- psu[psu$name == 'C9',]
plot(c9, add = TRUE)

zones <- c9
plot(zones)
zones


# Proportions for each class ----
# classes
xprop <- 1/10
class.targetProps <- rep(xprop, 10)
#class.targetProps <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125) # sum equals 1


# Get inclusion probablities according to class.targetProps ----

catB <- rc
rcdf <- raster::as.data.frame(rc, xy=T)
rcdf[,3] <- as.factor(rcdf[,3])
str(rcdf)
rcdf2 <- na.omit(rcdf) 
str(rcdf2)
test <- raster::freq(rc)
class(test)
test
test2 <- as.data.frame(test)
test2
test2 <- na.omit(test2) # remove NAs
test2
propsOfClass <- test2$count/sum(test2$count) # sum equals 1


#test3 <- t(test2)

# test <- table(rcdf)
# xy.list <- split(rcdf, seq(nrow(rcdf)))
# head(xy.list)
# test <- table(xy.list)

plot(catB)
inclProbs <- catB

zoneID <- extract(x=catB, y=zones, cellnumbers=TRUE)
#zoneID <- rcdf
# propsOfbathy <- table(catB@data@values[zoneID[[1]][,"cell"]])
# 
# propsOfbathy <- propsOfbathy / sum( propsOfbathy)

tmp <- class.targetProps / propsOfClass #the desired inclusion probs (unstandardised)

for( ii in 1:length( propsOfClass)){
  inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
}


inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum(inclProbs[zoneID[[1]][,"cell"]])

#inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8,9,10)] <- NA  #cheats way to crop
inclProbs@data@values[inclProbs@data@values %in% string.no.classes] <- NA # string has to be the same lenght as no. of classes
plot( inclProbs)



# to check if sum of cells (probs) = 1
sumr <- cellStats(inclProbs, 'sum')
sumr


# Save inclusion probabilities ----
writeRaster(inclProbs, paste(o.dir, paste('inclProbs', site, paste(classby, 'tif', sep = '.'), sep='_'), sep ='/'), overwrite=TRUE)

# Check that there are 8 inclusion probs --
test<- inclProbs
test <- as.data.frame(test, xy=T)
head(test)
test$Sen2_c1_unsup <- as.factor(test$Sen2_c1_unsup)
str(test)



####    MBH design   ####

numby <- 50 # how many samples

newSites <- quasiSamp( n=numby, potential.sites=coordinates(inclProbs), 
                       inclusion.probs=values(inclProbs), nSampsToConsider=5000)
newSites
# plot design
plot(inclProbs)
points( newSites[,c("x","y")], pch=20, col='black')

# Make points sp ----
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#writeOGR(newSites, o.dir, "AB_NPZ6_BOSS", driver = "ESRI Shapefile")


### Make sure the drops are ~ 50 m apart ----

#inclProbs <- raster(paste(d.dir, "inclProbs_forBOSS.tif", sep='/'))


#s <- readOGR(paste(o.dir, "GB_BOSS_d1.shp", sep='/'))
# make sure is in UTM
s <- newSites
#proj4string(s) <- proj4string(rc)
s2 <- spTransform(s, CRSobj = proj4string(ps))
s2

## calculate if 2 points fall within 25m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance


p1_matrix <- gWithinDistance(s2, dist = 50, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:
p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
s2[v1, ] # 40 features left


## add unique id to points ----
IDnumber <- paste0(1:50)
s$IDnumber <- IDnumber
s


writeOGR(s, o.dir, paste(site, classby, "MBH", sep='_'), driver = "ESRI Shapefile", overwrite = TRUE)
