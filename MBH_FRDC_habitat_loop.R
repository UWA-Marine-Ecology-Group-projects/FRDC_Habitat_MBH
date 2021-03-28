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




# Read raster of classes ----
# Classes based on sentinel 2 unsupervised classification of sentinel --


#first import all files in a single folder as a list --
rast.list <- list.files(path = r.dir, pattern='.tif$', all.files=TRUE, full.names=TRUE)

#import all raster files in folder using lapply --
allrasters <- lapply(rast.list, raster) # list of raster files
allrasters

# rc <- raster(paste(r.dir, "Sen2_c1_unsup.tif", sep='/'))
# rc
# plot(rc)

# Get the crs in lat long ----
rcrs <- proj4string(allrasters[[1]])

# Read FRDC polygons ----
ps <- readOGR(paste(s.dir, "icoast_frdc_polygons.shp", sep='/'))

# transform to lat long --
psu <- spTransform(ps, CRSobj = rcrs)
length(levels(psu@data$name))
polynameslist <- as.list(levels(psu@data$name))



for( ii in c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")){
  print (ii)
  nam <- paste(ii)
  assign(nam, psu[psu$name == nam,])
}

# List the site polygons, make sure same order as rasters ----
#zones <- do.call(as.list, lapply(polynameslist, readOGR))
zones <- list() # make a loop or something for this!
zones$C1 <- C1
zones$C10 <- C10
zones$C11 <- C11
zones$C12 <- C12
zones$C2 <- C2
zones$C3 <- C3
zones$C4 <- C4
zones$C5 <- C5
zones$C6 <- C6
zones$C7 <- C7
zones$C8 <- C8
zones$C9 <- C9
zones



####    Get Inclusion Probabilities ####


for(i in 1:length(allrasters)){
  
  # Get proportion based on number of classes --
  rast <- allrasters[[i]]
  rdf <- as.data.frame(rast, na.rm = TRUE, stringsAsFactors = TRUE)
  rdf[,1] <- as.factor(as.character(rdf[,1]))
  no.classes <- length(levels(rdf[,1]))
  xprop <- 1/no.classes
  class.targetProps <- rep(xprop, no.classes) # all classes same proportion - sum equals 1
  #class.targetProps <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125) # all classes same proportion - sum equals 1
  
  # Get proportions of incl probs for each class for each location according to class.targetProps --
  class.freq <- raster::freq(rast)
  class.freq.df <- as.data.frame(class.freq)
  class.freq.df <- na.omit(class.freq.df)
  propsOfClass <- class.freq.df$count/sum(class.freq.df$count)
  
  # Make raster of inclusion probs for each location --
  inclProbs <- rast
  zoneID <- extract( x=inclProbs, y=zones[[i]], cellnumbers=TRUE)
  tmp <- class.targetProps / propsOfClass #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfClass)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
    inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
    inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
    string.no.classes <- seq(1, no.classes)
    inclProbs@data@values[inclProbs@data@values %in% string.no.classes] <- NA  # string has to be the same lenght as no. of classes - cheats way to crop
    #inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  # cheats way to crop
    plot( inclProbs)
    
  # Name and save raster of incluusion probs --
    nam <- paste(paste('c', i, sep =''), "InclProbs_SentFine.tif", sep ='_')
    writeRaster(inclProbs, paste(o.dir, nam, sep='/'))
    #assign(nam, inclProbs)
}





####    MBH design   ####

# import all incl prob in a single folder as a list --
ip.list <- list.files(path = o.dir, pattern='.tif$', all.files=TRUE, full.names=TRUE)

# import all raster files in folder using lapply --
all.ips <- lapply(ip.list , raster) # list of raster files
all.ips

# Set number of samples per location ----
numby <- 50 # how many samples

for(ii in 1:length(all.ips)){
  inclProbs <- all.ips[[ii]] # for each location
  newSites <- quasiSamp( n=numby, potential.sites=coordinates(inclProbs), 
                       inclusion.probs=values(inclProbs), nSampsToConsider=5000) # run MBH
  newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs))) # make sp
  get.name <- names(inclProbs) # get name of raster with site code
  site.code <- stringr::str_extract(first.part, "^.{3}") # extract fist two letters of name
  namsp <- paste(site.code, "MBH_SentFine", sep ='_') # make name of sp
  writeOGR(newSites, o.dir, namsp, driver = 'ESRI Shapefile', overwrite = T)
}





