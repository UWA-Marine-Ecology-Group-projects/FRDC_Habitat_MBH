

### Join bathy and point coordinates ###

library( rgdal)
library( sp)
library( raster)

## make sure the spatial data used here is in the MEG-gis folder? so anyone can use this.. hmm think about this


d <- "G:/My Drive/meg_projects/Projects_WRL/Project_WRL_habitat and heatwaves/Data/Locations"
d2 <- "G:/My Drive/meg_projects/Projects_WRL/Project_WRL_habitat and heatwaves/Data/Bathy"

wa <- readOGR(paste(d, "WA_utm.shp", sep ='/'))
plot(wa)

fw <- readOGR(paste(d, "FreshWater_coords_utm_wgs84.shp", sep ='/'))
plot(fw, pch= 20, add=T)

b <- raster(paste(d2, "Dongara_Lidar_EPSG32750_FP2016.tif", sep='/'))
plot(b, add=T)

e <- drawExtent()

wa2 <- crop(wa, e)

plot(wa2)

e2 <- drawExtent()

wa3 <- crop(wa2, e2)

plot(wa3)

e3 <- drawExtent()
wa4 <- crop(wa3, e3)

plot(wa4, add=T)

# remove more 0m in bathy

b2 <- b[b > 0] <- NA
plot(b2, add=T)
plot(b)

c <- rasterToContour(b)
plot(c, add=T)

c2 <- rasterToContour(b, maxpixels = 500000)
plot(c2, add=T)
