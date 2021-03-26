###   ###     Join the lidar for locations  ###   ###

# Libraries ----

library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(sp)
library(sf)
library(raster)
library(rgdal)
library(spastat)


# set working directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')
d.dir <- paste(w.dir, "data", sep='/')
x.dir <- "D:/Bathymetry/MidwestLidar"


# Read lidar files ----

tr <- raster(paste(x.dir, "TR2016_Mean_Lidar", "TR2016_Mean_Lidar.tif", sep='/'))
se <- raster(paste(x.dir, "SE2016_Mean_Lidar", "SE2016_Mean_Lidar.tif", sep='/'))
pd <- raster(paste(x.dir, "PD2016_Mean_Lidar", "PD2016_Mean_Lidar.tif", sep='/'))
lp <- raster(paste(x.dir, "LP2016_Mean_Lidar", "LP2016_Mean_Lidar.tif", sep='/'))
la <- raster(paste(x.dir, "LA2016_Mean_Lidar", "LA2016_Mean_Lidar.tif", sep='/'))
ju <- raster(paste(x.dir, "JU2016_Mean_Lidar", "JU2016_Mean_Lidar.tif", sep='/'))
gn <- raster(paste(x.dir, "GN2016_Mean_Lidar", "GN2016_Mean_Lidar.tif", sep='/'))
fp <- raster(paste(x.dir, "FP2016_Mean_Lidar", "FP2016_Mean_Lidar.tif", sep='/'))
ce <- raster(paste(x.dir, "CE2016_Mean_Lidar", "CE2016_Mean_Lidar.tif", sep='/'))
ak <- raster(paste(x.dir, "AK2016_Mean_Lidar", "AK2016_Mean_Lidar.tif", sep='/'))
ab <- raster(paste(x.dir, "AB2016_Mean_Lidar", "AB2016_Mean_Lidar.tif", sep='/'))


# Get extent ----

wa <- raster("D:/Bathymetry/wcbbathy/depth/w001001.adf")

crs1 <- proj4string(tr)
crs2 <- proj4string(wa)

wa <- projectRaster(wa, crs = crs1)
ext1 <- wa@extent


# list rasters --
li.list <- list(tr, se, pd, lp, la, ju, gn, fp, ce, ak, ab)

tr2 <- raster::extend(tr, ext1)

li.list <- list(tr, se)

li.list$tolerance <- 1
li.list$filename <- 'test.tif'
li.list$overwrite <- TRUE

sw.lidar <- do.call(raster::merge, li.list)
plot(sw.lidar)

ext2 <- sw.lidar@extent

tr2 <- raster::extend(tr, ext2)
plot(tr2)
origin(tr2)

sw.lidar <- raster::mosaic(tr$TR2016_Mean_Lidar, se$SE2016_Mean_Lidar, tolerance = 0.5, fun = 'mean')
plot(sw.lidar)

plot(se)


se2 <- projectRaster(tr, se)
origin(se)
se3 <- raster::merge(se, se2)
plot(se2)
origin(se2)

test <- projectRaster(se, template)
test2 <- raster::merge(tr2, se)
plot(test2)

plot(tr)
plot(se, add=T)
test2 <- mosaic(tr, se, tolerance = 1, fun=sum)
plot(test2)



##########################################

# PARAMS ----
waorigin <- origin(wa)
wa
wa2 <- raster::disaggregate(wa, fact = c(47.7788, 55.3754), method = 'bilinear')
trres <- res(tr)
trcrs <- proj4string(tr)

reproject_align_raster<- function(rast, ref_rast=NULL, desired_origin, desired_res, desired_crs, method= "bilinear"){
  
  if (!is.null(ref_rast)) {
    desired_origin<- origin(ref_rast) #Desired origin
    desired_res<- res(ref_rast) #Desired resolution
    desired_crs<- crs(ref_rast) #Desired crs
  } #Set parameters based on ref rast if it was supplied
  if(length(desired_res)==1){
    desired_res<- rep(desired_res,2)}
  
  if(identical(crs(rast), desired_crs) & identical(origin(rast), desired_origin) & identical(desired_res, res(rast))){
    message("raster was already aligned")
    return(rast)} #Raster already aligned
  
  if(identical(crs(rast), desired_crs)){
    rast_orig_extent<- extent(rast)} else{
      rast_orig_extent<- extent(projectExtent(object = rast, crs = desired_crs))} #reproject extent if crs is not the same
  var1<- floor((rast_orig_extent@xmin - desired_origin[1])/desired_res[1])
  new_xmin<-desired_origin[1]+ desired_res[1]*var1 #Calculate new minimum x value for extent
  var2<- floor((rast_orig_extent@ymin - desired_origin[2])/desired_res[2])
  new_ymin<-desired_origin[2]+ desired_res[2]*var2 #Calculate new minimum y value for extent
  n_cols<- ceiling((rast_orig_extent@xmax-new_xmin)/desired_res[1]) #number of cols to be in output raster
  n_rows<- ceiling((rast_orig_extent@ymax-new_ymin)/desired_res[2]) #number of rows to be in output raster
  new_xmax<- new_xmin+(n_cols*desired_res[1]) #Calculate new max x value for extent
  new_ymax<- new_ymin+(n_rows*desired_res[2]) #Calculate new max y value for extent
  rast_new_template<- raster(xmn=new_xmin, xmx =new_xmax,  ymn=new_ymin, ymx= new_ymax, res=desired_res, crs= desired_crs) #Create a blank template raster to fill with desired properties
  if(!identical(desired_origin,origin(rast_new_template))){
    message("desired origin does not match output origin")
    stop()} #Throw error if origin doesn't match
  if(identical(crs(rast),desired_crs)){
    rast_new<- raster::resample(x=rast, y=rast_new_template, method = method)} else{
      rast_new<- projectRaster(from=rast, to=rast_new_template, method = method)} #Use projectRaster if crs doesn't match and resample if they do
  if(!identical(desired_origin,origin(rast_new))){
    message("desired origin does not match output origin")
    stop()} #Throw error if origin doesn't match
  return(rast_new)
}

#' Combines multiple rasters into one
#'
#' Combines multiple rasters into one by using merge after first reprojecting or resampling and aligning rasters by matching them up with a a specified origin, resolution, and coordinate reference system, or that of a reference raster.
#' @param raster_list a list of rasters
#' @param rast raster to be reprojected or resampled
#' @param ref_rast reference raster with desired properties (Alternatively can supply desired_origin, desired_res, and desired_crs)
#' @param desired_origin desired origin of output raster as a vector with length 2 (x,y)
#' @param desired_res  desired resolution of output raster. Either an integer or a vector of length 2 (x,y)
#' @param desired_crs desired coordinate reference system of output raster (CRS class)
#' @param method resampling method. Either "bilinear" for bilinear interpolation (the default), or "ngb" for using the nearest neighbor
#' @param display_progress Logical specifying whether or not to indicate progress
#' @importFrom  raster merge
#' @export

combine_rasters<- function(raster_list, ref_rast=NULL, desired_origin, desired_res, desired_crs, method= "bilinear", display_progress=TRUE){
  raster_list2<- vector("list", length = length(raster_list))
  for (i in 1:length(raster_list)) {
    if(display_progress){
      print(paste("Reprojecting", as.character(i), "of", as.character(length(raster_list))))}
    raster_list2[[i]]<- reproject_align_raster(raster_list[[i]], ref_rast=ref_rast, desired_origin=desired_origin, desired_res=desired_res, desired_crs=desired_crs, method= method)}
  for (j in 1:length(raster_list2)) {
    if(display_progress){
      print(paste("Combining raster", as.character(j), "of", as.character(length(raster_list2))))}
    if(j==1){
      combined_raster<- raster_list2[[j]]} else{
        combined_raster<- raster::merge(combined_raster, raster_list2[[j]])}}
  return(combined_raster)
}

xtest <- combine_rasters(li.list, 0, trres, trcrs, method ='bilinear', display_progress=TRUE)
xtest 
plot(xtest)
