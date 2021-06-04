# Purpose: Example area GRTS survey designs for WAMSI II cenral site 
#          within Northern Kimberly
# Author: Ben Radford (b.radford@aims.gov.au)
# Date: June 26, 2012
# Last Modified: June 03, 2015 (orginal)

# Load the spsurvey library

library(spsurvey)

# set the random number seed can be what ever you want but has to be consistent
set.seed(3)

#Set working directory which contains your input shapefiles
setwd("G:/Icoasts/FRDC_sample_design_out")
# list the shapefiles dbf tables (optional) to find the file name you want
dir(pattern=".dbf$")

# Read the attribute table from the shapefile
shapeattr <- read.dbf("c1_unsupdiss.shp")


# Create the design list so in this example there are 20 points the PanelOne value devided between
# 8 levels which are derived from the GRIDCODE attribute in ArcGIS which is a polygonised ISOClASS Raster 
#postprocessed with a attribute disolve 
sampledesign <- list(None=list(panel=c(PanelOne=80),
                               seltype="Unequal",
                               caty.n=c("1"=16,
                                        "2"=16,
                                        "3"=16,
                                        "4"=16,
                                        "5"=16)))



# makes the GRIDCODE interger into a R factor varaible which is require for GRTS
shapeattr$gcf <- factor(shapeattr$gridcode)

# Create the GRTS survey design
# set design to to your list sample design see "sampledesign" above 
# set in.shape to input shapefile
# set att.frame to your dbf table name minus extention
# set DesignID to a meaninful name
# set prj to the .prj file which will be the same name as you input shapfile and dbf file 
#(this file has to exist in your working directory)
# set out.shape to your output shapfile sample design name
Unequalsites <- grts(design=sampledesign,
                     src.frame="shapefile",
                     in.shape="c1_unsupdiss.shp",  
                     att.frame=shapeattr,
                     type.frame="area",
                     mdcaty="gcf",									
                     DesignID="example1", 
                     shapefile=TRUE,
                     prj="c1_unsupdiss.shp", 
                     out.shape="example_GRTS") 

# END 