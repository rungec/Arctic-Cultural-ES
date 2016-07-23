#Arctic Cultural Ecosystem Service Modelling
#code by Claire Runge runge@nceas.ucsb.edu
#Project in collaboration with Vera Hausner Associate Professor, Department of Arctic and Marine Biology, The Arctic University of Tromso vera.hausner@uit.no
#This script takes PPGIS data and generates predictive spatial maps of cultural ecosystem services across Norway

##########################################
#Setup

options(java.parameters = "-Xmx2g" ) #to give java more memory (2G)
options(stringsAsFactors=FALSE) # turn off automatic factor coersion
# options(scipen=9999)            # turn off plotting axis lables in scientific notation

library(dismo) #for maxent interface
library(ENMtools) #for AICc
library(raster)

wd <- "D:/Arctic Cultural ES/"
setwd(wd)

#load ces observations (PPGIS data)
markersN <- read.csv(paste0(wd, "Cultural PPGIS Data/PPGIS_Markers_north"), header=TRUE)
markersS <- read.csv(paste0(wd, "Cultural PPGIS Data/PPGIS_Markers_south"), header=TRUE)

#CES we are interested in
cultESlist <- c("", "")

#load environmental data
rastDir <- paste0(wd, "Spatial data/")
varFolder <- c("CORINE2006", 

need
template raster (mask raster)
Categorical		corrine data as .tif clipped to mask
Continuous		distance to houses as raster or housing density
Continuous		distance to roads as raster
??				accessibility - how is this different to distance to roads
Categorical 	land ownership /protected areas 


#use 10000 background points
#use multiple backgrounds?

#fit using S, test using N and vice versa
#http://www.inside-r.org/packages/cran/dismo/docs/maxent

