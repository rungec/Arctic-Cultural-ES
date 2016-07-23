#Arctic Cultural Ecosystem Service Modelling
#code by Claire Runge runge@nceas.ucsb.edu
#Project in collaboration with Vera Hausner Associate Professor, Department of Arctic and Marine Biology, The Arctic University of Tromso vera.hausner@uit.no
#This script takes PPGIS data and generates predictive spatial maps of cultural ecosystem services across Norway

##########################################
#Setup

options(stringsAsFactors=FALSE) # turn off automatic factor coersion
# options(scipen=9999)            # turn off plotting axis lables in scientific notation

library(dismo) #for maxent interface
library(ENMtools) #for AICc
library(raster)
library(foreign) #read dbfs

wd <- "D:/Arctic Cultural ES/"
setwd(wd)

#load ces observations (PPGIS data)
markersN <- read.dbf(paste0(wd, "Cultural PPGIS Data/PPGIS_Markers_north_UTM33N.dbf"))
markersS <- read.dbf(paste0(wd, "Cultural PPGIS Data/PPGIS_Markers_south_UTM33N.dbf"))

#load environmental data
rastDir <- paste0(wd, "Spatial data/")
varFolder <- c("CORINE2006", 

# need
# template raster (mask raster)
# Categorical		corrine data as .tif clipped to mask
# Continuous		distance to houses as raster or housing density
# Continuous		distance to roads as raster
# ??				accessibility - how is this different to distance to roads
# Categorical 	land ownership /protected areas 


##combine north & south ppgis data
names(markersS)[which(names(markersS)=="LoginID")] <- "LogID"
markersAll <- rbind(markersN, markersS[,!names(markersS) %in% c("zoom")])
markersAll$region = rep(c("north", "south"), times=c(nrow(markersN), nrow(markersS)))

#table of how many points for each type of CES in each region
counttable <- with(markersAll, table(category, region))
write.csv(counttable, paste0(wd, "Cultural PPGIS Data/NumberofPointsbyCEScategory.csv"))