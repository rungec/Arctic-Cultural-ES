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
library(rgdal) #read shps
library(rgeos) #for buffer

wd <- "C:/Claire/Arctic_Cultural_ES/"
setwd(wd)


########################
## Setting system preferences for better work with raster
# Define your temp folder
my_tmpdir=paste0(wd, "tmpRaster")

# Create it (handles the case where the folder already exists)
dir.create(my_tmpdir, showWarnings=F)

# Set the raster option to this folder
rasterOptions(tmpdir= my_tmpdir)

## Maximum Memory
# Often server has more RAM than desktop. Therefore the mamximum amount of memory that can be allocated 
# to process a raster can be increased. Aurora has a lot of RAM (384GB), so we recommend to use the following 
# settings:
rasterOptions(maxmemory =1e+09)
rasterOptions(chunksize=1e+08)

############################
#Process PPGIS data
#load ces observations (PPGIS data)
markersN <- read.dbf(paste0(wd, "Cultural PPGIS Data/PPGIS_Markers_north_UTM33N.dbf"))
markersS <- read.dbf(paste0(wd, "Cultural PPGIS Data/PPGIS_Markers_south_UTM33N.dbf"))

##combine north & south ppgis data
names(markersS)[which(names(markersS)=="LoginID")] <- "LogID"
markersAll <- rbind(markersN, markersS[,!names(markersS) %in% c("zoom")])
markersAll$region = rep(c("north", "south"), times=c(nrow(markersN), nrow(markersS)))

#table of how many points for each type of CES in each region
counttable <- with(markersAll, table(category, region))
write.csv(counttable, paste0(wd, "Cultural PPGIS Data/NumberofPointsbyCEScategory.csv"))

############################
#Processing shps to rasters
rastDir <- paste0(wd, "Spatial data/")
rastTemplate = raster(paste0(rastDir, "Processed/Templates and boundaries/Norway_template_raster.tif"))

#Protected areas
PAshp <- readOGR(paste0(rastDir, "Original/Norway protected areas"), "Protected area Norway Svalbard")
PArast <- rasterize(PAshp, rastTemplate, field="IUCNCat", fun='min', background=0, filename=paste0(rastDir, "Processed/Protected_areas_norway.tif"), format = "GTiff", datatype="INT2S")
PArastTF <- reclassify(PArast, rcl=matrix(c(0:7, 0, rep(1, 7)), ncol=2), filename=paste0(rastDir, "Processed/Protected_areas_norway_truefalse.tif"), format = "GTiff", datatype="INT2S")

#Important ecological areas
Ecolshp <- readOGR(paste0(rastDir, "Original/NATURBASE/Naturtyper"), "Naturtyper_flater")
Ecolbuff <- gBuffer(Ecolshp, byid=TRUE, id=Ecolshp$FID, width=49) #because the polygons in this dataset are smaller than the raster I buffer them by 49m which is just smaller than half the raster resolution of 100m
Ecolrast <- rasterize(Ecolbuff, rastTemplate, field="VERDIint", background=0, fun='max', filename=paste0(rastDir, "Processed/Ecological_areas_norway.tif"), format = "GTiff", datatype="INT2S")

#State commons
stateshp <- readOGR(paste0(rastDir, "Original/Statskog eiendom"), "Statskog eiendom 2014")
statesub <- stateshp[stateshp@data$EKAT %in% c(2,3),]
staterast <- rasterize(statesub, rastTemplate, field="EKAT", background=0, filename=paste0(rastDir, "Processed/State_commons_norway.tif"), format = "GTiff", datatype="INT2S")

##############################
#Changing resolution of rasters
#All potential rasters
rastList <- list(
Distance_to_Town=raster(paste0(rastDir, "Original/SSB/Tettsted2015/tet_dist2015")),
Distance_to_River=raster(paste0(rastDir, "Original/N50 Data/Rivers50/river_dist")),
Distance_to_Coast= raster(paste0(rastDir, "Original/N250 Data/Coastline250/coastdistance")),
Distance_to_Road=raster(paste0(rastDir, "Original/N250 Data/Roads250/roaddist")),
Distance_to_House=raster(paste0(rastDir, "Original/N250 Data/Buildings250/house_dist"))
)

a <- lapply(1:length(rastList), function(x) {
		newRast <- resample(rastList[[x]], rastTemplate, method='ngb')
		writeRaster(newRast,  filename=paste0(rastDir, "Processed/", names(rastList)[[x]], "_norway.tif"), format = "GTiff", datatype="INT2S")
		})
		

		
		
		
		
		
		
###########################	
#remove temp raster folder	
unlink(my_tmpdir, recursive = TRUE)
#########################
#END
#########################

#load environmental data

varFolder <- c("CORINE2006", 

# need
# template raster (mask raster)
# Categorical		corrine data as .tif clipped to mask
# Continuous		distance to houses as raster or housing density
# Continuous		distance to roads as raster
# ??				accessibility - how is this different to distance to roads
# Categorical 	land ownership /protected areas 


		
lapply(1:length(rastList), function(x) {
		#newRast <- resample(x, Norway_template_raster, method='ngb')
		print(paste0(rastDir, "Processed/", names(rastList)[[x]], ".tif"))
		}
		#crop(x, extent(Norway_template_raster)))
		)
		
		crop(x, extent(Norway_template_raster)))