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
library(dplyr) #for %>%

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
rastDir <- paste0(wd, "Spatial data/")
rastTemplate = raster(paste0(rastDir, "Processed/Templates and boundaries/Norway_template_raster.tif"))

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

#clip south data to municipality boundaries
markersSshp <- readOGR(paste0(wd, "PPGIS data/Original/shps"), "PPGIS_Markers_south_UTM33N")
Smuns <- readOGR(paste0(rastDir, "Processed/Templates and boundaries"), "South_municipalities")
markersSsub <- markersSshp[Smuns,] #clip to municipal boundaries
writeOGR(markersSsub, paste0(wd, "PPGIS data/Original/shps"), "PPGIS_Markers_south_UTM33N_municipal", driver="ESRI Shapefile") 
write.csv(markersSsub@data, paste0(wd, "PPGIS data/Original/shps/PPGIS_Markers_south_UTM33N_municipal.csv"), row.names=FALSE)

############################
#Processing shps to rasters
#Protected areas
#PAshp <- readOGR(paste0(rastDir, "Original/Norway protected areas"), "Protected area Norway Svalbard")
#PArast <- rasterize(PAshp, rastTemplate, field="IUCNCat", fun='min', background=0, filename=paste0(rastDir, "Processed/Protected_areas_norway.tif"), format = "GTiff", datatype="INT2S")
#PArastTF <- reclassify(PArast, rcl=matrix(c(0:7, 0, rep(1, 7)), ncol=2), filename=paste0(rastDir, "Processed/Protected_areas_norway_truefalse.tif"), format = "GTiff", datatype="INT2S")
PArast <- raster(paste0(rastDir, "Processed/Protected_areas_norway.tif"))
PArastTF <- reclassify(PArast, rcl=matrix(c(0:7, 0,1,1,1,1,2,2,0), ncol=2), filename=paste0(rastDir, "Processed/Protected_areas_norway_forothervalues.tif"), format = "GTiff", datatype="INT2S") #nature protection (iucn1-4) = 1; managed landscapes (iucn 5,6) = 2; not protected = 0
PArastTF <- reclassify(PArast, rcl=matrix(c(0:7, 0,1,2,3,3,4,4,0), ncol=2), filename=paste0(rastDir, "Processed/Protected_areas_norway_forbiological.tif"), format = "GTiff", datatype="INT2S") #group 1: nature reserve; national park, iucn(2,3,4) = 3 ; (IUCN class 5,6) = 4; not protected = 0

#Important ecological areas
Ecolshp <- readOGR(paste0(rastDir, "Original/NATURBASE/Naturtyper"), "Naturtyper_flater")
Ecolbuff <- gBuffer(Ecolshp, byid=TRUE, id=Ecolshp$FID, width=49) #because the polygons in this dataset are smaller than the raster I buffer them by 49m which is just smaller than half the raster resolution of 100m
Ecolrast <- rasterize(Ecolbuff, rastTemplate, field="VERDIint", background=0, fun='max', filename=paste0(rastDir, "Processed/Ecological_areas_norway.tif"), format = "GTiff", datatype="INT2S")

#State commons
stateshp <- readOGR(paste0(rastDir, "Original/Statskog eiendom"), "Statskog eiendom 2014")
statesub <- stateshp[stateshp@data$EKAT %in% c(2,3,4,6),]
staterast <- rasterize(statesub, rastTemplate, field="EKAT", background=0, filename=paste0(rastDir, "Processed/State_commons_norway_all.tif"), format = "GTiff", datatype="INT2S")
staterast2 <- reclassify(staterast, rcl=matrix(c(0:6, 0,rep(1, 6)), ncol=2), filename=paste0(rastDir, "Processed/State_commons_norway_binary.tif"), format = "GTiff", datatype="INT2S") #private or private commons=0; state (state or municipal commons) =1

#water features
# lakes <- readOGR(paste0(rastDir, "Original/N50 Data/Lakes50"), "Innsjo_Innsjo")
# rivers <- readOGR(paste0(rastDir, "Original/N50 Data/Rivers50"), "Elv_Elvenett")
# riverssub <- rivers[rivers@data$objType=='ElvBekkMidtlinje',]
# allwater <- gUnion(lakes, riverssub)
watershp <- readOGR(paste0(rastDir, "Original/Statskog eiendom"), "Statskog eiendom 2014")
statesub <- stateshp[stateshp@data$EKAT %in% c(2,3,4,6),]
staterast <- rasterize(statesub, rastTemplate, field="EKAT", background=0, filename=paste0(rastDir, "Processed/State_commons_norway_all.tif"), format = "GTiff", datatype="INT2S")

##############################
#CORRINE2012 moving window 
maskTemplate <- readOGR(paste0(rastDir, "Processed/Templates and boundaries"), "Norway_border_10kmbuffer")
corrineshp <- readOGR(paste0(rastDir, "Original/CORINE2012"), "CORINE2012_Norge_ab21a_UTM33N")
corrineshp@data$newcode <- 1
corrineclass <- list(broadleafforest=c(311, 313), coniferforest=c(312), heathshrub=c(321:324), sparselyvegetated=c(331:335), cropland=c(211, 212, 213, 221, 222, 223, 231, 241:244), wetland=c(411, 412, 422, 423))

#function to make a circular weights matrix of given radius and resolution
#NB radius must me an even multiple of res!
make_circ_filter<-function(radius, res){
  circ_filter<-matrix(NA, nrow=1+(2*radius/res), ncol=1+(2*radius/res))
  dimnames(circ_filter)[[1]]<-seq(-radius, radius, by=res)
  dimnames(circ_filter)[[2]]<-seq(-radius, radius, by=res)
  sweeper<-function(mat){
    for(row in 1:nrow(mat)){
      for(col in 1:ncol(mat)){
        dist<-sqrt((as.numeric(dimnames(mat)[[1]])[row])^2 +
          (as.numeric(dimnames(mat)[[1]])[col])^2)
        if(dist<=radius) {mat[row, col]<-1}
      }
    }
    return(mat)
  }
out<-sweeper(circ_filter)
out<- out*100/sum(out, na.rm=TRUE)
return(out)
}

#make moving window matrix
currRadius <- 1000 #diameter of moving window
mwm <- make_circ_filter(currRadius/2, 100)
dimnames(mwm)<-NULL

#create a raster of the amount of each land class in radius around central cell
a <- lapply(c(1:length(corrineclass)), function(x) {
		currclass <- corrineclass[[x]]
		currshp <- corrineshp[corrineshp@data$CLC12_KODE %in% currclass, ]
		#corrRast <- rasterize(currshp, rastTemplate, field="newcode", background=0, filename=paste0(rastDir, "Processed/Corrine2012_norway_", names(corrineclass)[[x]], ".tif"), format = "GTiff", datatype="INT2S")		#focRast <- focal(corrRast, w=matrix(100/(29*29), nrow=29, ncol=29), filename=paste0(rastDir, "Processed/Corrine2012_norway_", names(corrineclass)[[x]], "_3km.tif"), format = "GTiff", datatype="INT2S")
		corrRast <- raster(paste0(rastDir, "Processed/Corrine2012_norway_", names(corrineclass)[[x]], ".tif"))
		focRast <- focal(corrRast, w=mwm, filename=paste0(rastDir, "Processed/Corrine2012_norway_", names(corrineclass)[[x]], "_",as.character(currRadius/1000), "km.tif"), format = "GTiff", datatype="INT2S")
		maskedRast <- mask(focRast, maskTemplate, filename=paste0(rastDir, "Processed/masked/Corrine2012_norway_", names(corrineclass)[[x]], "_",as.character(currRadius/1000), "km.tif"), format = "GTiff", datatype="INT2S")
		writeRaster(maskedRast, filename=paste0(rastDir, "Processed/forMaxent/Corrine2012_norway_", names(corrineclass)[[x]], "_",as.character(currRadius/1000), "km.asc"), format = "ascii")
		return()
		})




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

#resample to match resultion and extent, and extend as these rasters are smaller than the template raster
a <- lapply(1:length(rastList), function(x) {
		newRast <- raster::resample(rastList[[x]], rastTemplate, method='ngb')
		newRast2 <- raster::extend(newRast, rastTemplate, value=NA)
		writeRaster(newRast2,  filename=paste0(rastDir, "Processed/", names(rastList)[[x]], "_norway.tif"), format = "GTiff", datatype="INT4S")
		})

#Extend the corrine rasters
corrineRastList <- list.files(paste0(rastDir, "Processed/"), "Corrine.*tif$", full.names=TRUE)
b <- lapply(1:length(corrineRastList), function(x) {
		newRast <- raster(corrineRastList[x])
		newRast2 <- raster::extend(newRast, rastTemplate, value=NA)
		writeRaster(newRast2,  filename=corrineRastList[x], format = "GTiff", datatype="INT2S", overwrite=TRUE)
		})

###########################
#Mask out areas outside study region as NA
maskTemplate <- readOGR(paste0(rastDir, "Processed/Templates and boundaries"), "Norway_border_10kmbuffer")

envStack <- raster::stack(list.files(paste0(rastDir, "Processed/"), "*.tif$", full.names=TRUE))
maskedStack <- mask(envStack, maskTemplate)

setwd(paste0(rastDir, "Processed/masked/"))
writeRaster(maskedStack, filename=names(maskedStack), bylayer=TRUE, format = "GTiff", datatype="INT4S")

###########################
#Convert rasters to .asc
setwd(paste0(rastDir, "Processed/forMaxent/"))
writeRaster(maskedStack, filename=names(maskedStack), bylayer=TRUE, format = "ascii")

###########################
#make a mask raster from alpine climate regions to limit Maxent background to these regions
alpineshp <- readOGR(paste0(rastDir, "Processed/Templates and boundaries"), "Norway_alpine")
alpinerast <- mask(rastTemplate, alpineshp)
writeRaster(alpinerast, filename=paste0(rastDir, "Processed/Templates and boundaries/Norway_alpine.tif"), format = "GTiff", datatype="INT2S")
writeRaster(alpinerast, filename=paste0(rastDir,"Processed/forMaxent/Norway_alpine.asc"), format = "ascii")

##########
#Make mask raster from municipality boundaries
muns <- readOGR(paste0(rastDir, "Original/BasemapN500"), "n500_municipalsf")
#dput(muns@data$NAVN)
#north study region
northmuns <- muns[muns@data$NAVN %in% c("Sørfold", "Bodø", "Fauske", 
"Meløy", "Gildeskål", "Saltdal", "Rødøy", "Rana", "Beiarn"),]
writeOGR(northmuns, paste0(rastDir, "Processed/Templates and boundaries"), "North_municipalities", driver="ESRI Shapefile")
northrast <- mask(rastTemplate, northmuns)
writeRaster(northrast, filename=paste0(rastDir, "Processed/Templates and boundaries/North_municipalities.tif"), format = "GTiff", datatype="INT2S")
writeRaster(northrast, filename=paste0(rastDir,"Processed/forMaxent/North_municipalities.asc"), format = "ascii")
#south study region
southmuns <- muns[muns@data$NAVN %in%c("Skjåk","Lom", "Vågå", "Vik","Balestrand","Luster","Årdal", "Vang","Voss","Sogndal","Leikanger", "Høyanger", "Skjåk","Lærdal", "Aurland"),]
writeOGR(southmuns, paste0(rastDir, "Processed/Templates and boundaries"), "South_municipalities", driver="ESRI Shapefile")
southrast <- mask(rastTemplate, southmuns)
writeRaster(southrast, filename=paste0(rastDir, "Processed/Templates and boundaries/South_municipalities.tif"), format = "GTiff", datatype="INT2S")
writeRaster(southrast, filename=paste0(rastDir,"Processed/forMaxent/South_municipalities.asc"), format = "ascii")

###########################
#Clip rasters

setwd(paste0(rastDir, "Processed/forMaxent_prediction/"))
maxentRastList <- list.files(paste0(rastDir, "Processed/forMaxent"), ".asc$", full.names=TRUE)

b <- lapply(1:length(maxentRastList), function(x) {
		currRast <- raster(maxentRastList[x])
		newRast2 <- raster::crop(currRast, extent(alpineshp), filename=basename(maxentRastList[x]), format="ascii")
		return()
		})


		
###########################	
#remove temp raster folder	
unlink(my_tmpdir, recursive = TRUE)
#########################
#END
#########################