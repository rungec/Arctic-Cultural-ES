#Arctic Cultural Ecosystem Service Modelling
#code by Claire Runge runge@nceas.ucsb.edu
#Project in collaboration with Vera Hausner Associate Professor, Department of Arctic and Marine Biology, The Arctic University of Tromso vera.hausner@uit.no
#This script takes PPGIS data and generates predictive spatial maps of cultural ecosystem services across Norway

##########################################
#Setup

options(java.parameters = "-Xmx2g" ) #to give java more memory (2G)
options(stringsAsFactors=FALSE) # turn off automatic factor coersion
# options(scipen=9999)            # turn off scientific notation

library(dismo) #for maxent interface
library(ENMeval) #for AICc
library(raster)
library(foreign) #to read dbfs
#library(rgdal) #to read shapefiles

wd <- "C:/Claire/Arctic_Cultural_ES/"
setwd(wd)
outDir <- paste0(wd, "Maxent runs/")

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

#number of core to use
ncore=3

###########################
#PRELIMINARY PROCESSING
#load ces observations (PPGIS data)
markersN <- read.csv(paste0(wd, "PPGIS Data/Processed/", "PPGIS_Markers_north_UTM33N_alpine.csv"), header=TRUE)
markersS <- read.csv(paste0(wd, "PPGIS Data/Processed/", "PPGIS_Markers_south_UTM33N_alpine.csv"), header=TRUE)

#CES we are interested in
cultESlist <- c("biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt/fish", "income", "pasture", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")
#All c("-boats", "-development", "-energy", "-fishing", "-grazing", "-helicopter", "-hunting", "-logging", "-predator", "-roads_atv", "-snowmobiles", "-tourism", "+boats", "+development", "+energy", "+fishing", "+grazing", "+helicopter", "+hunting", "+logging", "+predator", "+roads_atv", "+snowmobiles", "+tourism", "biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt/fish", "income", "otherchange", "pasture", "recnofac", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")

markersNsub <- subset(markersN, species %in% cultESlist)
markersSsub <- subset(markersS, species %in% cultESlist)

#load environmental data
rastDir <- paste0(wd, "Spatial data/Processed/forMaxent/")
varnames <- c("Corrine2006_norway", "Corrine2006_norway_noSea","Distance_to_Coast_norway.tif", "Distance_to_House_norway", "Distance_to_River_norway.tif", "Distance_to_Road_norway", "Distance_to_Town_norway.tif", "Ecological_areas_norway.", "Protected_areas_norway.tif", "State_commons_norway") # dput(list.files(rastDir, "*.tif$"))

envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))


#############################
#check correlation of environmental variables
cors <- layerStats(envStack, 'pearson', na.rm=TRUE)

#############################
#MAXENT RUNS
#############################

args <-  outputdirectory=, samplesfile)
#maxent(envStack, occurencedata, nbg=number of background points, factors=names of layers considered as categorical, args=arguments passed to maxent, removeDuplicates=FALSE, path=place to store output files)

#Conduct sensitivity analysis to evaluate optimal set of regularization multipliers and evaluate model fit #hinge only
#occ=two column matrix or data.fram of lon and lat in that order

NmodelEval <- lapply(1:length(cultESlist), function(x) {
		currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]), c("lon", "lat")]
		currEval <- ENMevaluate(occ=currOcc, env=envStack[[2:length(envStack@layers)]], 
			RMvalue=seq(0.5, 4, 0.5), #regularization parameters
			fc=c("H", "L"), #hinge & linear
			categoricals=c("Corrine2006_norway_noSea","Ecological_areas_norway.", "Protected_areas_norway.tif", "State_commons_norway"), 
			n.bg=10000,
			method='randomkfold',
			overlap=FALSE,
			kfolds=10,
			bin.output=TRUE,
			clamp=TRUE,
			rasterPreds=FALSE, 
			parallel=FALSE)
			return(currEval)
			})

SmodelEval <- lapply(1:length(cultESlist), function(x) {
		currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]), c("lon", "lat")]
		currEval <- ENMevaluate(occ=currOcc, env=envStack[[2:length(envStack@layers)]], 
			RMvalue=seq(0.5, 4, 0.5), #regularization parameters
			fc=c("H", "L"), #hinge & linear
			categoricals=c("Corrine2006_norway_noSea","Ecological_areas_norway.", "Protected_areas_norway.tif", "State_commons_norway"), 
			n.bg=10000,
			method='randomkfold',
			overlap=FALSE,
			kfolds=10,
			bin.output=TRUE,
			clamp=TRUE,
			rasterPreds=FALSE, 
			parallel=FALSE)
			return(currEval)
			})
			
#Test how well the different datasets predict each otherchange
ENMevaluate(occ, envStack[2:length(envStack@layers)], 
			bg.coords=, 
			occ.grp=, 
			bg.grp=, 
			RMvalue=seq(0.5, 4, 0.5), #use the values found above
			fc=c("H"), 
			categoricals=c("Corrine2006_norway_noSea","Ecological_areas_norway.", "Protected_areas_norway.tif", "State_commons_norway"), 
			n.bg=10000,
			method='user',
			overlap=TRUE,
			bin.output=TRUE,
			clamp=TRUE,
			rasterPreds=TRUE, 
			parallel=TRUE,
			numCores=3)

#Base run of N model using all data
Nmodel <- lapply(1:length(cultESlist), function(x) {
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),]
			currOutPath <- paste0(outDir, "North model/Base run/", as.character(cultESlist[x]))
			#run the model
			currMod <- maxent(envStack, currOcc, nbg=10000, factors=c(""), path=currOutPath, args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'threads=3'), path=currOutPath)
			#make predictive maps
			currMap <- predict(currMod, envStack, args=c('outputformat=cumulative', 'threads=3')), 
			writeRaster(currMap, filename=paste0(currOutPath,"/", cultESlist[x], "_basemodel.tif"), format="GTiff")
			
			names(currMod) <- cultESlist[x]
			return(currMod)
			})
			
me <- maxent(predictors, occtrain, factors='biome', path="C:/Claire/Arctic_Cultural_ES/Maxent runs/test", args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE'))
pmap <- predict(me, predictors, path="C:/Claire/Arctic_Cultural_ES/Maxent runs/test/test.tif", progress="text", format="GTiff", args=c('outputformat=cumulative', 'perspeciesresults=true', 'threads=3'))



make.args(RMvalues = seq(0.5, 4, 0.5), fc = c("H"), labels = FALSE) #hinge only

#Base run of S model using all data
Smodel <- lapply(1:)

#Check N model 10-fold cross validation, with response curves & jackknife 
currArgs <- c('responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false')

#S model Response curves & jackknife with k-fold cross validation

#sensitivity to background - use randomseed=TRUE
currArgs <- c('randomseed=TRUE', 'replicates=100', 'outputgrids=false') #check how long each run takes

#compare N-S predictions
#run through the different CultES 'species'
NmodelSdata <- lapply(1:length(cultESlist), function(x) {
		currCultES <- cultESlist[x]
		evaluate(Nmodel[x], p=markersSsub[markersSsub$species==currCultES,], x=envStack) 
		})

SmodelNdata <- lapply(1:length(cultESlist), function(x) {
		currCultES <- cultESlist[x]
		evaluate(Smodel[x], p=markersNsub[markersNsub$species==currCultES,], x=envStack) 
		})

a <- lapply(list.files(rastDir, "*.tif$", full.names=TRUE), function(x) raster(x))
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

		
###########################	
#remove temp raster folder	
unlink(my_tmpdir, recursive = TRUE)
#########################
#END
#########################
