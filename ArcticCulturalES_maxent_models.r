#Arctic Cultural Ecosystem Service Modelling
#code by Claire Runge runge@nceas.ucsb.edu
#Project in collaboration with Vera Hausner Associate Professor, Department of Arctic and Marine Biology, The Arctic University of Tromso vera.hausner@uit.no
#This script takes PPGIS data and generates predictive spatial maps of cultural ecosystem services across Norway

##########################################
#Setup

#options(java.parameters = "-Xmx16g" ) #to give java more memory (16G)
options(java.parameters = "-Xmx2g" ) #to give java more memory (2G)
options(stringsAsFactors=FALSE) # turn off automatic factor coersion
# options(scipen=9999)            # turn off scientific notation

library(dismo) #for maxent interface
library(ENMeval) #for AICc
library(raster)
library(foreign) #to read dbfs
library(rgdal) #to read shapefiles
#library(corrplot) #to plot correlation matrix
#library(parallel)

wd <- "C:/Claire/Arctic_Cultural_ES/"
#wd <- "/home/runge/Data/Arctic_Cultural_ES/"
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
#ncore=8

###########################
#PRELIMINARY PROCESSING
#load ces observations (PPGIS data)
markersN <- read.csv(paste0(wd, "PPGIS data/Processed/", "PPGIS_Markers_north_UTM33N_alpine.csv"), header=TRUE)
markersS <- read.csv(paste0(wd, "PPGIS data/Processed/", "PPGIS_Markers_south_UTM33N_alpine.csv"), header=TRUE)
#markersNshp <- readOGR(paste0(wd, "PPGIS data/Original/shps"), "PPGIS_Markers_north_UTM33N_alpine")
#markersSshp <- readOGR(paste0(wd, "PPGIS data/Original/shps"), "PPGIS_Markers_south_UTM33N_alpine")


#CES we are interested in
markersN$species[markersN$species=="hunt/fish"] <- "hunt_fish"
markersS$species[markersS$species=="hunt/fish"] <- "hunt_fish"

cultESlist <- c("biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt_fish", "income", "pasture", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")
#All c("-boats", "-development", "-energy", "-fishing", "-grazing", "-helicopter", "-hunting", "-logging", "-predator", "-roads_atv", "-snowmobiles", "-tourism", "+boats", "+development", "+energy", "+fishing", "+grazing", "+helicopter", "+hunting", "+logging", "+predator", "+roads_atv", "+snowmobiles", "+tourism", "biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt/fish", "income", "otherchange", "pasture", "recnofac", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")

markersNsub <- subset(markersN, species %in% cultESlist)
markersSsub <- subset(markersS, species %in% cultESlist)

#load environmental data
rastDir <- paste0(wd, "Spatial data/Processed/forMaxent/")

envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])


#############################
#remove unwanted variables
varnames <- c("Corrine2012_norway_broadleafforest_3km", "Corrine2012_norway_coniferforest_3km", "Corrine2012_norway_heathshrub_3km",  
"Corrine2012_norway_sparselyvegetated_3km", "Corrine2012_norway_wetland_3km", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbody", "Distance_to_industrialdevelopment_norway","State_commons_norway_binary", "Protected_areas_norway_forothervalues", "Protected_areas_norway_forbiological")# dput(list.files(rastDir, "*.asc$"))
#varnames <- c("Corrine2006_norway_noSea", "Distance_to_River_norway", "Distance_to_Road_norway", "Distance_to_Town_norway", "Distance_to_Coast_norway2","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway") # dput(list.files(rastDir, "*.tif$"))

envStack <- envStack[[varnames]]


#############################
#check correlation of environmental variables
envStackTIF <- raster::stack(list.files(paste0(wd, "Spatial data/Processed/masked/"), "*.tif$", full.names=TRUE))
names(envStackTIF) <- sapply(list.files(paste0(wd, "Spatial data/Processed/masked/"), "*.tif$"), function(x) strsplit(x, "\\.")[[1]][1])
envStackTIF <- envStackTIF[[varnames]]
cors <- layerStats(envStackTIF, 'pearson', na.rm=TRUE)
saveRDS(cors, file=paste0(outDir, "/Correlation of variables/CorrelationofEnvironmentalVariables2.rds")) #readRDS to open
#plot correlation matrix
outPath <-  
png(filename=paste0(outDir, "/Correlation of variables/CorrelationPlotofEnvironmentalVariables2.png"), width=620, height=480)
corrplot(cors[[1]], method='number', type='lower')
dev.off()

#############################
#set up background points
alpineMask <- raster(paste0(rastDir, "Norway_alpine.asc"))
bg <- randomPoints(alpineMask, n=10000)
write.csv(bg, paste0(outDir, "backgroundpoints.csv"), row.names=FALSE)

#plot the background points
bgshp <- SpatialPoints(bg, proj4string=crs(markersNshp))
alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "Norway_alpine")
png(paste0(outDir, "Map_of_backgroundpoints.png"), width=480, height=620)
	plot(alpineshp, border="grey70")
	plot(markersSshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
	plot(markersNshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
	plot(bgshp, col="darkslategray2", add=TRUE, pch=20, cex=0.5, alpha=0.5)
dev.off()
	
#############################
#MAXENT RUNS
#############################
#Conduct sensitivity analysis to evaluate optimal set of regularization multipliers and evaluate model fit #hinge and linear features only
#occ=two column matrix or data.fram of lon and lat in that order

#Sensitivity analysis of features & regularization of North dataset
cl <- makeCluster(ncore)
clusterExport(cl=cl, varlist=c("markersNsub", "bg", "envStack", "rastDir", "outDir", "cultESlist"))
NmodelEval <- parLapply(cl, 1:length(cultESlist), function(x) {
		currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]), c("lon", "lat")]
		currEval <- ENMeval::ENMevaluate(occ=currOcc, bg.coords=bg, env=envStack, 
			RMvalue=seq(0.5, 4, 0.5), #regularization parameters
			fc=c("H", "L"), #hinge & linear
			categoricals=c("Corrine2006_norway_noSea","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway"), 
			n.bg=10000,
			method='randomkfold',
			kfolds=10,
			bin.output=TRUE,
			rasterPreds=TRUE, 
			overlap=TRUE,
			parallel=FALSE)
			return(currEval)
		#save each model
		saveRDS(currEval, file=paste0(outDir, "North model/ENM eval/ENMevalofNmodel_", as.character(cultESlist[x]), ".rds"))
			})
saveRDS(NmodelEval, file=paste0(outDir, "North model/ENM eval/ENMevalofNmodel.rds"))
stopCluster(cl)

bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))
#Sensitivity analysis of features & regularization of South dataset
cl <- makeCluster(ncore)
clusterExport(cl=cl, varlist=c("markersSsub", "bg", "envStack", "rastDir", "outDir", "cultESlist"))
SmodelEval <- parLapply(cl, 1:length(cultESlist), function(x) {
		currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]), c("lon", "lat")]
		currEval <- ENMeval::ENMevaluate(occ=currOcc, bg.coords=bg, env=envStack, 
			RMvalue=seq(0.5, 4, 0.5), #regularization parameters
			fc=c("H", "L"), #hinge & linear
			categoricals=c("Corrine2006_norway_noSea","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway"), 
			n.bg=10000,
			method='randomkfold',
			kfolds=10,
			bin.output=TRUE,
			rasterPreds=TRUE, 
			overlap=TRUE,
			parallel=FALSE)
			return(currEval)
		saveRDS(currEval, file=paste0(outDir, "South model/ENM eval/ENMevalofSmodel_", as.character(cultESlist[x]), ".rds"))
			})
saveRDS(SmodelEval, file=paste0(outDir, "South model/ENM eval/ENMevalofSmodel.rds"))
stopCluster(cl)


#which.min(enmeval_results@results$AICc)
#which.max(enmeval_results@results$Mean.AUC)

#test how well the datasets predict each other using the predictions in the eval files	Warren's D similarity statistic		
NScompare <- lapply(1:length(cultESlist), function(x) {
				currEvalN <- readRDS(paste0(outDir, "North model/ENM eval/ENMevalofNmodel_", as.character(cultESlist[x]), ".rds"))
				currEvalS <- readRDS(paste0(outDir, "South model/ENM eval/ENMevalofSmodel_", as.character(cultESlist[x]), ".rds"))	
				bestN <- which.min(currEvalN@results$AICc) #pick model with lowest AICc
				bestS <- which.min(currEvalS@results$AICc) #pick model with lowest AICc
				NmodelRast <- currEvalN@predictions[[bestN]]
				SmodelRast <- currEvalS@predictions[[bestS]]
				return(nicheOverlap(NmodelRast, SmodelRast, mask=FALSE, stat='D'))
			}

########################
#Response curves & jackknife of N model
# to do:
# fill in regularisation parameters 'beta_hinge=4'

#set up cluster
# cl <- makeCluster(ncore)
# clusterExport(cl=cl, varlist=c("markersNsub", "markersNsub", "bg", "envStack", "rastDir", "outDir", "cultESlist"))

#run maxent model across all cultES
#Nmodel <- parLapply(cl, 1:length(cultESlist), function(x) {

#set up background points
#bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
Nmask <- raster(paste0(rastDir, "North_municipalities.asc"))
bg <- randomPoints(Nmask, n=10000)
write.csv(bg, paste0(outDir, "backgroundpoints_north.csv"), row.names=FALSE)

#set up variables
varnames <- c("North_municipalities", "Corrine2012_norway_broadleafforest_3km", "Corrine2012_norway_coniferforest_3km", "Corrine2012_norway_heathshrub_3km",  
"Corrine2012_norway_sparselyvegetated_3km", "Corrine2012_norway_wetland_3km", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbody", "Distance_to_industrialdisturbance_norway","State_commons_norway_binary", "Protected_areas_norway_forothervalues")
envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
envStack <- envStack[[varnames]]

#run model
Nmodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Response curves and jacknife2/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Protected_areas_norway_forothervalues", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})

###run model for biological & undisturbnature values
cultESlist <- c("biological", "undisturbnature")

#set up variables
varnames <- c("North_municipalities", "Corrine2012_norway_broadleafforest_3km", "Corrine2012_norway_coniferforest_3km", "Corrine2012_norway_heathshrub_3km",  
"Corrine2012_norway_sparselyvegetated_3km", "Corrine2012_norway_wetland_3km", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbody", "Distance_to_industrialdisturbance_norway","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
envStack <- envStack[[varnames]]

bg <- read.csv(paste0(outDir, "backgroundpoints_north.csv"))
Nmodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Response curves and jacknife2/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
#stopCluster(cl)

#'beta_hinge=4'

########################
#Response curves & jackknife of S model
# to do:
# fill in regularisation parameters 'beta_hinge=4'

#set up cluster
#cl <- makeCluster(ncore)
#clusterExport(cl=cl, varlist=c("markersNsub", "markersNsub", "bg", "envStack", "rastDir", "outDir", "cultESlist"))

#run maxent model for all cultES
#Smodel <- parLapply(cl, 1:length(cultESlist), function(x) {

#bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
#Smodel <- lapply(6:length(cultESlist), function(x) {
	

#set up background points
#bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
Smask <- raster(paste0(rastDir, "South_municipalities.asc"))
bg <- randomPoints(Smask, n=10000)
write.csv(bg, paste0(outDir, "backgroundpoints_south.csv"), row.names=FALSE)

#set up variables
varnames <- c("South_municipalities", "Corrine2012_norway_broadleafforest_3km", "Corrine2012_norway_coniferforest_3km", "Corrine2012_norway_heathshrub_3km",  
"Corrine2012_norway_sparselyvegetated_3km", "Corrine2012_norway_wetland_3km", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbody", "Distance_to_industrialdisturbance_norway","State_commons_norway_binary", "Protected_areas_norway_forothervalues")
envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
envStack <- envStack[[varnames]]

#run model
Smodel <- lapply(2:(length(cultESlist)-1), function(x) {	
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Response curves and jacknife2/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Protected_areas_norway_forothervalues", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})

			
###run model for biological & undisturbnature values
cultESlist <- c("biological", "undisturbnature")

#set up variables
varnames <- c("South_municipalities", "Corrine2012_norway_broadleafforest_3km", "Corrine2012_norway_coniferforest_3km", "Corrine2012_norway_heathshrub_3km",  
"Corrine2012_norway_sparselyvegetated_3km", "Corrine2012_norway_wetland_3km", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbody", "Distance_to_industrialdisturbance_norway","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
envStack <- envStack[[varnames]]

bg <- read.csv(paste0(outDir, "backgroundpoints_south.csv"))
#run model
Smodel <- lapply(1:length(cultESlist), function(x) {	
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Response curves and jacknife2/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
#stopCluster(cl)

########################
#Response curves & jackknife of model created with both north and south data 
# to do:
# fill in regularisation parameters 'beta_hinge=4'
markersCombined <- rbind(markersNsub, markersSsub)

#set up cluster
# cl <- makeCluster(ncore)
# clusterExport(cl=cl, varlist=c("markersCombined", "bg", "envStack", "rastDir", "outDir", "cultESlist"))

# #run maxent model for all cultES
# ComboModel <- parLapply(l, 1:length(cultESlist), function(x) {

bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
ComboModel <- lapply(1:length(cultESlist), function(x) {
ComboModel <- lapply(6:8, function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Response curves and jacknife/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Corrine2006_norway_noSea","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc'), path=currOutPath)
			
				
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
				
#stopCluster(cl)

########################
#Assuming these models all look ok, we now create a model with all the data
########################
#Base run of N model

#set up cluster
# cl <- makeCluster(ncore)
# clusterExport(cl=cl, varlist=c("markersNsub", "markersNsub", "bg", "envStack", "rastDir", "outDir", "cultESlist"))

#run maxent model across all cultES
#Nbasemodel <- parLapply(cl, 1:length(cultESlist), function(x) {

bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
Nbasemodel <- lapply(6:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Base run/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Corrine2006_norway_noSea","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
#stopCluster(cl)

#'beta_hinge=4'

########################
#Base run of S model

#set up cluster
#cl <- makeCluster(ncore)
#clusterExport(cl=cl, varlist=c("markersNsub", "markersNsub", "bg", "envStack", "rastDir", "outDir", "cultESlist"))

#run maxent model for all cultES
#Sbasemodel <- parLapply(cl, 1:length(cultESlist), function(x) {

bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
Sbasemodel <- lapply(6:length(cultESlist), function(x) {
			
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Base run/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Corrine2006_norway_noSea","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return()
			})
			
#stopCluster(cl)

########################
#Base run of model created with both north and south data 

markersCombined <- rbind(markersNsub, markersSsub)

#set up cluster
# cl <- makeCluster(ncore)
# clusterExport(cl=cl, varlist=c("markersCombined", "bg", "envStack", "rastDir", "outDir", "cultESlist"))

# #run maxent model for all cultES
# CombobaseModel <- parLapply(l, 1:length(cultESlist), function(x) {

bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
CombobaseModel <- lapply(1:length(cultESlist), function(x) {
CombobaseModel <- lapply(6:8, function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Base run/", as.character(cultESlist[x]))
			dir.create(currOutPath)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("Corrine2006_norway_noSea","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25'), path=currOutPath)
			
				
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
				
#stopCluster(cl)


######################
#Predict models
alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "Norway_alpine")

rastDir <- paste0(wd, "Spatial data/Processed/forMaxent_prediction/")
varnames <- c("Corrine2006_norway_noSea", "Distance_to_River_norway", "Distance_to_Road_norway", "Distance_to_Town_norway", "Distance_to_Coast_norway2","Ecological_areas_norway", "Protected_areas_norway", "State_commons_norway") # dput(list.files(rastDir, "*.tif$"))

cultESlist <- c("biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt_fish", "income", "pasture", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")

envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
envStack <- envStack[[varnames]] #remove unwanted variables

#using a for loop because of memory allocataion issues
for (y in c("Combined model", "North model", "South model"){
	for (x in seq_along(cultESlist)){
	
		currOutPath <- paste0(outDir, y, "/Base run/", as.character(cultESlist[x]))	
		#dir.create(currOutPath)
		currMod <- readRDS(paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
		print(paste0("starting prediction ", y, cultESlist[x]))
		print(Sys.time())
		#make predictive maps
		currMap <- dismo::predict(currMod, envStack, ext=extent(alpineshp), args=c('outputformat=logistic', "doclamp=TRUE"), progress='text', filename=paste0(currOutPath,"/", cultESlist[x], "_basemodel.tif"), format="GTiff") 
		#currMap <- dismo::predict(currMod, envStack, ext=extent(alpineshp), args=c('outputformat=logistic', "doclamp=TRUE", "writeclampgrid=TRUE", "writemess=TRUE"), progress='text', filename=paste0(currOutPath,"/", cultESlist[x], "_basemodel.tif"), format="GTiff") 
		print(paste0("finished prediction ", y, cultESlist[x]))
		print(Sys.time())
		rm(currMap) #trying to clear memory issues
		
		}
}

######################
#test N model against S data
bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
for (x in 1:length(cultESlist)){
	currOutPath <- paste0(outDir, "/North model/Test NS data/", as.character(cultESlist[x]))
	dir.create(currOutPath)
	currInp <- 	paste0(outDir, "/North model/Base run/", as.character(cultESlist[x]))
	currTestData <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
	currMod <- readRDS(paste0(currInp, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
	
	e1 <- evaluate(currMod, p=currTestData, a=bg, x=envStack)
}
#'testsamplesfile=xxx'

######################
#test S model against N data
bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
for (x in 1:length(cultESlist)){
	currOutPath <- paste0(outDir, "/South model/Test NS data/", as.character(cultESlist[x]))
	dir.create(currOutPath)
	currInp <- 	paste0(outDir, "/South model/Base run/", as.character(cultESlist[x]))
	currTestData <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
	currMod <- readRDS(paste0(currInp, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
	
	e1 <- evaluate(currMod, p=currTestData, a=bg, x=envStack)
}


#sensitivity to background - use randomseed=TRUE
currArgs <- c('randomseed=TRUE', 'replicates=100', 'outputgrids=false') #check how long each run takes
alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "Norway_alpine")



###LAMBDAS
#helpful discussion on extracting lambdas
# https://groups.google.com/forum/#!topic/maxent/iwd8kR_J0tk
#only variables with value that is not zero in the second column were used in the model
#first column: variable and feature type
#second column: lambda value
#third column: min value
#fourth column: max value
#
  # # read lambdas file
  # rf <- read.table(file.path(curpath, 'species.lambdas'), sep=',', fill=TRUE)
  # # get variables used in model - no 0 in 2nd column)
  # rf[!is.na(rf[3]) & rf[2] != 0,]







#maxent(envStack, occurencedata, nbg=number of background points, factors=names of layers considered as categorical, args=arguments passed to maxent, removeDuplicates=FALSE, path=place to store output files)			
me <- maxent(predictors, occtrain, factors='biome', path="C:/Claire/Arctic_Cultural_ES/Maxent runs/test", args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE'))
pmap <- predict(me, predictors, path="C:/Claire/Arctic_Cultural_ES/Maxent runs/test/test.tif", progress="text", format="GTiff", args=c('outputformat=cumulative', 'perspeciesresults=true', 'threads=3'))




#Check N model 10-fold cross validation, with response curves & jackknife 
currArgs <- c('responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false')

#S model Response curves & jackknife with k-fold cross validation



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
