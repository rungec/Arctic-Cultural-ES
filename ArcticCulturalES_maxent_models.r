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
library(corrplot) #to plot correlation matrix
#library(parallel)

wd <- "C:/Claire/Arctic_Cultural_ES/"
#wd <- "/home/runge/Data/Arctic_Cultural_ES/"
setwd(wd)

rastDir <- paste0(wd, "Spatial data/Processed/3_forMaxent_asc/")
outDir <- paste0(wd, "Maxent runs/")

currDate <- paste(strsplit(as.character(Sys.Date()), "-")[[1]], collapse="") #formatting the date for folder naming

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
markersS <- read.csv(paste0(wd, "PPGIS data/Processed/", "PPGIS_Markers_south_UTM33N_municipal.csv"), header=TRUE)
markersNshp <- readOGR(paste0(wd, "PPGIS data/Processed/shps"), "PPGIS_Markers_north_UTM33N_alpine")
markersSshp <- readOGR(paste0(wd, "PPGIS data/Processed/shps"), "PPGIS_Markers_south_UTM33N_municipal")


#CES we are interested in
markersN$species[markersN$species=="hunt/fish"] <- "hunt_fish"
markersS$species[markersS$species=="hunt/fish"] <- "hunt_fish"

cultESlist <- c("biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt_fish", "income", "pasture", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")
#All c("-boats", "-development", "-energy", "-fishing", "-grazing", "-helicopter", "-hunting", "-logging", "-predator", "-roads_atv", "-snowmobiles", "-tourism", "+boats", "+development", "+energy", "+fishing", "+grazing", "+helicopter", "+hunting", "+logging", "+predator", "+roads_atv", "+snowmobiles", "+tourism", "biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt/fish", "income", "otherchange", "pasture", "recnofac", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")

markersNsub <- subset(markersN, species %in% cultESlist)
markersSsub <- subset(markersS, species %in% cultESlist)

#load environmental data
envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])

#load bias grid
biasGrid <- raster(paste0(basename(rastDir), "/Bias_grids/BiasGrid_distancetoroad_nowater.asc"))

# #############################
# #remove unwanted variables
#included 
#exclude distance to road from model with bias grid
#replace "Governance_plus_protectedareas_norway" with "State_commons_norway_binary" & "Protected_areas_norway_forbiological" for biological & undisturbnature
varnames <- c("Norway_alpine", "North_municipalities_alpine", "South_municipalities, AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway","Distance_to_Road_norway", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway", "State_commons_norway_binary",  "Protected_areas_norway_forbiological")
envStack <- envStack[[varnames]]
length(varnames)

#not included
# c("Corrine2012_norway_wetland_1km","Distance_to_waterbodies_norway2", , "State_commons_norway_binary", "Distance_to_Coast_norway2" "Percent_Waterbodies_majorriversandlakesbiggerthan2ha_1_norway",
# "Protected_areas_norway_forbiological.tif")
# dput(list.files(rastDir, "*.asc$"))
# dput(list.files(rastDir, "*.tif$"))

#############################
###CHECK CORRELATION OF VARIABLES
#############################

#check correlation of environmental variables
envStackTIF <- raster::stack(list.files(paste0(wd, "Spatial data/Processed/2b_masked_no_water/"), "*.tif$", full.names=TRUE))
names(envStackTIF) <- sapply(list.files(paste0(wd, "Spatial data/Processed/2b_masked_no_water/"), "*.tif$"), function(x) strsplit(x, "\\.")[[1]][1])
#envStackTIF <- envStackTIF[[varnames]]
cors <- layerStats(envStackTIF, 'pearson', na.rm=TRUE)
saveRDS(cors, file=paste0(outDir, "/Correlation of variables/CorrelationofEnvironmentalVariables_nowater2.rds")) #readRDS to open

#plot correlation matrix
png(filename=paste0(outDir, "/Correlation of variables/CorrelationPlotofEnvironmentalVariables_nowater2.png"), width=1860, height=1440)
corrplot::corrplot(cors[[1]], method='number', type='lower')
dev.off()

#############################
###BACKGROUND
#############################

#set up background points for the combined NS models
	alpineMask <- raster(paste0(rastDir, "Norway_alpine.asc"))
	alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "Norway_alpine")
	bg <- randomPoints(alpineMask, n=10000, ext=extent(alpineshp), tryf=100)
	write.csv(bg, paste0(dirname(rastDir), "/Background_points/backgroundpoints_wholeregion.csv"), row.names=FALSE)

	#set up background points for the combined NS models with bias
	alpineMask <- raster(paste0(dirname(rastDir), "/Bias_grids/BiasGrid_distancetoroad_Norway_alpine.tif"))
	alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "Norway_alpine")
	bg <- randomPoints(alpineMask, n=10000, ext=extent(alpineshp), tryf=100, prob=TRUE)
	write.csv(bg, paste0(dirname(rastDir), "/Background_points/backgroundpoints_wholeregion_bias.csv"), row.names=FALSE)

	#plot the background points
		bgshp <- SpatialPoints(bg, proj4string=crs(markersNshp))
		png(paste0(dirname(rastDir), "/Background_points/Map_of_backgroundpoints_wholeregion.png"), width=480, height=620)
			plot(alpineshp, border="grey70")
			plot(markersSshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
			plot(markersNshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
			plot(bgshp, col="darkslategray2", add=TRUE, pch=20, cex=0.5, alpha=0.5)
		dev.off()

#set up background points for North models
	Nmask <- raster(paste0(dirname(rastDir), "/Templates and boundaries/North_municipalities_alpine.tif"))	
	alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "North_municipalities_alpine")
	bgN <- randomPoints(Nmask, n=10000, ext=extent(alpineshp), tryf=100)
	write.csv(bgN, paste0(dirname(rastDir), "/Background_points/backgroundpoints_north.csv"), row.names=FALSE)
	
#set up background points for North models with bias
	Nmask <- raster(paste0(dirname(rastDir), "/Bias_grids/BiasGrid_distancetoroad_North_municipalities_alpine.tif"))	
	alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "North_municipalities_alpine")
	bgN <- randomPoints(Nmask, n=10000, ext=extent(alpineshp), tryf=100, prob=TRUE)
	write.csv(bgN, paste0(dirname(rastDir), "/Background_points/backgroundpoints_north_bias.csv"), row.names=FALSE)
	
	#plot the background points
	bgshp <- SpatialPoints(bgN, proj4string=crs(markersNshp))
	png(paste0(dirname(rastDir), "/Background_points/Map_of_backgroundpoints_north.png"), width=480, height=620)
		plot(alpineshp, border="grey70")
		plot(markersSshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
		plot(markersNshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
		plot(bgshp, col="darkslategray2", add=TRUE, pch=20, cex=0.5, alpha=0.5)
	dev.off()

#set up background points for South models
	Smask <- raster(paste0(dirname(rastDir), "/Templates and boundaries/South_municipalities.tif"))
	alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "South_municipalities")
	bgS <- randomPoints(Smask, n=10000, ext=extent(alpineshp), tryf=100)
	write.csv(bgS, paste0(dirname(rastDir), "/Background_points/backgroundpoints_south.csv"), row.names=FALSE)	
	
	#set up background points for South models with bias
	Smask <- raster(paste0(dirname(rastDir), "/Bias_grids/BiasGrid_distancetoroad_South_municipalities.tif"))
	alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "South_municipalities")
	bgS <- randomPoints(Smask, n=10000, ext=extent(alpineshp), tryf=100, prob=TRUE)
	write.csv(bgS, paste0(dirname(rastDir), "/Background_points/backgroundpoints_south_bias.csv"), row.names=FALSE)	
	
	#plot the background points
	bgshp <- SpatialPoints(bgS, proj4string=crs(markersNshp))
	png(paste0(dirname(rastDir), "/Background_points/Map_of_backgroundpoints_south.png"), width=480, height=620)
		plot(alpineshp, border="grey70")
		plot(markersSshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
		plot(markersNshp, col="darkslateblue", add=TRUE, pch=20, cex=0.5, alpha=0.5)
		plot(bgshp, col="darkslategray2", add=TRUE, pch=20, cex=0.5, alpha=0.5)
	dev.off()
	
#############################
#CHOOSE REGULARIZATION MULTIPLIER
#############################
#Conduct sensitivity analysis to evaluate optimal set of regularization multipliers and evaluate model fit #hinge and linear features only
#occ=two column matrix or data.fram of lon and lat in that order
#Sensitivity analysis of features & regularization of recreation with North dataset

#set up variables
	x=9 #recreation #x=10 #scenic
	currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]), c("lon", "lat")]
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_north.csv"))
	currEnvStack <- envStack[[varnames[1:(length(varnames)-2)]]]

#run evaluation
	currEval <- ENMeval::ENMevaluate(occ=currOcc, bg.coords=bg, env=currEnvStack, 
			RMvalue=seq(0.5, 2.5, 0.5), #regularization parameters
			#fc=c("H", "L"), #hinge & linear
			fc=c("H"), #hinge only
			categoricals=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Governance_plus_protectedareas_norway"), 
			n.bg=10000,
			method='randomkfold',
			kfolds=10,
			bin.output=TRUE, #appends evaluations metrics for each evaluation bin to results table
			rasterPreds=FALSE, #if TRUE predict each model across input variables (for AICc)
			#overlap=TRUE, #pairwise metric of niche overlap
			parallel=TRUE, numCores=6)
			
		#save evaluation
		saveRDS(currEval, file=paste0(outDir, "North model/ENM eval/ENMevalofNmodel_", as.character(cultESlist[x]), ".rds"))
		write.csv(currEval@results, paste0(outDir, "North model/ENM eval/ENMevalofNmodel_", as.character(cultESlist[x]), ".csv"), row.names=FALSE)
			
#look at it
	data(currEval)
	currEval@results
	which.min(currEval@results$AICc)
	which.max(currEval@results$Mean.AUC)

### Plot prediction AICc
	png(paste0(outDir, "North model/ENM eval/PlotofAICc_forregularizationselection_", as.character(cultESlist[x]), ".png"))
		plot(currEval@predictions[[which (currEval@results$delta.AICc == 0) ]])
		points(currEval@occ.pts, pch=21, bg=currEval@occ.grp)
	dev.off()
	
### Plot model AUC
#Mean.AUC.DIFF gives difference between test & training AUCs averaged across 10 reps - expected to be high for models overfit to the training data
#full.AUC is AUC of full model (all data, no bins) - evaluation of degree to which a model can successfully discriminate occurrence from background localities
	png(paste0(outDir, "North model/ENM eval/PlotofAUC_forregularizationselection_", as.character(cultESlist[x]), ".png"))
		par(mfrow=c(1,2))
		plot(currEval@results$full.AUC, xaxt='n', xlab="Settings", ylab="AUC", main="Test AUC full model")
		axis(1, at=1:10, labels=currEval@results$settings, las=2)
				
		plot(currEval@results$Mean.AUC.DIFF, xaxt='n', xlab="Settings", ylab="Difference", main="Test vs Train AUC")
		axis(1, at=1:10, labels=currEval@results$settings, las=2)
	dev.off()


########################
###MAXENT RUNS NORTH MODEL
########################
#TEST MODEL
#Response curves, jackknife & 10-fold cross-validation of N model
#set up variables

	varnames <- c("North_municipalities_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_north.csv"))

#run model
Nmodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Response curves and jacknife ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
###RUN BASE MODEL of N model
#Assuming these models all look ok, we now create a model with all the data which we then use for prediction

#run maxent model across all cultES
Nbasemodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Base_run_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})

###RUN N MODELS biological & undisturbnature values
#set up variables for biological & undisturbnature values
	cultESlist <- c("biological", "undisturbnature")
	varnames <- c("North_municipalities_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_north.csv"))

#run model for biological & undisturbnature values
Nmodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Response curves and jacknife ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})

###RUN BASE MODEL of N model for biological & undisturbnature
#Assuming these models all look ok, we now create a model with all the data which we then use for prediction

#run maxent model across all cultES
Nbasemodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Base_run_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			

########################
###MAXENT RUNS NORTH MODEL WITH BIAS GRID
########################

###SET UP N DATA WITH BIAS GRID
#set up variables
	varnames <- c("North_municipalities_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway",  "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_north.csv"))

#TEST N MODEL WITH BIAS GRID
#Response curves, jackknife & 10-fold cross-validation of N model with bias grid
Nmodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Response curves and jacknife bias grid ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
###RUN BASE MODEL of N model WITH BIAS GRID
#Assuming these models all look ok, we now create a model with all the data which we then use for prediction

#run maxent model across all cultES
Nbasemodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Base_run_bias_grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})

#SET UP DATA N model WITH BIAS GRID for biological & undisturbnature values
	cultESlist <- c("biological", "undisturbnature")
	varnames <- c("North_municipalities_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_north.csv"))

##TEST N MODEL WITH BIAS GRID for biological & undisturbnature values
Nmodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Response curves and jacknife bias grid ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})

###RUN BASE MODEL of N model WITH BIAS GRID for biological & undisturbnature
#Assuming these models all look ok, we now create a model with all the data which we then use for prediction
#run maxent model across biological & undisturbnature
Nbasemodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "North model/Base_run_bias_grid_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			

########################
###MAXENT RUNS SOUTH MODEL
########################

#set up variables
	varnames <- c("South_municipalities", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","Governance_plus_protectedareas_norway")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_south.csv"))

###TEST RUN S MODEL
#Response curves, jackknife & 10-fold cross-validation of S model
Smodel <- lapply(2:(length(cultESlist)-1), function(x) {
	
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Response curves and jacknife ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
###RUN BASE MODEL of S model
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
Sbasemodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Base_run_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return()
			})
			#run maxent model across all cultES

###SET UP DATA S MODEL for biological & undisturbnature values
#set up variables for biological & undisturbnature values
	cultESlist <- c("biological", "undisturbnature")
	varnames <- c("South_municipalities", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_south.csv"))

###TEST RUN S MODEL for biological & undisturbnature values
#Response curves, jackknife & 10-fold cross-validation of S model
Smodel <- lapply(1:length(cultESlist), function(x) {	
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Response curves and jacknife ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "State_commons_norway_binary",  "Protected_areas_norway_forbiological"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
###RUN BASE MODEL of S model for biological & undisturbnature
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
Sbasemodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Base_run_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			

########################
###MAXENT RUNS SOUTH MODEL WITH BIAS GRID
########################
###SET UP DATA
#Response curves, jackknife & 10-fold cross-validation of S model with bias grid
#set up variables
	varnames <- c("South_municipalities", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway",  "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","Governance_plus_protectedareas_norway")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_south.csv"))

###TEST MODEL WITH BIAS GRID
#run model with k-fold cross validation & jackknife
Smodel <- lapply(2:(length(cultESlist)-1), function(x) {
	
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Response curves and jacknife bias grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
			
###RUN BASE MODEL of S model WITH BIAS GRID
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
Sbasemodel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Base_run_bias_grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return()
			})
			#run maxent model across all cultES

###SET UP DATA FOR S MODEL biological & undisturbnature values WITH BIAS GRID
#set up variables for biological & undisturbnature values
	cultESlist <- c("biological", "undisturbnature")
	varnames <- c("South_municipalities", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway",  "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_south.csv"))

###TEST MODEL FOR S MODEL biological & undisturbnature values WITH BIAS GRID
#run model for biological & undisturbnature values
Smodel <- lapply(1:length(cultESlist), function(x) {	
			currOcc <- markersSsub[markersSsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Response curves and jacknife bias grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "State_commons_norway_binary",  "Protected_areas_norway_forbiological"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
###RUN BASE MODEL WITH BIAS GRID of S model for biological & undisturbnature
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
Sbasemodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersNsub[markersNsub$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "South model/Base_run_bias_grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
########################
###MAXENT RUNS NORTHSOUTH COMBINED MODEL
########################
###SET UP DATA for NORTHSOUTH COMBINED MODEL
	varnames <- c("Norway_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	markersCombined <- rbind(markersNsub, markersSsub)
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_wholeregion.csv"))

###TEST MODEL	
#Response curves, jackknife & 10-fold cross-validation of model created with both north and south data 
ComboModel <- lapply(2:length(cultESlist-1), function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Response curves and jacknife ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1'), path=currOutPath)
				
			#save model
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
###RUN BASE MODEL of model created with both north and south data
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
CombobaseModel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Base_run_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)

			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return()
			})

###SET UP DATA FOR biological & undisturbnature values NORTHSOUTH COMBINED MODEL
	cultESlist <- c("biological", "undisturbnature")
	varnames <- c("Norway_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_wholeregion.csv"))

###TEST MODEL FOR biological & undisturbnature values NORTHSOUTH COMBINED MODEL
#run model for biological & undisturbnature values
Combomodel <- lapply(1:length(cultESlist), function(x) {	
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Response curves and jacknife ", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "State_commons_norway_binary",  "Protected_areas_norway_forbiological"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
###RUN BASE MODEL of Combo model for biological & undisturbnature NORTHSOUTH COMBINED MODEL
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
Combobasemodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Base_run_", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
			
########################
###MAXENT RUNS NORTHSOUTH COMBINED MODEL WITH BIAS GRID
########################
###SET UP DATA WITH BIAS GRID
	varnames <- c("Norway_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway",  "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	markersCombined <- rbind(markersNsub, markersSsub)
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_wholeregion.csv"))

###TEST NORTHSOUTH COMBINED MODEL WITH BIAS GRID
#Response curves, jackknife & 10-fold cross-validation of model created with both north and south data and bias grid
ComboModel <- lapply(2:length(cultESlist-1), function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Response curves and jacknife bias grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
				
			#save model
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
				
###RUN BASE NORTHSOUTH COMBINED MODEL WITH BIAS GRID
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
CombobaseModel <- lapply(2:(length(cultESlist)-1), function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Base_run_bias_grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)

			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Governance_plus_protectedareas_norway"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return()
			})

###SET UP DATA FOR biological & undisturbnature values NORTHSOUTH COMBINED MODEL WITH BIAS GRID
	cultESlist <- c("biological", "undisturbnature")
	varnames <- c("Norway_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway",  "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater","State_commons_norway_binary",  "Protected_areas_norway_forbiological")
	envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
	names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
	envStack <- envStack[[varnames]]
	names(envStack)[[1]] <- "mask"
	bg <- read.csv(paste0(dirname(rastDir), "/Background_points/backgroundpoints_wholeregion.csv"))

###TEST MODEL FOR biological & undisturbnature values NORTHSOUTH COMBINED MODEL WITH BIAS GRID
#run model for biological & undisturbnature values
Combomodel <- lapply(1:length(cultESlist), function(x) {	
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Response curves and jacknife bias grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c("AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "State_commons_norway_binary",  "Protected_areas_norway_forbiological"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE',  'responsecurves=TRUE', 'jackknife=TRUE', 'replicates=10', 'outputgrids=false', 'outputfiletype=asc', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	
###RUN BASE MODEL of Combo model for biological & undisturbnature NORTHSOUTH COMBINED MODEL WITH BIAS GRID
#Assuming these models all look ok, we now create a model with all the data which we later use for prediction
Combobasemodel <- lapply(1:length(cultESlist), function(x) {
			
			currOcc <- markersCombined[markersCombined$species==as.character(cultESlist[x]),c("lon", "lat")]
			currOutPath <- paste0(outDir, "Combined model/Base_run_bias_grid", currDate, "/", as.character(cultESlist[x]))
			dir.create(currOutPath, recursive=TRUE)
			
			#run the model
			currMod <- dismo::maxent(envStack, currOcc, a=bg, factors=c( "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water","Protected_areas_norway_forbiological", "State_commons_norway_binary"), args=c('hinge=TRUE', 'linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'autofeature=FALSE', 'writeplotdata=TRUE', 'responsecurves=TRUE', 'jackknife=TRUE', 'randomtestpoints=25', 'beta_hinge=1', 'biasfile=biasGrid'), path=currOutPath)
			
			#save model
			#names(currMod) <- cultESlist[x]
			saveRDS(currMod, file=paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
			return(currMod)
			})
	

########################
###PREDICT MODELS ONTO ENVIRONMENTAL VARIABLES
########################
#Predict models
alpineshp <- readOGR(paste0(wd, "Spatial data/Processed/Templates and boundaries"), "Norway_alpine")

rastDir <- paste0(wd, "Spatial data/Processed/4_forMaxent_prediction/")
varnames <- c("Norway_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway", "Distance_to_Road_norway", "AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway", "State_commons_norway_binary",  "Protected_areas_norway_forbiological") # dput(list.files(rastDir, "*.asc$"))

envStack <- raster::stack(list.files(rastDir, "*.asc$", full.names=TRUE))
names(envStack) <- sapply(list.files(rastDir, "*.asc$"), function(x) strsplit(x, "\\.")[[1]][1])
envStack <- envStack[[varnames]] #remove unwanted variables
names(envStack)[[1]] <- "mask"

###PREDICT across all models 
#using a for loop because of memory allocation issues
for (y in c("Combined model", "North model", "South model"){
	for (x in seq_along(cultESlist)){
	
		currOutPath <- paste0(outDir, y, "/Base_run_", currDate, "/", as.character(cultESlist[x]))	
		#dir.create(currOutPath, recursive=TRUE)
		currMod <- readRDS(paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
		print(paste0("starting prediction ", y, cultESlist[x]))
		print(Sys.time())
		#make predictive maps
		currMap <- dismo::predict(currMod, envStack, ext=extent(alpineshp), args=c('outputformat=logistic', "doclamp=TRUE"), progress='text', filename=paste0(currOutPath,"/", cultESlist[x], "_prediction_frombasemodel.tif"), format="GTiff") 
		#currMap <- dismo::predict(currMod, envStack, ext=extent(alpineshp), args=c('outputformat=logistic', "doclamp=TRUE", "writeclampgrid=TRUE", "writemess=TRUE"), progress='text', filename=paste0(currOutPath,"/", cultESlist[x], "_basemodel.tif"), format="GTiff") 
		print(paste0("finished prediction ", y, cultESlist[x]))
		print(Sys.time())
		rm(currMap) #trying to clear memory issues
		
		}
}

###PREDICT across all models WITH BIAS GRID
#using a for loop because of memory allocation issues
for (y in c("Combined model", "North model", "South model"){
	for (x in seq_along(cultESlist)){
	
		currOutPath <- paste0(outDir, y, "/Base_run_bias_grid", currDate, "/", as.character(cultESlist[x]))	
		#dir.create(currOutPath, recursive=TRUE)
		currMod <- readRDS(paste0(currOutPath, "/Maxentmodel_rds_", as.character(cultESlist[x]), ".rds"))
		print(paste0("starting prediction ", y, cultESlist[x]))
		print(Sys.time())
		#make predictive maps
		currMap <- dismo::predict(currMod, envStack, ext=extent(alpineshp), args=c('outputformat=logistic', "doclamp=TRUE"), progress='text', filename=paste0(currOutPath,"/", cultESlist[x], "_prediction_frombasemodelbiasgrid.tif"), format="GTiff") 
		#currMap <- dismo::predict(currMod, envStack, ext=extent(alpineshp), args=c('outputformat=logistic', "doclamp=TRUE", "writeclampgrid=TRUE", "writemess=TRUE"), progress='text', filename=paste0(currOutPath,"/", cultESlist[x], "_basemodel.tif"), format="GTiff") 
		print(paste0("finished prediction ", y, cultESlist[x]))
		print(Sys.time())
		rm(currMap) #trying to clear memory issues
		
		}
}






######################
###COMPARE MODELS
######################
#test N model against S data
bg <- read.csv(paste0(outDir, "backgroundpoints.csv"))[,2:3]
for (x in 1:length(cultESlist)){
	currOutPath <- paste0(outDir, "/North model/Test NS data/", as.character(cultESlist[x]))
	dir.create(currOutPath, recursive=TRUE)
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
	dir.create(currOutPath, recursive=TRUE)
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




#fit using S, test using N and vice versa
#http://www.inside-r.org/packages/cran/dismo/docs/maxent
# #test how well the datasets predict each other using the predictions in the eval files	Warren's D similarity statistic		
# NScompare <- lapply(1:length(cultESlist), function(x) {
				# currEvalN <- readRDS(paste0(outDir, "North model/ENM eval/ENMevalofNmodel_", as.character(cultESlist[x]), ".rds"))
				# currEvalS <- readRDS(paste0(outDir, "South model/ENM eval/ENMevalofSmodel_", as.character(cultESlist[x]), ".rds"))	
				# bestN <- which.min(currEvalN@results$AICc) #pick model with lowest AICc
				# bestS <- which.min(currEvalS@results$AICc) #pick model with lowest AICc
				# NmodelRast <- currEvalN@predictions[[bestN]]
				# SmodelRast <- currEvalS@predictions[[bestS]]
				# return(nicheOverlap(NmodelRast, SmodelRast, mask=FALSE, stat='D'))
			# }
		
###########################	
#remove temp raster folder	
unlink(my_tmpdir, recursive = TRUE)
#########################
#END
#########################
