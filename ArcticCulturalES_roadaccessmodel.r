#This script models the relationship between ppgis points and distance from roads
#and creates a raster of 'bias' across the landscape

library(raster)
library(rgdal)
library(ggplot2)
library(nlme)

##########################################
#Setup

options(stringsAsFactors=FALSE) # turn off automatic factor coersion
# options(scipen=9999)            # turn off plotting axis lables in scientific notation


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
rastDir <- paste0(wd, "Spatial data/Processed/")
rastTemplate = raster(paste0(rastDir, "Templates and boundaries/Norway_template_raster.tif"))

############################
#PPGIS data
southPPGIS <- readOGR(paste0(wd, "PPGIS data/Processed/shps"), "PPGIS_Markers_south_UTM33N_municipal") 
northPPGIS <- readOGR(paste0(wd, "PPGIS data/Processed/shps"), "PPGIS_Markers_north_UTM33N_alpine") 
#distance to road rasterOptions
disttoroadRast <- raster(paste0(rastDir, "2_masked_to_norway/Distance_to_Road_norway.tif"))

#############################
#extract distance to road for each ppgis points
disttoroadS <- extract(disttoroadRast, southPPGIS, df=TRUE, method='simple')
disttoroadN <- extract(disttoroadRast, northPPGIS, df=TRUE, method='simple')

allDat <- data.frame(disttoroad=append(disttoroadS[,2], disttoroadN[,2]), region=rep(c("South", "North"), times =c(nrow(disttoroadS), nrow(disttoroadN))))

ggplot(allDat, aes(disttoroad, fill=region)) +#or try group
 	geom_histogram() +
	theme_light(17) + #get rid of grey bkg and gridlines
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	labs(x="Distance to road", y="Number of records")+
	scale_x_continuous(expand = c(0.0, 0.0)) + scale_y_continuous(expand = c(0.01, 0.001))+ #set x and y limits
	theme(axis.title.x = element_text(vjust=-0.6),axis.title.y = element_text(vjust=1))+	#move xylabels away from graph
	theme(legend.position="bottom")#gets rid of legend

outPath <- paste0(wd, "/figures/Histogram_DistancetoRoad.png")
	ggsave(filename=outPath)

###############################
#Model relationship between number of records & distance to road
counttable <- table(allDat$disttoroad)
countdf <- data.frame(disttoroad=as.numeric(dimnames(counttable)[[1]]), frequ=as.vector(counttable))

with(countdf, plot(disttoroad, frequ))
with(countdf, plot(log(frequ)~log(disttoroad)))

mod1 <- with(countdf, lm(log(frequ)~-log(disttoroad))) #straight line, residuals are not normal
mod2 <- nls(frequ~exp(a+b*disttoroad), data=countdf, start=list(a=0, b=0)) #exponential line

with(countdf, plot(disttoroad, frequ))
lines(countdf$disttoroad, predict(mod2, list(disttoroad = countdf$disttoroad)), col="red") #mod2 is ok
#lines(countdf$disttoroad, exp(predict(mod1, list(disttoroad = countdf$disttoroad))), col="blue")#mod1 is a poor fit 

#Try a subset of the data
countdfsub <- countdf[countdf$disttoroad<2000,]
summary(countdfsub)
mod3 <- nls(frequ~exp(a+b*disttoroad), data=countdfsub, start=list(a=0, b=0))
summary(mod3)
plot(mod3)

with(countdfsub, plot(disttoroad, frequ))
lines(countdfsub$disttoroad, predict(mod3, list(disttoroad = countdfsub$disttoroad)), col="red")
#lines(countdfsub$disttoroad, predict(mod2, list(disttoroad = countdfsub$disttoroad)), col="blue")

summary(mod2)
summary(mod3)#models are essentially the same

#Check that north and south data have same relationship to distance to road
counttable2 <- as.dataframe(table(allDat$disttoroad, allDat$region))
t.test(counttable2[counttable2$Var2=="South", "Freq"],counttable2[counttable2$Var2=="North", "Freq"])


###############################
#Predict model aross raster of distance to road
#need to normalise values between 0 and 1

#biasrast <- predict(disttoroadRast, mod2, fun=predict.nls, newdata=list(disttoroad=values(disttoroadRast)), filename=outPath, format="GTiff")#couldnt get this to work so did it the long way
newpreds <- getValues(disttoroadRast)
maxpreds <- predict(mod2, list(disttoroad=c(0)))
RastPredictions <- predict(mod2, list(disttoroad = newpreds))
RastPredictionsNorm <- round(predict(mod2, list(disttoroad = newpreds))/maxpreds, 2) #normalise #if outside max distance seen in data, then 0, if close to road then 1

#set the values of these rasters to model predictions
outPath <- paste0(rastDir, "Bias_grids/BiasGrid_distancetoroad.tif")
biasRast <- disttoroadRast 
biasRast <- setValues(biasRast, RastPredictions)
writeRaster(biasRast, filename=outPath, format="GTiff",overwrite=TRUE)

outPath <- paste0(rastDir, "Bias_grids/BiasGrid_distancetoroad_normalised.tif") 
biasRastNorm <- disttoroadRast
biasRastNorm <- setValues(biasRastNorm, RastPredictionsNorm)
writeRaster(biasRastNorm, filename=outPath, format="GTiff", overwrite=TRUE)

###############################
#Plot the models and the bias raster

pdf(paste0(wd, "/Bias_grids/BiasGridModel_DistancetoRoad.pdf"))
par(mfrow=c(2,2), mar = c(4,4,1,1) + 0.1, mgp=c(2,1,0))

#plot panel 1 raw data
with(countdf, plot(disttoroad, frequ,
xlab="Distance to road (metres)", ylab="Frequency of PPGIS data", fg="grey70"))

#plot panel 2 log data
with(countdf, plot(log(frequ)~log(disttoroad),
xlab="Log distance to road (metres)", ylab="Log frequency of PPGIS data", fg="grey70"))

#plot panel 3 zoomed data plus model fit
with(countdfsub, plot(disttoroad, frequ,
xlab="Distance to road (metres)", ylab="Frequency of PPGIS data", fg="grey70"))
lines(countdfsub$disttoroad, predict(mod3, list(disttoroad = countdfsub$disttoroad)), col="red")
text(x=1500, y=2500, labels=expression(y == e^{8.1-0.006*x}))

#plot panel 4 normalised raster
plot(biasRastNorm)
title("Bias grid for Maxent")

dev.off()


#run maxent models with default beta value or test aic across different beta values
#run maxent models with and without bias grid
#write up model selection - why use maxent as opposed to random forest? (answer - because random forest needs absenses, which we don't have)
#look at where people value biodiversity
#compare with expert predictions
#methods paper
#N model & predict across whole region
#S model & predict across whole region
#combined all data & predict across whole region
#for all values
#then calculate niche overlap between the different models