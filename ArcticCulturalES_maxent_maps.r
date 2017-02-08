#Makes plots of the spatial predictions from ArcticCulturalES_maxent_models.r
#Not optimised for speed!

library(raster)
library(rasterVis)
library(rgdal)

wd <- "C:/Claire/Arctic_Cultural_ES/"
#wd <- "/home/runge/Data/Arctic_Cultural_ES/"
setwd(wd)

modelDir <- paste0(wd, "Maxent runs/")
plotDir <- paste0(wd, "figures/")

currDate <- paste(strsplit(as.character(Sys.Date()), "-")[[1]], collapse="") #formatting the date for folder naming

########################
## Setting system preferences for better work with raster
# Define your temp folder
# my_tmpdir=paste0(wd, "tmpRaster")

# # Create it (handles the case where the folder already exists)
# dir.create(my_tmpdir, showWarnings=F)

# # Set the raster option to this folder
# rasterOptions(tmpdir= my_tmpdir)

# ## Maximum Memory
# # Often server has more RAM than desktop. Therefore the mamximum amount of memory that can be allocated 
# # to process a raster can be increased. Aurora has a lot of RAM (384GB), so we recommend to use the following 
# # settings:
# rasterOptions(maxmemory =1e+09)
# rasterOptions(chunksize=1e+08)

######################
#a function to capitalise first letter of a word
capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


######################

cultESlist <- c("biological", "cabin", "cleanwater", "cultureident", "gathering", "hunt_fish", "income", "pasture", "recreation","scenic", "social", "specialplace", "spiritual", "therapuetic", "undisturbnature")


folderlist <- paste0(rep(c("Combined model", "North model", "South model"), each=2), paste0(c("/Base run ", "/Base run bias grid "), c("20170116/")))
plotnames <- c("Norway", "", "North", "", "South")
plotfilename <- c("_LatentvsRealised_Norway_", "", "_LatentvsRealised_North_", "", "_LatentvsRealised_South_")

#Set up masks
mtemprast <- raster(paste0(modelDir, folderlist[1], cultESlist[1], "/", cultESlist[1], "_prediction_alpine.tif"))
northshp <- readOGR("Spatial data/Processed/Templates and boundaries", "North_municipalities_alpine")
northmask <- raster("Spatial data/Processed/Templates and boundaries/North_municipalities_alpine.tif")
northmask <- crop(northmask, mtemprast)

southshp <- readOGR("Spatial data/Processed/Templates and boundaries", "South_municipalities")
southmask <- raster("Spatial data/Processed/Templates and boundaries/South_municipalities.tif")
southmask <- crop(southmask, mtemprast)			


for (x in cultESlist){
	for (i in c(1,3,5)){
			rreal <- raster(paste0(modelDir, folderlist[i], x, "/", x, "_prediction_alpine.tif"))
			rlat <- raster(paste0(modelDir, folderlist[i+1], x, "/", x, "_prediction_alpine.tif"))
			plotStack <- stack(rreal, rlat)
			names(plotStack) <- c("Realised", "Latent")
		if (i==1) {
			plotStack <- plotStack
		} else if (i==3){
			plotStack <- mask(plotStack, northmask)
			plotStack <- crop(plotStack, extent(northshp), snap='out')
		} else if (i==5) {
			plotStack <- mask(plotStack, southmask)
			plotStack <- crop(plotStack, extent(southshp), snap='out')
		}
		#make plots
		print(paste0("Starting plot ", x))
		png(filename=paste0(plotDir, capwords(x), plotfilename[i], currDate, ".png"),width=1024*0.707*2, height=1024, res=300)
		print(levelplot(plotStack, main=paste0(capwords(x), " - ", plotnames[i]), scales=list(draw=FALSE)) )#remove latlon
		dev.off()

	}		
			
}


#latent vs realised for combined model
#latent vs realised for N region
#latent vs realilsed for S region
