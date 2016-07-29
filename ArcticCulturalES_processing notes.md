---
title: "Arctic Cultural ES Data processing notes"
author: "Claire Runge - NCEAS"
output: html_document
---

This is a summary of the data processing and analysis performed in relation to the Arctic CultES project led by Vera Hausner and undertaken by Claire Runge in 2016.


***
##Preparation of PPGIS data
*ArcGIS*  
originally this data is in WGS84 (Google Earth)  
loaded to arcgis, converted to   
> PPGIS_Markers_south.shp & PPGIS_Markers_north.shp

projected to ETRS_1989_UTM_Zone_33N saved as   
> PPGIS_Markers_south_UTM33N.shp & PPGIS_Markers_north_UTM33N.shp

clip to alpine climate zone  
> PPGIS_Markers_south_UTM33N_alpine.shp & PPGIS_Markers_north_UTM33N_alpine.shp

this excludes 8726/9149= 4.6% and 9203/10222=10.0% of the data

***
##Preparation of Spatial data
###templates & boundaries
*ArcGIS*  
merged Basemap500::n500_nationalborderseal & n500_nationalborderlandl then feature to polygon  
> Norway_border.shp

buffered to 10km  
> Norway_border_10kmbuffer.shp

###Create template raster
*ArcGIS*  
Create Constant Raster with  
constant value: 0  
output data type: INTEGER  
output cell size: D:\Arctic Cultural ES\Spatial data\CORINE2006\Norway.tif  
output extent: same as Norway_border_10kmbuffer.shp  
environments::snap to CORINE2006\Norway.tif  
environments::mask Norway_border_10kmbuffer.shp  
>Norway_template_raster.tif

###Clip CORRINE2006
*ArcGIS*  
Spatial Analyst::Extract by Mask  
inputraster: D:\Arctic Cultural ES\Spatial data\CORINE2006\Norway.tif  
mask feature: Norway_template_raster.tif  
*In r ArcticCulturalES_createInputs.r* 
extend raster to match Norway_template_raster.tif   
> Corrine2006_norway.tif

###set sea as NA
*ArcGIS*  
reclassify Corrine2006_norway.tif value=44 -> NoData  
*In r ArcticCulturalES_createInputs.r* 
extend raster to match Norway_template_raster.tif  
> Corrine2006_norway_noSea.tif

###*27/07/16*
###Alpine climate regions within norway
*ArcGIS*  
Alpine= T_CLIM=='Z' in LANMAP2_LEV1.shp LANMAP-2 2007 Landscape map for Europe  
Clip to Norway_border_10kmbuffer.shp & dissolve  
> Norway_alpine.shp
*In r ArcticCulturalES_createInputs.r*
Made a mask to limit Maxent background to alpine climate regions (0 if alpine, else NA)  
> Norway_alpine.tif & Norway_alpine.asc

note that this shp has been oversimplified and many coastal islands are missing

###Protected areas 
Norway protected areas\Protected area Norway Svalbard.shp (use all IUCN categories)  
*ArcGIS*   
Added integer field "IUCNCat" (1-6=IUCN Ia to IUCN VI, 7=no IUCN category)  
*In r ArcticCulturalES_createInputs.r*  
Converted to raster using Norway_template_raster as template. Value =1-7  
> Protected_areas_norway.tif

###Ecologically important regions
NATURBASE/Naturtyper_flater.shp (use all categories; "VERDI" Localt viktig=locally important; Viktig=Important; Svaert viktig=very important)  
*ArcGIS*   
Added integer field VERDIint Localt viktig=10; Viktig=20; Svaert viktig=30  
*In r ArcticCulturalES_createInputs.r*  
as many of the polygons are very small, buffered by 49m. This should be enough so the polygon overlaps the raster cell centre if the polygon falls within the raster cell, but not so much that polygons that don't overlap a raster cell are assigned to that cell.
Converted to raster using Norway_template_raster as template. Value =10,20,30  
> Ecological_areas_norway.tif

###State land
Statskog eiendom/Statskog eiendom 2014 (use only categories EKAT=2 or 3)  
*In r ArcticCulturalES_createInputs.r*  
Converted to raster using Norway_template_raster as template. Value =2,3  
> State_commons_norway.tif

###Data that is already rasterized
*In r ArcticCulturalES_createInputs.r*  
resampled to Norway_template_raster. Rasters were already at 100m resolution and correct projection, just needed cropping and resampling to overlap. 
extend raster to match Norway_template_raster.tif

> distance_to_Town_norway.tif = "SSB/Tettsted2015/tet_dist2015"  
> distance_to_River_norway.tif = "N50 Data/Rivers50/river_dist"  
> distance_to_Coast_norway.tif = "N250 Data/Coastline250/coastdistance"  
> distance_to_Road_norway.tif = "N250 Data/Roads250/roaddist"  
> distance_to_House_norway.tif = "N250 Data/Buildings250/house_dist"  

###All data
applied mask based on Norway_border_10kmbuffer.shp, areas outside shp = NA  
> .tif rasters in folder "masked"

output as asciis for Maxent  
> .asc rasters in folder "forMaxent"

###
***
###*28/07/16*
##Maxent processing
*in R using ArcticCulturalES_maxent_models.r*

###Correlation of environmental variables
Checked correlation
Most variables show low correlation, with the exception of distance to road and distance to house which show pearson correlation coefficient of 0.76. Given road is a variable previously identified as important, we retain it in the models and remove distance to house to allow the response curves to be more easily interpreted. The next highest correlation was the corrine (land cover) dataset which showed correlation of -0.54 with distance to coast, -0.56 with distance to town.  
Saved as .rds which can be loaded to R using readRDS  

> folder: Maxent runs\Correlation of variables\
> CorrelationPlotofEnvironmentalVariables.png
> CorrelationofEnvironmentalVariables.rds

***
ArcGIS version 10.3.1
R version 3.2.3
***
