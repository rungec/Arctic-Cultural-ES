---
title: "Arctic Cultural ES Data processing notes"
author: "Claire Runge - NCEAS"
output: html_document
---

This is a summary of the data processing and analysis performed in relation to the Arctic CultES project led by Vera Hausner and undertaken by Claire Runge in 2016/2017.


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

###in ArcGIS from CORRINE2012 select CLC12_KODE <> 523 then dissolve
> Norway_nowater.shp

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

###CORRINE2012
*In r ArcticCulturalES_createInputs.r*
Calculated percentage of land cover in 1km (900m) circle around central cell

> Corrine2012_norway_xxx_1km.tif  
> broadleafforest = corine classes (311, 313)  
> coniferforest = corine classes (312)  
> heathshrub = corine classes (321:324)  
> sparselyvegetated = corine classes (331:335)  
> cropland = corine classes (211, 212, 213, 221, 222, 223, 231, 241:244)
> wetland = corine classes (411, 412, 422, 423)

###Industrial areas
*ArcGIS*  
Merged point files NVE/Vankraft_vannkraftverk (dams) & Vindkraft_Utbygd_Vindkraftverk (wind farms) and buffered by 5m to convert to polygon  
Merged line files Vankraft_vannvei (water pipelines) & Kraftnett_Kraftlinje (powerlines) and buffered by 5m to convert to polygon
Merged the point and line features from above with Corrine2012_minesdumpsconstruction.shp : Corrine 2012 where CLC12_KODE = c(131, 132, 133)  
> Industrial_disturbance_norway.shp  
Using Euclidean Distance in the Spatial Analyst toolbox, no maximum distance set, environments set to Norway_template_raster.tif  
*Note: this only transfers features that overlap the cell centre, so is incorrect. I didn't fix it in this layer, as we decided to exclude vannvei (water pipelines etc) as people don't seem so bothered by them i.e people will happily hike next to pipelines, and fish in dams*
> Distance_to_industrialdisturbance_norway.tif  
As above, without Vankraft_vannvei (water pipelines etc) 
> Industrial_disturbance_withoutvannvei_norway.shp
Dissolved
> Industrial_disturbance_withoutvannvei_norway_dissolve.shp

*in ArcGIS 30/12/16*
To create rasters first rasterized each layer individually, then combined.
Polyline to raster Kraftnett_Kraftlinje.shp, environments set to Norway_template_raster.tif
> Kraftnett_Kraftlinje.tif
Point to raster Kraftnett_Transformatorstasj.shp, environments set to Norway_template_raster.tif
>Kraftnett_Transformatorstasj.tif
Polygon to raster Corrine2012_minesdumpsconstruction.shp, environments set to Norway_template_raster.tif
> Corrine2012_minesdumpsconstr.tif
then...

*In r ArcticCulturalES_createInputs.r* 
Combined rasters using raster calculation
> Industrial_disturbance_withoutvannvei_norway.tif
> Percent_industrial_withoutvannvei_1km_norway.tif
> Percent_industrial_withoutvannvei_1km_norway_no_water.tif
Using ArcGIS, tabululate area of Percent_industrial_withoutvannvei_1km_norway_no_water.tif  
> Area_industry.csv
80% of land in Norway is further than 500m from industry.


###*27/07/16*
###Alpine climate regions within norway
*ArcGIS*  
Alpine= T_CLIM=='Z' in LANMAP2_LEV1.shp LANMAP-2 2007 Landscape map for Europe  
Clip to Norway_border_10kmbuffer.shp & dissolve  
> Norway_alpine.shp
Clipped to Norway_border.shp
> Norway_alpine_clip.shp
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
*06/08/16*  
Reclassified into five classes: 1= nature reserve; 2= national park; 3= iucn(3,4); 4= IUCN class 5,6); 0= not protected  
> Protected_areas_norway_forbiological.tif  
Reclassified into three classes: 1 = nature protection (iucn1-4); 2= managed landscapes (iucn 5,6); 0= not protected  
> Protected_areas_norway_forothervalues.tif  

###Ecologically important regions
NATURBASE/Naturtyper_flater.shp (use all categories; "VERDI" Localt viktig=locally important; Viktig=Important; Svaert viktig=very important)  
*ArcGIS*   
Added integer field VERDIint Localt viktig=10; Viktig=20; Svaert viktig=30  
*In r ArcticCulturalES_createInputs.r*  
as many of the polygons are very small, buffered by 49m. This should be enough so the polygon overlaps the raster cell centre if the polygon falls within the raster cell, but not so much that polygons that don't overlap a raster cell are assigned to that cell.
Converted to raster using Norway_template_raster as template. Value =10,20,30  
> Ecological_areas_norway.tif
*note: we decided not to include these in the models, and instead use them to compare whether and how well people see the same biological values as experts*

###State land
Statskog eiendom/Statskog eiendom 2014 (use only categories EKAT=2 or 3)  
*In r ArcticCulturalES_createInputs.r*  
Converted to raster using Norway_template_raster as template. Value =2,3  
> State_commons_norway.tif
07/08/16
Converted to raster using Norway_template_raster as template. private land or private commons = 0; State commons 2,3,4,6 = 1
> State_commons_norway_all.tif
> State_commons_norway_binary.tif

###State commons & protected areas combined  
00 private not protected; 01 government not protected; 10 nature protection, private; 11 nature protection government; 20 managed landsape private; 21 managed landscape government  
> Governance_plus_protectedareas_norway.tif  

###New Distance to Waterbodies
*in ArcGIS 06/08/16*
Combined Lakes & Rivers from N50 Data
lakes =Original/N50 Data/Lakes50/Innsjo_Innsjo.shp
rivers =Original/N50 Data/Rivers50/Elv_Elvenett.shp, objType=='ElvBekkMidtlinje'
> waterbodies.shp *note this dataset mistakenly only contains rivers - so don't use this layer, I didn't fix it as after discussion we decided to make another layer and only include major rivers & lakes >2ha, as below*
Using Euclidean Distance in the Spatial Analyst toolbox, no maximum distance set, environments set to Norway_template_raster.tif 
> Distance_to_waterbodies.tif
*Note: this only transfers features that overlap the cell centre, so is incorrect. I didn't fix it in this layer, as we decided to make a new layer including only large rivers & lakes*

*in ArcGIS 10/08/16*
Combined Lakes & Rivers from N50 Data
lakes =Original/N50 Data/Lakes50/Innsjo_Innsjo.shp select by attributes > 2ha
> Lakesbiggerthan2ha.shp
rivers =Original/N50 Data/Rivers50/Elv_Elvenett.shp, objType=='ElvBekkMidtlinje'
> MajorRiver.shp
Merged (with MajoRiver buffered to 5m)
> Waterbodies_majorriversandlakesbiggerthan2ha.shp
> Waterbodies_majorriversandlakesbiggerthan2ha_dissolve.shp

*in ArcGIS 30/12/16*
Converted to raster
Polyline to raster MajorRiver.shp
> MajorRiver.tif
Polygon to raster Lakesbiggerthan2ha.shp
> Lakesbiggerthan2ha.tif
then...

*In r ArcticCulturalES_createInputs.r* 
Combined rasters with Norway_template_raster.tif (raster calculation)
> Waterbodies_majorriversandlakesbiggerthan2ha.tif

#Percent water within 500m
> Percent_waterbodies_majorriversandlakesbiggerthan2ha_1km_norway.tif
#masked out water (sea, fjords)
> Percent_waterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water.tif
Raster with 1 if any water in radius, 0 if not
> AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway.tif
> AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water.tif
Using ArcGIS, tabululate area of AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water.tif  
> Area_ofland_within500mofwater.csv
43% of land in Norway is within 500m of water.

#Distance to road
*in ArcGIS*
The raster "N250 Data/Roads250/roaddist/distance_to_Road" had issues (1) values that are directly underneath roads were given value = 0, but so were areas that were water/sea. Also, areas that were further than 10km were given value of NA, meaning part of the study region was NA. For this reason, I remade this raster.

From shp of roads & ferry routes select by attributes
"OBJTYPE" 'Barmarkslgype' 'Sti' 'Traktorveg' 'VegSenterlinje' 'Bane' 'Barmarksllype' (roads, paths, tractor paths):
> input: N250 Data/Roads250/roaddist/no_n250_samferdsell", select by attributes
"OBJTYPE" 'Barmarkslgype' 'Sti' 'Traktorveg' 'VegSenterlinje' 'Bane' 'Barmarksllype' (roads, paths, tractor paths)
> output: Processed/shps/roads/no_n250_samferdsell_terrestrialonly.shp

Euclidean Distance, setting environments to Norway_template_raster resolution, extent, snap to. Maximum distance set to 40000 (40km) - values outside this are NA. Resolution = 100m
> Distance_to_Road_CR.tif
  


###Data that is already rasterized
*In r ArcticCulturalES_createInputs.r*  
resampled to Norway_template_raster. Rasters were already at 100m resolution and correct projection, just needed cropping and resampling to overlap. 
extend raster to match Norway_template_raster.tif

> distance_to_Town_norway.tif = "SSB/Tettsted2015/tet_dist2015"  
> distance_to_River_norway.tif = "N50 Data/Rivers50/river_dist"  
> *distance_to_Coast_norway.tif = "N250 Data/Coastline250/coastdistance" see edits below* 
> distance_to_Road_norway.tif = "N250 Data/Roads250/roaddist"  
> distance_to_House_norway.tif = "N250 Data/Buildings250/house_dist"  *not used - has issues*

###All data
applied mask based on Norway_border_10kmbuffer.shp, areas outside shp = NA  
> .tif rasters in folder "masked_to_norway"

applied mask based on Norway_nowater.shp, areas outside shp = NA
> .tif rasters in folder "masked_no_water"

output as asciis for Maxent  
> .asc rasters in folder "forMaxent_asc"

clipped extent to norway_alpine_clip.shp  
> .asc rasters in folder "forMaxent_predictions_alpine"

###Reprocessing of Distance_to_Coast_norway raster
This raster has NA values across much of the study region. This raster was subsetquently remade from N250 Data\Coastline250\coastlinel2.shp 
*ArcGIS*
Using Euclidean Distance in the Spatial Analyst toolbox, no maximum distance set, environments set to Norway_template_raster.tif  
*In r ArcticCulturalES_createInputs.r*  
Applied mask of Norway_alpine, saved as .asc  

> Distance_to_Coast_norway
*Note this was not included in the final modelling as its interpretation is based on supposition - it is an indirect metric of the things people percieve, like whether landscape is coastal or mountainous*  

###Bias grids for maxent  
Created bias grids for maxent using script *ArcticCulturalES_roadaccessmodel.r*
Modelled the frequency of PPGIS data by proximity to road using a non-linear least-squares model, then predicted that to the distance_to_Road_norway.tif raster. 

value is predicted frequency of ppgis data  (0-3000)
> BiasGrid_distancetoroad.tif  
> BiasGrid_distancetoroad_no_water.tif  #masked out sea & fjords
> BiasGrid_distancetoroad_no_water.asc  #used in Maxent
value is normalised predicted frequency between 0 & 1: 1=more sampled
> BiasGrid_distancetoroad_normalised.tif  
> BiasGrid_distancetoroad_normalised_no_water.tif  #masked out sea & fjords


###Background points
for NS combined models: 10000 points drawn from whole of study region as defined by Norway_alpine.shp
> backgroundpoints_wholeregion.csv
for N only models: 10000 points drawn from region defined by North_municipalities_alpine.shp
> backgroundpoints_north.csv
for S only models: 10000 points drawn from region defined by South_municipalities.shp
> backgroundpoints_south.csv

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

###Maxent model testing runs
####Response curves and jacknife
First run. Automatic regularization, background from whole study region.
After discussion, have come up with some different variables that are hav stronger justification, and are more readily explained by realistic hypotheses about people's behaviour and preferences. Primarily: Corine will be converted from a categorical layer to a series of layers describing the percentage of landcover type in a given radius around each cell. Will also include a layer describing the closeness of industrial development, and simplify the categorization of the protected area and state commons layers.  


####Response curves and jacknife2
Ran north & south models for all features. *regularization automatic*, background from north (municipalities unclipped) and south (municipality) study region boundaries respectively. 
Noted that the north municipalites extend into the sea, so there are a lot of points very far from rivers in the sea. 
Noticed that the distance to waterbodies only includes rivers; that the municipal boundaries mask was applied incorrectly in a couple of the runs; that I forgot to include cropland. 
Otherwise, models are looking good and make sense. Possibly slight overfitting.
Percentage of x in landscape rasters in *3km square* around central cell.
State commons - two categories (private, government)
Protected areas - 3 to 5 categories, more categories for bioligical and undisturbnature values; (nature protection, managed landscape, no protection).  

####Response curves and jacknife3  
Ran north & south models for biological and undisturbnature, and a few of the other values. *regularization (beta_hinge) =2*, background from north (municipalities clipped to alpine) and south (municipality) study region boundaries respectively.  
Distance to industrial disturbance is highly correlated with Distance to waterbodies (pearson correlation 0.8) and the biological models show strong relationship to industrial (high probability close to industrial) and no relationship with waterbodies. For this reason, this run was cancelled, and distance to industrial disturbance was replaced with percent industrial in 1km.  
#####variables used: 
Percentage of x in landscape rasters in *1km circle* around central cell.
#####For other values:  
State commons & protected areas combined into a single variable - Governance_plus_protectedareas_norway (00 private not protected; 01 government not protected; 10 nature protection, private; 11 nature protection government; 20 managed landsape private; 21 managed landscape government)  
#####For biological & undisturbnature:  
State commons - two categories (private, government)  
Protected areas - 5 categories, (nature protection, national park, iucn cat 3-4, managed landscape, no protection).  

c("North_municipalities_alpine","South_municipalities", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_heathshrub_1km",  
"Corrine2012_norway_sparselyvegetated_1km", "Corrine2012_norway_wetland_1km", "Corrine2012_norway_cropland_1km", "Distance_to_Coast_norway2", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbodies_norway2", "Distance_to_industrialdisturbance_norway2","State_commons_norway_binary",  "Protected_areas_norway_forbiological")  


####Response curves and jacknife4
Ran north & south models for biological and undisturbnature. regularization (beta_hinge) automatic (=1). Models show some overfitting.

variables used: c("North_municipalities_alpine", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_heathshrub_1km",  
"Corrine2012_norway_sparselyvegetated_1km", "Corrine2012_norway_wetland_1km", "Corrine2012_norway_cropland_1km", "Distance_to_Coast_norway2", "Distance_to_Road_norway", "Distance_to_Town2_norway", "Distance_to_waterbodies_norway2", "Percent_industrial_1km","State_commons_norway_binary",  "Protected_areas_norway_forbiological")

### 23/12/16
###Correlation of environmental variables
Checked correlation of rasters that had been masked to land (areas that were not classified as water in Corrine)
Saved as .rds which can be loaded to R using readRDS  

> folder: Maxent runs\Correlation of variables\
> CorrelationPlotofEnvironmentalVariables_nowater.png
> CorrelationofEnvironmentalVariables_nowater.rds

### 07/01/17
###Correlation of environmental variables
Checked correlation of new waterbodies & industrial rasters
 Most variables show pearson correlation of <0.4. Exceptions are Distance to road and distance to house which show pearson correlation of 0.73; Distance to Town2 shows correlation of 0.54 with distance to house and 0.47 with distance to road;
 Percent industrial 1km shows correlation of 0.46 with Corrine_2012 cropland; conifer forest shows negative correlation of -0.43 with sparse vegetated.  
Saved as .rds which can be loaded to R using readRDS  
> folder: Maxent runs\Correlation of variables\
> CorrelationPlotofEnvironmentalVariables_nowater2.png
> CorrelationofEnvironmentalVariables_nowater2.rds


###16/01/17
Test of regularisation parameters (=1, lower predictive power (more overfitted) but better fit to training data (higher AUC))
Ran N, S, Combined models with and without bias grids  
regularization (beta_hinge=1)  
> Response curves and jacknife 20170116: 10-fold cross validation, with jacknife  
> Base models 20170116:  25% of data set aside for testing  
> variables used: c("Norway_alpine", "North_municipalities_alpine", "South_municipalities, AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway","Distance_to_Road_CR_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway", "State_commons_norway_binary",  "Protected_areas_norway_forbiological")  

###17/01/17
Test of regularisation parameters (=1.5, higher predictive power (less overfitted) but lower fit to training data (lower AUC))
Ran N, S, Combined models with and without bias grids  
regularization (beta_hinge=1.5)  
> Response curves and jacknife 20170117: 10-fold cross validation, with jacknife  
> Base models 20170117:  25% of data set aside for testing  
> variables used: c("Norway_alpine", "North_municipalities_alpine", "South_municipalities, AnyWaterbodies_majorriversandlakesbiggerthan2ha_1km_norway_no_water", "Corrine2012_norway_broadleafforest_1km", "Corrine2012_norway_coniferforest_1km", "Corrine2012_norway_cropland_1km", "Corrine2012_norway_heathshrub_1km", "Corrine2012_norway_sparselyvegetated_1km", "Distance_to_Town2_norway","Distance_to_Road_CR_no_water", "Percent_industrial_1km_norway_nowater", "Governance_plus_protectedareas_norway", "State_commons_norway_binary",  "Protected_areas_norway_forbiological")  

###21/01/17
Decision process around predictions:
Options are to do clamping and whether MESS/clamp grids are needed.
The models were created in a given study region (north, south, combined north-south) with background drawn from the same study region. We would expect that there are areas outside these regions with qualities of the landscape that are not represented within these study regions (ie there might be places further from roads, or with more forest) that were not captured in the background sampling of the original (base or base + bias grid) model.
The question is, when we spatially predict these models, do we want to limit the 'habitat suitability' value to for these places to match that at the maximum (or minimum) limits for that particular environemntal variable (do clamping), or let the model prediction extrapolate from the modelled response (don't do clamping)? 
Answer: for the initial purposes, I will ignore clamping. (1) I suspect that the environemnt is well represented within the landscape, at least for the qualities that matter to these models (the main variables). This could be confirmed by exploring the ranges of the environmental variables in each of the study region background samples. (make swd of background, summary of each column), though of course there will be some combinations of environmental variables that will not be present. (2) most of these models are hinge, and thus the response curve is self-clamped to represent the average prediction outside response curve limits - look at the response curves and you will see they are flat horizontal lines outside curve limits. 

At some stage, it will be useful to generate a MESS grid for the landscape, to see where values exist that were not present during the model parameterization, and to generate clampgrids to determine a measure of uncertainty in the predictions.



***
ArcGIS version 10.3.1
R version 3.2.3
***
