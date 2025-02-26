
####
# This is a code prepared for the POWER ANALYSIS FOR OPTIMAL DESIGN OF A 
# PASSIVE ACOUSTIC MONITORING NETWORK FOR US EAST COAST OFFSHORE WIND 
# project

# This code can be used to calculate power in regional analysis (several windfarms at the same time)
# Prepared by Magda Chudzinska (SMRUC, St Andrews) January 2025

# loading packages
rm(list=ls())
#gc()
#.libPaths()
require(raster)
library(ggplot2)
library(viridis)
# url <- "https://download.r-forge.r-project.org/bin/windows/contrib/4.4/rgdal_1.6-7.zip"
# install.packages(url, type="source", repos=NULL)
# require(rgdal)
# library(rgeos)
require(sf)
require(sp)
library(tibble)
library(patchwork)
library(this.path)
require(mgcv)
require(exactextractr)
library(parallel)
library(doParallel)
library(foreach)
library(plot.matrix)
require(rlist)
require(dplyr)
require(mgcViz)
require(geosphere)
require(terra)
require(ggplot2)
require(RColorBrewer)

BigStart <- Sys.time()
#set.seed(2179155) # for reproducibility

# loading functions
source("PAM_Power_analysis_functions.R")

##### Step 0 - all the parameters #######
#########################################

# how many years of post-construction monitoring/operation
# note that >5 years will require generating more realization maps 
# delivering code to do that is not part of thos project
NmonitoringYears <- 5

#how many realization maps we use

# utm projection - final projection I use for all layers
# Roberts et al models have all covariates in lat and long so there will be a lot reshuffling between projections
newproj <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"
newprojLatLong <-  "+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84
+units=m +no_defs" # taken from Jason's maps

# which species (FW - fin whale, MW - minke whale, SW - sei whale, NA - north Altantic)

species <- "NA"

# which windfarm, for regional analysis it is only one
# so dont change it

windfarm <- "all"

# which PAM grid design: 

PAM <- "Original" #as given by Julia
PAM <- "Original_Large"
PAM <- "Original_Large_POWmoved"
# parameters related to making N realisation maps of the baseline abundance from Jason's model

if (NmonitoringYears == 5) N <- 250 # this is final number of realization maps

# parameters related to cue rate and detection range
# species - specific 

# cue production (in cues/hour/individual)
if (species == "MW") cues <- 6.04 # minke (Martin et al. (2013))
if (species == "FW") cues <- 45.08 # fin whales (Stimpert et al. (2015))
if (species == "SW") cues <- 10 
if (species == "NA") cues <- 6.4

# detection range (in km) (new based on EDF)
if (species == "MW") dr <- 8.6 
if (species == "SW") dr <- 21.1
if (species == "FW") dr <- 99.9 
if (species == "NA") dr <- 9.2


species_cue_rate <- cues #  HAS TO BE IN CUE RATES/HOUR
det_range <- dr #in km
cue_rate <-  species_cue_rate * 24 # cue rates/hour/individual, to get by month, I multiply it by 24h; I assume monthly densities are mean per day as they come from aerial surveys done on daily basis

# parameters related to operation
# we assume operation takes place 12 months/year

if (windfarm=="all"){

  monthsOperation <- c(1:12) # we assume all year
  monthsBaselineOp <- c(1:12) # we assume all year
}


nmonthsOperation<- length(monthsOperation) 
nmonthsBaselineOp <- length(monthsBaselineOp)

# parameters related to whale movement and re-distribution procedure
maxDist <- 100 * 1000 # maximum displacement distance for whales over 1 month (has to be in meters)

# parameters related to H8: global decline
# set what is population decline [ per year]
pop_decline <- 0.05

# parameters related to parallel computing
#parallel::stopCluster(cl = mycluster)
# make sure you dont use all your cores
parallel::detectCores()
numCores <- parallel::detectCores() - 6

mycluster <- parallel::makeCluster(
  numCores, 
  type = "PSOCK",
  outfile=""
)

# parameters related to power analysis, what is your significance threshold
alpha <- 0.05


##################################################################################################################
### Step 1: Preparing GIS layers 
##################################################################################################################


# first loading a raster - only for getting extend of study site
# so it does not matter for which species and month it is

fw4 <- raster("GIS/Density_04.img")

# windfarm footprint and buffer

# loading windfarm outlines, different for each year of monitoring operation
# but we agreed to run one simulation for 2029
# I leave the other one in the code

# wind_foot24 <- read_sf("GIS/WELAs_24_projected.shp")
# wind_foot25 <- read_sf("GIS/WELAs_24-25_projected.shp")
# wind_foot26 <- read_sf("GIS/WELAs_24-26_projected.shp")
# wind_foot27 <- read_sf("GIS/WELAs_24-27_projected.shp")
# wind_foot28 <- read_sf("GIS/WELAs_24-28_projected.shp")
wind_foot29 <- read_sf("GIS/WELAs_24-29_projected.shp")
# wind_foot30 <- read_sf("GIS/WELAs_24-30_projected.shp")
# wind_foot31 <- read_sf("GIS/WELAs_24-31_projected.shp")


if (PAM == "Original") {
  pamAllLL <-  read_sf("GIS/Devices_for_power_analysis_MainPower.shp")
  pamAllLL$Type <- "Original"
}

if (PAM == "Original_Large") {
  pamAllLL <-  read_sf("GIS/Devices_for_power_analysis_MainPower.shp")
  pamAllLL$Type <- "Original"
  
  pamAllLL <- pamAllLL %>%
    select(c("LATITUDE","LONGITUDE","geometry","Type"))
  
  pamAllLLLarge <-  read_sf("GIS/Extra_PAM_largeScale.shp")
  pamAllLLLarge$Type <- "Large"
  pamAllLLLarge <- pamAllLLLarge %>%
    select(c("Latitude","Longitude","geometry","Type"))
  
  colnames(pamAllLLLarge) <- colnames(pamAllLL)
  pamAllLLLarge <- st_transform(pamAllLLLarge, crs(pamAllLL))
  pamAllLL <- rbind.data.frame(pamAllLL,pamAllLLLarge)
}

if (PAM == "Original_Large_POWmoved") {
  pamAllLL <-  read_sf("GIS/Extra_PAM_largeScale_minusPOWERON.shp")
  pamAllLL$Type <- "PowMoved"
  
  pamAllLL <- pamAllLL %>%
    select(c("LATITUDE","LONGITUDE","geometry","Type"))
  
  pamAllLLLarge <-  read_sf("GIS/Extra_PAM_largeScale.shp")
  pamAllLLLarge$Type <- "Large"
  pamAllLLLarge <- pamAllLLLarge %>%
    select(c("Latitude","Longitude","geometry","Type"))
  
  pamAllLLLarge <- st_transform(pamAllLLLarge, crs(pamAllLL))
  colnames(pamAllLLLarge) <- colnames(pamAllLL)
  pamAllLL <- rbind.data.frame(pamAllLL,pamAllLLLarge)
}


# it is better if all data are in utm. Polygon already is so I will only convert the raster
# but because Jason Roberts's covariates are in lat long and I make map realisations in lat long, I have to stick to lat long for now and only later convert into utm
# windfarms are already in lat long but I transform it anyway

# wind_foot24LL <- st_transform(wind_foot24, newprojLatLong)
# wind_foot25LL <- st_transform(wind_foot25, newprojLatLong)
# wind_foot26LL <- st_transform(wind_foot26, newprojLatLong)
# wind_foot27LL <- st_transform(wind_foot27, newprojLatLong)
# wind_foot28LL <- st_transform(wind_foot28, newprojLatLong)
wind_foot29LL <- st_transform(wind_foot29, newprojLatLong)
# wind_foot30LL <- st_transform(wind_foot30, newprojLatLong)
# wind_foot31LL <- st_transform(wind_foot31, newprojLatLong)

# changing to sp objects as many functions do not work for sf
# wind_foot24LLsp <- sf:::as_Spatial(wind_foot24LL)
# wind_foot25LLsp <- sf:::as_Spatial(wind_foot25LL)
# wind_foot26LLsp <- sf:::as_Spatial(wind_foot26LL)
# wind_foot27LLsp <- sf:::as_Spatial(wind_foot27LL)
# wind_foot28LLsp <- sf:::as_Spatial(wind_foot28LL)
wind_foot29LLsp <- sf:::as_Spatial(wind_foot29LL)
# wind_foot30LLsp <- sf:::as_Spatial(wind_foot30LL)
# wind_foot31LLsp <- sf:::as_Spatial(wind_foot31LL)

# cropping big Jason's raster map to the study site
ee <- extent(wind_foot29LL)
ee[3] <- ee[3]-60000
ee[4] <- ee[4]+60000
ee[2] <- ee[2]+100000
studySite <- crop(fw4,ee)

# cropping and transforming PAM to the study site
pamAllLLSS <- st_transform(pamAllLL, newprojLatLong)
pamAllLLSS <- st_crop(pamAllLLSS,ee)

# quick check in lat long -  ALL LOOK FINE

plot(studySite)
plot(wind_foot29LL,add=T)
plot(pamAllLLSS, add=T, pch=16, col="red")


###############################################################################################################
### Step 2 Generating n realization of density map
###############################################################################################################

# Generation of 1250 (NmonitoringYears * N) realization maps from Jason Roberts's models is all done in a separate, 
# species - specific code. Generating more realizations would require access to all model covariates
# and model specifications and both require additional permission from Jason Roberts.
# Delivering these input data is not part of the project and here we provide 250 realizations maps per year for each species 
# as an example. 
# these are realizations for each month of the year cropped to the study site and turned in data frames
# these maps are already in # whales/25 km2

if (species == "FW") {load("BaselineRealisationMaps/FinWhale/Dens_SS_operationDF1250.RData")}
if (species == "NA") {load("BaselineRealisationMaps/NARW/Dens_SS_operationDF1250.RData")}
if (species == "SW") {load("BaselineRealisationMaps/SeiWhale/Dens_SS_operationDF1250.RData")}
if (species == "MW") {load("BaselineRealisationMaps/MinkeWhale/Dens_SS_operationDF1250.RData")}


#################################################################################################################################
## Step 5: Generating change in distribution - H7 - all are gone from footprint
#          Generating change in distribution - H6 - attraction to windfarm
#          Generating change in distribution - H8 - global decline
#################################################################################################################################

# because the method from Saana not only takes ages but also generates a massive file, I am going to do it differently
# as we are analyzing a large area, we assume that for H7, all whales which are at the footprints of the lease sites, are just gone
# if I try to redistribute them, they would all go very far from the sites as all sites are as one polygon

# testing some simpler option
# first finding all grids which overlap with 100km buffer from all sites but not overlapping with footprint
# and identify the one which are >= mean.
# I then add the whales which overlap with footprint (the one which should be displayed) to these grids which are >=mean
# to it, I calculate number of displayed whales, divide it by number of grids >=mean and add to these grids 

# first creating 100km buffer around each group of wind farms and then finding grids which overlap with the buffer
# which overlap with wind farm and the one which overlap with buffer but not wind farm
# I do it for each year of operation

# wind_foot24LL100 <- terra::buffer(wind_foot24LLsp, width=maxDist)
# wind_foot25LL100 <- terra::buffer(wind_foot25LLsp, width=maxDist)
# wind_foot26LL100 <- terra::buffer(wind_foot26LLsp, width=maxDist)
# wind_foot27LL100 <- terra::buffer(wind_foot27LLsp, width=maxDist)
# wind_foot28LL100 <- terra::buffer(wind_foot28LLsp, width=maxDist)
wind_foot29LL100 <- terra::buffer(wind_foot29LLsp, width=maxDist)
# wind_foot30LL100 <- terra::buffer(wind_foot30LLsp, width=maxDist)
# wind_foot31LL100 <- terra::buffer(wind_foot31LLsp, width=maxDist)

plot(studySite)
plot(wind_foot29LL100, add=T)

# putting all these footprint objects to a list so I can loop over them later
# it made sense when we were planning to run different windfarm operating at each year 
# I leave the commented code for now 

# wind_footLLBuffer <- ls(pattern = "LL100")
# wind_footLLBuffer <- lapply(wind_footLLBuffer, function(x) {get(x)})
# 
# wind_footLL <- ls(pattern = "LLsp")
# wind_footLL <- lapply(wind_footLL, function(x) {get(x)})

wind_footLLBuffer <- wind_foot29LL100
wind_footLL <- wind_foot29LLsp

# overlapping with footprint

temp <- as.data.frame(rasterize(wind_footLL, studySite, field=1), xy=T)
intersFP <- which(temp$layer==1)

# overlapping with buffer

temp <- as.data.frame(rasterize(wind_footLLBuffer, studySite, field=1), xy=T)
intersBuffer <- which(temp$layer==1)

# overlapping with buffer but not footprint

intersBuffer_noFP <- setdiff(intersBuffer, intersFP)

# assigning years to each realization map, needed for later in the code


for (i in 1:length(Dens_SS_operationDF)){
  for (k in 1:NmonitoringYears){
    rr <- k*N - N
    for (l in c((rr+1):(k*N))){
    Dens_SS_operationDF[[i]][[l]]$Year <- k
    }
  }
}

# For H6 first I find all grids which overlap with windfarm footprint (all sites) and then I multiple the density there by 2
# this is not adjusted to a loop if different farms operate at different year


for (i in 1:nmonthsOperation){		
  for (j in 1:length(Dens_SS_operationDF[[i]])){ #
    
    ## H6
    # it is easy peasy for H6, densities of whales at the footprint are multiplied by 2
    # i loop it over each year of operation and add each year as a new column
    
    Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H6 <- Dens_SS_operationDF[[i]][[j]]$DensOp_base
    Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H6[intersFP] <- Dens_SS_operationDF[[i]][[j]]$DensOp_base[intersFP]*2

    ## H7
    
    Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H7 <- Dens_SS_operationDF[[i]][[j]]$DensOp_base
    temp <- Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H7[intersFP]
    temp100 <- Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H7[intersBuffer_noFP]
    mean100 <- mean(temp100, na.rm=T)
    temp100[temp100>=mean100 & !is.na(temp100)] <- temp100[temp100>=mean100 & !is.na(temp100)] + sum(temp, na.rm=T)/length(temp100[temp100>=mean100 & !is.na(temp100)])
    Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H7[intersBuffer_noFP] <- temp100
    Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H7[intersFP] <- 0 # removing whales from the footprint
    
    ## H8, rest below
    Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H8 <- Dens_SS_operationDF[[i]][[j]]$DensOp_base
  }
  print(i)
}

##H8
# H8 was never tested for regional/operation scenarios so I have to define it here
# I will test pop_decline global decline in the population by reducing density in each grid by that number
# but this has to be incremental decline, so the start of each year is the end of the next one
# so I have to do it per year as well

for (i in 1:nmonthsOperation){
  for (l in (N+1):length(Dens_SS_operationDF[[i]])){

    Dens_SS_operationDF[[i]][[l]]$DensOp_AfterShifting_H8 <- Dens_SS_operationDF[[i]][[l-N]]$DensOp_AfterShifting_H8 - (Dens_SS_operationDF[[i]][[l-N]]$DensOp_AfterShifting_H8 * pop_decline)
    
  }
  print(i)
}


###############################################################################################################
# Step 6:	Transforming whale densities into cue densities
# using distance sampling equation
# D = n(1-c)/(cue rate * pi * detection_range^2). Assuming c=0 (false positive rate) so n would be
# n = D*(cue_rate * pi * detection_range^2) # units would be = #animals/km2 * km 2= #animals
#############################################################################################################

for (i in 1:nmonthsOperation){
  for (j in 1:length(Dens_SS_operationDF[[i]])){
    Dens_SS_operationDF[[i]][[j]]$CueDens_Base <- Dens_SS_operationDF[[i]][[j]]$DensOp_base * cue_rate 
    Dens_SS_operationDF[[i]][[j]]$CueDens_H6 <- Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H6 * cue_rate 
    Dens_SS_operationDF[[i]][[j]]$CueDens_H7 <- Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H7 * cue_rate 
    Dens_SS_operationDF[[i]][[j]]$CueDens_H8 <- Dens_SS_operationDF[[i]][[j]]$DensOp_AfterShifting_H8 * cue_rate 
    
  }
}

# if (species == "FW") {save(Dens_SS_operationDF,file="BaselineRealisationMaps/FinWhale/Dens_SS_operationDF1250.RData")}
# if (species == "NA") {save(Dens_SS_operationDF,file="BaselineRealisationMaps/NARW/Dens_SS_operationDF1250.RData")}
# if (species == "SW") {save(Dens_SS_operationDF,file="BaselineRealisationMaps/SeiWhale/Dens_SS_operationDF1250.RData")}
# if (species == "MW") {save(Dens_SS_operationDF,file="BaselineRealisationMaps/MinkeWhale/Dens_SS_operationDF1250.RData")}

##############################################################################################################################################################################################################################
# Step 6: Calculating cue density detected by each pam in the piling 
###########################################################################################################################################################################################################################

# here I have to calculate mean density of cues (per monthly map) at detection range for each pam
# because detection range is a circle and we are working on grid, for animals with small detection range, it will matter
# whether we calculate mean and mean of all the overlapping grids or as mean weighted by proportion of this grid being covered
# we, therefore, have to use the latter approach
# we ignore the fact that some pam detection areas overlap 

# first, creating buffer around each pam = detection range (has to be in meters hence * by 1000)
pamAllLLSSsp <- as(pamAllLLSS, "Spatial")
v <- vect(cbind(pamAllLLSSsp@coords[,1], pamAllLLSSsp@coords[,2]), crs=crs(pamAllLLSS))
pamBufferAll <- buffer(v, 1000*det_range)
plot(pamBufferAll)
points(v)
pamBufferAll <- as(pamBufferAll, "Spatial")
pamBufferAllsf <- sf::st_as_sf(pamBufferAll)
# pamAllLLSSsp <- sf::as_Spatial(pamAllLLSS)
# pamBufferAll <- terra::buffer(pamAllLLSSsp, width=1000*det_range)


plot(studySite)
plot(pamAllLLSS, add=T, pch=16, col="red")
plot(pamBufferAll, add=T)

pamAllLLSSDF <- data.frame(pamAllLLSS)
pamAllLLSSDF$id <- c(1:nrow(pamAllLLSSDF))

# making final plot for the report
# with species in the title and an example of distribution in January
# and pam devices colour coded by Poweron and no poweron

#janEx <- Dens_SS_operationDF[[1]][[15]]
#rm(Dens_SS_operationDF)

# getting land
land <- read_sf("GIS/ne_10m_land.shp")
land <- st_transform(land, newprojLatLong)
land <- st_crop(land,ee)

# buffPlot <- ggplot(janEx)+
#   geom_tile(data=janEx,aes(x=x, y=y,fill=DensOp_base))+
#   scale_fill_viridis(alpha=0.8,option="rocket", direction=-1,na.value="white")+
#   geom_sf(data=land, col="grey")+
#   geom_sf(data=wind_foot29LL, col="black", fill="orange", alpha=0.7)+
#   geom_sf(data=pamAllLLSS, col="darkblue", size=0.4)+
#   geom_sf(data=pamBufferAllsf, col="black", fill=NA)+
#   labs(fill=expression(paste("ind/25 ",km^{2})))+
#   ggtitle(species)+
#   theme_bw()+
#   theme(legend.position="bottom")

# ggsave(buffPlot, file=paste("Analysis/",species,"/",PAM,"PAM_withBuffers.png", sep=""),
#        width=15,
#        height=15,
#        units = "cm")
### I have to first make rasters with shifted cues for all hypothesis

# for Base, we only use N simulations, no matter number of monitoring months
CueDens_SS_OpRaster_H1 <- list()

for (i in 1:nmonthsOperation){
  CueDens_SS_OpRaster_H1[[i]] <- list()
  for (j in 1:N){
    
    temp <- Dens_SS_operationDF[[i]][[j]][,c("x","y","CueDens_Base")]
    CueDens_SS_OpRaster_H1[[i]][[j]] <- rasterFromXYZ(temp)
  }
  print(i)
}

# I will save it in case I just want to re-do power
save(CueDens_SS_OpRaster_H1,file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H1.RData", sep=""))

CueDens_SS_OpRaster_H6 <- list()
CueDens_SS_OpRaster_H7 <- list()
CueDens_SS_OpRaster_H8 <- list()

for (i in 1:nmonthsOperation){
  CueDens_SS_OpRaster_H6[[i]] <- list()
  CueDens_SS_OpRaster_H7[[i]] <- list()
  CueDens_SS_OpRaster_H8[[i]] <- list()
  
  for (j in 1:length(Dens_SS_operationDF[[i]])){
    
    temp <- Dens_SS_operationDF[[i]][[j]][,c("x","y","CueDens_H6")]
    CueDens_SS_OpRaster_H6[[i]][[j]] <- rasterFromXYZ(temp)
    
    temp <- Dens_SS_operationDF[[i]][[j]][,c("x","y","CueDens_H7")]
    CueDens_SS_OpRaster_H7[[i]][[j]] <- rasterFromXYZ(temp)
    
    temp <- Dens_SS_operationDF[[i]][[j]][,c("x","y","CueDens_H8")]
    CueDens_SS_OpRaster_H8[[i]][[j]] <- rasterFromXYZ(temp)
    
  }
  print(i)
}

# save your results for re-use
save(CueDens_SS_OpRaster_H6,file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H6.RData", sep=""))
save(CueDens_SS_OpRaster_H7,file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H7.RData", sep=""))
save(CueDens_SS_OpRaster_H8,file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H8.RData", sep=""))

# extracting density of cue overlaying with each pam as well as calculating proportion covered of each overlaying grid
# and then calculating the mean for each PAM

##### for testing alternative PAM designs, the code does not have to be rerun all the way till here but the rasters from
# above lines can be reloaded and code re-run from here
# you have to reload all settings and gis layers but no need to redistribute animals again
# but no need to rerun redistribution and changing to raster

load(file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H1.RData", sep=""))
load(file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H6.RData", sep=""))
load(file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H7.RData", sep=""))
load(file=paste("Analysis/",species,"/CueDens_SS_OpRaster_H8.RData", sep=""))



CalculateMeanCuePerPam <- function(df){
  
  ee <- exactextractr::exact_extract(
    df,
    pamBufferAll,
    coverage_area = F
  )
  ee2 <- do.call("rbind.data.frame",
                 lapply(ee, function(x) {
                   mean(x[, 1]*x[,2], na.rm = T)
                 }))
  colnames(ee2) <- c("Mean_Cue")
  ee2
}


# because this computer tents to crash from time to time, I will make a separate loop for each hypothesis, it takes about 4-5 min per version of hypothesis (set of rasters) for 
# piling hypothesis

parallel::stopCluster(cl = mycluster)

mycluster <- parallel::makeCluster(
  numCores, 
  type = "PSOCK",
  outfile=""
)
doParallel::registerDoParallel(cl = mycluster)	


# H1-operation
CueDens_Op_H1 <- list()

for (i in 1:nmonthsOperation){
  CueDens_Op_H1[[i]] <- list()

  CueDens_Op_H1[[i]] <-
    foreach (j = 1:length(CueDens_SS_OpRaster_H1[[i]])) %dopar% {
      CalculateMeanCuePerPam(CueDens_SS_OpRaster_H1[[i]][[j]])
    }
  print(i)
}



# H6 - H8
st <- Sys.time()
CueDens_H6 <- list()
CueDens_H7 <- list()
CueDens_H8 <- list()

for (i in 1:nmonthsOperation){
  CueDens_H6[[i]] <- list()
  CueDens_H7[[i]] <- list()
  CueDens_H8[[i]] <- list()

  CueDens_H6[[i]] <-
    foreach (j = 1:length(CueDens_SS_OpRaster_H6[[i]])) %dopar% {
      CalculateMeanCuePerPam(CueDens_SS_OpRaster_H6[[i]][[j]])
    }
  
  CueDens_H7[[i]] <-
    foreach (j = 1:length(CueDens_SS_OpRaster_H7[[i]])) %dopar% {
      CalculateMeanCuePerPam(CueDens_SS_OpRaster_H7[[i]][[j]])
    }
  
  CueDens_H8[[i]] <-
    foreach (j = 1:length(CueDens_SS_OpRaster_H8[[i]])) %dopar% {
      CalculateMeanCuePerPam(CueDens_SS_OpRaster_H8[[i]][[j]])
    }
  
  print(i)
}

save(CueDens_Op_H1,file=paste("Analysis/",species,"/",PAM,"CueDens_Op_H1.RData", sep=""))
save(CueDens_H6,file=paste("Analysis/",species,"/",PAM,"CueDens_H6.RData", sep=""))
save(CueDens_H7,file=paste("Analysis/",species,"/",PAM,"CueDens_H7.RData", sep=""))
save(CueDens_H8,file=paste("Analysis/",species,"/",PAM,"CueDens_H8.RData", sep=""))

load(file=paste("Analysis/",species,"/",PAM,"CueDens_Op_H1.RData", sep=""))
load(file=paste("Analysis/",species,"/",PAM,"CueDens_H6.RData", sep=""))
load(file=paste("Analysis/",species,"/",PAM,"CueDens_H7.RData", sep=""))
load(file=paste("Analysis/",species,"/",PAM,"CueDens_H8.RData", sep=""))

en <- Sys.time()
en-st

###################################################
# Step 7a: Calculating power using PG analysis
###################################################

# calculating distance of each PAM do the source
# for operation, source is just one and it is footprint

plot(studySite)
plot(pamAllLLSS, add=T, pch=16, col="red")
nrow(pamAllLLSS) #90 devices

t <- st_distance(pamAllLLSS, wind_foot29LL)
t <- apply(t,1,min)
pamAllLLSSDF$OpDist <- as.numeric(t)
pamAllLLSS$OpDist <- as.numeric(t)

# I dont want distance to source to be zero as I cant use log then so I cheat and add tiny number, data are in meters anyway
pamAllLLSSDF$OpDist[pamAllLLSSDF$OpDist==0] <- 0.01
pamAllLLSSDF$id <- as.numeric(rownames(pamAllLLSSDF))
  
# plotting

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

ggplot()+
  geom_sf(data=wind_foot29LL)+
  geom_sf(data=pamAllLLSS, aes(col=OpDist))+
  scale_colour_gradientn(colours = myPalette(100), limits=c(1, 20000))

# adding the close_sensor and close_sensor_andON
# close_sensors are the one within 100km from a windfarm
pamAllLLSSDF$close_sensor <- 1
pamAllLLSSDF$close_sensor[pamAllLLSSDF$OpDist>maxDist] <- 0
table(pamAllLLSSDF$close_sensor) #still non is further away than 100km


# PREPARING DATA FOR ALL HYPOTHESIS

# baseline (it only has 1 year)

for (i in 1:length(CueDens_Op_H1[[1]])){
  for (j in 1:nmonthsOperation){
    
    pamAllLLSSDF$OpDist
    CueDens_Op_H1[[j]][[i]]$Dist <- pamAllLLSSDF$OpDist
    CueDens_Op_H1[[j]][[i]]$pam <- as.numeric(rownames(CueDens_Op_H1[[j]][[i]]))
    CueDens_Op_H1[[j]][[i]]$Phase <- "Base"
    CueDens_Op_H1[[j]][[i]]$Month <- j
    CueDens_Op_H1[[j]][[i]]$Year <- 0
    CueDens_Op_H1[[j]][[i]]$close_sensor <- pamAllLLSSDF$close_sensor
    CueDens_Op_H1[[j]][[i]]$close_sensor_andON <- 0
    CueDens_Op_H1[[j]][[i]]$close_sensor_andON[CueDens_Op_H1[[j]][[i]]$close_sensor_andON==0] <- 0
  }
}

# Operation 

for (i in 1:length(CueDens_H6[[1]])){
  for (j in 1:nmonthsOperation){
    CueDens_H6[[j]][[i]]$Dist <- pamAllLLSSDF$OpDist
    CueDens_H6[[j]][[i]]$pam <- as.numeric(rownames(CueDens_H6[[j]][[i]]))
    CueDens_H6[[j]][[i]]$Phase <- "Operation"
    CueDens_H6[[j]][[i]]$Month <- j
    CueDens_H6[[j]][[i]]$Year <- ceiling(i/N)
    CueDens_H6[[j]][[i]]$close_sensor <- pamAllLLSSDF$close_sensor
    CueDens_H6[[j]][[i]]$close_sensor_andON <- 1
    CueDens_H6[[j]][[i]]$close_sensor_andON[CueDens_H6[[j]][[i]]$close_sensor_andON==0] <- 0
    
    CueDens_H7[[j]][[i]]$Dist <- pamAllLLSSDF$OpDist
    CueDens_H7[[j]][[i]]$pam <- as.numeric(rownames(CueDens_H7[[j]][[i]]))
    CueDens_H7[[j]][[i]]$Phase <- "Operation"
    CueDens_H7[[j]][[i]]$Month <- j
    CueDens_H7[[j]][[i]]$Year <- ceiling(i/N)
    CueDens_H7[[j]][[i]]$close_sensor <- pamAllLLSSDF$close_sensor
    CueDens_H7[[j]][[i]]$close_sensor_andON <- 1
    CueDens_H7[[j]][[i]]$close_sensor_andON[CueDens_H7[[j]][[i]]$close_sensor_andON==0] <- 0
    
    
    CueDens_H8[[j]][[i]]$Dist <- pamAllLLSSDF$OpDist
    CueDens_H8[[j]][[i]]$pam <- as.numeric(rownames(CueDens_H8[[j]][[i]]))
    CueDens_H8[[j]][[i]]$Phase <- "Operation"
    CueDens_H8[[j]][[i]]$Month <- j
    CueDens_H8[[j]][[i]]$Year <- ceiling(i/N)
    CueDens_H8[[j]][[i]]$close_sensor <- pamAllLLSSDF$close_sensor
    CueDens_H8[[j]][[i]]$close_sensor_andON <- 1
    CueDens_H8[[j]][[i]]$close_sensor_andON[CueDens_H8[[j]][[i]]$close_sensor_andON==0] <- 0
    
    
  }
}

save(CueDens_Op_H1,file=paste("Analysis/",species,"/",PAM,"CueDens_Op_H1.RData", sep=""))
save(CueDens_H6,file=paste("Analysis/",species,"/",PAM,"CueDens_H6.RData", sep=""))
save(CueDens_H7,file=paste("Analysis/",species,"/",PAM,"CueDens_H7.RData", sep=""))
save(CueDens_H8,file=paste("Analysis/",species,"/",PAM,"CueDens_H8.RData", sep=""))

# merging all months for each of the N realizations

BaseOp <- list()

# Note that we only have one year of Base
# so as input to the model, I need 12 month of Base and 12 X Nmonitoring years months for Operation
# the list below will have as many rows as number of months * number of PAMs
j=1
for (i in 1:length(CueDens_Op_H1[[1]])){
  j=1
  BaseOp[[i]] <- CueDens_Op_H1[[j]][[1]]

  for (j in 2:nmonthsOperation){
    
    BaseOp[[i]] <- rbind.data.frame(BaseOp[[i]],CueDens_Op_H1[[j]][[i]])
  }
}

# this should have the size of number of months * number of PAMs * NmonitoringYears
if (NmonitoringYears==5){
PG_H6 <- list()
PG_H7 <- list()
PG_H8 <- list()

for (k in 1:NmonitoringYears){
 rr <- k*N - N

j=1
for (i in 1:N){ #length(CueDens_H6[[1]])/NmonitoringYears)
  j=1
  PG_H6[[i+rr]] <- CueDens_H6[[j]][[i+rr]]
  PG_H7[[i+rr]] <- CueDens_H7[[j]][[i+rr]]
  PG_H8[[i+rr]] <- CueDens_H8[[j]][[i+rr]]
  
  for (j in 2:nmonthsOperation){
    
    PG_H6[[i+rr]] <- rbind.data.frame(PG_H6[[i+rr]],CueDens_H6[[j]][[i+rr]])
    PG_H7[[i+rr]] <- rbind.data.frame(PG_H7[[i+rr]],CueDens_H7[[j]][[i+rr]])
    PG_H8[[i+rr]] <- rbind.data.frame(PG_H8[[i+rr]],CueDens_H8[[j]][[i+rr]])
  }
}
}


# As I need all years of Operation in one file, the loop before should do it
# at the end I should have a list of N elements
# and finally merging Base and Operation
PG_H6_full <- list()
PG_H7_full <- list()
PG_H8_full <- list()

for (i in 1:N){
  PG_H6_full[[i]]<- rbind.data.frame(BaseOp[[i]],PG_H6[[i]],PG_H6[[i+1*N]],PG_H6[[i+2*N]],PG_H6[[i+3*N]],PG_H6[[i+4*N]])
  PG_H7_full[[i]]<- rbind.data.frame(BaseOp[[i]],PG_H7[[i]],PG_H7[[i+1*N]],PG_H7[[i+2*N]],PG_H7[[i+3*N]],PG_H7[[i+4*N]])
  PG_H8_full[[i]]<- rbind.data.frame(BaseOp[[i]],PG_H8[[i]],PG_H8[[i+1*N]],PG_H8[[i+2*N]],PG_H8[[i+3*N]],PG_H8[[i+4*N]])
  
}

PG_H6 <- PG_H6_full
PG_H7 <- PG_H7_full
PG_H8 <- PG_H8_full

}

# final step, merge it with x and y coordinates of each pam

pamAllLLSSDFred <- pamAllLLSSDF %>%
  select(c(id, LATITUDE, LONGITUDE))

colnames(pamAllLLSSDFred) <- c("pam","x","y")

for (m in 1:length(PG_H6)){
  PG_H6[[m]] <- PG_H6[[m]] %>%
    left_join(pamAllLLSSDFred, by=c("pam"))
  
  PG_H7[[m]] <- PG_H7[[m]] %>%
    left_join(pamAllLLSSDFred, by=c("pam"))
  
  PG_H8[[m]] <- PG_H8[[m]] %>%
    left_join(pamAllLLSSDFred, by=c("pam"))
  #print(m)
}

# finally, producing data frame for H1 - no effect.
# here I shuffle baseline NmonitoringYears times for each month

# THE BELOW HAS TO BE ADJUSTED IF WE HAVE SOME SENSORS>100KM
PG_H1 <- PG_H6
for (i in 1:length(PG_H1)){
  PG_H1[[i]] <- PG_H1[[i]][PG_H1[[i]]$Phase=="Base",]
}

PG_H1B <- PG_H1

for (i in 1:length(PG_H1)){
  for (k in 1:NmonitoringYears){
    randi <- sample(c(1:N),1)
    temp <- PG_H1[[randi]]
    temp$Phase <- "Operation"
    PG_H1B[[i]] <- rbind(PG_H1B[[i]],temp)
  }

  npams <- nrow(pamAllLLSSDF)
  PG_H1B[[i]]$Year <- c(rep(0,12*npams),
                       rep(1,12*npams),
                       rep(2,12*npams),
                       rep(3,12*npams),
                       rep(4,12*npams),
                       rep(5,12*npams))
  
  PG_H1B[[i]]$close_sensor_andON <- c(rep(0,12*npams),
                                     rep(1,12*npams*NmonitoringYears))
}

PG_H1 <- PG_H1B
save(PG_H6,file=paste("Analysis/",species,"/",PAM,"PG_H6.RData", sep=""))
save(PG_H7,file=paste("Analysis/",species,"/",PAM,"PG_H7.RData", sep=""))
save(PG_H8,file=paste("Analysis/",species,"/",PAM,"PG_H8.RData", sep=""))
save(PG_H1,file=paste("Analysis/",species,"/",PAM,"PG_H1.RData", sep=""))

load(file=paste("Analysis/",species,"/",PAM,"PG_H6.RData", sep=""))
load(file=paste("Analysis/",species,"/",PAM,"PG_H7.RData", sep=""))
load(file=paste("Analysis/",species,"/",PAM,"PG_H8.RData", sep=""))
load(file=paste("Analysis/",species,"/",PAM,"PG_H1.RData", sep=""))

# function specifying model for H8 where year has to be a continuous variable
# below are functions using bam and gam function. Bam is much faster
extractingFromGamH8 <- function(InputData){
  ndf <- InputData
  
  # baseline would assume no windfarm operating
  # what we want is the significance of s(Dist, by=close_sensor_andON) in the power analysis and year
  
  tryCatch({
    suppressWarnings(mod <- mgcv::gam(as.integer(Mean_Cue) ~  te(x,y) + 
                                        s(Dist, k=5) +
                                        s(Dist, by=close_sensor, k=5) +
                                        s(Dist, by=close_sensor_andON, m=1, pc=maxDist, k=5) + 
                                        as.factor(Month) + 
                                        s(Year, k=3), data=ndf,family = quasipoisson(link = "log"), method="REML"));
    
    
     msum <- summary(mod)
  }
  , error = function(e) {an.error.occured <<- TRUE}
  )
  ifelse (exists("msum")==TRUE,
          {
            pDist_CloseOn <- msum$s.table[19]
            pYear<- msum$s.table[20]
          },
          {
            pDist_CloseOn <- NA
            pYear<- NA
          })
  
  return(c(pDist_CloseOn,pYear))
}

extractingFromBamH8 <- function(InputData){
  ndf <- InputData
  
  # baseline would assume no windfarm operating
  # what we want is the significance of s(Dist, by=close_sensor_andON) in the power analysis and year
  
  tryCatch({
    suppressWarnings(mod <- mgcv::bam(as.integer(Mean_Cue) ~  te(x,y) + 
                                        s(Dist, k=5) +
                                        s(Dist, by=close_sensor, k=5) +
                                        s(Dist, by=close_sensor_andON, m=1, pc=maxDist, k=5) + 
                                        as.factor(Month) + 
                                        s(Year, k=3), data=ndf,family = quasipoisson(link = "log"), method="fREML",
                                      discrete=TRUE));
    
    
    msum <- summary(mod)
  }
  , error = function(e) {an.error.occured <<- TRUE}
  )
  ifelse (exists("msum")==TRUE,
          {
            pDist_CloseOn <- msum$s.table[19]
            pYear<- msum$s.table[20]
          },
          {
            pDist_CloseOn <- NA
            pYear<- NA
          })
  
  return(c(pDist_CloseOn,pYear))
}

# for H1,6 and 7 year should be a factor but in H1 we want phase to be a factor too
extractingFromGamH67 <- function(InputData){
  ndf <- InputData
  
  tryCatch({
    suppressWarnings(mod <- mgcv::gam(as.integer(Mean_Cue) ~  te(x,y) + 
                                        s(Dist, k=5) +
                                        s(Dist, by=close_sensor, k=5) +
                                        s(Dist, by=close_sensor_andON, m=1, pc=maxDist, k=5) + 
                                        as.factor(Month) + 
                                        as.factor(Year), data=ndf,family = quasipoisson(link = "log"), method="REML"));
    
    
    msum <- summary(mod)
  }
  , error = function(e) {an.error.occured <<- TRUE}
  )
  ifelse (exists("msum")==TRUE,
          {
            pDist_CloseOn <- msum$s.table[16]
          },
          {
            pDist_CloseOn <- NA
          })
  
  return(c(pDist_CloseOn))
}

extractingFromBamH67 <- function(InputData){
  ndf <- InputData
  
  tryCatch({
    suppressWarnings(mod <- mgcv::bam(as.integer(Mean_Cue) ~  te(x,y) + 
                                        s(Dist, k=5) +
                                        s(Dist, by=close_sensor, k=5) +
                                        s(Dist, by=close_sensor_andON, m=1, pc=maxDist, k=5) + 
                                        as.factor(Month) + 
                                        as.factor(Year), data=ndf,family = quasipoisson(link = "log"), method="fREML",
                                      discrete=TRUE));
    
    
    msum <- summary(mod)
  }
  , error = function(e) {an.error.occured <<- TRUE}
  )
  ifelse (exists("msum")==TRUE,
          {
            pDist_CloseOn <- msum$s.table[16]
          },
          {
            pDist_CloseOn <- NA
          })
  
  return(c(pDist_CloseOn))
}

extractingFromBamH1 <- function(InputData){
  ndf <- InputData
  
  tryCatch({
    suppressWarnings(mod <- mgcv::bam(as.integer(Mean_Cue) ~  te(x,y) + 
                                        s(Dist, k=5) +
                                        s(Dist, by=close_sensor, k=5) +
                                        s(Dist, by=factor(close_sensor_andON), m=1, pc=maxDist, k=5) + 
                                        as.factor(Month) + 
                                        as.factor(Year), data=ndf,family = quasipoisson(link = "log"), method="fREML",
                                      discrete=TRUE));
    
    
    msum <- summary(mod)
  }
  , error = function(e) {an.error.occured <<- TRUE}
  )
  ifelse (exists("msum")==TRUE,
          { 
            pDist_CloseOn0 <- msum$s.table[19]
            pDist_CloseOn1 <- msum$s.table[20]
          },
          {
            pDist_CloseOn0 <- NA
            pDist_CloseOn1 <- NA
          })
  
  return(c(pDist_CloseOn0,pDist_CloseOn1))
}
# parallel::stopCluster(cl = mycluster)
# 
# mycluster <- parallel::makeCluster(
#   numCores, 
#   type = "PSOCK",
#   outfile=""
# )
# doParallel::registerDoParallel(cl = mycluster)	

# ok, the below does not work in parallel system and I dont know why
# I think it has something to do with using te(x,y)
# ss <- Sys.time()
# pss_H6_all<- foreach (i = 1:N)  %dopar% { extractingFromGam(PG_H6[[i]])}
# pss_H7_all<- foreach (i = 1:N)  %dopar% { extractingFromGam(PG_H7[[i]])}
# 
# e <- Sys.time()
# e-s
pss_H1_all <- list()
pss_H6_all <- list()
pss_H7_all <- list()
pss_H8_all <- list()
# doing all three at the same time crushes the machine

for (i in 1:length(PG_H1)){#length(PG_H1)
  pss_H1_all[[i]] <- extractingFromBamH1(PG_H1[[i]])
  print(i)
}
save(pss_H1_all, file=paste("Analysis/",species,"/",PAM,"pss_H1_all.RData", sep=""))


for (i in 1:length(PG_H6)){
  pss_H6_all[[i]] <- extractingFromBamH67(PG_H6[[i]])
  print(i)
}

save(pss_H6_all, file=paste("Analysis/",species,"/",PAM,"pss_H6_all.RData", sep=""))

for (i in 1:length(PG_H7)){
  pss_H7_all[[i]] <- extractingFromBamH67(PG_H7[[i]])
  print(i)
}

save(pss_H7_all, file=paste("Analysis/",species,"/",PAM,"pss_H7_all.RData", sep=""))

for (i in 1:length(PG_H8)){
  pss_H8_all[[i]] <- extractingFromBamH8(PG_H8[[i]])
  print(i)
}

save(pss_H8_all, file=paste("Analysis/",species,"/",PAM,"pss_H8_all.RData", sep=""))

# calculating power: proportion of models where:
# For H6 and H7 interaction (s(Dist, by=close_sensor_andON) is significant, we dont care about the rest
# For H8 interaction is not significant, but year is
# For H1 interaction is significant but should not 

PG_power67 <- function(x){
  x1 <- do.call("rbind.data.frame",x)
  colnames(x1) <- c("pDist_CloseOn")
  x1$signDist_CloseOn <- 0
  x1$signDist_CloseOn[x1$pDist_CloseOn<=alpha] <- 1
  return(length(x1$signDist_CloseOn[x1$signDist_CloseOn==1])/length(x))
  
}

PG_power1 <- function(x){
  x1 <- do.call("rbind.data.frame",x)
  colnames(x1) <- c("pDist_CloseOn0","pDist_CloseOn1")
  x1$sign <- 0
  x1$sign[x1$pDist_CloseOn0>=alpha & x1$pDist_CloseOn1<alpha] <- 1
  return(length(x1$sign[x1$sign==1])/length(x))
  
}

PG_power8 <- function(x){
  x1 <- do.call("rbind.data.frame",x)
  colnames(x1) <- c("pDist_CloseOn","pYear")
  x1$sign <- 0
  x1$sign[x1$pDist_CloseOn>alpha & x1$pYear<=alpha] <- 1
  return(length(x1$sign[x1$sign==1])/length(x))
  
}

pss_H1_all2 <- PG_power1(pss_H1_all)
pss_H6_all2 <- PG_power67(pss_H6_all)
pss_H7_all2 <- PG_power67(pss_H7_all)
pss_H8_all2 <- PG_power8(pss_H8_all)

HS <- c("H1","H6","H7","H8")
resAll <- c(pss_H1_all2, 
            pss_H6_all2,
            pss_H7_all2, 
            pss_H8_all2)

PG_results <- cbind.data.frame(HS,resAll)
PG_results$species <- species
PG_results$PAM <- PAM
PG_results

save(PG_results, file=paste("Analysis/",species,"/",PAM,"PG_results_allPAM.RData", sep=""))
Sys.time()
