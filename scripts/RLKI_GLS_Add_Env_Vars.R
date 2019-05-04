# Script to add environmental variables to GLS tracking data
# from Red-legged Kittiwakes in the Bering Sea and western North Pacific Ocean
# The oceanographic datasets must be in raster form stored locally.
#
# Masters Research
#
# Abram Fleishman
# Updated 23 Dec 2018
rm(list= ls())

# Load Packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(raster)
library(velox)
library(sp)
library(rgeos)
library(ncdf4)
library(stringr)
library(maptools)
library(marmap)
library(sf)

# Sett point buffer in meters
bufferM<-185000

# Set main dir: Sys.info()[6] is the username for the computer.
# fill in the "" with your user name
if(Sys.info()[["user"]]=="abramfleishman") {
  Migration<-"~/Dropbox (CMI)/RLKI_Migrations"
  BS_ENV<-"~/Dropbox (CMI)/RLKI_Large_Data/Bering_Environmental_Data"
  LandMaskPath<-"~/Dropbox (CMI)/RLKI_Large_Data/Bering_Environmental_Data/LandMask"
  OutDataPath<-paste0("~/Dropbox (CMI)/RLKI_Large_Data/tracks_buldir_vars_",bufferM/1000,"km_25.rds")
  DataInPath<-paste0(Migration, "/Data/processedLocs/probGLS/RLKI_GLS_Buldir_allyearsalldata.rds")
  DataInPath_STG<-paste0(Migration, "/Data/processedLocs/probGLS/RLKI_GLS_201617_forAMNWR_2018-10-30.rda")
}


# OutDataPath<-paste0("//NAS1/Public/all_tracks_11-17_vars_",bufferM/1000,"km_25.rds")
# load tracks file
# tracks<-readRDS(paste0(Migration,"/Data/processedLocs/probGLS/RLKI_GLS_alltracks.rda"))
tracks<-readRDS(DataInPath) %>% bind_rows(readRDS(DataInPath_STG))
#
# # make date col
tracks$date<-as.Date(tracks$datetime)
tracks$week<-week(tracks$datetime)
# n_distinct(tracks$loggerID)
# # tracks<-tracks %>% filter(month%in%c(9,10,11,12,1,2,3,4,5,6))
# # Add ENV vars ------------------------------------------------------------
# # Data from the CMEMS Global Analysis Forecast ----------------------------
#
# # File dates
#
# # 024_1494782985210.nc 2013-14
# # 024_1494711054915.nc 2015-16
# # 024_1494776551825.nc 2014-15
# # 024_1494798426931.nc 2010-06-01 to 2010-12-31
# # 024_1494798789193.nc 2010-12-31 2011-06-01
# # 024_1501874976336.nc 2015-16 noz uo vo
# # 024_1501877327442.nc 2016-17 noz uo vo
# # 024_1545798301435.nc 2016-17 noz uo vo
# # 024_1546115422194.nc 2016-17 noz uo vo
# # 024_1546116043746.nc 2016-17 noz uo vo
# # 024_1546116288864.nc 2016-17 noz uo vo
# # 024_1546116465020.nc 2016-17 noz uo vo
#
# # Bring in the current data uo,vo, and ssh
# uo11a<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494798426931.nc") ,varname="uo")
# uo11b<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494798789193.nc") ,varname="uo")
# uo14<- stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494782985210.nc") ,varname="uo")
# uo15<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494776551825.nc") ,varname="uo")
# uo16<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1501874976336.nc") ,varname="uo")
# uo17<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1501877327442.nc") ,varname="uo")
# extent(uo17)
# shift_crop<-function(x,y){
#   x1 <- crop(x, extent(-180, 0, 34, 67))
#   x1<-shift(x1,x=360)
#   x2 <- crop(x, extent(0, 180, 34, 67))
#   m <- merge(x1, x2)
#   m<-crop(m,y = extent(y))
#   names(m)<-names(x)
#   return(m)
# }
#
# uo17a<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1545798301435.nc") ,varname="uo")
# uo17a<-shift_crop(uo17a,uo17)
# uo17b<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546115422194.nc") ,varname="uo")
# uo17b<-shift_crop(uo17b,uo17)
#
# uo17c<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116043746.nc") ,varname="uo")
# uo17c<-shift_crop(uo17c,uo17)
#
# uo17d<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116288864.nc") ,varname="uo")
# uo17d<-shift_crop(uo17d,uo17)
#
# uo17e<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116465020.nc") ,varname="uo")
# names(uo17e)
# uo17e<-shift_crop(uo17e,uo17)
#
# uo18<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1556952426071.nc") ,varname="uo")
# names(uo18)
# uo18<-shift_crop(uo18,uo17)
#
# extent(uo17)
# extent(uo17a)
# uo<-stack(uo11a,uo11b,uo14,uo15,uo16,uo17,uo17a,uo17b,uo17c,uo17d,uo17e,uo18)
# # writeRaster(uo,filename =paste0(BS_ENV,'/global-analysis-forecast-phy-001-024/uo_compiled.nc' ))
# saveRDS(uo,paste0(BS_ENV,'/global-analysis-forecast-phy-001-024/uo_compiled.rds' ))
# tail(names(uo))
# vo11a<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494798426931.nc") ,varname="vo")
# vo11b<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494798789193.nc") ,varname="vo")
# vo14<- stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494782985210.nc") ,varname="vo")
# vo15<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494776551825.nc") ,varname="vo")
# vo16<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1501874976336.nc") ,varname="vo")
# vo17<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1501877327442.nc") ,varname="vo")
#
# vo17a<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1545798301435.nc") ,varname="vo")
# vo17a<-shift_crop(vo17a,vo17)
#
# vo17b<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546115422194.nc") ,varname="vo")
# vo17b<-shift_crop(vo17b,vo17)
#
# vo17c<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116043746.nc") ,varname="vo")
# vo17c<-shift_crop(vo17c,vo17)
#
# vo17d<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116288864.nc") ,varname="vo")
# vo17d<-shift_crop(vo17d,vo17)
#
# vo17e<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116465020.nc") ,varname="vo")
# vo17e<-shift_crop(vo17e,vo17)
#
# vo18<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1556952426071.nc") ,varname="vo")
# vo18<-shift_crop(vo18,vo17)
#
# vo<-stack(vo11a,vo11b,vo14,vo15,vo16,vo17,vo17a,vo17b,vo17c,vo17d,vo17e,vo18)
# saveRDS(vo,paste0(BS_ENV,'/global-analysis-forecast-phy-001-024/vo_compiled.rds' ))
#
# ssh11a<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494798426931.nc") ,varname="zos")
# ssh11b<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494798789193.nc") ,varname="zos")
# ssh14<- stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494782985210.nc") ,varname="zos")
# ssh15<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1494776551825.nc") ,varname="zos")
# ssh16<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1501874976336.nc") ,varname="zos")
# ssh17<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1501877327442.nc") ,varname="zos")
#
# ssh17a<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1545798301435.nc") ,varname="zos")
# ssh17a<-shift_crop(ssh17a,ssh17)
#
# ssh17b<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546115422194.nc") ,varname="zos")
# ssh17b<-shift_crop(ssh17b,ssh17)
#
# ssh17c<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116043746.nc") ,varname="zos")
# ssh17c<-shift_crop(ssh17c,ssh17)
#
# ssh17d<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116288864.nc") ,varname="zos")
# ssh17d<-shift_crop(ssh17d,ssh17)
#
# ssh17e<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1546116465020.nc") ,varname="zos")
# ssh17e<-shift_crop(ssh17e,ssh17)
#
# ssh18<-stack(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/global-analysis-forecast-phy-001-024_1556952426071.nc") ,varname="zos")
# ssh18<-shift_crop(ssh18,ssh17)
#
# ssh<-stack(ssh11a,ssh11b,ssh14,ssh15,ssh16,ssh17,ssh17a,ssh17b,ssh17c,ssh17d,ssh17e,ssh18)
# saveRDS(ssh,paste0(BS_ENV,'/global-analysis-forecast-phy-001-024/ssh_compiled.rds' ))
#
#
# # # To resample vars --------------------------------------------------------
# #
# # To aggragate raster cells to downsample resulution in chunks of 100 layers
# agg <- function( stack_s = ssh,  fun = "mean", fact = 3 ){
#   if(dim(stack_s)[3]>100) {
#     idx<-seq(1,dim(stack_s)[3],100)
#   } else {
#     idx <- c(1,dim(stack_s)[3])
#   }
#
#   for(i in 1:(length(idx) -1) ){
#     print(paste("index:", i," of ", length(idx) -1 ))
#     stack1 <- stack_s[[ idx[i]:lead( idx )[i] ]]
#     if( i==1) stack1_down <- aggregate( stack1, fact = fact, fun = fun, na.rm = T )
#     if( i>1) stack1_down <- stack(stack1_down,aggregate( stack1, fact = fact, fun = fun, na.rm = T ))
#   }
#   return(stack1_down)
# }
#
# ssh<-agg(stack_s = ssh,fun = "mean",fact = 3)
# uo<-agg(stack_s = uo,fun = "mean",fact = 3)
# vo<-agg(stack_s = vo,fun = "mean",fact = 3)
#
# saveRDS(ssh,paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/ssh_0.25.rds"))
# saveRDS(uo,paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/uo_0.25.rds"))
# saveRDS(vo,paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/vo_0.25.rds"))
names(tracks)
# Load the resampled vars (all 0.25 degre resolution)
ssh<-readRDS(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/ssh_0.25.rds"))
uo<-readRDS(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/uo_0.25.rds"))
vo<-readRDS(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/vo_0.25.rds"))

# Create a date table for the global-analysis-forecast-phy vars
DateTable<-NULL

for(i in 1:length(ssh@layers)) DateTable[i]<-ssh[[i]]@data@names

DateTable<-data.frame(names=DateTable,  DateTime=ymd_hms(gsub("X","",DateTable))) %>%
  mutate(
    year=year(DateTime),month=month(DateTime),day=day(DateTime),
    hour=hour(DateTime),minute=minute(DateTime),second=second(DateTime),
    Date=as.Date(DateTime)
  )

# Loop for global-analysis-forecast-phy vars ---------------------------------

# Create empty columns to populate
tracks$eke<-NA
tracks$eke_g<-NA

tracks$uo<-NA
tracks$vo<-NA

tracks$ssh<-NA
tracks$ssh_g<-NA

# Get the unique dates in the tracks data
dates<-sort(unique(tracks$date))

# go through each date in the tracks data, extract the vars at each of the
# points with a buffer, aggragate the buffered data to end up with a single
# value for each point.
#
# Currently using the scale of error 925km to buffer and taking the mean.
# install.packages("velox")

for(i in 1:length(dates)){
  print(dates[i])

  datex<-as.character(DateTable$names[DateTable$Date==dates[i]])[1]
  if(length(datex)<1|is.na(datex)) next

  pt<-data.frame(lon360=tracks$lon360[which(tracks$date==dates[i]&is.na(tracks$ssh))],
                 lat=tracks$lat[which(tracks$date==dates[i]&is.na(tracks$ssh))])

  coordinates(pt)<-cbind(tracks$lon360[which(tracks$date==dates[i]&is.na(tracks$ssh))],
                         tracks$lat[which(tracks$date==dates[i]&is.na(tracks$ssh))])
  crs(pt)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180")

  pt_sp<-spTransform(pt,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))

  spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)

  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
  spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))
  # plot(spdf,add=T)
  # plot(ssh[[datex]])
  # ssh
  vx<- velox(ssh[[datex]])
  tracks$ssh[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  if(sum(is.na(tracks$ssh[which(tracks$date==dates[i])]))>0) stop("boom")
  vx<- velox(terrain(ssh[[datex]],opt = "slope",neighbors = 8))
  tracks$ssh_g[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  # eke
  vx<-velox(1/2*(vo[[datex]]^2+uo[[datex]]^2))
  tracks$eke[which(tracks$date==dates[i])]<-vx$extract(spdf, fun = function(x) mean(x, na.rm = T))

  vx<-velox(terrain(( 1/2 * (vo[[datex]]^2 + uo[[datex]]^2)),opt = "slope",neighbors = 8))
  tracks$eke_g[which(tracks$date == dates[i])]<-vx$extract(spdf, fun = function(x) mean(x, na.rm = T))

  # uo
  vx<- velox(uo[[datex]])
  tracks$uo[which(tracks$date==dates[i])]<-vx$extract(spdf, fun = function(x) mean(x, na.rm = T))

  # vo
  vx<- velox(vo[[datex]])
  tracks$vo[which(tracks$date==dates[i])]<-vx$extract(spdf, fun = function(x) mean(x, na.rm = T))

}

saveRDS(tracks,OutDataPath)

# tracks<-readRDS(OutDataPath)
head(tracks)

# SST from the NOAA 0.25° x 0.25° blended data set ----------------------------
# NOAA High-resolution Blended Analysis: Daily Values using AVHRR only

# read in all the yearly files and stack them
#
# sst2010<-stack(file.path(LandMaskPath,'sst.day.mean.2010.v2.nc'))
# sst2011<-stack(file.path(LandMaskPath,'sst.day.mean.2011.v2.nc'))
# sst2012<-stack(file.path(LandMaskPath,'sst.day.mean.2012.v2.nc'))
# sst2013<-stack(file.path(LandMaskPath,'sst.day.mean.2013.v2.nc'))
# sst2014<-stack(file.path(LandMaskPath,'sst.day.mean.2014.v2.nc'))
# sst2015<-stack(file.path(LandMaskPath,'sst.day.mean.2015.v2.nc'))
sst2016<-stack(file.path(LandMaskPath,'sst.day.mean.2016.v2.nc'))
sst2017<-stack(file.path(LandMaskPath,'sst.day.mean.2017.v2.nc'))
sst2018<-stack(file.path(LandMaskPath,'sst.day.mean.2018.v2.nc'))
# sst2017
# sst2017a
# nc_open(file.path(LandMaskPath,'sst.day.mean.2017.nc'))
# nc_open(file.path(LandMaskPath,'sst.day.mean.2017.v2.nc'))


# Combine stacks
sst<-stack(sst2016,sst2017,sst2018)
tracks %>%
  ungroup %>%
  filter(is.na(ssh)) %>% pull(datetime) %>% as.Date %>% unique %>% sort

ggplot(tracks %>%
         ungroup %>%
         filter(!is.na(ssh)) %>%
         ungroup #%>% sample_n(30000)
)+geom_point(aes(lon360,lat,col=ssh))
# Loop to add sst and sst_g
tracks$sst<-NA
tracks$sst_g<-NA
for(i in 1:length(dates)){
  print(dates[i])
  datex<-paste0("X",year(dates[i]),".",
                str_pad(string = month(dates[i]),width = 2,pad = "0",side = "left"),".",
                str_pad(string = day(dates[i]),width = 2,pad = "0",side = "left"))

  pt<-data.frame(lon360=tracks$lon360[which(tracks$date==dates[i]&is.na(tracks$sst))],
                 lat=tracks$lat[which(tracks$date==dates[i]&is.na(tracks$sst))])

  coordinates(pt)<-cbind(tracks$lon360[which(tracks$date==dates[i]&is.na(tracks$sst))],
                         tracks$lat[which(tracks$date==dates[i]&is.na(tracks$sst))])
  crs(pt)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180")


  pt_sp<-spTransform(pt,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))

  spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)

  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
  spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))

  # sst
  vx<- velox(sst[[datex]])
  tracks$sst[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  vx<- velox(terrain(sst[[datex]],opt = "slope",neighbors = 8))
  tracks$sst_g[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))
}

saveRDS(tracks,OutDataPath)
# tracks<-readRDS(OutDataPath)
rm(list=ls(pattern = "sst"))
rm(list=ls(pattern = "uo"))
rm(list=ls(pattern = "vo"))
rm(list=ls(pattern = "ssh"))
rm(list=ls(pattern = "uo"))

# ice shp predata prep ----------------------------------------------------
# ftp://sidads.colorado.edu/DATASETS/NOAA/G02186/shapefiles/4km/
# http://nsidc.org/data/docs/noaa/g02186_masie/index.html

# For bulk unzipping
zips<-list.files(paste0(BS_ENV,"/Sea_Ice_shps/zip"),pattern = "zip",recursive = T,full.names = T)
# head(zips)
for(i in zips) unzip(zipfile = i,exdir = paste0(BS_ENV,"/Sea_Ice_shps/unzip"))


shps<-list.files(paste0(BS_ENV,"/Sea_Ice_shps/unzip"),pattern="shp$",recursive = T,full.names = T)

shpdate<-str_extract(shps,"[0-9]{7}")
shpyear<-str_extract(shpdate,"[0-9]{4}")
shpdoy<-str_extract(shpdate,"[0-9]{3}$")
missingice<-data.frame(year=c(2010,2010,2010,2011,2011,2012,2012,2014,2014,2014,2014,2015,2016,2017,2017),
                       doy =c( 010, 031, 181, 024, 040, 250, 252, 241, 290, 293, 294, 108, 195,307,308))


# Function to open ice for a year and doy
openIce<-function(year=2017,yearT=2017,doy=1,path=paste0(BS_ENV,"/Sea_Ice_shps/unzipped")){
  # There are some missing ice files and to make sure it doesnt break I made this
  # while statement
  while(TRUE%in%(year==missingice$year&doy==missingice$doy)) { doy=doy-1 }
  # path to read from
  Path<-paste0(path,"/masie_ice_r00_v01_",paste0(year,str_pad(doy,3,pad="0")),"_4km.shp")

  # read Shp
  # ice<-rgdal::readOGR(dsn = dirname(Path),layer = tools::file_path_sans_ext(basename(Path)))
  # ice<-readShapePoly(fn =tools::file_path_sans_ext(Path) )
  # plot(ice)
  # crs(ice)<- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84
  # +towgs84=0,0,0"
  #
  # # Project data into laea
  # ice <- spTransform(ice, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))

  if(!file.exists(Path)) next
  ice<-sf::read_sf(dirname(Path),layer =tools::file_path_sans_ext(basename(Path))) %>%
    sf::st_transform( "+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs") %>%
    as("Spatial")

  # add the coords with propernaes
  # ice$x_laea<-ice$long
  # ice$y_laea<-ice$lat
  ice$doy<-doy
  ice$year=year
  ice$yearT<-yearT
  return(ice)
}

# Make a spatial Points DF with the tracks
tracks_sp <- SpatialPointsDataFrame(coords=cbind(tracks$lon,tracks$lat),
                                    data=as.data.frame(tracks),
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

tracks_sp<-spTransform(tracks_sp,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                   +ellps=WGS84 +units=m +no_defs"))

# Loop through and open ice each day then calculate distance from ice to point


tracks$dist2ice<-NA
Dates<-unique(tracks$date)
for(j in 1:length(Dates)){
  # make subset of data for each day
  tracks1<-subset(tracks_sp,date==Dates[j])
  # this opens all the ice files for the date
  ice<-openIce(year=unique(year(unique(tracks1$date))),yearT=unique(tracks1$yearT),doy = yday(unique(tracks1$date)),
               path=paste0(BS_ENV,"/Sea_Ice_shps/unzip"))

  # calculate the distance from every point to every edge
  tempdist<-  gDistance(tracks1,ice,byid = T)

  # replace dist2ice in the ORIGINAL tracks dataframe with the minimum distance
  tracks$dist2ice[tracks$date==Dates[j]]<-apply(tempdist, 2, function(x) min(x, na.rm = TRUE))

  print(Dates[j])

}

saveRDS(tracks,OutDataPath)

# Bathymetry     --------------------------------------------------------------
#
# Download NOAA ETOPO1 data
# res 15 = 15 minutes = .25 degrees
aleu <- getNOAA.bathy(130, -140, 30, 75, resolution = 15,
                      antimeridian = T,keep = T) #resolution 1.85km*27=49.9km2
# Make it a raster
bathy<-as.raster(aleu)
plot(bathy)
bathy[bathy>=0]<-NA
head(tracks)
bathyvx<- velox(bathy)
aspectvx<- velox(terrain(bathy,opt="aspect",neighbors = 8))
slopevx<- velox(terrain(bathy,opt="slope",neighbors = 8))
TPIvx<- velox(terrain(bathy,opt="TPI",neighbors = 8))


# Download NOAA ETOPO1 data
# res 15 = 15 minutes = .25 degrees
tracks$Depth<-NA
tracks$Aspect<-NA
tracks$Slope<-NA
tracks$BPI<-NA


for(i in unique(tracks$doy)){
  print(paste("DOY:",i))


  tracks_sp <- SpatialPointsDataFrame(coords=cbind(tracks$lon360[tracks$doy==i&is.na(tracks$Depth)],
                                                   tracks$lat[tracks$doy==i&is.na(tracks$Depth)]),
                                      data=data.frame(tracks$loggerID[tracks$doy==i&is.na(tracks$Depth)]),
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_wrap=180"))

  pt_sp<-spTransform(tracks_sp,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                   +ellps=WGS84 +units=m +no_defs"))

  spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)

  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
  spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))

  # Attach Depth ------------------------------------------------------------
  tracks$Depth[tracks$doy==i] <- bathyvx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  tracks$Aspect[tracks$doy==i] <- aspectvx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  tracks$Slope[tracks$doy==i] <- slopevx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  tracks$BPI[tracks$doy==i] <- TPIvx$extract(spdf, fun=function(x)mean(x,na.rm=T))
}

saveRDS(tracks,OutDataPath)
#
ggplot()+
  geom_polygon(data=map_data(map = 'world2'),aes(long,lat,group=group),fill=8,color="black")+
  geom_point(data=tracks %>%  filter(!is.na(Depth)) ,aes(x=lon360,y=lat,col=Depth))+
  coord_fixed(xlim = c(140,220), ylim=c(30,70))+
  theme(legend.position = "bottom") +
  facet_wrap(~yearT)

ggplot()+
  geom_polygon(data=map_data(map = 'world2'),aes(long,lat,group=group),fill=8,color="black")+
  geom_point(data=tracks %>%  filter(!is.na(sst)) ,aes(x=lon360,y=lat,col=sst))+
  coord_fixed(xlim = c(140,220), ylim=c(30,70))+
  theme(legend.position = "bottom") +
  facet_wrap(~yearT)+scale_colour_viridis_c(option = "B")
ggplot()+
  geom_polygon(data=map_data(map = 'world2'),aes(long,lat,group=group),fill=8,color="black")+
  geom_point(data=tracks %>%  filter(!is.na(dist2ice)) ,aes(x=lon360,y=lat,col=dist2ice))+
  coord_fixed(xlim = c(140,220), ylim=c(30,70))+
  theme(legend.position = "bottom") +
  facet_wrap(~yearT)+scale_colour_viridis_c(option = "B")
ggplot()+
  geom_polygon(data=map_data(map = 'world2'),aes(long,lat,group=group),fill=8,color="black")+
  geom_point(data=tracks %>%  filter(!is.na(ssh)) ,aes(x=lon360,y=lat,col=ssh))+
  coord_fixed(xlim = c(140,220), ylim=c(30,70))+
  theme(legend.position = "bottom") +
  facet_wrap(~yearT)+scale_colour_viridis_c(option = "B")
# # Modis CHL monthly means data --------------------------------------------
# # monthly 30 day
# #  # code to resample to 0.25 degrees
# modisCHL<-read_rds(paste0(BS_ENV,"/NCAR_Reanalysis1/erdMBchlamday10-17_raster.rds"))
# modisCHL<-agg(stack_s = modisCHL,fun = "mean",fact = 10)
#
# # Make date table
# datesCHL<-data.frame(xdate=names(modisCHL),date=ymd(gsub("X","",names(modisCHL)))) %>%
#   mutate(year=year(date),month=month(date)) %>%
#   dplyr::select(-date)
#
# dates<-data.frame(date=sort(unique(tracks$date))) %>%
#   mutate(year=year(date),month=month(date)) %>%
#   left_join(datesCHL)
# head(dates)
#
# # Loop to add modis chl
# tracks$chl<-NA
# tracks$chl_g<-NA
# dim(modisCHL)
#
# for(i in 1:length(dates$date)){
#   print(dates[i,])
#   if(is.na(dates$xdate[i])) next
#   pt<-data.frame(lon360=tracks$lon360[which(tracks$date==dates$date[i])],
#                  lat=tracks$lat[which(tracks$date==dates$date[i])])
#
#   coordinates(pt)<-cbind(tracks$lon360[which(tracks$date==dates$date[i])],
#                          tracks$lat[which(tracks$date==dates$date[i])])
#   crs(pt)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180")
#
#   pt_sp<-spTransform(pt,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
#                              +ellps=WGS84 +units=m +no_defs"))
#
#   spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)
#
#   spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
#   spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))
#
#   # chl
#   vx<- velox(modisCHL[[dates$xdate[i]]])
#   tracks$chl[tracks$date==dates$date[i]] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))
#
#   vx<- velox(terrain(modisCHL[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$chl_g[tracks$date==dates$date[i]] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))
#
# }
#
# saveRDS(tracks,OutDataPath)
#
# # Chem/Bio Vars from dataset-global-analysis-forecast-bio model ---------------
# # currently I only have the 2012 - 2017 data 7day res 0.5 degrees
# # The first date in the dataset is 2011-12-30
# #
# # nc_open(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc")
# # dataset-global-analysis-forecast-bio-001-014_1509132488705.nc 12-14
# # dataset-global-analysis-forecast-bio-001-014_1494903464030.nc 15-17
#
# # 2012-14
# Fe12<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="Fe")
# O212<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="O2")
# PO412<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="PO4")
# PHYC12<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="PHYC")
# CHL12<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="CHL")
# NO312<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="NO3")
# PP12<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="PP")
# SI12<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1509132488705.nc") ,varname="Si")
#
# # 2015-17
# Fe15<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="Fe")
# O215<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="O2")
# PO415<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="PO4")
# PHYC15<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="PHYC")
# CHL15<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="CHL")
# NO315<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="NO3")
# PP15<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="PP")
# SI15<-stack(paste0(BS_ENV,"/global-analysis-forecast-bio/dataset-global-analysis-forecast-bio-001-014_1494903464030.nc") ,varname="Si")
#
# # STACK THE VARS
# Fe<-stack(Fe12,Fe15)
# O2<-stack(O212,O215)
# PO4<-stack(PO412,PO415)
# PHYC<-stack(PHYC12,PHYC15)
# CHL<-stack(CHL12,CHL15)
# NO3<-stack(NO312,NO315)
# PP<-stack(PP12,PP15)
# SI<-stack(SI12,SI15)
#
# # Date Table
# datesMod<-data.frame(xdate=names(O2),date=as.Date("1950-01-01 00:00:00")+hours(gsub("X","",names(O2)))) %>%
#   mutate(year=year(date),
#          month=month(date),
#          week=week(date)) %>%
#   dplyr::select(-date)
#
# dates<-data.frame(date=sort(unique(tracks$date))) %>%
#   mutate(year=year(date),month=month(date),week=week(date)) %>%
#   dplyr::select(-date) %>%
#   left_join(datesMod) %>% distinct() %>%
#   mutate(xdate=if_else(is.na(xdate),lead(xdate),xdate)) %>% filter(year>2011)
#
# # Loop to add bio vars
# tracks$Fe<-NA
# tracks$Fe_g<-NA
# tracks$O2<-NA
# tracks$O2_g<-NA
# tracks$PO4<-NA
# tracks$PO4_g<-NA
# tracks$PHYC<-NA
# tracks$PHYC_g<-NA
# tracks$CHL<-NA
# tracks$CHL_g<-NA
# tracks$NO3<-NA
# tracks$NO3_g<-NA
# tracks$PP<-NA
# tracks$PP_g<-NA
# tracks$SI<-NA
# tracks$SI_g<-NA
#
# for(i in 1:length(dates$xdate)){
#   print(dates[i,])
#   if(is.na(dates$xdate[i])) next
#   weeks<-tracks$year==dates$year[i]&
#     tracks$month==dates$month[i]&
#     tracks$week==dates$week[i]
#   pt<-data.frame(lon360=tracks$lon360[weeks],
#                  lat=tracks$lat[weeks])
#
#   coordinates(pt)<-cbind(tracks$lon360[weeks],
#                          tracks$lat[weeks])
#   crs(pt)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180")
#
#   pt_sp<-spTransform(pt,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
#                              +ellps=WGS84 +units=m +no_defs"))
#
#   spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)
#
#   spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
#   spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))
#
#   # Fe
#   vx<- velox(Fe[[dates$xdate[i]]])
#   tracks$Fe[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(Fe[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$Fe_g[weeks] <- vx$extract(spdf, fun=mean)
#
#   #O2
#   vx<- velox(O2[[dates$xdate[i]]])
#   tracks$O2[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(O2[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$O2_g[weeks] <- vx$extract(spdf, fun=mean)
#
#   #PO4
#   vx<- velox(PO4[[dates$xdate[i]]])
#   tracks$PO4[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(PO4[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$PO4_g[weeks] <- vx$extract(spdf, fun=mean)
#
#   #PHYC
#   vx<- velox(PHYC[[dates$xdate[i]]])
#   tracks$PHYC[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(PHYC[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$PHYC_g[weeks] <- vx$extract(spdf, fun=mean)
#
#
#   #CHL
#   vx<- velox(CHL[[dates$xdate[i]]])
#   tracks$CHL[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(CHL[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$CHL_g[weeks] <- vx$extract(spdf, fun=mean)
#
#   #NO3
#   vx<- velox(NO3[[dates$xdate[i]]])
#   tracks$NO3[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(NO3[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$NO3_g[weeks] <- vx$extract(spdf, fun=mean)
#
#   #PP
#   vx<- velox(PP[[dates$xdate[i]]])
#   tracks$PP[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(PP[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$PP_g[weeks] <- vx$extract(spdf, fun=mean)
#
#   #SI
#   vx<- velox(SI[[dates$xdate[i]]])
#   tracks$SI[weeks] <- vx$extract(spdf, fun=mean)
#
#   vx<- velox(terrain(SI[[dates$xdate[i]]],opt = "slope",neighbors = 8))
#   tracks$SI_g[weeks] <- vx$extract(spdf, fun=mean)
#
# }

# head(as.data.frame(tracks))
# table(is.na(tracks$O2))

# saveRDS(tracks,OutDataPath)


# Distance to Shore and Shelf in m -------------------------------------------

# Download NOAA ETOPO1 data
# Highest resolution or the 0.25 degree res?
aleu <- marmap::getNOAA.bathy(130, -140, 30, 75, resolution = 1,
                              antimeridian = T,keep = T) #resolution 1.85km*27=49.9km2

# aleu <-raster('marmap_coord_-140;30;130;75_res_1_anti.csv')
# Make it a raster
bathy<-as.raster(aleu)
bathy[bathy>=0]<-NA
plot(aleu)
# get the 200m isobath
Shelf<-rasterToContour(bathy,level=-200)

# reproject 200m isobath in to LAEA
Shelf<-spTransform(Shelf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))


# Make a spatial Points DF with the tracks
tracks_sp <- SpatialPointsDataFrame(coords=cbind(tracks$lon,tracks$lat),
                                    data=as.data.frame(tracks),
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

tracks_sp_laea<-spTransform(tracks_sp,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                          +ellps=WGS84 +units=m +no_defs"))

# use rgeos::gDistance() to calculate the distance from every point to every edge
Shelfdist<-  gDistance(tracks_sp_laea,Shelf,byid = T)

# take the minumum distance from each point to the edge
tracks$dist2shelf<-apply(Shelfdist, 2, function(x) min(x, na.rm = TRUE))

# Distance to Shore in m-------------------------------------------------------
#
# Extract the shoreline
Shore<-rasterToContour(bathy,level=-1)

# reproject
Shore<-spTransform(Shore,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))

# calculate distances
Shoredist<-  gDistance(tracks_sp_laea,Shore,byid = T)

# find minimum distance
tracks$dist2shore<-apply(Shoredist, 2, function(x) min(x, na.rm = TRUE))

saveRDS(tracks,OutDataPath)
#
tracks<-readRDS(OutDataPath)

# # Inside bering sea? ------------------------------------------------------
# # use Bering_sea.shp to calculate if a point is in the Bering sea
# # http://www.marineregions.org/gazetteer.php?p=details&id=4310
#
# tracksTemp<-data.frame(tracks)
# coordinates(tracksTemp)<-cbind(tracksTemp$lon,tracksTemp$lat)
# proj4string(tracksTemp)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# ber<-read_sf(paste0(BS_ENV,"/iho"),layer="iho",crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% as("Spatial")
# proj4string(ber)
# plot(ber)
# InBer<-over(tracksTemp,ber)
#
# tracks$bs<-as.character(InBer$name)
# tracks$bs[is.na(tracks$bs)]<-"N.Pacific"
#
#
# # Longhurst Provinces -----------------------------------------------------
#
# longh<-read_sf(paste0(BS_ENV,"/longhurst_v4_2010"),layer="Longhurst_world_v4_2010",crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% as("Spatial")
# longhT<-over(tracksTemp,longh)
#
# tracks$longhurst<-as.character(longhT$ProvCode)
#
# saveRDS(tracks,OutDataPath)
#
# length(unique(tracks$loggerID))
# table(tracks$loggerID)
# names(tracks)
#
# # IHO Provinces -----------------------------------------------------
#
# IHO<-read_sf(paste0(BS_ENV,"/World_Seas_IHO_v2"),layer="World_Seas_IHO_v2",crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% as("Spatial")
# IHOT<-over(tracksTemp,IHO)
# table(IHOT$NAME)
# tracks$IHO<-as.character(IHOT$NAME)
#
# saveRDS(tracks,OutDataPath)
#
# length(unique(tracks$loggerID))
# table(tracks$loggerID)
# names(tracks)
#
# # MEOW Ecoregions -----------------------------------------------------
#
# MEOW<-read_sf(paste0(BS_ENV,"/MEOW-TNC"),layer="meow_ecos",crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% as("Spatial")
# plot(IHO)
# MEOWT<-over(tracksTemp,MEOW)
# table(MEOWT$Lat_Zone  )
# tracks$MEOWeco<-as.character(MEOWT$ECOREGION)
#
# ggplot()+
#   geom_point(data=tracks,aes(x=lon360,y=lat,col=MEOWeco),alpha=1,shape=16,size=1)+
#   # coord_fixed(xlim=c(-1000,3500),ylim=c(-1400,3100))+
#   # scale_color_manual(values = c("#F8DF4F","#1C366B","#C4CFD0","#1DACE8","#F24D29","#E5C4A1"))+
#   # geom_polygon(data = w2,aes(x=long,y=lat,group=group),fill="grey")+
#   labs(col="")+theme(legend.position = "top")
#
#
# saveRDS(tracks,OutDataPath)
#
# length(unique(tracks$loggerID))
# table(tracks$loggerID)
# names(tracks)
