#RAO Jan 9 2019: needs to be moditied to work with Buldir datafiles
# also needs modifications to make saving directories cross compatable

#general non-spatial packages
library(readr) #for read_csv that inturprets dates for you
library(lubridate) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

wrap360 = function(lon) {lon360<-ifelse(lon<0,lon+360,lon); return(lon360)}
w2hr<-map_data('world')
w2hr_sub<-w2hr[w2hr$region%in%c("USA","Russia","China","Japan","Canada"),]

source("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK/Functions/KDOverlap_Its.R")

# Read in RLKIs -----------------------------------------------------------
#if(Sys.info()[6]=="rachaelorben") pdir<-"/Users/rachaelorben/Dropbox/Research/BLKI_MigrationsTelomeres" ##RAO

trj<-readRDS("/Users/rachaelorben/Dropbox/Research/RLKI_Migrations/Data/processedLocs/probGLS/RLKI_GLS_allyearsalldata.rda")
trj<-trj%>%filter(trip==1)
trj<-trj[trj$band!=127307204,]
#trj$band[is.na(trj$band)==TRUE]<-173300873

trj$day<-day(trj$datetime)
trj$doySept<-as.numeric((as.Date(trj$datetime))-ymd(paste(trj$yearT-1,"09","01")))
summary(trj$doySept)
summary(trj)

A<-unique(trj$band)

N<-trj%>%group_by(band)%>%summarize(n=n())
N

#allbirds
quartz(height=6, width=6)
ggplot(trj, aes(x=lon360, y=lat, group=band))+
  coord_map("lambert",45,60,xlim = c(130,240),ylim=c(35,70))+
  geom_polygon(data=w2hr,aes(long,lat,group=group),fill="grey90",color="grey",show.legend=T) + 
  geom_path(aes(colour = factor(band)))+facet_wrap(~yearT)+
  theme(legend.position = "none")


colnames(trj)
library(adehabitatLT)
library(adehabitatHR)

# Bathymetry for plotting -------------------------------------------------
#Landmask in Adehabitat
library(marmap)
library(sp)
library(rgdal)
library(raster)

library(SDMTools)
library(maptools)
library(ggmap)
library(mapdata)
Migration<-"~/Dropbox/Research/RLKI_Migrations"
bathy2<-readRDS(paste0(Migration,"/LandOceanGrid_rasterplot_resolution5.rda"))
if(Sys.info()[7]=="rachaelorben") userdir<-'/Users/rachaelorben/Dropbox/Research/GlobalFishingWatch'
bathy2<-readRDS(paste0(userdir,"/Analysis/compileddata/Bathymetryforggplot.rda"))


# egg groups -----------------------------------------------------------
# regeime mean deviations - ALL Birds:
if(Sys.info()[["user"]]=="rachaelorben") {SGRK_Data<-"~/Dropbox/Research/RLKIinc/SGRK_DATA"}
Binfo<-readRDS(paste0(SGRK_Data, "/Physiology/DeployMatrix_Physiology_Compiled_year_deploy_p1_18Jan18.rds"))
Binfo<-Binfo%>%dplyr::select(band_999noband,island,year_deploy_p1,deployID,Sex)
data.table::setnames(Binfo, "band_999noband", "band")

trj<-left_join(trj,Binfo,by="band")
trj$day<-day(trj$datetime)

quartz(height=6, width=6)
ggplot()+
  geom_path(data=trj, aes(x=wrap360(lon), y=lat, group=band,color=Sex))+
  coord_map("lambert",45,60,xlim = c(130,210),ylim=c(35,70))+
  geom_path(data=w2hr,aes(long,lat,group=group),show.legend=T) 


# Spatial Pts Df ---------------------------------------------------
trjI<-trj[is.na(trj$Sex)==FALSE,]
summary(trjI)

mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjI$lon,trjI$lat)),
                                         data=trjI, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Project data into laea 
trj.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                             +ellps=WGS84 +units=km +no_defs"))



# add the laea coordinates to the dataframe
head(coordinates(trj.spdf.t))
trj.spdf.t$x_laea<-round(coordinates(trj.spdf.t)[,1],0)
trj.spdf.t$y_laea<-round(coordinates(trj.spdf.t)[,2],0)

mpt<-data.frame(id=(trj.spdf.t$Sex))#pools all individuals for each year
coordinates(mpt)<-cbind(trj.spdf.t$x_laea,trj.spdf.t$y_laea)
proj4string(mpt)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"

# plot the projected data (laea)
plot(trj.spdf.t)
map.axes()
# Regular grid with no land mask ---------------------------
# define the coordinates (all the detections fall inside the polygon)
# 2. Compute a grid around the fixes.
buffer_x <- as.integer((max(trj.spdf.t$x_laea) - min(trj.spdf.t$x_laea)) * 0.5/100) * 100
buffer_y <- as.integer((max(trj.spdf.t$y_laea) - min(trj.spdf.t$y_laea)) * 0.5/100) * 100
buffer <- max(buffer_x, buffer_y)
xy_sp <- SpatialPoints(data.frame(x = c((as.integer((max(trj.spdf.t$x_laea) + 100)/100) * 100 + buffer),
                                        (as.integer((min(trj.spdf.t$x_laea) - 100)/100) * 100 - buffer)),
                                  y = c((as.integer((max(trj.spdf.t$y_laea) + 100)/100) * 100 + buffer),
                                        (as.integer((min(trj.spdf.t$y_laea) - 100)/100) * 100 - buffer))))
customGrid <- ascgen(xy_sp, cellsize = 50)

# Kernals --------------------------------------------------------
library(adehabitatHR)
library(sp)
library(adehabitatMA)

ud <- kernelUD(mpt, h = 80, grid=customGrid)
BA_over50<-kerneloverlaphr(ud , method="BA", percent=50, conditional=TRUE)
BA_over50

###saveRDS(ud,paste0(Migration,"/Data/KDestimates_stress.rda"))
image(ud)

# Plot Sex Plot kernels by changing id value 1&2 -----------------------------
id=2
ud25 <- getverticeshr(ud[[id]], percent=25, standardize=TRUE)
ud50 <- getverticeshr(ud[[id]], percent=50, standardize=TRUE)
ud75 <- getverticeshr(ud[[id]], percent=75, standardize=TRUE)
ud95 <- getverticeshr(ud[[id]], percent=95, standardize=TRUE)

proj4string(ud25) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud25 <- spTransform(ud25, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud50) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50 <- spTransform(ud50, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud75) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud75 <- spTransform(ud75, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95 <- spTransform(ud95, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

p2<-ggplot()+ geom_tile(data=bathy2,aes(x=V1,y=V2,fill=Depth))+
  scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen") + 
  geom_polygon(data=fortify(ud95),aes(x=wrap360(long),y=lat,group=group),alpha=0.75,fill="green",show.legend=T)+
  geom_polygon(data=fortify(ud75),aes(x=wrap360(long),y=lat,group=group),fill="yellow",show.legend=T)+
  geom_polygon(data=fortify(ud50),aes(x=wrap360(long),y=lat,group=group),fill="orange",show.legend=T)+
  geom_polygon(data=fortify(ud25),aes(x=wrap360(long),y=lat,group=group),fill="red",show.legend=T)+
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),color="black",size=0.1, alpha=.05)+
  coord_fixed(ratio=1.7,xlim = c(142,240),ylim=c(30,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold",size=32,hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white"))
library(gridExtra)
quartz()
grid.arrange(p1,p2,nrow=2,right="Male = 2 (bottom), Female = 1(top)")
quartz.save("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/UDs/SexKD_1114151617_alltime.png",type = "png",dpi = 300)


# Plot 50%:  -----------------------------
ud50H <- getverticeshr(ud[[1]], percent=50, standardize=TRUE)
ud50L <- getverticeshr(ud[[2]], percent=50, standardize=TRUE)
ud95H <- getverticeshr(ud[[1]], percent=95, standardize=TRUE)
ud95L <- getverticeshr(ud[[2]], percent=95, standardize=TRUE)

proj4string(ud50H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50H <- spTransform(ud50H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud50L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50L <- spTransform(ud50L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95H <- spTransform(ud95H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95L <- spTransform(ud95L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

bathy2$Depth[bathy2$Depth>0]<-0
quartz()
ggplot()+ geom_tile(data=bathy2,aes(x=V1,y=V2,fill=Depth))+
  scale_fill_gradientn(colours = c("gray15","gray25", "gray50", "gray75","gray95"),name="Depth") +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_polygon(data=fortify(ud95L),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="turquoise",show.legend=T)+
  geom_polygon(data=fortify(ud95H),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="red",show.legend=T)+
  geom_polygon(data=fortify(ud50L),aes(x=wrap360(long),y=lat,group=group),alpha=0.65,fill="turquoise",col="turquoise",show.legend=T)+
  geom_polygon(data=fortify(ud50H),aes(x=wrap360(long),y=lat,group=group),alpha=0.40,fill="red",col="red",show.legend=T)+
  #annotate("text",x=143,y=64,label="B)",color="white",size=7)+
  coord_fixed(ratio=1.7,xlim = c(142,240),ylim=c(35,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold",size=32,hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white"))
quartz.save("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/UDs/SexKD_50UD_1114151617_alltime.png",type = "png",dpi = 200)

# sliding intervals - pvalues -------------------------------------------------
trj1<-trj[is.na(trj$Sex)==FALSE,]
trip.id<-trj1$band
groupid=trj1$Sex
days<-seq(0,255,5)
intervals<-c(10,20,30)

KDits_datalist<-OverlapBA_its(trj1, 
                              tripid=bands, 
                              groupid=trj1$Sex, 
                              its=1000, 
                              h=80, 
                              gridcell=50, 
                              days, 
                              intervals)
saveRDS(KDits_datalist,"/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/SexKD_50UD95UD_1114151617_RandomOverlap.rda")
KDits_datalist<-readRDS("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/SexKD_50UD95UD_1114151617_RandomOverlap.rda")


# Overlap tile plot -------------------------------------------------------

RIfeb_50a<-KDits_datalist[[1]]
RIfeb_95a<-KDits_datalist[[2]]

colnames(RIfeb_50a)<-c("d1","interval","pval","BA","meanBA","maxBA","minBA","sdBA")
colnames(RIfeb_95a)<-c("d1","interval","pval","BA","meanBA","maxBA","minBA","sdBA")

RIfeb_50a$dd<-ymd(paste("2000","09","01"))+days(RIfeb_50a$d1)
RIfeb_95a$dd<-ymd(paste("2000","09","01"))+days(RIfeb_95a$d1)

quartz(width=4,height=3)
ggplot()+
  geom_tile(data=RIfeb_50a, aes(x=dd,y=interval, fill=BA),size=1)+
  geom_point(data=RIfeb_50a[RIfeb_50a$pval<0.05,], aes(x=dd,y=interval), color="red",size=1.5)+
  scale_fill_gradient2(low="white", high="black") +
  theme_classic()+ylab("Interval")+
  xlab(NULL)+
  scale_x_date(date_breaks = "month",date_labels =  "%b")
quartz(width=4,height=3)
ggplot()+
  geom_tile(data=RIfeb_95a, aes(x=dd,y=interval, fill=BA),size=1)+
  geom_point(data=RIfeb_95a[RIfeb_95a$pval<0.05,], aes(x=dd,y=interval), color="red",size=1.5)+
  scale_fill_gradient2(low="white", high="black") +
  theme_classic()+ylab("Interval")+
  xlab(NULL)
#quartz.save(paste0(pdir,"Stress50_pvals.png"),type = "png",dpi = 600)

ggplot()+
  geom_point(data=RIfeb, aes(x=d1,interval, color=pval),size=10)+
  scale_color_viridis(discrete = F,direction = -1) 

scale_color_gradientn(colors=rev(rainbow(50,end=.7)))


# Spatial Pts Df ---------------------------------------------------
#fall egg seperation period
data=RIfeb_50a[RIfeb_50a$pval<0.05 & RIfeb_50a$interval==20,]
ggplot()+geom_point(data=data,aes(x=d1,y=pval))
#data2=RIfeb_50a[RIfeb_95a$pval<0.05 & RIfeb_95a$interval==20,]
#ggplot()+geom_point(data=data2,aes(x=d1,y=pval))
trjI_split1<-trjI[trjI$doySept>39 & trjI$doySept<(40+20),]
min(trjI_split1$datetime); max(trjI_split1$datetime)

trjI_split2<-trjI[trjI$doySept>59 & trjI$doySept<(60+20),]
min(trjI_split2$datetime); max(trjI_split2$datetime)

trjI_split3<-trjI[trjI$doySept>79 & trjI$doySept<(80+20),]
min(trjI_split3$datetime); max(trjI_split3$datetime)

trjI_splitmin<-trjI[trjI$doySept>54 & trjI$doySept<(55+10),]
min(trjI_splitmin$datetime); max(trjI_splitmin$datetime)


# 40-60 -------------------------------------------------------------------
mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjI_split1$lon,trjI_split1$lat)),
                                         data=trjI_split1, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Project data into laea 
trj.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                             +ellps=WGS84 +units=km +no_defs"))



# add the laea coordinates to the dataframe
head(coordinates(trj.spdf.t))
trj.spdf.t$x_laea<-round(coordinates(trj.spdf.t)[,1],0)
trj.spdf.t$y_laea<-round(coordinates(trj.spdf.t)[,2],0)

mptstress<-data.frame(id=(trj.spdf.t$egg))#pools all individuals for each year
coordinates(mptstress)<-cbind(trj.spdf.t$x_laea,trj.spdf.t$y_laea)
proj4string(mptstress)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"

# plot the projected data (laea)
plot(trj.spdf.t)
map.axes()

ud <- kernelUD(mptstress, h = 80, grid=customGrid)

# Plot 50%:  -----------------------------
image(ud)
ud50H <- getverticeshr(ud[[1]], percent=50, standardize=TRUE)
ud50L <- getverticeshr(ud[[2]], percent=50, standardize=TRUE)
ud95H <- getverticeshr(ud[[1]], percent=95, standardize=TRUE)
ud95L <- getverticeshr(ud[[2]], percent=95, standardize=TRUE)

proj4string(ud50H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50H <- spTransform(ud50H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud50L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50L <- spTransform(ud50L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95H <- spTransform(ud95H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95L <- spTransform(ud95L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

bathy2$Depth[bathy2$Depth>0]<-0
#quartz()
p1<-ggplot()+ geom_tile(data=bathy2,aes(x=V1,y=V2,fill=Depth))+
  scale_fill_gradientn(colours = c("gray15","gray25", "gray50", "gray75","gray95"),name="Depth") +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_polygon(data=fortify(ud95L),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud95H),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="yellow",show.legend=T)+
  geom_polygon(data=fortify(ud50L),aes(x=wrap360(long),y=lat,group=group),alpha=0.65,fill="blue",col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud50H),aes(x=wrap360(long),y=lat,group=group),alpha=0.40,fill="yellow",col="yellow",show.legend=T)+
  annotate("text",x=143,y=64,label="A)",color="white",size=5)+
  annotate("text",x=201,y=64,label="Oct 11-30",color="white",size=3,hjust="left")+
  coord_fixed(ratio=1.7,xlim = c(142,215),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold",size=32,hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white"),
        (legend.position="none"))
#quartz.save("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/UDs/EggNoEggKD_50UD_1114151617_Oct.png",type = "png",dpi = 200)


# Nov ---------------------------------------------------------------------

mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjI_split2$lon,trjI_split2$lat)),
                                         data=trjI_split2, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Project data into laea 
trj.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                             +ellps=WGS84 +units=km +no_defs"))



# add the laea coordinates to the dataframe
head(coordinates(trj.spdf.t))
trj.spdf.t$x_laea<-round(coordinates(trj.spdf.t)[,1],0)
trj.spdf.t$y_laea<-round(coordinates(trj.spdf.t)[,2],0)

mptstress<-data.frame(id=(trj.spdf.t$egg))#pools all individuals for each year
coordinates(mptstress)<-cbind(trj.spdf.t$x_laea,trj.spdf.t$y_laea)
proj4string(mptstress)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"

# plot the projected data (laea)
plot(trj.spdf.t)
map.axes()

ud <- kernelUD(mptstress, h = 80, grid=customGrid)

# Plot 50%:  -----------------------------
image(ud)
ud50H <- getverticeshr(ud[[1]], percent=50, standardize=TRUE)
ud50L <- getverticeshr(ud[[2]], percent=50, standardize=TRUE)
ud95H <- getverticeshr(ud[[1]], percent=95, standardize=TRUE)
ud95L <- getverticeshr(ud[[2]], percent=95, standardize=TRUE)

proj4string(ud50H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50H <- spTransform(ud50H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud50L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50L <- spTransform(ud50L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95H <- spTransform(ud95H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95L <- spTransform(ud95L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

bathy2$Depth[bathy2$Depth>0]<-0
quartz()
p2<-ggplot()+ geom_tile(data=bathy2,aes(x=V1,y=V2,fill=Depth))+
  scale_fill_gradientn(colours = c("gray15","gray25", "gray50", "gray75","gray95"),name="Depth") +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_polygon(data=fortify(ud95L),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud95H),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="yellow",show.legend=T)+
  geom_polygon(data=fortify(ud50L),aes(x=wrap360(long),y=lat,group=group),alpha=0.65,fill="blue",col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud50H),aes(x=wrap360(long),y=lat,group=group),alpha=0.40,fill="yellow",col="yellow",show.legend=T)+
  annotate("text",x=143,y=64,label="B)",color="white",size=5)+
  annotate("text",x=201,y=64,label="Oct 31-Nov 19",color="white",size=3,hjust="left")+
  coord_fixed(ratio=1.7,xlim = c(142,215),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold",size=32,hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white",
                                        legend.position="none"))
#quartz.save("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/UDs/EggNoEggKD_50UD_1114151617_Nov.png",type = "png",dpi = 200)


# Dec ---------------------------------------------------------------------

mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjI_split3$lon,trjI_split3$lat)),
                                         data=trjI_split3, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Project data into laea 
trj.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                             +ellps=WGS84 +units=km +no_defs"))



# add the laea coordinates to the dataframe
head(coordinates(trj.spdf.t))
trj.spdf.t$x_laea<-round(coordinates(trj.spdf.t)[,1],0)
trj.spdf.t$y_laea<-round(coordinates(trj.spdf.t)[,2],0)

mptstress<-data.frame(id=(trj.spdf.t$egg))#pools all individuals for each year
coordinates(mptstress)<-cbind(trj.spdf.t$x_laea,trj.spdf.t$y_laea)
proj4string(mptstress)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"

# plot the projected data (laea)
plot(trj.spdf.t)
map.axes()

ud <- kernelUD(mptstress, h = 80, grid=customGrid)

# Plot 50%:  -----------------------------
image(ud)
ud50H <- getverticeshr(ud[[1]], percent=50, standardize=TRUE)
ud50L <- getverticeshr(ud[[2]], percent=50, standardize=TRUE)
ud95H <- getverticeshr(ud[[1]], percent=95, standardize=TRUE)
ud95L <- getverticeshr(ud[[2]], percent=95, standardize=TRUE)

proj4string(ud50H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50H <- spTransform(ud50H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud50L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50L <- spTransform(ud50L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95H <- spTransform(ud95H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95L <- spTransform(ud95L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

bathy2$Depth[bathy2$Depth>0]<-0
#quartz()
p3<-ggplot()+ geom_tile(data=bathy2,aes(x=V1,y=V2,fill=Depth))+
  scale_fill_gradientn(colours = c("gray15","gray25", "gray50", "gray75","gray95"),name="Depth") +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_polygon(data=fortify(ud95L),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud95H),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="yellow",show.legend=T)+
  geom_polygon(data=fortify(ud50L),aes(x=wrap360(long),y=lat,group=group),alpha=0.65,fill="blue",col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud50H),aes(x=wrap360(long),y=lat,group=group),alpha=0.40,fill="yellow",col="yellow",show.legend=T)+
  annotate("text",x=143,y=64,label="C)",color="white",size=5)+
  annotate("text",x=201,y=64,label="Nov 20-Dec 9",color="white",size=3,hjust="left")+
  coord_fixed(ratio=1.7,xlim = c(142,215),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold",size=32,hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white",
                                        legend.position="none"))
quartz.save("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/UDs/EggNoEggKD_50UD_1114151617_Dec.png",type = "png",dpi = 200)

grid.arrange(p1,p2,p3,ncol=1)

# min BA ---------------------------------------------------------------------

mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjI_splitmin$lon,trjI_splitmin$lat)),
                                         data=trjI_splitmin, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Project data into laea 
trj.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                             +ellps=WGS84 +units=km +no_defs"))



# add the laea coordinates to the dataframe
head(coordinates(trj.spdf.t))
trj.spdf.t$x_laea<-round(coordinates(trj.spdf.t)[,1],0)
trj.spdf.t$y_laea<-round(coordinates(trj.spdf.t)[,2],0)

mptstress<-data.frame(id=(trj.spdf.t$egg))#pools all individuals for each year
coordinates(mptstress)<-cbind(trj.spdf.t$x_laea,trj.spdf.t$y_laea)
proj4string(mptstress)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"

# plot the projected data (laea)
plot(trj.spdf.t)
map.axes()

ud <- kernelUD(mptstress, h = 80, grid=customGrid)

# Plot 50%:  -----------------------------
image(ud)
ud50H <- getverticeshr(ud[[1]], percent=50, standardize=TRUE)
ud50L <- getverticeshr(ud[[2]], percent=50, standardize=TRUE)
ud95H <- getverticeshr(ud[[1]], percent=95, standardize=TRUE)
ud95L <- getverticeshr(ud[[2]], percent=95, standardize=TRUE)

proj4string(ud50H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50H <- spTransform(ud50H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud50L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud50L <- spTransform(ud50L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95H) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95H <- spTransform(ud95H, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

proj4string(ud95L) <- CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs") 
ud95L <- spTransform(ud95L, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 

#quartz()
ggplot()+ geom_tile(data=bathy2,aes(x=V1,y=V2,fill=Depth))+
  scale_fill_gradientn(colours = c("gray15","gray25", "gray50", "gray75","gray95"),name="Depth") +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_polygon(data=fortify(ud95L),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud95H),aes(x=wrap360(long),y=lat,group=group),alpha=0.10,col="yellow",show.legend=T)+
  geom_polygon(data=fortify(ud50L),aes(x=wrap360(long),y=lat,group=group),alpha=0.65,fill="blue",col="blue",show.legend=T)+
  geom_polygon(data=fortify(ud50H),aes(x=wrap360(long),y=lat,group=group),alpha=0.40,fill="yellow",col="yellow",show.legend=T)+
  #annotate("text",x=143,y=64,label="C)",color="white",size=5)+
  annotate("text",x=201,y=64,label="Oct 26 - Nov 4",color="white",size=3,hjust="left")+
  coord_fixed(ratio=1.7,xlim = c(142,215),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+
  theme(strip.text = element_text(face = "bold",size=32,hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white"),
        legend.position="none")
quartz.save("/Users/rachaelorben/Dropbox/Research/RLKIinc/SGRK_Results/UDs/EggNoEggKD_50UD_1114151617_Dec.png",type = "png",dpi = 200)

