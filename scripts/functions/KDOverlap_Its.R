# trip.id<-trjI$trip.id
# stress<-trjI$stress


# group ID should be a 2 factor variable (e.g. sex, island)
# interval is the interval overlap should be calculated for in days (e.g. 14, 30, 90)
# intervals<-c(20,25,30)
# days<-seq(10,260,10) #jdate
# days<-seq(10,10,10)
# intervals<-c(20)

#function is not cleanly independent and may rely on particular unspecified columns 
#days since Sept1
#trj$day<-day(trj$datetime)
#trj$doySept<-as.numeric((as.Date(trj$datetime))-ymd(paste(trj$yearT-1,"09","01")))

#lat and lon are column names for coordinates

OverlapBA_its<-function (trj, tripid, groupid, its=1000, h=80, gridcell=50, days, intervals){
  bands<-unique(trip.id) #unique trips
  require(sp)
  ungroupid<-unique(groupid)
  a<-length(unique(trip.id[groupid==ungroupid[1]]))
  b<-length(unique(trip.id[groupid==ungroupid[2]]))
  c<-a+b
  
  trj$groupid<-groupid
  trj$trip.id<-trj$band
  trj<- trj[!is.na(trj$lat),]
  mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trj$lon,trj$lat)),
                                           data=trj, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  # Project data into laea 
  trj.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                               +ellps=WGS84 +units=km +no_defs"))
  
  
  
  # add the laea coordinates to the dataframe
  trj.spdf.t$x_laea<-round(coordinates(trj.spdf.t)[,1],0)
  trj.spdf.t$y_laea<-round(coordinates(trj.spdf.t)[,2],0)
  
  buffer_x <- as.integer((max(trj.spdf.t$x_laea) - min(trj.spdf.t$x_laea)) * 0.5/100) * 100
  buffer_y <- as.integer((max(trj.spdf.t$y_laea) - min(trj.spdf.t$y_laea)) * 0.5/100) * 100
  buffer <- max(buffer_x, buffer_y)
  xy_sp <- SpatialPoints(data.frame(x = c((as.integer((max(trj.spdf.t$x_laea) + 100)/100) * 100 + buffer),
                                          (as.integer((min(trj.spdf.t$x_laea) - 100)/100) * 100 - buffer)),
                                    y = c((as.integer((max(trj.spdf.t$y_laea) + 100)/100) * 100 + buffer),
                                          (as.integer((min(trj.spdf.t$y_laea) - 100)/100) * 100 - buffer))))
  customGrid <- ascgen(xy_sp, cellsize = gridcell)
  
  
  
  RIfeb_50<-data.frame()
  RIfeb_95<-data.frame()
  for (k in days){#110:152
    print(k)
    for (j in intervals){
      trjint<-trj[trj$doySept>k & trj$doySept<(k+j),]
      trjint<- trjint[!is.na(trjint$lat),]
      mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjint$lon,trjint$lat)),
                                               data=trjint, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
      # Project data into laea 
      mptstress.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                                         +ellps=WGS84 +units=km +no_defs"))
      # add the laea coordinates to the dataframe
      mptstress.spdf.t$x_laea<-round(coordinates(mptstress.spdf.t)[,1],0)
      mptstress.spdf.t$y_laea<-round(coordinates(mptstress.spdf.t)[,2],0)
      
      mptstress<-data.frame(id=(mptstress.spdf.t$groupid))#pools all individuals for each year
      coordinates(mptstress)<-cbind(mptstress.spdf.t$x_laea,mptstress.spdf.t$y_laea)
      proj4string(mptstress)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"
      
      udstress <- kernelUD(mptstress, h = h,grid=customGrid)
      BA_o50<-kerneloverlaphr(udstress , method="BA", percent=50, conditional=TRUE)
      BA_o50<-BA_o50[1,2]
      BA_o95<-kerneloverlaphr(udstress , method="BA", percent=95, conditional=TRUE)
      BA_o95<-BA_o95[1,2]
      
      RandomIndices_50<-NULL
      RandomIndices_95<-NULL
      for (i in 1:its){
        Idx<-sample(1:c, c, replace = FALSE, prob = NULL)
        bandsH<-(bands[Idx<(a+1)])
        bandsL<-(bands[Idx>a])
        
        trjL<-trjint[trjint$trip.id%in%bandsL,]
        trjL$groupRan<-"low"
        #unique(trjL$trip.id)
        trjH<-trjint[trjint$trip.id%in%bandsH,]
        trjH$groupRan<-"high"
        #unique(trjH$trip.id)
        trjR<-rbind(trjL,trjH)
        
        mptstress.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(trjR$lon,trjR$lat)),
                                                 data=trjR, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
        # Project data into laea 
        mptstress.spdf.t <- spTransform(mptstress.spdf,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                                           +ellps=WGS84 +units=km +no_defs"))
        
        # add the laea coordinates to the dataframe
        mptstress.spdf.t$x_laea<-round(coordinates(mptstress.spdf.t)[,1],0)
        mptstress.spdf.t$y_laea<-round(coordinates(mptstress.spdf.t)[,2],0)
        
        mptstress<-data.frame(id=(mptstress.spdf.t$groupRan))#pools all individuals for each year
        coordinates(mptstress)<-cbind(mptstress.spdf.t$x_laea,mptstress.spdf.t$y_laea)
        proj4string(mptstress)<-"+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000 +ellps=WGS84 +units=km +no_defs"
        
        udstress1 <- kernelUD(mptstress, h = h,grid=customGrid)
        BA_50<-kerneloverlaphr(udstress1 , method="BA", percent=50, conditional=TRUE)
        BA_95<-kerneloverlaphr(udstress1 , method="BA", percent=95, conditional=TRUE)
        
        RandomIndices_50<-c(RandomIndices_50,BA_50[1,2])
        RandomIndices_95<-c(RandomIndices_95,BA_95[1,2])
      }
      
      pval_50<-length(RandomIndices_50[RandomIndices_50<BA_o50])/its
      print(pval_50)
      pval_95<-length(RandomIndices_95[RandomIndices_95<BA_o95])/its
      print(pval_95)
      
      data5<-c(k,j,pval_50,BA_o50,mean(RandomIndices_50),min(RandomIndices_50),max(RandomIndices_50),sd(RandomIndices_50))
      RIfeb_50<-rbind(RIfeb_50,data5)
      
      data9<-c(k,j,pval_95,BA_o95,mean(RandomIndices_95),min(RandomIndices_95),max(RandomIndices_95),sd(RandomIndices_95))
      RIfeb_95<-rbind(RIfeb_95,data9)
    }
  }
  data<-list(RIfeb_50,RIfeb_95)
  return(data)
}

# A<-OverlapBA_its(trj=trjI, tripid=trip.id, 
#                  groupid=stress, its=1000, 
#                  h=80, gridcell=50, 
#                  days=days,intervals=intervals)
