# This script reads in all the pr files (output from probGLS) and compiles the
# "most probable track" and adds band number from the handling data.  It also
# adds columns that are useful across analyses
#
# Abram Fleishman
# 23 Dec 2018

#Load all the packages
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(dplyr)
library(probGLS)
library(readxl)

# Set main dir: Sys.info()[6] is the username for the computer.
# fill in the "" with your user name
if(Sys.info()[["user"]]=="abramfleishman") {
  Migration<-"~/Dropbox/RLKI/RLKI_Migrations";
  prPaths<-"~/Dropbox (CMI)/RLKI_Large_Data/Buldir RLKI Geolocator HPT";
}

if(Sys.info()[["user"]]=="rachaelorben") {
  Migration<-"~/Dropbox/Research/RLKI_Migrations"
  prPaths<-"~/Desktop/RKLI_FinalGLSTracks_2010-16/Final";
}


# Read in the pr files from probGLS
Files<-list.files(prPaths, pattern = "\\.rds",full.names = T,recursive = T)

# Compile 2016 ------------------------------------------------------------
outData<-NULL
for(i in 1:length(Files)){
  twlPath<-Files[i]

  pr<-readRDS(file = twlPath)
  Meta<-separate(data.frame(File=basename(dirname(twlPath))),File,
                 into=c("Tag","Band"),sep="_")
  year_deploy<-str_extract(basename(twlPath),"[0-9]{2}[JunSepAug1817]{5}") %>% dmy(.) %>% year(.) -1
  mpt<-data.frame(pr$`most probable track`)
  mpt$lon360<-ifelse(mpt$lon<0,mpt$lon+360,mpt$lon)
  mpt$loggerID<-Meta$Tag
  mpt$band<-Meta$Band
  mpt$year_deploy<-year_deploy
  mpt$island<-"Buldir"
  mpt$species<-"RLKI"
  outData<-bind_rows(outData,mpt)
}
head(outData)
mpt<-outData %>%
  group_by(loggerID,year_deploy) %>%
  mutate(yearT=year(dtime%m+%months(6)),
         trip=yearT-min(year)) %>%
  group_by(loggerID,trip) %>%
  mutate(tripDur=as.numeric(round(max(dtime)-min(dtime))),
         tripDay=as.numeric(floor(time_length(dtime-min(dtime),unit = "days")))+1,
         species="RLKI") %>%
  rename(datetime=dtime)
head(mpt)

saveRDS(mpt,paste0(Migration,"/Data/processedLocs/probGLS/RLKI_GLS_Buldir_allyearsalldata.rds"))

mpt<-readRDS(paste0(Migration,"/Data/processedLocs/probGLS/RLKI_GLS_Buldir_allyearsalldata.rda"))
