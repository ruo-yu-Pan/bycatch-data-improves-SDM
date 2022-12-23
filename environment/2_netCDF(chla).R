wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data/"
setwd(wd)

library(ncdf4)
library(plyr)
library(reshape2)
#source("D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_code/environment/function_DMwR_utils.R")


# read data and merge
chl_YMD <- data.frame()
Chla_monthly_filename = list.files("./raw_data/Environmental Data/Download/Chla_Monthly_4km/",full.names=T)

for(i in 1:length(Chla_monthly_filename)) {
  
  nc <- nc_open(Chla_monthly_filename[i])
  
  #Get coordinate (including time) variables
  lat <- ncvar_get( nc = nc, varid = 'lat')
  lon <- ncvar_get( nc = nc, varid = 'lon')
  
  # get 
  dname <- "chlor_a"  
  chla_array <- ncvar_get(nc, dname)
  fillvalue <- ncatt_get(nc, dname,"_FillValue")
  
  #Replace netCDF fillvalues with R NAs
  chla_array[chla_array==fillvalue$value] <- NA
  
  # Get data within study area
  lat_in_ground = lat[which(lat >= 23 & lat <= 30)]
  lon_in_ground = lon[which(lon>=119 & lon <= 126)]
  chla_array_in_ground = chla_array[which(lon>=119 & lon <= 126),which(lat >= 23 & lat <= 30)]
  
  chla_array_in_ground = cbind(lat_in_ground,chla_array_in_ground)
  colnames(chla_array_in_ground) = c("lat",lon_in_ground)
  chla_array_in_ground = as.data.frame(chla_array_in_ground)
  
  #get time
  start_time <- ncatt_get(nc, 0, "time_coverage_start")$value
  as.character(as.Date(as.POSIXlt(start_time)))
  
  nc_close(nc)
  
  
  # create a dataframe 
  
  chla <- melt(chla_array_in_ground, id.vars=c("lat"),
               variable.name = "lon", 
               value.name = "chla")
  chla$lon = as.numeric(as.character(chla$lon))
  chla$`Y-M-D` <- start_time
  chla$Year <- substring(chla$`Y-M-D`, 1, 4)
  chla$Month <- substring(chla$`Y-M-D`, 6,7)
  chla$Day <- substring(chla$`Y-M-D`, 9,10)
  
  # Computing by groups within data.frames
  chla$Lon <- floor(as.numeric((as.character(chla$lon)))/0.1)*0.1
  chla$Lat <- floor(as.numeric((as.character(chla$lat)))/0.1)*0.1
  chla$chla <- as.numeric(as.character(chla$chla))
  
  chl_YMD <- rbind(chl_YMD,chla)
  
  print(i)
}

chla_YM <- ddply(chl_YMD,c("Year","Month","Day","Lon","Lat"),summarise,chla=mean(chla,na.rm=TRUE))
chla_M <- ddply(chl_YMD,c("Month","Day","Lon","Lat"),summarise,chla=mean(chla,na.rm=TRUE))

chla_YM$chla[which(is.nan(chla_YM$chla))]=NA
chla_M$chla[which(is.nan(chla_M$chla))]=NA

#--------------------------------------------------------------------------------------------------------
write.table(chla_YM , file="./raw_data/Environmental Data/original_data/chla_YM.csv",sep=",",na = "NA",
              append=FALSE,col.names=TRUE,row.names=FALSE )
  
write.table(chla_M   ,file="./raw_data/Environmental Data/original_data/chla_M.csv",sep=",",na = "NA",
              append=FALSE,col.names=TRUE,row.names=FALSE )

#--------------------------------------------------------------------------------------------------------
  
library("ggplot2")
library("maps")
library("colorRamps")

# R color palettes
op <- par(mfrow=c(6,1),mar=c(0,0,2,0))
n <- 16


mp <- map_data("world")


plot <-  ggplot(data=chla_YM, aes(x=Lon, y=Lat, fill=chla)) + 
  geom_raster(interpolate=TRUE) +
  geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
  scale_fill_gradientn(colours = blue2green2red(n) ,na.value = NA) +
  theme_bw() + facet_wrap(Month~Year, ncol=11) +  
  coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))+
  ggtitle("Chlorophyll Concentration")

ggsave(plot, file="./raw_data/Environmental Data/original_data/Chla.png" ,width =25, height =40, dpi=96, units = "in", device='png',limitsize = FALSE) 

--------------------------------------------------------------------------------------------------------

chla <- read.csv("D://農委會計劃案//Environmental Data//EV//FillNA//EV_Chla_Knn.csv")

colnames(chla) <- c("Year","Month","lon","lat", "Value")


ggplot(data=chla, aes(x=lon, y=lat, fill=Value)) + 
  geom_raster(interpolate=TRUE) +
  geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
  scale_fill_gradientn(colours = blue2green2red(n) ,na.value = NA, limits=c(0, 10)) +
  theme_bw() + facet_wrap(Month~Year, ncol=10) +  
  coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))

ggsave(plot, file="Chla.png" ,width =25, height =40, dpi=96, units = "in", device='png',limitsize = FALSE) 

