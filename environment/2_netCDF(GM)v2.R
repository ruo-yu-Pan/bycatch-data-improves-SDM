
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

Root_dir = "D:/Ruo_data/2020_Assis_FishTactic/"

library(ncdf4)
library(plyr)
library(reshape2)

library(wvtool) #for Sobel operator
source(paste0(Root_dir,"Fishtactics_SDM_code_check/environment/function_DMwR_utils.R"))
source(paste0(Root_dir,"Fishtactics_SDM_code_check/environment/function_convK.R"))
YMD <- data.frame()


#-------------------------------------------
  
nc = nc_open("./raw_data/Environmental Data/Download/global-reanalysis-phy-001-030-monthly.nc")
print(nc)


#Get a variable
#Get the the variable and its attributes, and verify the size of the array.
#Get coordinate (including time) variables
lat <- ncvar_get( nc = nc, varid = 'latitude')
lon <- ncvar_get( nc = nc, varid = 'longitude')
time <- ncvar_get( nc = nc, varid = 'time')

i <- 1

for(i in 1:125) {
  
  
  # get 
  array <- ncvar_get( nc = nc, varid = "thetao", start = c(1,1,1,i), count = c(85,85,1,1))
  fillvalue <- ncatt_get(nc, "thetao" ,"_FillValue")
  dim(array)
  
  
  #get time
  t <- ncvar_get( nc = nc, varid = 'time', start = c(i), count = c(1))
  
  
  #Replace netCDF fillvalues with R NAs
  array[array==fillvalue$value] <- NA
  
  
  #Gradient computation (creat function first)
  array01 <- knnImputation(array,k=10)
  tmp_conv <- edge.detect(x = array01, thresh1=1, thresh2=15, noise="gaussian", noise.s=3, method="Sobel")
  
  # Add dummy values around the matrix
  X <- NA
  tmp_conv_1 <- rbind(X, tmp_conv )
  tmp_conv_1 <- rbind(tmp_conv_1, X )
  tmp_conv_1 <- cbind(X, tmp_conv_1 )
  tmp_conv_1 <- cbind(tmp_conv_1, X )
  
  # Create a data frame
  lonlat <- as.matrix(expand.grid(lon,lat))
  lonlat_df <- data.frame (lonlat)
  head(lonlat, 10)
  dim(lonlat)
  
  
  # vector of `EV` values
  GM_vec_long <- as.vector(tmp_conv_1)
  length(GM_vec_long)
  
  
  # reshape the vector into a matrix
  GM_mat <- matrix(GM_vec_long)
  dim(GM_mat)
  head(na.omit(GM_mat))
  
  
  # create a dataframe 
  GM <- data.frame(cbind(lonlat,GM_mat))
  names(GM) <- c("lon","lat","GM")
  GM$`Y-M-D` <- as.character(as.Date(as.POSIXlt(t*60*60, origin="1950-01-01")))
  GM$Year <- substring(GM$`Y-M-D`, 1, 4)
  GM$Month <- substring(GM$`Y-M-D`, 6,7)
  GM$Day <- substring(GM$`Y-M-D`, 9,10)
  
  
  
  # Computing by groups within data.frames
  GM$Lon <- floor(as.numeric((as.character(GM$lon)))/0.1)*0.1
  GM$Lat <- floor(as.numeric((as.character(GM$lat)))/0.1)*0.1
  GM$GM <- as.numeric(as.character(GM$GM))
  
  YMD <- rbind(YMD,GM)
  
  print(i)
  
  i <- i+1
  
}

nc_close(nc)


GM_YM <- ddply(YMD,c("Year","Month","Lon","Lat"),summarise, GM=mean(GM,na.rm=TRUE))
GM_M <- ddply(YMD,c("Month","Lon","Lat"),summarise, GM=mean(GM,na.rm=TRUE))

write.csv(GM_YM,"./raw_data/Environmental data/original_data/EV(GM)_YM.csv", row.names = FALSE)
write.csv(GM_M,"./raw_data/Environmental data/original_data/EV(GM)_M.csv", row.names = FALSE)


----------------------------------------------------------------
  library("ggplot2")
library("maps")
library("colorRamps")

# R color palettes
op <- par(mfrow=c(6,1),mar=c(0,0,2,0))
n <- 16

mp <- map_data("world")



plot <-  ggplot(data=GM_YM, aes(x=Lon, y=Lat, fill=GM)) + 
  geom_raster(interpolate=TRUE) +
  geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
  scale_fill_gradientn(colours = blue2green2red(n) ,na.value = NA) +
  theme_bw() + facet_wrap(Month~Year, ncol=12) +  
  coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))+
  ggtitle("Temperature Front")

ggsave(plot, file="/home/kuo/?ୱ/Niza/2021/GMv2.png" ,width =25, height =40, dpi=96, units = "in", device='png',limitsize = FALSE) 

dev.off()