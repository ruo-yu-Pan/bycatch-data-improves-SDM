
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library(ncdf4)
library(plyr)
library(reshape2)


# Temperature,Salinity,Sea floor temperature, Sea surface height, Mixed layer thickness,
#Northward velocity, Eastward velocity in 0.5m deep

EVshortname = c("thetao","so","bottomT","zos","mlotst","vo","uo")
EVname = c("SST","SSS","SBT","SSH","MLD","vo","uo")

#-------------------------------------------

nc = nc_open("./raw_data/Environmental Data/Download/global-reanalysis-phy-001-030-monthly.nc")
print(nc)

YM = data.frame()
M = data.frame()

for (i in 1:7) {
  
  #Get a variable
  #Get the the variable and its attributes, and verify the size of the array.
  #Get coordinate (including time) variables
  lat = ncvar_get( nc = nc, varid = 'latitude')
  lon = ncvar_get( nc = nc, varid = 'longitude')
  time = ncvar_get( nc = nc, varid = 'time')
  
  
  # get 
  dname = EVshortname[i]  
  array = ncvar_get(nc, dname)
  dlname = ncatt_get(nc, dname,"long_name")
  dunits = ncatt_get(nc, dname,"units")
  fillvalue = ncatt_get(nc, dname,"_FillValue")
  dadd_offset = ncatt_get(nc, dname,"add_offset")$value
  dscale_factor = ncatt_get(nc, dname,"scale_factor")$value
  dim(array)
  
  # convert time -- split the time units string into fields
  
  mdy = as.character(as.Date(as.POSIXlt(time*60*60, origin="1950-01-01")))
  
  #Replace netCDF fillvalues with R NAs
  array[array==fillvalue$value] = NA
  
  
  # Create a data frame
  lonlatmdy = as.matrix(expand.grid(lon,lat,mdy))
  lonlatmdy_df = data.frame (lonlatmdy)
  lonlatmdy_df$Year = substring(lonlatmdy_df$Var3, 1, 4)
  lonlatmdy_df$Month = substring(lonlatmdy_df$Var3, 6,7)
  head(lonlatmdy, 10)
  dim(lonlatmdy)
  
  
  # vector of `EV` values
  EV_vec_long = as.vector(array)
  length(EV_vec_long)
  
  # transform EV value back to original value
  EV_vec_long = EV_vec_long*dscale_factor+dadd_offset
  
  # reshape the vector into a matrix
  EV_mat = matrix(EV_vec_long)
  dim(EV_mat)
  head(na.omit(EV_mat))
  
  
  # create a dataframe with spatial and temporal information 
  EV = data.frame(cbind(lonlatmdy_df,EV_mat))
  colnames(EV)[c(1:3,6)] = c("Lon","Lat","Y-M-D","EV")
  
  
  # Computing by groups within data.frames
  EV$lon = floor(as.numeric((as.character(EV$Lon)))/0.1)*0.1
  EV$lat = floor(as.numeric((as.character(EV$Lat)))/0.1)*0.1
  EV$var = EVname[i]
  EV$EV = as.numeric(as.character(EV$EV))
  
  EV_YM = plyr::ddply(EV,c("Year","Month","lon","lat","var"),summarise, EV=mean(EV,na.rm=TRUE))
  EV_M = plyr::ddply(EV,c("Month","lon","lat","var"),summarise, EV=mean(EV,na.rm=TRUE))
  
  YM = rbind(YM,EV_YM)
  M = rbind(M,EV_M)

}

nc_close(nc)


grid_YM_wide = dcast(YM, Year+Month+lon+lat ~ var, value.var="EV")
grid_M_wide = dcast(M, Month+lon+lat ~ var, value.var="EV")

write.csv(grid_YM_wide,"./raw_data/Environmental data/original_data/EV(Re)_YM.csv", row.names = FALSE)
write.csv(grid_M_wide,"./raw_data/Environmental data/original_data/EV(Re)_M.csv", row.names = FALSE)

## calculate EKE
EV = read.csv("./raw_data/Environmental data/original_data/EV(Re)_YM.csv")
EV$EKE = sqrt((EV$uo)^2+(EV$vo)^2)
write.csv(EV,"./raw_data/Environmental data/original_data/EV(Re)_YM.csv", row.names = FALSE)



# Plot ----------------------------------------------------------------
library("ggplot2")
library("maps")
library("colorRamps")

# R color palettes
n <- 16
mp <- map_data("world")

ReEV_YM <- read.csv("./raw_data/Environmental data/original_data/EV(Re)_YM.csv")

EV_fullname <- c("Year","Month","lon","lat","Density ocean mixed layer thickness",
                 "Sea floor potential temperature","Sea surface height",
                 "Salinity(0.5m)","Temperature(0.5m)",
                 "Eastward velocity", "Northward velocity","Eddy kinetic energy")


for (i in 5:12) {
  reEV_YM <- ReEV_YM[,c(1,2,3,4,i)]
  colnames(reEV_YM) <- c("Year","Month","lon","lat", "Value")
  
  plot <-  ggplot(data=reEV_YM , aes(x=lon, y=lat, fill=Value)) + 
    geom_raster(interpolate=TRUE) +
    geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
    scale_fill_gradientn(colours = blue2green2red(n) ,na.value = NA) +
    theme_bw() + facet_wrap(Month~Year, ncol=11) +  
    coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))+
    ggtitle(EV_fullname[i])

  ggsave(plot, file = paste0("./raw_data/Environmental data/original_data/",colnames(ReEV_YM)[i],".png") ,width =28, height =40, dpi=96, units = "in", device='png',limitsize = FALSE) 
}
