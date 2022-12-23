wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

#Reshaping from raster to rectangular
library(ncdf4)
library(plyr)
library(reshape2)


nc = nc_open("./raw_data/Environmental Data/Download/gebco_2020.nc")
print(nc)

#-------------------------------------------
#Get a variable
#Get the the variable and its attributes, and verify the size of the array.
  
#Get coordinate (including time) variables
lat <- ncvar_get( nc = nc, varid = 'lat')
lon <- ncvar_get( nc = nc, varid = 'lon')


# get 
array <- ncvar_get(nc, "elevation")
fillvalue <- ncatt_get(nc, "elevation", "_FillValue")
dim(array)


#Replace netCDF fillvalues with R NAs
array[array==fillvalue$value] <- NA

nc_close(nc)

# Create a data frame
lonlat <- as.matrix(expand.grid(lon,lat))
lonlat_df <- data.frame (lonlat)
dim(lonlat)


# vector of `tmp` values
Bathy_vec_long <- as.vector(array)
length(Bathy_vec_long)


# reshape the vector into a matrix
Bathy_mat <- matrix(Bathy_vec_long)
dim(Bathy_mat)
head(na.omit(Bathy_mat))

# create a dataframe 
Bathy <- data.frame(cbind(lonlat_df,Bathy_mat))
names(Bathy) <- c("lon","lat","Bathy")

#--------------------------------------------------------------
#Delete Bathy>0

ind = which(Bathy$Bathy > 0)
Bathy_Ocean <- Bathy[-ind,]
ind_Ocean <- which(Bathy_Ocean$Bathy > 0)

#--------------------------------------------------------------
# Computing by groups within data.frames (grid for 0.01)
Bathy_Ocean$Lon <- floor(as.numeric((as.character(Bathy_Ocean$lon)))/0.1)*0.1
Bathy_Ocean$Lat <- floor(as.numeric((as.character(Bathy_Ocean$lat)))/0.1)*0.1
Bathy_Ocean$Bathy<- as.numeric(as.character(Bathy_Ocean$Bathy))

Bathy_Ocean01 <- ddply(Bathy_Ocean,c("Lon","Lat"), summarise, Bathy=mean(Bathy,na.rm=TRUE))

write.table(Bathy_Ocean01 , file="./raw_data/Environmental data/original_data/BathyOcean.csv",sep="," , na = "NA",
            append=FALSE, col.names=TRUE, row.names=FALSE)

----------------------------------------------------------------

library("ggplot2")
library("maps")
library("colorRamps")
library("RColorBrewer")

# R color palettes
n <- 16

mp <- map_data("world")

  

plot <-  ggplot(data=Bathy_Ocean01, aes(x=Lon, y=Lat, fill=Bathy)) + 
  geom_raster(interpolate=TRUE) +
  geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
  scale_fill_gradientn(colours = c("darkblue","blue4","dodgerblue3","aliceblue") ,na.value = NA) +
  theme_bw() +
  coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))+
  ggtitle("Bathy")
  
ggsave(plot, file="./raw_data/Environmental data/original_data/Bathy.png" ,width =25, height =40, dpi=96, units = "in", device='png',limitsize = FALSE) 

