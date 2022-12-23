wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library("viridis") 
library("dplyr")
library("maps")
library("mapdata")
library("sp")

Cutlassfish = read.csv("./raw_data/Fishery_data/TWnorthern_TorchLightNet_fishery.csv")
Cutlassfish_fix = Cutlassfish


# arrange and unify the unit ###########################################################

## Lat and Lon exchange #######
Cutlassfish_fix$Lat[which(Cutlassfish$Lat>100 & Cutlassfish$Lat<200)] = 
  Cutlassfish$Lon[which(Cutlassfish$Lat>100 & Cutlassfish$Lat<200)]
Cutlassfish_fix$Lon[which(Cutlassfish$Lat>100 & Cutlassfish$Lat<200)] = 
  Cutlassfish$Lat[which(Cutlassfish$Lat>100 & Cutlassfish$Lat<200)]

Cutlassfish_fix[which(Cutlassfish$Lat>100 & Cutlassfish$Lat<200),]

## Transform unit #######

# Cutlassfish_fix$Lat[which(Cutlassfish$Lat>500)] = 
#   Cutlassfish$Lat[which(Cutlassfish$Lat>500)]%/%100+
#   Cutlassfish$Lat[which(Cutlassfish$Lat>500)]%%100/60

# Cutlassfish_fix$Lon[which(Cutlassfish$Lat>500)] = 
#   Cutlassfish$Lon[which(Cutlassfish$Lat>500)]%/%100+
#   Cutlassfish$Lon[which(Cutlassfish$Lat>500)]%%100/60

#Cutlassfish_fix[which(Cutlassfish$Lat>200),]


Cutlassfish_fix = Cutlassfish_fix[-which(Cutlassfish_fix$Lat==0|Cutlassfish_fix$Lon==0|Cutlassfish_fix$Lon>200),]
Cutlassfish_fix = Cutlassfish_fix[-which(is.na(Cutlassfish_fix$Lat)|is.na(Cutlassfish_fix$Lon)),]


plot(c(115,130),xlim=c(115,130),
     ylim=c(20,32),type="n",xlab = "Longtitude",ylab="Latitude")

Cutlassfish_fix %>% dplyr::filter(Fishing_method_label ==1) %>%
  with(points(Lon,Lat,pch=1,col=viridis(2)[1],cex=0.8))
Cutlassfish_fix %>% dplyr::filter(Fishing_method_label ==3) %>%
  with(points(Lon,Lat,pch=1,col=viridis(2)[2],cex=0.8))

maps::map('worldHires',xlim=c(115,130),
          ylim=c(20,32),fill=F,add=T)




# transform into grid data ###########################################################

## grid information ###########################################################

grid_size_lat = 0.1
grid_size_lon = 0.1

lat_vec =  seq(23,30,grid_size_lat)
lon_vec =  seq(119,126,grid_size_lon)

lat_grid = seq(23.05,29.95,grid_size_lat)
lon_grid = seq(119.05,125.95,grid_size_lat)
grid_info = expand.grid(lat_grid,lon_grid)

lat_grid_lab = c(1:length(lat_grid))
lon_grid_lab = c(1:length(lon_grid))
grid_info_lab = expand.grid(lat_grid_lab,lon_grid_lab)
grid_info$grid_lab = paste(grid_info_lab$Var1,grid_info_lab$Var2,sep="-")


## Total data (Seine+angle) ##########################################################
Cutlassfish_total = Cutlassfish_fix[-which(Cutlassfish_fix$Lon<119|Cutlassfish_fix$Lon>126|Cutlassfish_fix$Lat>30|Cutlassfish_fix$Lat<23),]

### test-1 ############################################

min_lat =min(Cutlassfish_total$Lat)
max_lat =max(Cutlassfish_total$Lat)
min_lon =min(Cutlassfish_total$Lon)
max_lon =max(Cutlassfish_total$Lon)


lat_grp = cut(Cutlassfish_total$Lat,lat_vec,lat_grid_lab)
lon_grp = cut(Cutlassfish_total$Lon,lon_vec,lon_grid_lab)
grid_lab = paste(lat_grp,lon_grp,sep="-")
grid_data = as.data.frame(table(grid_lab))

grid_data$grid_lab = as.character(grid_data$grid_lab)
grid_data_merge = merge(x = grid_data, y = grid_info, by = "grid_lab", all.x = TRUE)
grid_data_enoughdata = grid_data_merge[which(grid_data_merge$Freq>=3),]


# calculation for polygon
hull <- chull(as.matrix(grid_data_enoughdata[,c(3,4)]))
hull <- c(hull, hull[1])

plot(c(118,126),xlim=c(118,126),ylim=c(23,30),type="n",
     xlab = "Longtitude",ylab="Latitude",main="Total")

maps::map('worldHires',xlim=c(118,126),ylim=c(23,30),fill=F,add=T)

points(grid_data_enoughdata$Var2,grid_data_enoughdata$Var1,pch=19,cex=0.3,lwd=0.3)
polygon(grid_data_enoughdata[hull,4],
        grid_data_enoughdata[hull,3],col=rgb(0,0,0,alpha=0.3))

### final ###########################################################

Cutlassfish_total = Cutlassfish_total[-which(Cutlassfish_total$Lat<24),]

lat_lab_tot = cut(Cutlassfish_total$Lat,lat_vec,lat_grid_lab)
lon_lab_tot = cut(Cutlassfish_total$Lon,lon_vec,lon_grid_lab)
grid_lab_tot = paste(lat_lab_tot,lon_lab_tot,sep="-")
grid_data_tot = as.data.frame(table(grid_lab_tot))

colnames(grid_data_tot)[1] = "grid_lab"
grid_data_tot$grid_lab = as.character(grid_data_tot$grid_lab)

grid_data_tot_merge = merge(x = grid_data_tot, y = grid_info, by = "grid_lab", all.x = TRUE)
grid_data_tot_enoughdata = grid_data_tot_merge[which(grid_data_tot_merge$Freq>=3),]

# calculation for polygon
hull <- chull(as.matrix(grid_data_tot_enoughdata[,c(3,4)]))
hull <- c(hull, hull[1])

jpeg("./output/fig/fishing_ground.jpeg",width = 1400,height = 1600,units = "px",res = 300)
plot(c(119,126),xlim=c(119,126),ylim=c(23,30),type="n",
     xlab = "Longitude",ylab="Latitude",main="")

maps::map('worldHires',xlim=c(118,126),ylim=c(23,30),fill=F,add=T)

points(grid_data_tot_enoughdata$Var2,grid_data_tot_enoughdata$Var1,pch=19,cex=0.3,lwd=0.3)
polygon(grid_data_tot_enoughdata[hull,4],
        grid_data_tot_enoughdata[hull,3],col=rgb(0,0,0,alpha=0.3))

dev.off()

## total grid in polygon ###########################################
grid_in_poly = grid_info[which(point.in.polygon(grid_info$Var1,grid_info$Var2,
                                                grid_data_tot_enoughdata$Var1[hull],grid_data_tot_enoughdata$Var2[hull])!=0),]
plot(grid_in_poly$Var2,grid_in_poly$Var1,xlim=c(118,126),ylim=c(23,30),cex=0.5)

### remove the data on land
library(maptools)
data(wrld_simpl)

crs_wgs84 = CRS(SRS_string = "EPSG:4326")
slot(wrld_simpl, "proj4string") <- crs_wgs84

points <- grid_in_poly[,c(2,1)]  # Note that I reversed OP's ordering of lat/long
colnames(points) = c("lon","lat")
pts <- SpatialPoints(points, proj4string=CRS(SRS_string='EPSG:4326'))
ii <- which(!is.na(over(pts, wrld_simpl)$FIPS)==T) # find which is on land

grid_in_poly_insea = grid_in_poly[-ii,]
grid_in_poly_insea = grid_in_poly_insea[-which(grid_in_poly_insea$grid_lab=="22-26"),] # found no env data at this place


plot(grid_in_poly_insea$Var2,grid_in_poly_insea$Var1,xlim=c(118,126),ylim=c(23,30),cex=0.5)
colnames(grid_in_poly_insea)[1:2]=c("Lat","Lon")

write.csv(grid_in_poly_insea,"./compiled_data/map/grid_in_ground.csv")

par(mfrow=c(1,1))
plot(grid_in_poly_insea$Var2,grid_in_ground$Var1,xlim=c(118,126),ylim=c(23,30),cex=0.4)
map('worldHires',xlim=c(118,126),ylim=c(23,30),fill=T,add=T,col = "white")
polygon(grid_data_tot_enoughdata[hull,4],
        grid_data_tot_enoughdata[hull,3],col=rgb(0,0,0,alpha=0.3))



## data used (data in polygon) ####################################################

Cutlassfish_total$grid_lab = grid_lab_tot
Cutlassfish_total = Cutlassfish_total%>%dplyr::filter(is.na(Fishing_method_label)==F)
Cutlassfish_total$fishing_loc = 1

Cutlassfish_tot_in_polygon = merge(x = Cutlassfish_total, y = grid_in_poly_insea, by = "grid_lab", all.x = TRUE)
Cutlassfish_tot_in_polygon = Cutlassfish_tot_in_polygon[-which(is.na(Cutlassfish_tot_in_polygon$Lat.y)),]
Cutlassfish_tot_in_polygon$Lat.x = Cutlassfish_tot_in_polygon$Lat.y
Cutlassfish_tot_in_polygon$Lon.x = Cutlassfish_tot_in_polygon$Lon.y
Cutlassfish_tot_in_polygon = Cutlassfish_tot_in_polygon[,-c(14:15)]

write.csv(Cutlassfish_tot_in_polygon,".\\compiled_data\\fishery_data_in_ground\\Cutlassfish_total.csv")



