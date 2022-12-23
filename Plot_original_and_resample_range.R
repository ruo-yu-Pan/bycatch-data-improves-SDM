
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)


library("dplyr")
library("maps")
library("mapdata")
library("sp")
library("sf")

# import presence record data and number of occ data
Cutlassfish_bycatch_filename = list.files("./compiled_data/occ_data_in_ground/grid_bycatch/",full.names=T)
Cutlassfish_target_filename = list.files("./compiled_data/occ_data_in_ground/grid_target/",full.names=T)

grid_in_ground = read.csv("./compiled_data/map/grid_in_ground.csv")

jpeg("./output/fig/fishing_ground.jpeg",width = 1400,height = 1600,units = "px",res = 300)
plot(c(119,126),xlim=c(119,126),ylim=c(23,30),type="n",
     xlab = "Longitude",ylab="Latitude",main="")
maps::map('worldHires',xlim=c(118,126),ylim=c(23,30),fill=F,add=T)

points(grid_data_tot_enoughdata$Var2,grid_data_tot_enoughdata$Var1,pch=19,cex=0.3,lwd=0.3)
polygon(grid_data_tot_enoughdata[hull,4],
        grid_data_tot_enoughdata[hull,3],col=rgb(0,0,0,alpha=0.3))

dev.off()




jpeg("./output/fig/Two_gear_polygon.jpg",width = 800,height = 2200,res=300,units = "px")
par(mfrow=c(3,1),mar=c(1,1,3,1),oma=c(4,3,0,0))
for(m in 9:11){
  Cutlassfish_bycatch_in_polygon = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
  Cutlassfish_target_in_polygon = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)
  
  # calculate sampling range
  unique_lon_lat_b = unique(as.matrix(Cutlassfish_bycatch_in_polygon[which(Cutlassfish_bycatch_in_polygon$fishing_loc==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat_b)
  hull_b2 <- c(hull_b, hull_b[1])
  
  unique_lon_lat_t = unique(as.matrix(Cutlassfish_target_in_polygon[which(Cutlassfish_target_in_polygon$fishing_loc==1),c("Lon","Lat")]))
  # if(m==9){unique_lon_lat_t_rm = unique_lon_lat_t[-which(unique_lon_lat_t[,2]>26.5|
  #                                                          unique_lon_lat_t[,1]<121.2),]
  # }else if(m==10){unique_lon_lat_t_rm = unique_lon_lat_t[-which(unique_lon_lat_t[,1]>122.5),]
  # }else{unique_lon_lat_t_rm = unique_lon_lat_t}
  # 
  hull_t <- chull(unique_lon_lat_t)
  hull_t2 <- c(hull_t, hull_t[1])
  
  # hull_t_rm <- chull(unique_lon_lat_t_rm)
  # hull_t_rm2 <- c(hull_t_rm, hull_t_rm[1])
  
  # overlapp
  grid_overlap = grid_in_ground[
    which(point.in.polygon(grid_in_ground$Lon,
                           grid_in_ground$Lat,
                           unique_lon_lat_b[hull_b2,1],
                           unique_lon_lat_b[hull_b2,2])!=0&
            point.in.polygon(grid_in_ground$Lon,
                             grid_in_ground$Lat,
                             unique_lon_lat_t[hull_t2,1],
                             unique_lon_lat_t[hull_t2,2])!=0),
  ]
  hull_ovlap = chull(grid_overlap[,c("Lat","Lon")])
  hull_ovlap2 = c(hull_ovlap,hull_ovlap[1])
  
  plot(c(119,126),xlim=c(119,126),ylim=c(23,30),type="n",
       xlab = "",ylab="")
  title(paste0("(",letters[m-8],") ",month.name[m]), line = 0.5, adj=0)
  maps::map('worldHires',xlim=c(118,126),ylim=c(24,30),fill=F,add=T)
  
  polygon(grid_overlap$Lon[hull_ovlap2],
          grid_overlap$Lat[hull_ovlap2],
          col=rgb(0,0,0,alpha=0.3),
          border=NA)
  polygon(unique_lon_lat_b[hull_b2,1],
          unique_lon_lat_b[hull_b2,2],
          col=rgb(0,0,0,alpha=0),
          border="black")+
  polygon(unique_lon_lat_t[hull_t2,1],
          unique_lon_lat_t[hull_t2,2],
          col=rgb(0.5,0,0,alpha=0),
          border="red")
  

  
  unique_lon_lat_b = as.data.frame(unique_lon_lat_b)
  unique_lon_lat_t = as.data.frame(unique_lon_lat_t)
  
  # if(m<11){
  #   range_buffered = st_read(
  #     paste0("./compiled_data/constrained_sampling_range/Sb1_",
  #            formatC(m,width = 2,flag = 0),"/Sb1_",
  #            formatC(m,width = 2,flag = 0),".shp"))
  #   
  #   #jpeg(paste0("./compiled_data/constrained_sampling_range/Sb1_",
  #   #            formatC(m,width = 2,flag = 0),".jpeg"),width = 4,height = 6,res = 300,units = "in")
  #   p = ggplot(data=world) +
  #     geom_sf(fill=NA) +
  #     
  #     #geom_sf(data = range_sf_polygon_proj, alpha = 0.5, colour = "black") +
  #     geom_sf(data = range_buffered, fill = NA,colour="blue")+
  #     coord_sf(xlim = c(119, 126), ylim = c(23, 30), expand = FALSE)+
  #     geom_polygon(data=unique_lon_lat_b[hull_b2,],mapping = aes(x=Lon,y=Lat),
  #                  fill=NA,colour="black")+
  #     ggtitle(paste0(month.name[m],"_Sb1"))+
  #     theme_bw()
  #   #dev.off()
  #   ggsave(paste0("./compiled_data/constrained_sampling_range/Sb1_",
  #                               formatC(m,width = 2,flag = 0),".jpeg"))
  # }
  
}

mtext("Longitude",side=1,outer = T,line=2)
mtext("Latitude",side=2,outer=T,line=2)

dev.off()


# Plot fishing ground and Two gear sampling range together
# (also used the 1_data_mining.R)

jpeg("./output/fig/fishing_ground_AND_Twogear.jpeg",width = 1600,height = 1800,units = "px",res = 300)
par(mfrow=c(2,2),mar=c(2,2,2,1),oma=c(4,3,0,0))

plot(c(119,126),xlim=c(119,126),ylim=c(23,30),type="n",
     xlab = "Longitude",ylab="Latitude",main="")

maps::map('worldHires',xlim=c(118,126),ylim=c(23,30),fill=F,add=T)

points(grid_data_tot_enoughdata$Var2,grid_data_tot_enoughdata$Var1,pch=19,cex=0.3,lwd=0.3)
polygon(grid_data_tot_enoughdata[hull,4],
        grid_data_tot_enoughdata[hull,3],col=rgb(0,0,0,alpha=0.3))


for(m in 9:11){
  Cutlassfish_bycatch_in_polygon = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
  Cutlassfish_target_in_polygon = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)
  
  # calculate sampling range
  unique_lon_lat_b = unique(as.matrix(Cutlassfish_bycatch_in_polygon[which(Cutlassfish_bycatch_in_polygon$occ==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat_b)
  hull_b2 <- c(hull_b, hull_b[1])
  
  unique_lon_lat_t = unique(as.matrix(Cutlassfish_target_in_polygon[which(Cutlassfish_target_in_polygon$occ==1),c("Lon","Lat")]))
  # if(m==9){unique_lon_lat_t_rm = unique_lon_lat_t[-which(unique_lon_lat_t[,2]>26.5|
  #                                                          unique_lon_lat_t[,1]<121.2),]
  # }else if(m==10){unique_lon_lat_t_rm = unique_lon_lat_t[-which(unique_lon_lat_t[,1]>122.5),]
  # }else{unique_lon_lat_t_rm = unique_lon_lat_t}
  # 
  hull_t <- chull(unique_lon_lat_t)
  hull_t2 <- c(hull_t, hull_t[1])
  
  # hull_t_rm <- chull(unique_lon_lat_t_rm)
  # hull_t_rm2 <- c(hull_t_rm, hull_t_rm[1])
  
  # overlapp
  grid_overlap = grid_in_ground[
    which(point.in.polygon(grid_in_ground$Lon,
                           grid_in_ground$Lat,
                           unique_lon_lat_b[hull_b2,1],
                           unique_lon_lat_b[hull_b2,2])!=0&
            point.in.polygon(grid_in_ground$Lon,
                             grid_in_ground$Lat,
                             unique_lon_lat_t[hull_t2,1],
                             unique_lon_lat_t[hull_t2,2])!=0),
  ]
  hull_ovlap = chull(grid_overlap[,c("Lat","Lon")])
  hull_ovlap2 = c(hull_ovlap,hull_ovlap[1])
  
  plot(c(119,126),xlim=c(119,126),ylim=c(23,30),type="n",
       xlab = "",ylab="")
  #title(month.name[m], line = 0.5, adj=0)
  maps::map('worldHires',xlim=c(118,126),ylim=c(24,30),fill=F,add=T)
  
  polygon(grid_overlap$Lon[hull_ovlap2],
          grid_overlap$Lat[hull_ovlap2],
          col=rgb(0,0,0,alpha=0.3),
          border=NA)
  polygon(unique_lon_lat_b[hull_b2,1],
          unique_lon_lat_b[hull_b2,2],
          col=rgb(0,0,0,alpha=0),
          border="black")+
    polygon(unique_lon_lat_t[hull_t2,1],
            unique_lon_lat_t[hull_t2,2],
            col=rgb(0.5,0,0,alpha=0),
            border="red")
  
  unique_lon_lat_b = as.data.frame(unique_lon_lat_b)
  unique_lon_lat_t = as.data.frame(unique_lon_lat_t)
  
}

mtext("Longitude",side=1,outer = T,line=2)
mtext("Latitude",side=2,outer=T,line=2)
dev.off()
















# --------------------------------------------------------------
library(sp)
library(sf)
library(pracma)
library(ggplot2)

library("rnaturalearth")
library("rnaturalearthdata")
world = ne_countries(scale = "medium", returnclass = "sf")

m=9
#m=10
Cutlassfish_bycatch_occ = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)

# calculate sampling range
unique_lon_lat = unique(as.matrix(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$occ==1),c("Lon","Lat")]))
hull_b <- chull(unique_lon_lat)
hull_b2 <- c(hull_b, hull_b[1])
sampling_range = 
  grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                        grid_in_ground$Lat,
                                        unique_lon_lat[hull_b2,1],
                                        unique_lon_lat[hull_b2,2])!=0),]
# the shrinkage of the polygon --------------------------------
# - transform to sf polygon 
sampling_polygon = as.data.frame(unique_lon_lat[hull_b,])
colnames(sampling_polygon) = c("Lon","Lat")
range_sf = st_as_sf(sampling_polygon, coords = c("Lon", "Lat"), crs = 4326)
range_sf_polygon = sampling_polygon %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# - transform to meter unit
lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {utm + 32600
  }else{utm + 32700 }
}

coords_example = c(123,26.5)
EPSG_2_UTM <- lonlat2UTM(coords_example) # 32651, find the coordinate reference

range_sf_polygon_proj = st_transform(range_sf_polygon, EPSG_2_UTM) # transform

# - buffer size
range_sf_proj = st_transform(range_sf, EPSG_2_UTM) # transform

lon_dist = st_bbox(range_sf_proj)[3]-st_bbox(range_sf_proj)[1]
lat_dist = st_bbox(range_sf_proj)[4]-st_bbox(range_sf_proj)[2]
#poly_area = st_area(range_sf_polygon_proj)
buffer_level = 0.12
buffer_size = -0.5*(lon_dist+lat_dist)*buffer_level

buffer_level2 = 0.24
buffer_size2 = -0.5*(lon_dist+lat_dist)*buffer_level2

# - compute buffered polygons
range_buffered <- st_buffer(range_sf_polygon_proj, dist=buffer_size) 

range_buffered2 <- st_buffer(range_sf_polygon_proj, dist=buffer_size2) 


Cutlassfish_bycatch_occ = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
Cutlassfish_target_occ = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)

# interdiscipline ---------------------------------------------
# calculate target data sampling range 
Cutlassfish_target_occ = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)

# calculate target data sampling range
unique_lon_lat_t = unique(as.matrix(Cutlassfish_target_occ[which(Cutlassfish_target_occ$fishing_loc==1),c("Lon","Lat")]))
hull_t <- chull(unique_lon_lat_t)
hull_t2 <- c(hull_t, hull_t[1])
sampling_range_target = 
  grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                        grid_in_ground$Lat,
                                        unique_lon_lat_t[hull_t2,1],
                                        unique_lon_lat_t[hull_t2,2])!=0),]
sampling_polygon_t = as.data.frame(unique_lon_lat_t[hull_t,])
colnames(sampling_polygon_t) = c("Lon","Lat")
range_sf_t = st_as_sf(sampling_polygon_t, coords = c("Lon", "Lat"), crs = 4326)
range_sf_t_polygon = sampling_polygon_t %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# - transform to meter unit
lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {utm + 32600
  }else{utm + 32700 }
}

coords_example = c(123,26.5)
EPSG_2_UTM <- lonlat2UTM(coords_example) # 32651, find the coordinate reference

range_sf_t_polygon_proj = st_transform(range_sf_t_polygon, EPSG_2_UTM) # transform


# Plot ------------------------------------------------------



png(paste0("./compiled_data/constrained_sampling_range/ALL_",
           formatC(m,width = 2,flag = 0),".png"),width = 4,height = 6,res = 300,units = "in")
ggplot(data=world) +
  geom_sf(fill=NA) +
  geom_sf(data = range_sf_polygon_proj, alpha = 0.5, colour = "black") +
  geom_sf(data = range_buffered, fill = NA,colour="blue")+
  geom_sf(data = range_buffered2, fill = NA,colour="green")+
  geom_sf(data = range_sf_t_polygon_proj, fill = NA, colour = "red") +
  coord_sf(xlim = c(119, 126), ylim = c(23, 30), expand = FALSE)+
  ggtitle(paste0(month.name[m],"_Sb1"))+
  theme_bw()
dev.off()




