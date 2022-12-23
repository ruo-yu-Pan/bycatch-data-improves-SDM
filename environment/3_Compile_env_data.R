

wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

Root_dir = "D:/Ruo_data/2020_Assis_FishTactic/"

library("dplyr")
source(paste0(Root_dir,"Fishtactics_SDM_code_check/environment/function_DMwR_utils.R"))

EVdata_Re = read.csv("./raw_data/Environmental data/original_data/EV(Re)_YM.csv")
EVdata_Bathy = read.csv("./raw_data/Environmental data/original_data/BathyOcean.csv")
EVdata_Chla = read.csv("./raw_data/Environmental data/original_data/chla_YM.csv")
EVdata_GM = read.csv("./raw_data/Environmental data/original_data/EV(GM)_YM.csv")

# merge to data at 0.1*0.1 resolution
grid_size_lat = 0.1
grid_size_lon = 0.1

lat_vec =  seq(23,30,grid_size_lat) # the cut point of the grids
lon_vec =  seq(119,126,grid_size_lon)

lat_grid = seq(23.05,29.95,grid_size_lat) # the central point of the grid
lon_grid = seq(119.05,125.95,grid_size_lat)
grid_info = expand.grid(lat_grid,lon_grid)

lat_grid_lab = c(1:length(lat_grid))
lon_grid_lab = c(1:length(lon_grid))
grid_info_lab = expand.grid(lat_grid_lab,lon_grid_lab)
grid_info$grid_lab = paste(grid_info_lab$Var1,grid_info_lab$Var2,sep="-")

EVdata_YM = NULL
y=2009
for(m in 8:12){
  # Bathy data
  lat_Bathy = as.character(cut(EVdata_Bathy$Lat+0.05,lat_vec,lat_grid_lab))
  lon_Bathy = as.character(cut(EVdata_Bathy$Lon+0.05,lon_vec,lon_grid_lab))
  EVdata_Bathy$grid_lab = paste(lat_Bathy,lon_Bathy,sep="-")
  EVdata_Bathy_allgrid = merge(x=EVdata_Bathy,y=grid_info,by="grid_lab",all.x=T)
  
  # Reanalysis data
  EVdata_Re_unit = EVdata_Re %>% filter(Year==y&Month==m)
  
  lat_Re = as.character(cut(EVdata_Re_unit$lat+0.05,lat_vec,lat_grid_lab))
  lon_Re = as.character(cut(EVdata_Re_unit$lon+0.05,lon_vec,lon_grid_lab))
  EVdata_Re_unit$grid_lab = paste(lat_Re,lon_Re,sep="-")
  EVdata_BathyRe_allgrid = merge(x=EVdata_Bathy_allgrid,y=EVdata_Re_unit,by="grid_lab",all.x=T)
  
  # GM data
  EVdata_GM_unit = EVdata_GM %>% filter(Year==y&Month==m)
  lat_GM = as.character(cut(EVdata_GM_unit$Lat+0.05,lat_vec,lat_grid_lab))
  lon_GM = as.character(cut(EVdata_GM_unit$Lon+0.05,lon_vec,lon_grid_lab))
  EVdata_GM_unit$grid_lab = paste(lat_GM,lon_GM,sep="-")
  EVdata_BathyReGM_allgrid = merge(x=EVdata_BathyRe_allgrid,y=EVdata_GM_unit,by="grid_lab",all.x=T)
  
  
  # Chla
  EVdata_Chla_unit = EVdata_Chla %>% filter(Year==y&Month==m)
  
  lat_Chla = as.character(cut(EVdata_Chla_unit$Lat+0.05,lat_vec,lat_grid_lab))
  lon_Chla = as.character(cut(EVdata_Chla_unit$Lon+0.05,lon_vec,lon_grid_lab))
  EVdata_Chla_unit$grid_lab = paste(lat_Chla,lon_Chla,sep="-")
  EVdata_BathyReGMChla_allgrid = merge(x=EVdata_BathyReGM_allgrid,
                                       y=EVdata_Chla_unit,by="grid_lab",all.x=T)
  
  
  EVdata_YM = rbind(EVdata_YM,EVdata_BathyReGMChla_allgrid) 
}

for(y in 2010:2019){
  for(m in 1:12){
    
    # Bathy data
    lat_Bathy = as.character(cut(EVdata_Bathy$Lat+0.05,lat_vec,lat_grid_lab))
    lon_Bathy = as.character(cut(EVdata_Bathy$Lon+0.05,lon_vec,lon_grid_lab))
    EVdata_Bathy$grid_lab = paste(lat_Bathy,lon_Bathy,sep="-")
    EVdata_Bathy_allgrid = merge(x=EVdata_Bathy,y=grid_info,by="grid_lab",all.x=T)
    
    # Reanalysis data
    EVdata_Re_unit = EVdata_Re %>% filter(Year==y&Month==m)
    
    lat_Re = as.character(cut(EVdata_Re_unit$lat+0.05,lat_vec,lat_grid_lab))
    lon_Re = as.character(cut(EVdata_Re_unit$lon+0.05,lon_vec,lon_grid_lab))
    EVdata_Re_unit$grid_lab = paste(lat_Re,lon_Re,sep="-")
    EVdata_BathyRe_allgrid = merge(x=EVdata_Bathy_allgrid,y=EVdata_Re_unit,by="grid_lab",all.x=T)
    
    # GM data
    EVdata_GM_unit = EVdata_GM %>% filter(Year==y&Month==m)
    lat_GM = as.character(cut(EVdata_GM_unit$Lat+0.05,lat_vec,lat_grid_lab))
    lon_GM = as.character(cut(EVdata_GM_unit$Lon+0.05,lon_vec,lon_grid_lab))
    EVdata_GM_unit$grid_lab = paste(lat_GM,lon_GM,sep="-")
    EVdata_BathyReGM_allgrid = merge(x=EVdata_BathyRe_allgrid,y=EVdata_GM_unit,by="grid_lab",all.x=T)
    
    
    # Chla
    EVdata_Chla_unit = EVdata_Chla %>% filter(Year==y&Month==m)
    
    lat_Chla = as.character(cut(EVdata_Chla_unit$Lat+0.05,lat_vec,lat_grid_lab))
    lon_Chla = as.character(cut(EVdata_Chla_unit$Lon+0.05,lon_vec,lon_grid_lab))
    EVdata_Chla_unit$grid_lab = paste(lat_Chla,lon_Chla,sep="-")
    EVdata_BathyReGMChla_allgrid = merge(x=EVdata_BathyReGM_allgrid,
                                       y=EVdata_Chla_unit,by="grid_lab",all.x=T)
    
    
    EVdata_YM = rbind(EVdata_YM,EVdata_BathyReGMChla_allgrid) 
  }
}


EVdata_YM_clear = EVdata_YM[,c(1,4,11:18,23,27:29)]

# k near-neighbor imputation to fill NA
EVdata_YM_fillNA = knnImputation(EVdata_YM_clear[,-1])



EVdata_YM_fillNA = cbind(EVdata_YM$Year,
                         EVdata_YM$Month,
                         EVdata_YM$grid_lab,
                         EVdata_YM_fillNA)
EVdata_YM_fillNA = EVdata_YM_fillNA[,c(1:3,14,15,4:13,16)]
colnames(EVdata_YM_fillNA)[1:3] = c("Year","Month","grid_lab")
write.csv(EVdata_YM_fillNA,"./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv")

# Correlation examination
library(corrplot)
library(RColorBrewer)

EVdata_YM_fillNA = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv")
EVdata_YM_fillNA = EVdata_YM_fillNA[,7:17]

Env_Cor <-cor(EVdata_YM_fillNA)

jpeg("./raw_data/Environmental data/Env_cor.jpeg",width = 4,height = 4,units = "in",res = 300)
corrplot(Env_Cor, type="upper", order="hclust",method="square",
         col=brewer.pal(n=8, name="RdYlBu"),number.cex=0.6,
         addCoef.col = "black",diag = F)
         # Combine with significance
         #p.mat = p.mat, sig.level = 0.01, insig = "blank")

dev.off()

# Plot ----------------------------------------------------------------
library("ggplot2")
library("maps")
library("colorRamps")

# R color palettes
n <- 16
mp <- map_data("world")

EVdata_YM_fillNA <- read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv")
EVdata_YM_fillNA_9to11 = EVdata_YM_fillNA %>% filter(Month%in%c(9:11))

EV_fullname <- c(rep(NA,5),
                 "Bathymetric",
                 "Density ocean mixed layer thickness",
                 "Sea floor potential temperature",
                 "Sea surface height",
                 "Salinity(0.5m)","Temperature(0.5m)",
                 "uo","vo",
                 "Eddy kinetic energy",
                 "Temperature gratitude magnitude",
                 "Chlorophyll a concentration")

for (i in 6) {
  EV_YM <- EVdata_YM_fillNA_9to11[,c(1,2,4,5,i)]
  colnames(EV_YM) <- c("Year","Month","lon","lat", "Value")
  
  plot <-  ggplot(data=EV_YM , aes(x=lon, y=lat, fill=Value)) + 
    geom_raster(interpolate=TRUE) +
    geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
    scale_fill_gradientn(colours = blue2green2red(n) ,na.value = NA) +
    theme_bw() + facet_wrap(Month~Year, ncol=11) +  
    coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))+
    ggtitle(EV_fullname[i])
  
  ggsave(plot, file = paste0("./raw_data/Environmental data/FillNA/",colnames(EVdata_YM_fillNA)[i],".png") ,width =28, height =10, dpi=96, units = "in", device='png',limitsize = FALSE) 
}

for (i in c(7:11,14:16)) {
  EV_YM <- EVdata_YM_fillNA_9to11[,c(1,2,4,5,i)]
  colnames(EV_YM) <- c("Year","Month","lon","lat", "Value")
  
  plot <-  ggplot(data=EV_YM , aes(x=lon, y=lat, fill=Value)) + 
    geom_raster(interpolate=TRUE) +
    geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = "white") +
    scale_fill_gradientn(colours = blue2green2red(n) ,na.value = NA) +
    theme_bw() + facet_wrap(Month~Year, ncol=11) +  
    coord_fixed(1, xlim = c(119, 126), ylim = c(23, 30))+
    ggtitle(EV_fullname[i])
  
  ggsave(plot, file = paste0("./raw_data/Environmental data/FillNA/",colnames(EVdata_YM_fillNA)[i],".png") ,width =28, height =10, dpi=96, units = "in", device='png',limitsize = FALSE) 
}

