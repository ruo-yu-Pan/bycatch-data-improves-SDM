wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data/"
setwd(wd)

library("sp")
library("raster")
library("rgdal")
library("MASS")
library("dplyr")


Env_data = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv")
Env_data = Env_data %>% filter(Year ==2009,Month==9)
map_pt = t(xtabs(~Lon+Lat,data=Env_data))
map_pt = ifelse(map_pt==0,NA,1)



kernel_density_func = function(occ,effort_kd_filename,m,sz,r){
  effort_dens = kde2d(x=occ$Lon,y=occ$Lat,
                      n=c(70,70),
                      lims = c(119,126,23,30))
  dens_est = t(effort_dens[[3]])*map_pt
  dens_est = dens_est[70:1,]
  effort_dens.ras <- raster(dens_est)
  extent(effort_dens.ras) <- c(119.05,125.95,23.05,29.95)
  crs(effort_dens.ras) <- "+init=epsg:4326"
  
   writeRaster(effort_dens.ras, 
               paste0("./compiled_data/",effort_kd_filename,".tif"),
               overwrite=T)
  
  # transform to long dataframe ---------------------
  raspt <- as.data.frame(rasterToPoints(effort_dens.ras))
  colnames(raspt) = c("Lon","Lat","kd")
  raspt$Month = m
  raspt$Smp_size = sz
  raspt$resample = r
  
  return(raspt)
}

kernel_density_notwrite_func = function(occ,m,sz,r){
  effort_dens = kde2d(x=occ$Lon,y=occ$Lat,
                      n=c(70,70),
                      lims = c(119,126,23,30))
  dens_est = t(effort_dens[[3]])*map_pt
  dens_est = dens_est[70:1,]
  effort_dens.ras <- raster(dens_est)
  extent(effort_dens.ras) <- c(119.05,125.95,23.05,29.95)
  crs(effort_dens.ras) <- "+init=epsg:4326"
  
  # writeRaster(effort_dens.ras, 
  #             paste0("./compiled_data/",effort_kd_filename,".tif"),
  #             overwrite=T)
  
  # transform to long dataframe ---------------------
  raspt <- as.data.frame(rasterToPoints(effort_dens.ras))
  colnames(raspt) = c("Lon","Lat","kd")
  raspt$Month = m
  raspt$Smp_size = sz
  raspt$resample = r
  
  return(raspt)
}