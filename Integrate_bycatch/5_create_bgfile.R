
#########################################################
# This script create the bg point 
# without and with consideration of sampling bias
#########################################################

wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)


# reference for with sampling bias: 
# https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/
# and
# https://github.com/jamiemkass/ENMeval/issues/26

library("sp")
library("raster")
library("rgdal")
library("MASS")
library("dplyr")


# grid information
grid_size_lat = 0.1
grid_size_lon = 0.1
lat_vec =  seq(23,30,grid_size_lat)
lon_vec =  seq(119,126,grid_size_lon)

# map information
Env_data = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv",stringsAsFactors = F)
Env_data_oneyr = Env_data %>% dplyr::filter(Year ==2009,Month==9)
map_pt = t(xtabs(~Lon+Lat,data=Env_data_oneyr))
map_pt = ifelse(map_pt==0,NA,1)

# import env in ground data
Env = list()
for(m in 1:3){
  Env_mth = Env_data%>%filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}


# -------------------------------------------------------------------------
# with sampling bias --------------------------------------------------------
# - bg = 10000 #####

## both ####
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both/",full.names=T)

for(m in 9:11){
  
  Cutlassfish_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Cutlassfish_occ = Cutlassfish_occ[which(Cutlassfish_occ$fishing_loc==1),]
  #spatialpresence = unique(presence)
  effort_dens = kde2d(x=Cutlassfish_occ$Lon,y=Cutlassfish_occ$Lat,
                      n=c(70,70),
                      lims = c(119,126,23,30))
  dens_est = t(effort_dens[[3]])*map_pt
  dens_est = dens_est[70:1,]
  effort_dens.ras <- raster(dens_est)
  extent(effort_dens.ras) <- c(119.05,125.95,23.05,29.95)
  crs(effort_dens.ras) <- "+init=epsg:4326"
  
  # select background point based on bias
  biased_bg = NULL
  n_bg = 1000
  effort_dens.ras.nona = !is.na(effort_dens.ras)
  
  for (y in 2009:2019){
    biased_bg_yr <- xyFromCell(effort_dens.ras, 
                               sample(which(!is.na(values(effort_dens.ras))), 
                                      n_bg, 
                                      prob=values(effort_dens.ras)[which(!is.na(values(effort_dens.ras)))],
                                      replace = T))
    
    env_lat_grp = cut(biased_bg_yr[,2],lat_vec,c(1:70))
    env_lon_grp = cut(biased_bg_yr[,1],lon_vec,c(1:70))
    biased_bg_yr = as.data.frame(biased_bg_yr)
    biased_bg_yr$grid_lab=paste0(env_lat_grp,"-",env_lon_grp)
    
    Env_mth = Env[[m-8]]
    Env_mth_yr = Env_mth%>%filter(Year==y)
    
    biased_bg_env_yr = merge(x=biased_bg_yr,y=Env_mth_yr,by="grid_lab",all.x=T)
    biased_bg_env_yr = biased_bg_env_yr[,-c(2:4)]
    
    biased_bg = rbind(biased_bg,biased_bg_env_yr) 
  }
  
  biased_bg$n_fold = sample(c(1:5),nrow(biased_bg),replace = T)
  write.csv(biased_bg,paste0("./compiled_data/bg_file/grid_both/bg10e4Bgbias_",
                             formatC(m,width = 2,flag = 0),".csv"))
}



## bycatch #####
### - origin
Cutlassfish_bycatch_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_bycatch/",full.names=T)

for(m in 9:11){
  #N_lab = length(resample_size_bycatch[[m-8]])
  Cutlassfish_occ = read.csv(Cutlassfish_bycatch_occ_filename[m],row.names = 1)
  Cutlassfish_fishLoc = Cutlassfish_occ[which(Cutlassfish_occ$fishing_loc==1),]
  #spatialpresence = unique(presence)
  effort_dens = kde2d(x=Cutlassfish_fishLoc$Lon,y=Cutlassfish_fishLoc$Lat,
                      n=c(70,70),
                      lims = c(119,126,23,30))
  dens_est = t(effort_dens[[3]])*map_pt
  dens_est = dens_est[70:1,]
  effort_dens.ras <- raster(dens_est)
  extent(effort_dens.ras) <- c(119.05,125.95,23.05,29.95)
  crs(effort_dens.ras) <- "+init=epsg:4326"
  
  # select background point based on bias
  biased_bg = NULL
  n_bg = 1000
  effort_dens.ras.nona = !is.na(effort_dens.ras)
  
  for (y in 2009:2019){
    biased_bg_yr <- xyFromCell(effort_dens.ras, 
                               sample(which(!is.na(values(effort_dens.ras))), 
                                      n_bg, 
                                      prob=values(effort_dens.ras)[which(!is.na(values(effort_dens.ras)))],
                                      replace = T))
    
    env_lat_grp = cut(biased_bg_yr[,2],lat_vec,c(1:70))
    env_lon_grp = cut(biased_bg_yr[,1],lon_vec,c(1:70))
    biased_bg_yr = as.data.frame(biased_bg_yr)
    biased_bg_yr$grid_lab=paste0(env_lat_grp,"-",env_lon_grp)
    
    Env_mth = Env[[m-8]]
    Env_mth_yr = Env_mth%>%filter(Year==y)
    
    biased_bg_env_yr = merge(x=biased_bg_yr,y=Env_mth_yr,by="grid_lab",all.x=T)
    biased_bg_env_yr = biased_bg_env_yr[,-c(2:4)]
    
    biased_bg = rbind(biased_bg,biased_bg_env_yr) 
  }
  
  biased_bg$n_fold = sample(c(1:5),nrow(biased_bg),replace = T)
  write.csv(biased_bg,paste0("compiled_data/bg_file/grid_bycatch/bg10e4Bgbias_",
                             formatC(m,width = 2,flag = 0),".csv"))
}


## target #####
### - origin
Cutlassfish_target_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_target/",full.names=T)

for(m in 9:11){
  #N_lab = length(resample_size_target[[m-8]])
  Cutlassfish_occ = read.csv(Cutlassfish_target_occ_filename[m-3],row.names = 1)
  Cutlassfish_occ = Cutlassfish_occ[which(Cutlassfish_occ$fishing_loc==1),]
  #spatialpresence = unique(presence)
  effort_dens = kde2d(x=Cutlassfish_occ$Lon,y=Cutlassfish_occ$Lat,
                      n=c(70,70),
                      lims = c(119,126,23,30))
  dens_est = t(effort_dens[[3]])*map_pt
  dens_est = dens_est[70:1,]
  effort_dens.ras <- raster(dens_est)
  extent(effort_dens.ras) <- c(119.05,125.95,23.05,29.95)
  crs(effort_dens.ras) <- "+init=epsg:4326"
  
  # select background point based on bias
  biased_bg = NULL
  n_bg = 1000
  effort_dens.ras.nona = !is.na(effort_dens.ras)
  
  for (y in 2009:2019){
    biased_bg_yr <- xyFromCell(effort_dens.ras, 
                               sample(which(!is.na(values(effort_dens.ras))), 
                                      n_bg, 
                                      prob=values(effort_dens.ras)[which(!is.na(values(effort_dens.ras)))],
                                      replace = T))
    
    env_lat_grp = cut(biased_bg_yr[,2],lat_vec,c(1:70))
    env_lon_grp = cut(biased_bg_yr[,1],lon_vec,c(1:70))
    biased_bg_yr = as.data.frame(biased_bg_yr)
    biased_bg_yr$grid_lab=paste0(env_lat_grp,"-",env_lon_grp)
    
    Env_mth = Env[[m-8]]
    Env_mth_yr = Env_mth%>%filter(Year==y)
    
    biased_bg_env_yr = merge(x=biased_bg_yr,y=Env_mth_yr,by="grid_lab",all.x=T)
    biased_bg_env_yr = biased_bg_env_yr[,-c(2:4)]
    
    biased_bg = rbind(biased_bg,biased_bg_env_yr) 
  }
  
  biased_bg$n_fold = sample(c(1:5),nrow(biased_bg),replace = T)
  write.csv(biased_bg,paste0("compiled_data/bg_file/grid_target/bg10e4Bgbias_",
                             formatC(m,width = 2,flag = 0),".csv"))
}

