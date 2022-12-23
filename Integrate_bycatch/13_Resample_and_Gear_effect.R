

######## 
# This script is for trying to manipulate a fixed sampling range
########

wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

Root_dir = "D:/Ruo_data/2020_Assis_FishTactic/"

###################################################
# RESAMPLE DATA and DATA PROPERTY #################
library(sp)
library(sf)
library(pracma)
library(ggplot2)

library("rnaturalearth")
library("rnaturalearthdata")
world = ne_countries(scale = "medium", returnclass = "sf")

source(paste0(Root_dir,"Fishtactics_SDM_code_check/Function_CV_fold_generator.R"))
source(paste0(Root_dir,"Fishtactics_SDM_code_check/Function_examine_kernel_density.R"))

grid_in_ground = read.csv("./compiled_data/map/grid_in_ground.csv",row.names = 1)

# import 2D kernel density of sampling effort
#effort_dens_all_bycatch = read.csv("./compiled_data/effort_kernel_density/Whole_data/effort_kd2d_bycatch.csv")
#effort_dens_all_target = read.csv("./compiled_data/effort_kernel_density/Whole_data/effort_kd2d_target.csv")

# -----------------------------------------------------------------------------------
# import presence record data and number of occ data
Cutlassfish_bycatch_filename = list.files("./compiled_data/occ_data_in_ground/grid_bycatch/",full.names=T)
Cutlassfish_target_filename = list.files("./compiled_data/occ_data_in_ground/grid_target/",full.names=T)

num_occ_bycatch = NULL
for(m in 9:11){
  Cutlassfish_bycatch_occ = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
  num_occ_bycatch = c(num_occ_bycatch,length(which(Cutlassfish_bycatch_occ$occ==1)))
}

num_occ_target = NULL
for(m in 6:8){
  Cutlassfish_target_occ = read.csv(Cutlassfish_target_filename[m],row.names = 1)
  num_occ_target = c(num_occ_target,length(which(Cutlassfish_target_occ$occ==1)))
}

# number of data (1798 for each YearMth)
n_data = nrow(Cutlassfish_bycatch_occ)

# resample sample size
samesize = NULL

for(m in 1:3){
  if(num_occ_bycatch[m]>num_occ_target[m]){
    samesize_mth = num_occ_target[m]
    samesize =  c(samesize,samesize_mth)
  }else{
    samesize_mth = num_occ_bycatch[m]
    samesize =  c(samesize,samesize_mth)
  }
}
samesize = round(samesize/2)

# -------------------------------------------------------------------------
# (1) Sample the bycatch data with different sample size 
#     within the overlapped sampling range of targeted data 
# (2) And examine the data property of resample data #####

## -bycatch resample ####
data_property_occ_bycatch_Sbt = NULL
data_property_effort_bycatch_Sbt = NULL
kd_data_occ_bycatch_Sbt = NULL
kd_data_effort_bycatch_Sbt = NULL

for(m in 9:11){
  Cutlassfish_bycatch_occ = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
  Cutlassfish_target_occ = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)
  
  # calculate bycatch data sampling range
  unique_lon_lat_b = unique(as.matrix(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$fishing_loc==1),c("Lon","Lat")]))
  
  hull_b <- chull(unique_lon_lat_b)
  hull_b2 <- c(hull_b, hull_b[1])
  sampling_range_bycatch = 
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat_b[hull_b2,1],
                                          unique_lon_lat_b[hull_b2,2])!=0),]
  
  # calculate target data sampling range
  unique_lon_lat_t = unique(as.matrix(Cutlassfish_target_occ[which(Cutlassfish_target_occ$fishing_loc==1),c("Lon","Lat")]))
  
  #if(m==9){unique_lon_lat_t = unique_lon_lat_t[-which(unique_lon_lat_t[,2]>26.5),]}
  #if(m==10){unique_lon_lat_t = unique_lon_lat_t[-which(unique_lon_lat_t[,1]<121),]}
  
  hull_t <- chull(unique_lon_lat_t)
  hull_t2 <- c(hull_t, hull_t[1])
  sampling_range_target = 
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat_t[hull_t2,1],
                                          unique_lon_lat_t[hull_t2,2])!=0),]
  
  plot(unique_lon_lat_b[,1],unique_lon_lat_b[,2],xlim=c(118,124.5),ylim=c(24,29))
  points(unique_lon_lat_t[,1],unique_lon_lat_t[,2],col=2)
  
  # resample the data
  intrsectSmpRange_bycatch_pos = which(point.in.polygon(Cutlassfish_bycatch_occ$Lon,
                                                        Cutlassfish_bycatch_occ$Lat,
                                                        unique_lon_lat_b[hull_b2,1],
                                                        unique_lon_lat_b[hull_b2,2])!=0 &
                                         point.in.polygon(Cutlassfish_bycatch_occ$Lon,
                                                          Cutlassfish_bycatch_occ$Lat,
                                                          unique_lon_lat_t[hull_t2,1],
                                                          unique_lon_lat_t[hull_t2,2])!=0 )
  Cutlass_bycatch_occ_pos = which(Cutlassfish_bycatch_occ$occ==1)
  intrsectSmpRange_Cutlass_bycatch_occ_pos = which(c(1:n_data) %in% intrsectSmpRange_bycatch_pos &
                                                     c(1:n_data) %in% Cutlass_bycatch_occ_pos)
  
  if(length(intrsectSmpRange_Cutlass_bycatch_occ_pos)>=samesize[m-8]){
    for(r in 1:20){
      sample_unit = sample(x=intrsectSmpRange_Cutlass_bycatch_occ_pos,
                           samesize[m-8],replace=F)
      intrsec_Cutlass_resample_bycatch = Cutlassfish_bycatch_occ
      intrsec_Cutlass_resample_bycatch$fishing_loc_resample = intrsec_Cutlass_resample_bycatch$fishing_loc
      intrsec_Cutlass_resample_bycatch[which(c(1:n_data) %in% sample_unit ==F),c("occ","fishing_loc_resample")]=0
      
      intrsec_Cutlass_resample_bycatch = generate_CV_fold(intrsec_Cutlass_resample_bycatch,nfold = 5)
      
      
      write.csv(intrsec_Cutlass_resample_bycatch,
                paste0("./compiled_data/occ_data_in_ground/resample/",
                       formatC(m, width = 2, flag = 0),
                       "/bycatch_r",formatC(r, width = 2, flag = 0),".csv"))
      
      effort_smpsize = length(which(intrsec_Cutlass_resample_bycatch$fishing_loc_resample==1))
      
      # Data Property -------------------------------------------------
      
      #-----------------------------------------------------------------------------
      # for occ
      # - calculate sampling range
      unique_lon_lat = unique(as.matrix(intrsec_Cutlass_resample_bycatch[sample_unit,c("Lon","Lat")]))
      real_intrsec_hull_b <- chull(unique_lon_lat)
      real_intrsec_hull_b <- c(real_intrsec_hull_b, real_intrsec_hull_b[1])
      
      real_intrsec_sampling_range = 
        grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                              grid_in_ground$Lat,
                                              unique_lon_lat[real_intrsec_hull_b,1],
                                              unique_lon_lat[real_intrsec_hull_b,2])!=0),]
      extent_sampling_range = 100*(nrow(real_intrsec_sampling_range)/nrow(grid_in_ground))
      
      # - number of site
      n_site = nrow(unique_lon_lat)
      
      # - calculate the aggregation of sampling probability
      # -- CV of kernel density
      kd_resample = kernel_density_func(intrsec_Cutlass_resample_bycatch[which(intrsec_Cutlass_resample_bycatch$occ==1),],
                                        paste0("occ_kernel_density/Resample_data/tif/bycatch_",m,"_",r),
                                        m,NA,r)
      kd_range = 
        kd_resample[which(point.in.polygon(kd_resample$Lon,
                                           kd_resample$Lat,
                                           unique_lon_lat[real_intrsec_hull_b,1],
                                           unique_lon_lat[real_intrsec_hull_b,2])!=0),]
      kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
      
      # -- distribution of kernel density
      kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100, 
                           kd_val = sort(kd_range$kd/max(kd_range$kd)))
      #plot(kd_rad,xlab = "std_rank",ylab = "std_kd",pch = 19,type = "b",lwd = 0.5)
      kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
      
      # save
      kd_data_occ_bycatch_Sbt = rbind(kd_data_occ_bycatch_Sbt,kd_range)
      data_property = c(m, samesize[m-8],effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site,r)
      data_property_occ_bycatch_Sbt = rbind(data_property_occ_bycatch_Sbt,data_property)
      
      #-----------------------------------------------------------------------------
      # for effort
      # - calculate sampling range
      unique_lon_lat = unique(as.matrix(intrsec_Cutlass_resample_bycatch[
        which(intrsec_Cutlass_resample_bycatch$fishing_loc_resample==1),c("Lon","Lat")]))
      real_intrsec_hull_b <- chull(unique_lon_lat)
      real_intrsec_hull_b <- c(real_intrsec_hull_b, real_intrsec_hull_b[1])
      
      real_intrsec_sampling_range = 
        grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                              grid_in_ground$Lat,
                                              unique_lon_lat[real_intrsec_hull_b,1],
                                              unique_lon_lat[real_intrsec_hull_b,2])!=0),]
      extent_sampling_range = 100*(nrow(real_intrsec_sampling_range)/nrow(grid_in_ground))
      
      # - number of site
      n_site = nrow(unique_lon_lat)
      
      # - calculate the aggregation of sampling probability
      # -- CV of kernel density
      kd_resample = kernel_density_func(intrsec_Cutlass_resample_bycatch[which(intrsec_Cutlass_resample_bycatch$fishing_loc_resample==1),],
                                        paste0("effort_kernel_density/Resample_data/tif/bycatch_",m,"_",r),
                                        m,NA,r)
      kd_range = 
        kd_resample[which(point.in.polygon(kd_resample$Lon,
                                           kd_resample$Lat,
                                           unique_lon_lat[real_intrsec_hull_b,1],
                                           unique_lon_lat[real_intrsec_hull_b,2])!=0),]
      kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
      
      # -- distribution of kernel density
      kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100, 
                           kd_val = sort(kd_range$kd/max(kd_range$kd)))
      #plot(kd_rad,xlab = "std_rank",ylab = "std_kd",pch = 19,type = "b",lwd = 0.5)
      kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
      
      # save
      kd_data_effort_bycatch_Sbt = rbind(kd_data_effort_bycatch_Sbt,kd_range)
      data_property = c(m, samesize[m-8], effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site,r)
      data_property_effort_bycatch_Sbt = rbind(data_property_effort_bycatch_Sbt,data_property)
    }
  }  
}


data_property_occ_bycatch_Sbt = as.data.frame(data_property_occ_bycatch_Sbt)
data_property_effort_bycatch_Sbt = as.data.frame(data_property_effort_bycatch_Sbt)
colnames(data_property_occ_bycatch_Sbt) = c("Month","Smp_size","effort_smpsize","Extent_smp_range","kd_CV","kd_rad_auc","No_site","resample")
colnames(data_property_effort_bycatch_Sbt) = c("Month","Smp_size","effort_smpsize","Extent_smp_range","kd_CV","kd_rad_auc","No_site","resample")

write.csv(kd_data_occ_bycatch_Sbt,"./compiled_data/occ_kernel_density/Resample_data/kd_data_bycatch.csv")
write.csv(kd_data_effort_bycatch_Sbt,"./compiled_data/effort_kernel_density/Resample_data/kd_data_bycatch.csv")
write.csv(data_property_occ_bycatch_Sbt,"./compiled_data/data_property/resample_bycatch_occ.csv")
write.csv(data_property_effort_bycatch_Sbt,"./compiled_data/data_property/resample_bycatch_effort.csv")

## -target resample ####
data_property_occ_target_Stb = NULL
data_property_effort_target_Stb = NULL
kd_data_occ_target_Stb = NULL
kd_data_effort_target_Stb = NULL


for(m in 9:11){
  Cutlassfish_bycatch_occ = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
  Cutlassfish_target_occ = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)
  
  # calculate bycatch data sampling range
  unique_lon_lat_b = unique(as.matrix(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$fishing_loc==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat_b)
  hull_b2 <- c(hull_b, hull_b[1])
  sampling_range_bycatch = 
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat_b[hull_b2,1],
                                          unique_lon_lat_b[hull_b2,2])!=0),]
  
  # calculate target data sampling range
  unique_lon_lat_t = unique(as.matrix(Cutlassfish_target_occ[which(Cutlassfish_target_occ$fishing_loc==1),c("Lon","Lat")]))
  # if(m==9){unique_lon_lat_t = unique_lon_lat_t[-which(unique_lon_lat_t[,2]>26.5),]}
  # if(m==10){unique_lon_lat_t = unique_lon_lat_t[-which(unique_lon_lat_t[,1]<121),]}
  
  hull_t <- chull(unique_lon_lat_t)
  hull_t2 <- c(hull_t, hull_t[1])
  sampling_range_target = 
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat_t[hull_t2,1],
                                          unique_lon_lat_t[hull_t2,2])!=0),]
  
  plot(unique_lon_lat_b[,1],unique_lon_lat_b[,2],xlim=c(118,124.5),ylim=c(24,28))
  points(unique_lon_lat_t[,1],unique_lon_lat_t[,2],col=2)
  
  # resample the data
  intrsectSmpRange_target_pos = which(point.in.polygon(Cutlassfish_target_occ$Lon,
                                                       Cutlassfish_target_occ$Lat,
                                                       unique_lon_lat_b[hull_b2,1],
                                                       unique_lon_lat_b[hull_b2,2])!=0 &
                                        point.in.polygon(Cutlassfish_target_occ$Lon,
                                                         Cutlassfish_target_occ$Lat,
                                                         unique_lon_lat_t[hull_t2,1],
                                                         unique_lon_lat_t[hull_t2,2])!=0 )
  Cutlass_target_occ_pos = which(Cutlassfish_target_occ$occ==1)
  intrsectSmpRange_Cutlass_target_occ_pos = which(c(1:n_data) %in% intrsectSmpRange_target_pos &
                                                    c(1:n_data) %in% Cutlass_target_occ_pos)

    if(length(intrsectSmpRange_Cutlass_target_occ_pos)>=samesize[m-8]){
      for(r in 1:20){
        sample_unit = sample(x=intrsectSmpRange_Cutlass_target_occ_pos,
                             samesize[m-8],replace=F)
        intrsec_Cutlass_resample_target = Cutlassfish_target_occ
        intrsec_Cutlass_resample_target$fishing_loc_resample = intrsec_Cutlass_resample_target$fishing_loc
        intrsec_Cutlass_resample_target[which(c(1:n_data) %in% sample_unit ==F),c("occ","fishing_loc_resample")]=0
        
        intrsec_Cutlass_resample_target = generate_CV_fold(intrsec_Cutlass_resample_target,nfold = 5)
        
        write.csv(intrsec_Cutlass_resample_target,
                  paste0("./compiled_data/occ_data_in_ground/resample/",
                         formatC(m, width = 2, flag = 0),
                         "/target_r",formatC(r, width = 2, flag = 0),".csv"))
        
        effort_smpsize = length(which(intrsec_Cutlass_resample_target$fishing_loc_resample==1))
        # Data Property -------------------------------------------------
        
        #-----------------------------------------------------------------------------
        # for occ
        # - calculate sampling range
        unique_lon_lat = unique(as.matrix(intrsec_Cutlass_resample_target[sample_unit,c("Lon","Lat")]))
        real_intrsec_hull_t <- chull(unique_lon_lat)
        real_intrsec_hull_t <- c(real_intrsec_hull_t, real_intrsec_hull_t[1])
        
        real_intrsec_sampling_range = 
          grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                                grid_in_ground$Lat,
                                                unique_lon_lat[real_intrsec_hull_t,1],
                                                unique_lon_lat[real_intrsec_hull_t,2])!=0),]
        extent_sampling_range = 100*(nrow(real_intrsec_sampling_range)/nrow(grid_in_ground))
        
        # - number of site
        n_site = nrow(unique_lon_lat)
        
        # - calculate the aggregation of sampling probability
        # -- CV of kernel density
        kd_resample = kernel_density_func(intrsec_Cutlass_resample_target[which(intrsec_Cutlass_resample_target$occ==1),],
                                          paste0("occ_kernel_density/Resample_data/tif/target_",m,"_",r),
                                          m,NA,r)
        kd_range = 
          kd_resample[which(point.in.polygon(kd_resample$Lon,
                                             kd_resample$Lat,
                                             unique_lon_lat[real_intrsec_hull_t,1],
                                             unique_lon_lat[real_intrsec_hull_t,2])!=0),]
        kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
        
        # -- distribution of kernel density
        kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100, 
                             kd_val = sort(kd_range$kd/max(kd_range$kd)))
        #plot(kd_rad,xlab = "std_rank",ylab = "std_kd",pch = 19,type = "b",lwd = 0.5)
        kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
        
        # save
        kd_data_occ_target_Stb = rbind(kd_data_occ_target_Stb,kd_range)
        
        data_property = c(m, samesize[[m-8]],effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site,r)
        data_property_occ_target_Stb = rbind(data_property_occ_target_Stb,data_property)
        
        #-----------------------------------------------------------------------------
        # for effort
        # - calculate sampling range
        unique_lon_lat = unique(as.matrix(intrsec_Cutlass_resample_target[
          which(intrsec_Cutlass_resample_target$fishing_loc_resample==1),c("Lon","Lat")]))
        real_intrsec_hull_t <- chull(unique_lon_lat)
        real_intrsec_hull_t <- c(real_intrsec_hull_t, real_intrsec_hull_t[1])
        
        real_intrsec_sampling_range = 
          grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                                grid_in_ground$Lat,
                                                unique_lon_lat[real_intrsec_hull_t,1],
                                                unique_lon_lat[real_intrsec_hull_t,2])!=0),]
        extent_sampling_range = 100*(nrow(real_intrsec_sampling_range)/nrow(grid_in_ground))
        
        # - number of site
        n_site = nrow(unique_lon_lat)
        
        # - calculate the aggregation of sampling probability
        # -- CV of kernel density
        kd_resample = kernel_density_func(intrsec_Cutlass_resample_target[which(intrsec_Cutlass_resample_target$fishing_loc==1),],
                                          paste0("effort_kernel_density/Resample_data/tif/target_",m,"_",r),
                                          m,NA,r)
        kd_range = 
          kd_resample[which(point.in.polygon(kd_resample$Lon,
                                             kd_resample$Lat,
                                             unique_lon_lat[real_intrsec_hull_t,1],
                                             unique_lon_lat[real_intrsec_hull_t,2])!=0),]
        kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
        
        # -- distribution of kernel density
        kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100, 
                             kd_val = sort(kd_range$kd/max(kd_range$kd)))
        #plot(kd_rad,xlab = "std_rank",ylab = "std_kd",pch = 19,type = "b",lwd = 0.5)
        kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
        
        # save
        kd_data_effort_target_Stb = rbind(kd_data_effort_target_Stb,kd_range)
        
        data_property = c(m, samesize[m-8],effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site,r)
        data_property_effort_target_Stb = rbind(data_property_effort_target_Stb,data_property)
        
      }
    }
}


data_property_occ_target_Stb = as.data.frame(data_property_occ_target_Stb)
data_property_effort_target_Stb = as.data.frame(data_property_effort_target_Stb)
colnames(data_property_occ_target_Stb) = c("Month","Smp_size","effort_smpsize","Extent_smp_range","kd_CV","kd_rad_auc","No_site","resample")
colnames(data_property_effort_target_Stb) = c("Month","Smp_size","effort_smpsize","Extent_smp_range","kd_CV","kd_rad_auc","No_site","resample")

write.csv(kd_data_occ_target_Stb,"./compiled_data/occ_kernel_density/Resample_data/kd_data_target.csv")
write.csv(kd_data_effort_target_Stb,"./compiled_data/effort_kernel_density/Resample_data/kd_data_target.csv")
write.csv(data_property_occ_target_Stb,"./compiled_data/data_property/resample_target_occ.csv")
write.csv(data_property_effort_target_Stb,"./compiled_data/data_property/resample_target_effort.csv")

data_property_effort_target_Stb = read.csv("./compiled_data/data_property/resample_target_effort.csv")
data_property_effort_bycatch_Sbt = read.csv("./compiled_data/data_property/resample_bycatch_effort.csv")

tapply(data_property_effort_target_Stb$Extent_smp_range,
       data_property_effort_target_Stb$Month,
       summary)
tapply(data_property_effort_bycatch_Sbt$Extent_smp_range,
       data_property_effort_bycatch_Sbt$Month,
       summary)
tapply(data_property_occ_target_Stb$kd_rad_auc,
       data_property_occ_target_Stb$Month,
       summary)
tapply(data_property_occ_bycatch_Sbt$kd_rad_auc,
       data_property_occ_bycatch_Sbt$Month,
       summary)

###################################################
# CREATE BG POINTS ###############################

library("sp")
library("raster")
library("rgdal")
library("MASS")


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
  Env_mth = Env_data%>%dplyr::filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}

# import resample sample size 
for(m in 9:11){
  Cutlassfish_bycatch_resample_filename = 
    list.files(paste0("./compiled_data/occ_data_in_ground/resample/",formatC(m,width=2,flag = 0)),
               pattern = glob2rx("bycatch*.csv"),full.names=T)
  Cutlassfish_bycatch_resample_kd_filename = 
    list.files("./compiled_data/effort_kernel_density/Resample_data/tif",
               pattern = glob2rx(paste0("bycatch_",m,"*.tif")),full.names=T)
  
  for(f in 1:length(Cutlassfish_bycatch_resample_filename)){
    rX = substr(Cutlassfish_bycatch_resample_filename[f],56,58)
    
    # import effort kernel density
    effort_dens.ras = raster(Cutlassfish_bycatch_resample_kd_filename[f])
    
    # select background point based on bias
    biased_bg = NULL
    n_bg = 1000
    
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
      Env_mth_yr = Env_mth%>%dplyr::filter(Year==y)
      
      biased_bg_env_yr = merge(x=biased_bg_yr,y=Env_mth_yr,by="grid_lab",all.x=T)
      biased_bg_env_yr = biased_bg_env_yr[,-c(2:4)]
      
      biased_bg = rbind(biased_bg,biased_bg_env_yr) 
    }
    biased_bg$n_fold = sample(c(1:5),nrow(biased_bg),replace = T)
    
    write.csv(biased_bg,paste0("compiled_data/bg_file/resample/bycatch_",
                               formatC(m,width = 2,flag = 0),"_",rX,".csv"))
  }
}

for(m in 9:11){
  Cutlassfish_target_resample_filename = 
    list.files(paste0("./compiled_data/occ_data_in_ground/resample/",formatC(m,width=2,flag = 0)),
               pattern = glob2rx("target*.csv"),full.names=T)
  Cutlassfish_target_resample_kd_filename = 
    list.files("./compiled_data/effort_kernel_density/Resample_data/tif",
               pattern = glob2rx(paste0("target_",m,"*.tif")),full.names=T)
  
  for(f in 1:length(Cutlassfish_target_resample_filename)){
    rX = substr(Cutlassfish_target_resample_filename[f],55,57)
    
    # import effort kernel density
    effort_dens.ras = raster(Cutlassfish_target_resample_kd_filename[f])
    
    # select background point based on bias
    biased_bg = NULL
    n_bg = 1000
    
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
      Env_mth_yr = Env_mth%>%dplyr::filter(Year==y)
      
      biased_bg_env_yr = merge(x=biased_bg_yr,y=Env_mth_yr,by="grid_lab",all.x=T)
      biased_bg_env_yr = biased_bg_env_yr[,-c(2:4)]
      
      biased_bg = rbind(biased_bg,biased_bg_env_yr) 
    }
    biased_bg$n_fold = sample(c(1:5),nrow(biased_bg),replace = T)
    
    write.csv(biased_bg,paste0("./compiled_data/bg_file/resample/target_",
                               formatC(m,width = 2,flag = 0),"_",rX,".csv"))
  }
}

#########################################################################
# TRAIN MAXENT ####

library(dismo)
library(dplyr)
library("rJava")
library(ggplot2)
library(gridExtra)
library(viridis)

Sys.setenv(JAVA_Home="./Maxent_Java/jdk-11.0.2")

# function for MaxEnt model training and prediction
Maxent_train_pred_func = function(occ_data=Cutlassfish_occ,
                                  Env_mth_data=Env_mth,
                                  bg_env_data=bg_env_mth,
                                  save_dir){
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = occ_data %>%dplyr::filter(occ==1&Year==y)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  # MaxEnt for each fold
  Maxent_model_summary = NULL
  Maxent_model_list = list()
  pred_occ_test = NULL
  for(f in 1:5){
    mod_occ_train = c(occ_env$occ[which(occ_env$nth_fold!=f)],
                      rep(0,length(which(bg_env_data$n_fold!=f))))
    mod_env_train = rbind(occ_env[which(occ_env$nth_fold!=f),
                                  c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")],
                          bg_env_data[which(bg_env_data$n_fold!=f),
                                      c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")])
    
    # maxent
    
    ifelse(!dir.exists(paste0(save_dir,"/CV_",f)), 
           dir.create(file.path(paste0(save_dir,"/CV_",f))),FALSE)
    
    maxent_mod = maxent(x=mod_env_train,
                        p=mod_occ_train,
                        args=c(
                          'pictures=true',
                          'linear=true',
                          'quadratic=true',
                          'product=false',
                          'threshold=false',
                          'hinge=true',
                          'responsecurves=true',
                          'jackknife=true',
                          'askoverwrite=true'
                          #'replicatetype=',
                          #'replicates=5',
                          #'randomseed=true'
                        ),
                        path=paste0(paste0(save_dir,"/CV_",f)))
    
    mod_result = maxent_mod@results
    Maxent_model_summary = cbind(Maxent_model_summary,mod_result)
    Maxent_model_list[[f]] = maxent_mod
    
    # predict for test data
    
    mod_occ_test = c(occ_env$occ[which(occ_env$nth_fold==f)],
                     rep(0,length(which(bg_env_data$n_fold==f))))
    mod_env_test = rbind(occ_env[which(occ_env$nth_fold==f),
                                 c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")],
                         bg_env_data[which(bg_env_data$n_fold==f),
                                     c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")])
    
    pred_occ = data.frame(pred = predict(maxent_mod,mod_env_test))
    pred_occ = cbind(mod_occ_test,pred_occ)
    pred_occ$nfold = f
    pred_occ_test = rbind(pred_occ_test,pred_occ)
    
  }
  output_list = list(Maxent_model_summary,Maxent_model_list,pred_occ_test)
}

### -resample
Cutlassfish_resample_foldername=list.files("./compiled_data/occ_data_in_ground/resample/",full.names=T)


for (m in 9:11){
  # list the bycatch resample file
  Cutlassfish_bycatch_resample_filename = list.files(Cutlassfish_resample_foldername[m-8],pattern = glob2rx("bycatch*.csv"),full.names=T)
  # list the bg file of bycatch resample data 
  bg_env_bycatch_resample_filename=list.files("./compiled_data/bg_file/resample/",
                                              pattern =glob2rx(paste0("bycatch_",formatC(m,width = 2,flag = 0),"*.csv")) ,full.names=T)
  
  for(file in 1: length(Cutlassfish_bycatch_resample_filename)){
    
    rX = substr(Cutlassfish_bycatch_resample_filename[file],56,58)
    r = as.numeric(substr(rX,2,3))
    
    save_dir = paste0("./output/Maxent_5CV/Maxent_result_detail/resample/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias/bycatch_",rX)
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    Maxent_output = 
      Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_bycatch_resample_filename[file],row.names = 1,stringsAsFactors = F),
                             Env_mth_data=Env[[m-8]],
                             bg_env_data=read.csv(bg_env_bycatch_resample_filename[file],row.names = 1,stringsAsFactors = F),
                             save_dir=save_dir)
    # save
    if(r==1){
      # clean and empty data
      Maxent_model_summary_resample = NULL
      Maxent_model_resample = list()
      Maxent_pred_resample = NULL
      
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
    }else if(r==20){
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
      
      write.csv(Maxent_model_summary_resample,paste0("./output/Maxent_5CV/Maxent_result_summary/resample/Maxent_Summary_bg10e4Bgbias_bycatch_",formatC(m,width = 2,flag = 0),".csv"))
      save(Maxent_model_resample,file=paste0("./output/Maxent_5CV/Maxent_result_summary/resample/Maxent_mod_bg10e4Bgbias_bycatch_",formatC(m,width = 2,flag = 0),".RData"))
      write.csv(Maxent_pred_resample,paste0("./output/Maxent_5CV/Maxent_pred/resample/pred_bg10e4Bgbias_bycatch_",formatC(m,width = 2,flag = 0),".csv"))
    }else{
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
    }
  }
}
for (m in 9:11){
  # list the target resample file
  Cutlassfish_target_resample_filename = list.files(Cutlassfish_resample_foldername[m-8],pattern = glob2rx("target*.csv"),full.names=T)
  # list the bg file of target resample data 
  bg_env_target_resample_filename=list.files("./compiled_data/bg_file/resample/",
                                              pattern =glob2rx(paste0("target_",formatC(m,width = 2,flag = 0),"*.csv")) ,full.names=T)
  
  for(file in 1: length(Cutlassfish_target_resample_filename)){
    rX = substr(Cutlassfish_target_resample_filename[file],55,57)
    r = as.numeric(substr(rX,2,3))
    
    save_dir = paste0("./output/Maxent_5CV/Maxent_result_detail/resample/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias/target_",rX)
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    Maxent_output = 
      Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_target_resample_filename[file],row.names = 1,stringsAsFactors = F),
                             Env_mth_data=Env[[m-8]],
                             bg_env_data=read.csv(bg_env_target_resample_filename[file],row.names = 1,stringsAsFactors = F),
                             save_dir=save_dir)
    # save
    if(r==1){
      # clean and empty data
      Maxent_model_summary_resample = NULL
      Maxent_model_resample = list()
      Maxent_pred_resample = NULL
      
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
    }else if(r==20){
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
      
      write.csv(Maxent_model_summary_resample,paste0("./output/Maxent_5CV/Maxent_result_summary/resample/Maxent_Summary_bg10e4Bgbias_target_",formatC(m,width = 2,flag = 0),".csv"))
      save(Maxent_model_resample,file=paste0("./output/Maxent_5CV/Maxent_result_summary/resample/Maxent_mod_bg10e4Bgbias_target_",formatC(m,width = 2,flag = 0),".RData"))
      write.csv(Maxent_pred_resample,paste0("./output/Maxent_5CV/Maxent_pred/resample/pred_bg10e4Bgbias_target_",formatC(m,width = 2,flag = 0),".csv"))
    }else{
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
    }
  }
}

#########################################################################
# PREDICTION IN FISHING GROUND FOR AUC ####
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)

random_bg_filename=list.files("./compiled_data/bg_file/random_for_validation",full.names=T)

Maxent_mod_summary_filename = list.files("./output/Maxent_5CV/Maxent_result_summary/resample",
                                         pattern = glob2rx("*bg10e4Bgbias*.RData"), full.names = T)

for(file in 1:length(Maxent_mod_summary_filename)){
  
  dash_pos = gregexpr("_",Maxent_mod_summary_filename[file])[[1]]
  gear = substr(Maxent_mod_summary_filename[file],dash_pos[[6]]+1,dash_pos[[7]]-1)
  mm = as.numeric(substr(Maxent_mod_summary_filename[file],dash_pos[[7]]+1,dash_pos[[7]]+2))
  
  # for import resample cross validation group
  Cutlassfish_resample_filename=list.files(
    paste0("./compiled_data/occ_data_in_ground/resample/",formatC(mm,width = 2,flag = 0),"/"),
    pattern = glob2rx(paste0(gear,"*")),
    full.names=T) 
  
  # load model
  load(file=Maxent_mod_summary_filename[file])
  
  #
  obs_occ = read.csv(Cutlassfish_both_occ_filename[mm],row.names = 1)
  Env_mth_data = Env[[mm-8]]
  bg_env = read.csv(random_bg_filename[mm-8],row.names = 1)
  
  obs_occ$Month = mm
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y&occ==1)
    Env_mth_yr = Env_mth_data %>% dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  mod_env = rbind(occ_env,bg_env)
  # prediction
  fg_map_pred_resample = NULL
  
  for(r in 1:20){
    Cutlassfish_resample = read.csv(Cutlassfish_resample_filename[r],row.names = 1)
    
    # arrange
    resample_nthfold = NULL
    for(y in 2009:2019){
      resample_occ_yr = Cutlassfish_resample[
        which(obs_occ$Year==y&
                obs_occ$occ==1),]
      resample_nthfold=c(resample_nthfold,resample_occ_yr$nth_fold)
    }
    resample_nthfold = c(resample_nthfold,bg_env$nth_fold)
    mod_env$resample_nthfold = resample_nthfold
    
    for(f in 1:5){
      maxent_mod_indiv = Maxent_model_resample[[5*(r-1)+f]]
      mod_env$pred_fg = dismo::predict(maxent_mod_indiv,mod_env)
      mod_env$resample = r
      mod_env$mod_fold = f
      fg_map_pred_resample = rbind(fg_map_pred_resample,mod_env)
    }
  }
  write.csv(fg_map_pred_resample,
            paste0("./output/Maxent_5CV/Maxent_pred_fg_for_AUC/resample/pred_bg10e4Bgbias_",
                   gear, "_fg_", formatC(mm,width = 2,flag = 0),".csv"))
}


############################################################################
# CALCULATE AUC ############

library(pROC)
library(ecospat)
library(dplyr)


pred_fg_filename = list.files("./output/Maxent_5CV/Maxent_pred_fg_for_AUC/resample",
                              pattern = glob2rx("*bg10e4Bgbias*.csv"),full.names = T)

for(fil in 1:length(pred_fg_filename)){
  AUC_value = NULL
  
  pred_fg = read.csv(pred_fg_filename[fil])
  
  dash_pos = gregexpr("_",pred_fg_filename[fil])[[1]]
  gear = substr(pred_fg_filename[fil],dash_pos[[7]]+1,dash_pos[[8]]-1)
  mm = substr(pred_fg_filename[fil],dash_pos[[9]]+1,dash_pos[[9]]+2)
  
  for(r in 1:20){
    for(mod_f in 1:5){
      pred_fg_train = pred_fg %>% dplyr::filter(resample==r&mod_fold==mod_f&resample_nthfold!=mod_f)
      pred_fg_test = pred_fg %>% dplyr::filter(resample==r&mod_fold==mod_f&resample_nthfold==mod_f)
      
      # Train AUC
      train_auc_result = roc(response = pred_fg_train$occ, 
                             predictor = pred_fg_train$pred_fg,
                             levels = c(0,1), direction='<')$auc
      
      # Test AUC
      test_auc_result = roc(response = pred_fg_test$occ, 
                            predictor = pred_fg_test$pred_fg,
                            levels = c(0,1), direction='<')$auc
      
      auc_result = c(train_auc_result,test_auc_result,r,mod_f)
      
      
      # save
      AUC_value = rbind(AUC_value,auc_result)
    }
  }
  AUC_value = as.data.frame(AUC_value)
  colnames(AUC_value)=c("train_AUC","test_AUC","resample","mod_fold")
  
  write.csv(AUC_value,paste0("./output/Model_performance/AUC_fg_randomBg/resample/AUC_bg10e4Bgbias_",gear,"_",mm,".csv"))
}


##############################################################################
# ORGANIZE MODEL PERFORMANCE AND DATA PROPERTY #########################

library(dplyr)
## 5foldCV allmap training AUC, testing AUC ,  
## and self-mod Threshold data #######

training_AUC_resample = NULL
testing_AUC_resample = NULL
Thresh_resample = NULL

Cutlassfish_resample_summary_filename = list.files("./output/Maxent_5CV/Maxent_result_summary/resample",
                                                   pattern = glob2rx("*bg10e4Bgbias*.csv"),full.names = T)

AUC_resample_filename = list.files("./output/Model_performance/AUC_fg_randomBg/resample",
                                           pattern = glob2rx("*bg10e4Bgbias*.csv"),full.names = T)


for(file in 1:length(Cutlassfish_resample_summary_filename)){
  
  
  dash_pos = gregexpr("_",Cutlassfish_resample_summary_filename[file])[[1]]
  gear = substr(Cutlassfish_resample_summary_filename[file],dash_pos[[6]]+1,dash_pos[[7]]-1)
  mm = as.numeric(substr(Cutlassfish_resample_summary_filename[file],dash_pos[[7]]+1,dash_pos[[7]]+2))
  
  AUC_fg_resample = read.csv(AUC_resample_filename[file],row.names = 1)
  training_AUC_resample = c(training_AUC_resample,AUC_fg_resample$train_AUC)
  testing_AUC_resample = c(testing_AUC_resample,AUC_fg_resample$test_AUC)
  
  Cutlassfish_resample_5CV_summary=read.csv(Cutlassfish_resample_summary_filename[file],row.names = 1)
  thresh_10P = unlist(Cutlassfish_resample_5CV_summary["X10.percentile.training.presence.Cloglog.threshold",])
  Thresh_resample = rbind(Thresh_resample,cbind(as.numeric(as.character(thresh_10P)),rep(gear,100),rep(mm,100),AUC_fg_resample$resample))
}

Thresh_resample = as.data.frame(Thresh_resample)
Thresh_resample[,c(1,3,4)] = apply(Thresh_resample[,c(1,3,4)],2,function(x)as.numeric(as.character(x)))
colnames(Thresh_resample) = c("thresh_10P","gear","month","resample")
Thresh_resample$gear = as.character(Thresh_resample$gear)


## 5foldCV thresh-dependent model performance ########
source(paste0(Root_dir,"Fishtactics_SDM_code_check/Function_ThreshDepend_ModPerf.R"))

# import predicted data in fishing ground
fg_pred_resample_filename = list.files("./output/Maxent_5CV/Maxent_pred_fg_for_AUC/resample",
                                       glob2rx("*bg10e4Bgbias*.csv"),full.names = T)

# import true occ data
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)

# calculate ModPerf
ModPerf_gear_resample_pre = NULL
for(file in 1:length(fg_pred_resample_filename)){
  
  dash_pos = gregexpr("_",fg_pred_resample_filename[file])[[1]]
  gear_name = substr(fg_pred_resample_filename[file],dash_pos[[7]]+1,dash_pos[[8]]-1)
  mm = as.numeric(substr(fg_pred_resample_filename[file],dash_pos[[9]]+1,dash_pos[[9]]+2))
  
  fg_pred_resample = read.csv(fg_pred_resample_filename[file],row.names = 1)
  
  for(r in 1:20){
    n_data = nrow(fg_pred_resample)/20
    fg_pred_resample_onesample = fg_pred_resample[(((n_data*(r-1))+1):(n_data*r)),]
    thresh_r = Thresh_resample%>%dplyr::filter(month==mm &gear==gear_name & resample==r)
    ## Mod Perf
    ModPerf_gear = 
      Thresh_depend_ModPerf_resample_func(mth_pred_data = fg_pred_resample_onesample,
                                          thresh = thresh_r[,1])
    ModPerf_gear_df = cbind(as.data.frame(ModPerf_gear),mm,r,c(1:5),gear_name)
    ModPerf_gear_resample_pre = rbind(ModPerf_gear_resample_pre,
                                      ModPerf_gear_df)
  }
}

## Variable contribution #####################################
Cutlassfish_resample_summary_filename = list.files("./output/Maxent_5CV/Maxent_result_summary/resample",pattern = glob2rx("*bg10e4Bgbias*.csv"),full.names = T)

Var_contrib_resample = NULL
for(file in 1:length(Cutlassfish_resample_summary_filename)){
  Cutlassfish_resample_5CV_summary=read.csv(Cutlassfish_resample_summary_filename[file],row.names = 1)
  
  Var_contrib_resample = rbind(Var_contrib_resample,
                               t(Cutlassfish_resample_5CV_summary[15:22,]))
}
colnames(Var_contrib_resample) = c("Bathy","EKE","TGM","MLD","SSH","SSS","SST","Chla")


## Merge all model performance #######################
Mod_perform_resample = cbind(training_AUC_resample,
                             testing_AUC_resample,
                             ModPerf_gear_resample_pre[,1:10],
                             Var_contrib_resample,
                             ModPerf_gear_resample_pre[,11:14])
colnames(Mod_perform_resample) = c(
  "training_AUC","testing_AUC",
  colnames(ModPerf_gear_resample_pre)[1:10],
  "Bathy","EKE","TGM","MLD","SSH","SSS","SST","Chla",
  colnames(ModPerf_gear_resample_pre)[11:14])

# Mod_perform_resample = cbind(Mod_perform_resample,OR_GEAR_resample[,6:11])
# Mod_perform_resample$Sample_area = as.character(Mod_perform_resample$Sample_area)
# Mod_perform_resample$Sample_size = as.character(Mod_perform_resample$Sample_size)


# Import data property and merge with model performance ######
data_property_filename = list.files("./compiled_data/data_property/",pattern = glob2rx("resample_*effort*.csv"),full.names = T)

data_property_resample_b = read.csv(data_property_filename[1],row.names = 1,stringsAsFactors = F)
data_property_resample_b$gear = "bycatch"

data_property_resample_t = read.csv(data_property_filename[2],row.names = 1,stringsAsFactors = F)
data_property_resample_t$gear = "target"

data_property_resample = rbind(data_property_resample_b,
                               data_property_resample_t)

#data_property_resample_ord = data_property_resample[order(data_property_resample$gear,data_property_resample$Month,data_property_resample$Sample_area,data_property_resample$N_lab),]
data_property_resample_copy <- as.data.frame(lapply(data_property_resample, rep, each=5))


# Merge with ModPerf
Mod_perform_resample = 
Mod_perform_data_property_all = cbind(Mod_perform_resample[,1:20],
                                      data_property_resample_copy[,c(2,4,6:9,1)])                              

write.csv(Mod_perform_data_property_all, "./output/Model_performance/ModPerf_DataProp_resample_5CV.csv")


#################################################################################
# COMPARE resampled BYCATCH AND TARGET 

modperf_dataprop = read.csv("./output/Model_performance/ModPerf_DataProp_resample_5CV.csv",row.names = 1)

sum_range=tapply(modperf_dataprop$Extent_smp_range,
       list(modperf_dataprop$Month,modperf_dataprop$gear),
       function(x){mu = mean(x);s = min(x);l=max(x);return(c(mu,s,l))})
sum_even = tapply(modperf_dataprop$kd_rad_auc,
       list(modperf_dataprop$Month,modperf_dataprop$gear),
       function(x){mu = mean(x);s = min(x);l=max(x);return(c(mu,s,l))})

# calculate 5CV mean
ModPerf_DatProp_5CV_mean = unique(modperf_dataprop[,c(21:27)])

mod_perform_5CVmean = NULL
mod_perform_5CVvar = NULL
for(i in 1:20){
  CV_mean= NULL
  CV_var= NULL
  j=1
  while(j<121){
    CV_mean = c(CV_mean, mean(modperf_dataprop[((5*(j-1)+1):(5*j)),i],na.rm=T))
    CV_var = c(CV_var, var(modperf_dataprop[((5*(j-1)+1):(5*j)),i],na.rm=T))
    j=j+1
  }
  mod_perform_5CVmean = cbind(mod_perform_5CVmean,CV_mean)
  mod_perform_5CVvar = cbind(mod_perform_5CVvar,CV_var)
}
ModPerf_DatProp_5CV_mean = cbind(mod_perform_5CVmean,ModPerf_DatProp_5CV_mean)
ModPerf_DatProp_5CV_var = cbind(mod_perform_5CVvar,ModPerf_DatProp_5CV_mean)
colnames(ModPerf_DatProp_5CV_mean)[1:20] = 
  c("train_AUC.mu","test_AUC.mu",
    "TSS.mu","F2.mu","OR.mu",
    "sensitivity.mu","specificity.mu","Precision.mu",
    "N00.mu","N01.mu","N10.mu","N11.mu",
    "Bathy.mu","EKE.mu","TGM.mu","MLD.mu",
    "SSH.mu","SSS.mu","SST.mu","Chla.mu")


sum_range=tapply(ModPerf_DatProp_5CV_mean$Extent_smp_range,
                 list(ModPerf_DatProp_5CV_mean$Month,ModPerf_DatProp_5CV_mean$gear),
                 function(x){mu = mean(x);s = min(x);l=max(x);return(c(mu,s,l))})
sum_even = tapply(ModPerf_DatProp_5CV_mean$kd_rad_auc,
                  list(ModPerf_DatProp_5CV_mean$Month,ModPerf_DatProp_5CV_mean$gear),
                  function(x){mu = mean(x);s = min(x);l=max(x);return(c(mu,s,l))})

# transform function
transform_perc <- function(percentage_vec) {
  # See Cribari-Neto & Zeileis (2010)
  (percentage_vec * (length(percentage_vec) - 1) + 0.5)/length(percentage_vec)
}

# Compare bycatch and target
library(mgcv)
library(itsadug)
library(ggplot2)
library(gridExtra)
library(scales) 
library(dplyr)

mod_perf_name = c("train_AUC.mu","test_AUC.mu",
                  "TSS.mu","F2.mu","OR.mu",
                  "sensitivity.mu","specificity.mu","Precision.mu")


def_color <- hue_pal()(3)                             # Identify hex codes
def_color 

beta_model_5CV = list()
beta_model_5CV_summary = NULL
BandT_smallrange_pic_range = list()
BandT_smallrange_pic_even = list()
BandT_smallrange_pic_gear = list()

for(m in 9:11){
  jpeg(paste0("./output/Gear_effect/BandT_GearEffect_",month.abb[m],".jpeg"),
       width = 6,height = 6,res=300,units = "in")
  par(mfrow=c(3,3))
  for(mp in 1:8){  
    
    ModPerf_DatProp_clear = ModPerf_DatProp_5CV_mean[,c(mod_perf_name[mp],
                                                        "gear",
                                                        "Month",
                                                        "Smp_size",
                                                        "Extent_smp_range",
                                                        "kd_rad_auc",
                                                        "No_site"
    )]
    colnames(ModPerf_DatProp_clear) = c("ModPerf","gear", "month","Smp_Size",
                                        "Extent_smp_range","kd_rad_auc","No_site")
    
    
    
    ModPerf_DatProp_mth=ModPerf_DatProp_clear %>% dplyr::filter(month==m)
    if(mp==3){
      ModPerf_DatProp_mth$ModPerf_trans = transform_perc((ModPerf_DatProp_mth$ModPerf/2)+0.5)
    }else(
      ModPerf_DatProp_mth$ModPerf_trans = transform_perc(ModPerf_DatProp_mth$ModPerf)
    )
    
    mod_beta = gam(ModPerf_trans~
                     gear+s(Extent_smp_range,bs="ts")+s(kd_rad_auc,bs="ts"),
                   select = T,
                   ModPerf_DatProp_mth,family=betar(link="logit"))
    beta_model_5CV[[8*(m-9)+mp]] = mod_beta
    mod_beta_sum = summary(mod_beta)
    beta_sig_info = paste(
      as.character(round(mod_beta_sum$p.pv[1],2)),
      as.character(round(mod_beta_sum$p.pv[2],2)),
      as.character(round(mod_beta_sum$s.pv[1],2)),
      as.character(round(mod_beta_sum$s.pv[2],2)),
      sep=" ")
    
    beta_mod_info = c(
      round(mod_beta_sum$p.coeff[2],2),
      round(mod_beta_sum$p.pv[2],2),
      round(mod_beta_sum$edf[1],2),
      round(mod_beta_sum$s.pv[1],2),
      round(mod_beta_sum$edf[2],2),
      round(mod_beta_sum$s.pv[2],2),
      round(mod_beta_sum$r.sq,2),
      round(mod_beta_sum$dev.expl,2),
      round(mod_beta_sum$sp.criterion,2),
      m,mod_perf_name[mp]
      )
    beta_model_5CV_summary = rbind(beta_model_5CV_summary,
                                   beta_mod_info)
    
    
    print(beta_sig_info)
    ilink <- family(mod_beta)$linkinv
    newdat_range_1 =
      data.frame(
        Extent_smp_range=
          seq(min(ModPerf_DatProp_mth$Extent_smp_range[which(ModPerf_DatProp_mth$gear=="bycatch")]),
              max(ModPerf_DatProp_mth$Extent_smp_range[which(ModPerf_DatProp_mth$gear=="bycatch")]),
              length.out=10),
        kd_rad_auc = 
          rep(mean(ModPerf_DatProp_mth$kd_rad_auc),10))
    newdat_range_2 =
      data.frame(
        Extent_smp_range=
          seq(min(ModPerf_DatProp_mth$Extent_smp_range[which(ModPerf_DatProp_mth$gear=="target")]),
              max(ModPerf_DatProp_mth$Extent_smp_range[which(ModPerf_DatProp_mth$gear=="target")]),
              length.out=10),
        kd_rad_auc = 
          rep(mean(ModPerf_DatProp_mth$kd_rad_auc),10))
    newdat_even_1 =
      data.frame(
        Extent_smp_range=
          rep(mean(ModPerf_DatProp_mth$Extent_smp_range),10),
        kd_rad_auc = 
          c(seq(min(ModPerf_DatProp_mth$kd_rad_auc[which(ModPerf_DatProp_mth$gear=="bycatch")]),
                max(ModPerf_DatProp_mth$kd_rad_auc[which(ModPerf_DatProp_mth$gear=="bycatch")]),
                length.out=10)))
    newdat_even_2 =
      data.frame(
        Extent_smp_range=
          rep(mean(ModPerf_DatProp_mth$Extent_smp_range),10),
        kd_rad_auc = 
          seq(min(ModPerf_DatProp_mth$kd_rad_auc[which(ModPerf_DatProp_mth$gear=="target")]),
              max(ModPerf_DatProp_mth$kd_rad_auc[which(ModPerf_DatProp_mth$gear=="target")]),
              length.out=10))
    newdat_range = rbind(newdat_range_1,newdat_range_2)
    newdat_even =  rbind(newdat_even_1,newdat_even_2)
    newdat_range$gear = rep(c("bycatch","target"),each=10)
    newdat_even$gear = rep(c("bycatch","target"),each=10)
    
    Preddata_range <- bind_cols(newdat_range, setNames(as_tibble(predict(mod_beta, newdat_range, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))
    Preddata_even <- bind_cols(newdat_even, setNames(as_tibble(predict(mod_beta, newdat_even, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))
    
    ## create the interval and backtransform
    if(mp==3){
      Preddata_range <- mutate(Preddata_range,
                               fit_resp  = (2*ilink(fit_link)-1),
                               right_upr = (2*ilink(fit_link + (2 * se_link))-1),
                               right_lwr = (2*ilink(fit_link - (2 * se_link))-1))
      Preddata_even <- mutate(Preddata_even,
                              fit_resp  = (2*ilink(fit_link)-1),
                              right_upr = (2*ilink(fit_link + (2 * se_link))-1),
                              right_lwr = (2*ilink(fit_link - (2 * se_link))-1))
    }else{
      Preddata_range <- mutate(Preddata_range,
                               fit_resp  = ilink(fit_link),
                               right_upr = ilink(fit_link + (2 * se_link)),
                               right_lwr = ilink(fit_link - (2 * se_link)))
      Preddata_even <- mutate(Preddata_even,
                              fit_resp  = ilink(fit_link),
                              right_upr = ilink(fit_link + (2 * se_link)),
                              right_lwr = ilink(fit_link - (2 * se_link)))
    }
    
    # plot
    plot_parametric(mod_beta,pred=list(gear=c("target","bycatch")),parametricOnly =T)
    
    if(mp==3){ylowlimit=-1;yupplimit=1}else if(mp==5){
      ylowlimit=1;yupplimit=0}else{ylowlimit=0;yupplimit=1}
    
    p_range= ggplot()+
      geom_ribbon(data =Preddata_range, 
                  aes(ymin=right_lwr, ymax=right_upr, x=Extent_smp_range, fill = gear), alpha = 0.3,colour=NA)+
      
      geom_line(data=Preddata_range,mapping = aes(x=Extent_smp_range,y=fit_resp,colour = gear),lwd=1.5)+
      geom_point(data=ModPerf_DatProp_mth,aes(x=Extent_smp_range,
                                              y=ModPerf,
                                              colour=gear),alpha=0.4,cex=1.5)+
      # scale_color_brewer(palette="Dark2")+
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_manual(breaks = c("bycatch", "target"), 
                        values=def_color[c(2,3)])+
      scale_color_manual(breaks = c("bycatch", "target"), 
                        values=def_color[c(2,3)])+
      theme_bw()+ylab(mod_perf_name[mp])+
      ggtitle(paste0("p val = ",beta_sig_info))+
      xlim(min(ModPerf_DatProp_mth$Extent_smp_range)-0.5,max(ModPerf_DatProp_mth$Extent_smp_range)+0.5)+
      ylim(ylowlimit,yupplimit)+
      theme(legend.position = "none",
            plot.title = element_text(size = 8))
    # 
    p_even= ggplot()+
      geom_ribbon(data =Preddata_even, aes(ymin=right_lwr, ymax=right_upr, x=kd_rad_auc, fill = gear), alpha = 0.3,colour=NA)+
      
      geom_line(data=Preddata_even,mapping = aes(x=kd_rad_auc,y=fit_resp,colour = gear),lwd=1.5)+
      geom_point(data=ModPerf_DatProp_mth,aes(x=kd_rad_auc,
                                              y=ModPerf,
                                              colour=gear),alpha=0.4,cex=1.5)+
      # scale_color_brewer(palette="Dark2")+
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_manual(breaks = c("bycatch", "target"), 
                        values=def_color[c(2,3)])+
      scale_color_manual(breaks = c("bycatch", "target"), 
                         values=def_color[c(2,3)])+
      theme_bw()+ylab(mod_perf_name[mp])+
      ggtitle(paste0("p val = ",beta_sig_info))+
      xlim(min(ModPerf_DatProp_mth$kd_rad_auc)-0.5,max(ModPerf_DatProp_mth$kd_rad_auc)+0.5)+
      ylim(ylowlimit,yupplimit)+
      theme(legend.position = "none",
            plot.title = element_text(size = 8))

    BandT_smallrange_pic_range[[mp]]=p_range    #for Overlapped
    BandT_smallrange_pic_even[[mp]]=p_even     #for Overlapped
  }
  dev.off()
  
  jpeg(paste0("./output/Gear_effect/BandT_RangeEffect_",month.abb[m],".jpeg"),
       width = 6,height = 6,res=300,units = "in")   #for Overlapped
  grid.arrange(grobs=BandT_smallrange_pic_range,ncol=3)
  dev.off()
  
  jpeg(paste0("./output/Gear_effect/BandT_EvenEffect_",month.abb[m],".jpeg"),
       width = 6,height = 6,res=300,units = "in")   #for Overlapped
  grid.arrange(grobs=BandT_smallrange_pic_even,ncol=3)
  dev.off()
}

beta_model_5CV_summary = as.data.frame(beta_model_5CV_summary)
colnames(beta_model_5CV_summary) = c("p_estimate","p_pval","range_edf","range_pval",
                                     "kdrad_edf","kdrad_pval","R_sq","dev_expl",
                                     "neg_REML","month","ModPerf")

write.csv(beta_model_5CV_summary,"./output/Gear_effect/gear_effect_summary.csv")
save(beta_model_5CV,file="./output/Gear_effect/BandT_gear_model.RData")
load(file="./output/Gear_effect/BandT_gear_model.RData")
summary(beta_model_5CV[[1]])
