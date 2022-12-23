
####################################
# This script includes 2 parts
# (1) examine sample size, sampling range and sampling kernel density of original data
# (2) resample data - control sample size
# (3) examine the extent of sampling range and 
#     sampling kernel density of resampled data
####################################


wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

Root_dir = "D:/Ruo_data/2020_Assis_FishTactic/"

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
Cutlassfish_both_filename = list.files("./compiled_data/occ_data_in_ground/grid_both/",full.names=T)
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

# number of data (1797 for each YearMth)
n_data = nrow(Cutlassfish_bycatch_occ)

# -------------------------------------------------------------------------
# Examine the data property of the original data 
# -------------------------------------------------------------------------

kd_data_occ_both = NULL
kd_data_effort_both = NULL
data_property_occ_both = NULL
data_property_effort_both = NULL

## both gear ######
for(m in 9:11){
  Cutlassfish_both_occ = read.csv(Cutlassfish_both_filename[m],row.names = 1)
  effort_smpsize = length(which(Cutlassfish_both_occ$fishing_loc==1))
  sz = NA
  #-----------------------------------------------------------------------------
  # for occ
  # calculate sampling range
  unique_lon_lat = unique(as.matrix(Cutlassfish_both_occ[which(Cutlassfish_both_occ$occ==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat)
  hull_b <- c(hull_b, hull_b[1])
  sampling_range =
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat[hull_b,1],
                                          unique_lon_lat[hull_b,2])!=0),]
  extent_sampling_range = 100*(nrow(sampling_range)/nrow(grid_in_ground))
  
  # - number of site
  n_site = nrow(unique_lon_lat)
  
  # calculate the aggregation of sampling probability
  # -- CV of kernel density
  kd_resample = kernel_density_func(Cutlassfish_both_occ[which(Cutlassfish_both_occ$occ==1),],
                                    paste0("occ_kernel_density/Original_data/tif/both_",formatC(m,width = 2,flag = 0)),
                                    m,sz,NA)
  kd_range =
    kd_resample[which(point.in.polygon(kd_resample$Lon,
                                       kd_resample$Lat,
                                       unique_lon_lat[hull_b,1],
                                       unique_lon_lat[hull_b,2])!=0),]
  kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
  
  # -- distribution of kernel density
  kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100,
                       kd_val = sort(kd_range$kd/max(kd_range$kd)))
  kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
  
  #save
  kd_data_occ_both = rbind(kd_data_occ_both,kd_range)
  data_property_occ_both = rbind(data_property_occ_both,
                                 c(m, effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site))
  
  #-----------------------------------------------------------------------------
  # for effort
  # calculate sampling range
  unique_lon_lat = unique(as.matrix(Cutlassfish_both_occ[which(Cutlassfish_both_occ$fishing_loc==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat)
  hull_b <- c(hull_b, hull_b[1])
  sampling_range =
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat[hull_b,1],
                                          unique_lon_lat[hull_b,2])!=0),]
  extent_sampling_range = 100*(nrow(sampling_range)/nrow(grid_in_ground))
  
  # - number of site
  n_site = nrow(unique_lon_lat)
  
  # calculate the aggregation of sampling probability
  # -- CV of kernel density
  kd_resample = kernel_density_func(Cutlassfish_both_occ[which(Cutlassfish_both_occ$fishing_loc==1),],
                                    paste0("effort_kernel_density/Original_data/tif/both_",formatC(m,width = 2,flag = 0)),
                                    m,sz,NA)
  kd_range =
    kd_resample[which(point.in.polygon(kd_resample$Lon,
                                       kd_resample$Lat,
                                       unique_lon_lat[hull_b,1],
                                       unique_lon_lat[hull_b,2])!=0),]
  kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
  
  # -- distribution of kernel density
  kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100,
                       kd_val = sort(kd_range$kd/max(kd_range$kd)))
  kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
  
  #save
  kd_data_effort_both = rbind(kd_data_effort_both,kd_range)
  data_property_effort_both = rbind(data_property_effort_both,
                                    c(m, effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site))
  
}

write.csv(kd_data_occ_both,"./compiled_data/occ_kernel_density/Original_data/kd_data_both.csv")
write.csv(kd_data_effort_both,"./compiled_data/effort_kernel_density/Original_data/kd_data_both.csv")
write.csv(data_property_occ_both,"./compiled_data/data_property/Original_data/both_occ.csv")
write.csv(data_property_effort_both,"./compiled_data/data_property/Original_data/both_effort.csv")





## bycatch ######

kd_data_occ_bycatch = NULL
kd_data_effort_bycatch = NULL

data_property_occ_bycatch = NULL
data_property_effort_bycatch = NULL

for(m in 9:11){
  Cutlassfish_bycatch_occ = read.csv(Cutlassfish_bycatch_filename[m],row.names = 1)
  #N_lab = length(resample_size_bycatch[[m-8]])
  effort_smpsize = length(which(Cutlassfish_bycatch_occ$fishing_loc==1))
  
  #-----------------------------------------------------------------------------
  # for occ
  # calculate sampling range
  unique_lon_lat = unique(as.matrix(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$occ==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat)
  hull_b <- c(hull_b, hull_b[1])
  sampling_range =
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat[hull_b,1],
                                          unique_lon_lat[hull_b,2])!=0),]
  extent_sampling_range = 100*(nrow(sampling_range)/nrow(grid_in_ground))

  # - number of site
  n_site = nrow(unique_lon_lat)

  # calculate the aggregation of sampling probability
  # -- CV of kernel density
  kd_resample = kernel_density_func(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$occ==1),],
                                    paste0("occ_kernel_density/Original_data/tif/bycatch_",formatC(m,width = 2,flag = 0)),
                                    m,NA,NA)
  kd_range =
    kd_resample[which(point.in.polygon(kd_resample$Lon,
                                       kd_resample$Lat,
                                       unique_lon_lat[hull_b,1],
                                       unique_lon_lat[hull_b,2])!=0),]
  kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
  
  # -- distribution of kernel density
  kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100,
                       kd_val = sort(kd_range$kd/max(kd_range$kd)))
  kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)

  #save
  kd_data_occ_bycatch = rbind(kd_data_occ_bycatch,kd_range)
  data_property_occ_bycatch = rbind(data_property_occ_bycatch,
                                c(m,effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site))

  #-----------------------------------------------------------------------------
  # for effort
  # calculate sampling range
  unique_lon_lat = unique(as.matrix(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$fishing_loc==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat)
  hull_b <- c(hull_b, hull_b[1])
  sampling_range =
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat[hull_b,1],
                                          unique_lon_lat[hull_b,2])!=0),]
  extent_sampling_range = 100*(nrow(sampling_range)/nrow(grid_in_ground))
  
  # - number of site
  n_site = nrow(unique_lon_lat)
  
  # calculate the aggregation of sampling probability
  # -- CV of kernel density
  kd_resample = kernel_density_func(Cutlassfish_bycatch_occ[which(Cutlassfish_bycatch_occ$fishing_loc==1),],
                                    paste0("effort_kernel_density/Original_data/tif/bycatch_",formatC(m,width = 2,flag = 0)),
                                    m,NA,NA)
  kd_range =
    kd_resample[which(point.in.polygon(kd_resample$Lon,
                                       kd_resample$Lat,
                                       unique_lon_lat[hull_b,1],
                                       unique_lon_lat[hull_b,2])!=0),]
  kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
  
  # -- distribution of kernel density
  kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100,
                       kd_val = sort(kd_range$kd/max(kd_range$kd)))
  kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
  
  #save
  kd_data_effort_bycatch = rbind(kd_data_effort_bycatch,kd_range)
  data_property_effort_bycatch = rbind(data_property_effort_bycatch,
                                c(m,effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site))
  
}

write.csv(kd_data_occ_bycatch,"./compiled_data/occ_kernel_density/Original_data/kd_data_bycatch.csv")
write.csv(kd_data_effort_bycatch,"./compiled_data/effort_kernel_density/Original_data/kd_data_bycatch.csv")
write.csv(data_property_occ_bycatch,"./compiled_data/data_property/Original_data/bycatch_occ.csv")
write.csv(data_property_effort_bycatch,"./compiled_data/data_property/Original_data/bycatch_effort.csv")


## target ######

kd_data_occ_target = NULL
kd_data_effort_target = NULL
data_property_occ_target = NULL
data_property_effort_target = NULL


for(m in 9:11){
  Cutlassfish_target_occ = read.csv(Cutlassfish_target_filename[m-3],row.names = 1)
  #N_lab = length(resample_size_target[[m-8]])
  effort_smpsize = length(which(Cutlassfish_target_occ$fishing_loc==1))
  
  #-----------------------------------------------------------------------------
  # for occ
  # calculate sampling range
  unique_lon_lat = unique(as.matrix(Cutlassfish_target_occ[which(Cutlassfish_target_occ$occ==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat)
  hull_b <- c(hull_b, hull_b[1])
  sampling_range =
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat[hull_b,1],
                                          unique_lon_lat[hull_b,2])!=0),]
  extent_sampling_range = 100*(nrow(sampling_range)/nrow(grid_in_ground))

  # - number of site
  n_site = nrow(unique_lon_lat)

  # calculate the aggregation of sampling probability
  # -- CV of kernel density
  kd_resample = kernel_density_func(Cutlassfish_target_occ[which(Cutlassfish_target_occ$occ==1),],
                                    paste0("occ_kernel_density/Original_data/tif/target_",formatC(m,width = 2,flag = 0)),
                                    m,NA,NA)
  kd_range =
    kd_resample[which(point.in.polygon(kd_resample$Lon,
                                       kd_resample$Lat,
                                       unique_lon_lat[hull_b,1],
                                       unique_lon_lat[hull_b,2])!=0),]
  kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
  
  # -- distribution of kernel density
  kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100,
                       kd_val = sort(kd_range$kd/max(kd_range$kd)))
  kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)

  #save
  kd_data_occ_target = rbind(kd_data_occ_target,kd_range)
  data_property_occ_target = rbind(data_property_occ_target,
                                   c(m,effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site))
  
  #-----------------------------------------------------------------------------
  # for effort
  # calculate sampling range
  unique_lon_lat = unique(as.matrix(Cutlassfish_target_occ[which(Cutlassfish_target_occ$fishing_loc==1),c("Lon","Lat")]))
  hull_b <- chull(unique_lon_lat)
  hull_b <- c(hull_b, hull_b[1])
  sampling_range =
    grid_in_ground[which(point.in.polygon(grid_in_ground$Lon,
                                          grid_in_ground$Lat,
                                          unique_lon_lat[hull_b,1],
                                          unique_lon_lat[hull_b,2])!=0),]
  extent_sampling_range = 100*(nrow(sampling_range)/nrow(grid_in_ground))
  
  # - number of site
  n_site = nrow(unique_lon_lat)
  
  # calculate the aggregation of sampling probability
  # -- CV of kernel density
  kd_resample = kernel_density_func(Cutlassfish_target_occ[which(Cutlassfish_target_occ$fishing_loc==1),],
                                    paste0("effort_kernel_density/Original_data/tif/target_",formatC(m,width = 2,flag = 0)),
                                    m,NA,NA)
  kd_range =
    kd_resample[which(point.in.polygon(kd_resample$Lon,
                                       kd_resample$Lat,
                                       unique_lon_lat[hull_b,1],
                                       unique_lon_lat[hull_b,2])!=0),]
  kd_cv = sd(kd_range$kd)* 100 / mean(kd_range$kd)
  
  # -- distribution of kernel density
  kd_rad <- data.frame(std_rank = (1:nrow(kd_range)/nrow(kd_range))*100,
                       kd_val = sort(kd_range$kd/max(kd_range$kd)))
  kd_rad_auc = trapz(kd_rad$std_rank, kd_rad$kd_val)
  
  #save
  kd_data_effort_target = rbind(kd_data_effort_target,kd_range)
  data_property_effort_target = rbind(data_property_effort_target,
                                   c(m, effort_smpsize,extent_sampling_range,kd_cv,kd_rad_auc,n_site))
}

## save #####
write.csv(kd_data_occ_target,"./compiled_data/occ_kernel_density/Original_data/kd_data_target.csv")
write.csv(kd_data_effort_target,"./compiled_data/effort_kernel_density/Original_data/kd_data_target.csv")
write.csv(data_property_occ_target,"./compiled_data/data_property/Original_data/target_occ.csv")
write.csv(data_property_effort_target,"./compiled_data/data_property/Original_data/target_effort.csv")




