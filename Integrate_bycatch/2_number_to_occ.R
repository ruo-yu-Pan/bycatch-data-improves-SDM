
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

Root_dir = "D:/Ruo_data/2020_Assis_FishTactic/"

library(dplyr)

grid_in_ground = read.csv("./compiled_data/map/grid_in_ground.csv",row.names = 1,stringsAsFactors = F)
source(paste0(Root_dir,"Fishtactics_SDM_code_check/Function_CV_fold_generator.R"))

# Transform number to presence-absence ###################################
## total data ##########################################################

Cutlassfish_tot_in_polygon = read.csv("./compiled_data/fishery_data_in_ground/Cutlassfish_total.csv", row.names = 1,stringsAsFactors = F)
Cutlassfish_tot_in_polygon$No_Cutlassfish[which(is.na(Cutlassfish_tot_in_polygon$No_Cutlassfish))]=0
Cutlassfish_tot_in_polygon$No_Cutlassfish = as.numeric(Cutlassfish_tot_in_polygon$No_Cutlassfish)
#Cutlassfish_tot_in_polygon$fishing_loc = 1

for(m in 1:12){
  grid_in_ground_occ = NULL
  for(y in 2009:2019){
    Cutlassfish_tot_yr_mth = Cutlassfish_tot_in_polygon %>% filter(Year ==y,Month==m)
    Cutlassfish_tot_yr_mth_occ = tapply(Cutlassfish_tot_yr_mth$No_Cutlassfish,
                                       Cutlassfish_tot_yr_mth$grid_lab,
                                       sum,na.rm = T)
    Cutlassfish_tot_yr_mth_occ = as.data.frame(ifelse(Cutlassfish_tot_yr_mth_occ>0,1,0))
    colnames(Cutlassfish_tot_yr_mth_occ)[1] = "occ" 
    Cutlassfish_tot_yr_mth_occ$grid_lab = rownames(Cutlassfish_tot_yr_mth_occ)
    if(nrow(Cutlassfish_tot_yr_mth_occ)>0){
      Cutlassfish_tot_yr_mth_occ$fishing_loc = 1}
    
    grid_in_ground_yr_mth_occ = merge(x=grid_in_ground,y=Cutlassfish_tot_yr_mth_occ,by="grid_lab",all.x=T)
    grid_in_ground_yr_mth_occ$occ[which(is.na(grid_in_ground_yr_mth_occ$occ))]=0
    if(nrow(Cutlassfish_tot_yr_mth_occ)>0){
      grid_in_ground_yr_mth_occ$fishing_loc[which(is.na(grid_in_ground_yr_mth_occ$fishing_loc))]=0
    }else{grid_in_ground_yr_mth_occ$fishing_loc=0
    }
    
    grid_in_ground_yr_mth_occ$Year = y
    
    grid_in_ground_occ = rbind(grid_in_ground_occ,grid_in_ground_yr_mth_occ)
  }
  grid_in_ground_occ_CVfold = generate_CV_fold(grid_in_ground_occ,nfold = 5)
  print(length(which(grid_in_ground_occ$occ==1)))
  write.csv(grid_in_ground_occ_CVfold,paste0("./compiled_data/occ_data_in_ground/grid_both/occ_",ifelse(m %/%10 ==0,paste0("0",m),m),".csv"))
}

## bycatch data ##########################################################
Cutlassfish_bycatch_in_polygon = Cutlassfish_tot_in_polygon %>%dplyr::filter(Fishing_method_label==1)

for(m in 1:12){
  grid_in_ground_occ = NULL
  for(y in 2009:2019){
    Cutlassfish_bycatch_yr_mth = Cutlassfish_bycatch_in_polygon %>% dplyr::filter(Year ==y,Month==m)
    Cutlassfish_bycatch_yr_mth_occ = tapply(Cutlassfish_bycatch_yr_mth$No_Cutlassfish,
                                        Cutlassfish_bycatch_yr_mth$grid_lab,
                                        sum,na.rm = T)
    Cutlassfish_bycatch_yr_mth_occ = as.data.frame(ifelse(Cutlassfish_bycatch_yr_mth_occ>0,1,0))
    colnames(Cutlassfish_bycatch_yr_mth_occ)[1] = "occ" 
    Cutlassfish_bycatch_yr_mth_occ$grid_lab = rownames(Cutlassfish_bycatch_yr_mth_occ)
    if(nrow(Cutlassfish_bycatch_yr_mth_occ)>0){
      Cutlassfish_bycatch_yr_mth_occ$fishing_loc = 1}
    
    grid_in_ground_yr_mth_occ = merge(x=grid_in_ground,y=Cutlassfish_bycatch_yr_mth_occ,by="grid_lab",all.x=T)
    grid_in_ground_yr_mth_occ$occ[which(is.na(grid_in_ground_yr_mth_occ$occ))]=0
    if(nrow(Cutlassfish_bycatch_yr_mth_occ)>0){
      grid_in_ground_yr_mth_occ$fishing_loc[which(is.na(grid_in_ground_yr_mth_occ$fishing_loc))]=0
    }else{grid_in_ground_yr_mth_occ$fishing_loc=0
    }
    grid_in_ground_yr_mth_occ$Year = y
    
    grid_in_ground_occ = rbind(grid_in_ground_occ,grid_in_ground_yr_mth_occ)
  }
  grid_in_ground_occ_CVfold = generate_CV_fold(grid_in_ground_occ,nfold = 5)
  print(length(which(grid_in_ground_occ$occ==1)))
  write.csv(grid_in_ground_occ_CVfold,paste0("./compiled_data/occ_data_in_ground/grid_bycatch/occ_",ifelse(m %/%10 ==0,paste0("0",m),m),".csv"))
}

## target data ##########################################################
Cutlassfish_target_in_polygon = Cutlassfish_tot_in_polygon %>%dplyr::filter(Fishing_method_label==3)

for(m in c(1:4,8:12)){ # May have only 1 data and June to July have no data
  grid_in_ground_occ = NULL
  for(y in 2009:2019){
    Cutlassfish_target_yr_mth = Cutlassfish_target_in_polygon %>% dplyr::filter(Year ==y,Month==m)
    Cutlassfish_target_yr_mth_occ = tapply(Cutlassfish_target_yr_mth$No_Cutlassfish,
                                          Cutlassfish_target_yr_mth$grid_lab,
                                        sum,na.rm = T)
    Cutlassfish_target_yr_mth_occ = as.data.frame(ifelse(Cutlassfish_target_yr_mth_occ>0,1,0))
    colnames(Cutlassfish_target_yr_mth_occ)[1] = "occ" 
    Cutlassfish_target_yr_mth_occ$grid_lab = rownames(Cutlassfish_target_yr_mth_occ)
    if(nrow(Cutlassfish_target_yr_mth_occ)>0){
      Cutlassfish_target_yr_mth_occ$fishing_loc = 1}
    
    grid_in_ground_yr_mth_occ = merge(x=grid_in_ground,y=Cutlassfish_target_yr_mth_occ,by="grid_lab",all.x=T)
    grid_in_ground_yr_mth_occ$occ[which(is.na(grid_in_ground_yr_mth_occ$occ))]=0
    if(nrow(Cutlassfish_target_yr_mth_occ)>0){
      grid_in_ground_yr_mth_occ$fishing_loc[which(is.na(grid_in_ground_yr_mth_occ$fishing_loc))]=0
    }else{grid_in_ground_yr_mth_occ$fishing_loc=0
    }
    grid_in_ground_yr_mth_occ$Year = y
    
    grid_in_ground_occ = rbind(grid_in_ground_occ,grid_in_ground_yr_mth_occ)
  }
  grid_in_ground_occ_CVfold = generate_CV_fold(grid_in_ground_occ,nfold = 5)
  print(length(which(grid_in_ground_occ$occ==1)))
  write.csv(grid_in_ground_occ_CVfold,paste0("./compiled_data/occ_data_in_ground/grid_target/occ_",ifelse(m %/%10 ==0,paste0("0",m),m),".csv"))
}

