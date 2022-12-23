wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library(dismo)
library(dplyr)
library("rJava")
library(ggplot2)
library(gridExtra)
library(viridis)

Sys.setenv(JAVA_Home="./Maxent_Java/jdk-11.0.2")

# import environment data
Env_data = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv",stringsAsFactors = F)
Env = list()
for(m in 1:3){
  Env_mth = Env_data%>%filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}

# import true occ data name
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)

# create random bg for validation (only for first time, already create) ####
# for(m in 9:11){
#   obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
#   Env_mth_data = Env[[m-8]]
#   obs_occ$Month = m
# 
#   bg_env = NULL
#   for(y in 2009:2019){
#     occ_yr = obs_occ %>%filter(Year==y)
#     sample_random_bg = sample(1:nrow(occ_yr),1000)
#     random_bg_yr = occ_yr[sample_random_bg,]
# 
#     Env_mth_yr = Env_mth_data %>%filter(Year==y)
# 
#     bg_env_yr = merge(x=random_bg_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
#     bg_env=rbind(bg_env,bg_env_yr)
#   }
# 
#   bg_env$nth_fold = sample(c(1:5),nrow(bg_env),replace = T)
#   write.csv(bg_env,
#             paste0("./compiled_data/bg_file/random_for_validation/bg10e4_fg_",formatC(m,width = 2,flag = 0),".csv"))
# 
# }

# prediction for fishing ground -------------------------------------------------------
random_bg_filename=list.files("./compiled_data/bg_file/random_for_validation",full.names=T)


## both ####

fg_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_both/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Env_mth_data = Env[[m-8]]
  bg_env = read.csv(random_bg_filename[m-8],row.names = 1)
  
  obs_occ$Month = m
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y&occ==1)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  mod_env = rbind(occ_env,bg_env)
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    mod_env$pred_fg = dismo::predict(maxent_mod_indiv,mod_env)
    mod_env$mod_fold = f
    fg_map_pred = rbind(fg_map_pred,mod_env)
  }
}
write.csv(fg_map_pred,"./output/Maxent_5CV/Maxent_pred_fg_for_AUC/grid_both/pred_bg10e4Bgbias_both_fg.csv")


## byactch ####

fg_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_bycatch/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Env_mth_data = Env[[m-8]]
  bg_env = read.csv(random_bg_filename[m-8],row.names = 1)
  
  obs_occ$Month = m
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y&occ==1)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  mod_env = rbind(occ_env,bg_env)
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    mod_env$pred_fg = dismo::predict(maxent_mod_indiv,mod_env)
    mod_env$mod_fold = f
    fg_map_pred = rbind(fg_map_pred,mod_env)
  }
}
write.csv(fg_map_pred,"./output/Maxent_5CV/Maxent_pred_fg_for_AUC/grid_bycatch/pred_bg10e4Bgbias_bycatch_fg.csv")


## target ####

fg_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_target/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Env_mth_data = Env[[m-8]]
  bg_env = read.csv(random_bg_filename[m-8],row.names = 1)
  
  obs_occ$Month = m
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y&occ==1)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  mod_env = rbind(occ_env,bg_env)
  
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    mod_env$pred_fg = dismo::predict(maxent_mod_indiv,mod_env)
    mod_env$mod_fold = f
    fg_map_pred = rbind(fg_map_pred,mod_env)
  }
}
write.csv(fg_map_pred,"./output/Maxent_5CV/Maxent_pred_fg_for_AUC/grid_target/pred_bg10e4Bgbias_target_fg.csv")

