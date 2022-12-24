
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library(dismo)
library(dplyr)
library("rJava")
library(ggplot2)
library(gridExtra)
library(viridis)

Sys.setenv(JAVA_Home="./Maxent_Java/jdk-11.0.2")

library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
#load("./compiled_data/resample_size.RData")

# import environment data
Env_data = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv",stringsAsFactors = F)
Env = list()
for(m in 1:3){
  Env_mth = Env_data%>%filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}

# function for MaxEnt model training and prediction
Maxent_train_pred_func = function(occ_data=Cutlassfish_occ,
                                  Env_mth_data=Env_mth,
                                  bg_env_data=bg_env_mth,
                                  save_dir){
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = occ_data %>%filter(occ==1&Year==y)
    Env_mth_yr = Env_mth_data %>%filter(Year==y)
    
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

# Train Maxent model and predict ------------------------------------------------
# NOTE! This will take long time to compute and create a lot of files

Train_mod = F

## both ----------------------------------------------------------------------

if(Train_mod==T){
  Cutlassfish_both_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)
  bg_env_both_filename=list.files("./compiled_data/bg_file/grid_both",
                                  pattern =glob2rx("bg10e4Bgbias_*.csv") ,full.names=T)
  
  for (m in 9:11){
    
    save_dir = paste0("./output/Maxent_5CV/Maxent_result_detail/grid_both/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias")
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    Maxent_output = 
      Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_both_filename[m],row.names = 1,stringsAsFactors = F),
                             Env_mth_data=Env[[m-8]],
                             bg_env_data=read.csv(bg_env_both_filename[m-8],row.names = 1,stringsAsFactors = F),
                             save_dir=save_dir)
    
    Maxent_mod = Maxent_output[[2]]
    write.csv(Maxent_output[[1]],paste0("./output/Maxent_5CV/Maxent_result_summary/grid_both/Maxent_Summary_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".csv"))
    save(Maxent_mod,file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_both/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
    write.csv(Maxent_output[[3]],paste0("./output/Maxent_5CV/Maxent_pred/grid_both/pred_bg10e4Bgbias_",formatC(m,width = 2,flag = 0),".csv"))
  }
}



## bycatch ----------------------------------------------------------------------
### - origin

if(Train_mod==T){
  Cutlassfish_bycatch_filename=list.files("./compiled_data/occ_data_in_ground/grid_bycatch",full.names=T)
  bg_env_bycatch_filename=list.files("./compiled_data/bg_file/grid_bycatch",
                                     pattern =glob2rx("bg10e4Bgbias_*.csv") ,full.names=T)
  
  
  for (m in 9:11){
    
    save_dir = paste0("./output/Maxent_5CV/Maxent_result_detail/grid_bycatch/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias")
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    Maxent_output = 
      Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_bycatch_filename[m],row.names = 1,stringsAsFactors = F),
                             Env_mth_data=Env[[m-8]],
                             bg_env_data=read.csv(bg_env_bycatch_filename[m-8],row.names = 1,stringsAsFactors = F),
                             save_dir=save_dir)
    
    Maxent_mod = Maxent_output[[2]]
    write.csv(Maxent_output[[1]],paste0("./output/Maxent_5CV/Maxent_result_summary/grid_bycatch/Maxent_Summary_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".csv"))
    save(Maxent_mod,file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_bycatch/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
    write.csv(Maxent_output[[3]],paste0("./output/Maxent_5CV/Maxent_pred/grid_bycatch/pred_bg10e4Bgbias_",formatC(m,width = 2,flag = 0),".csv"))
  }
}


## target ----------------------------------------------------------------------
### - origin

if(Train_mod==T){
  Cutlassfish_target_filename=list.files("./compiled_data/occ_data_in_ground/grid_target",full.names=T)
  bg_env_target_filename=list.files("./compiled_data/bg_file/grid_target",
                                    pattern =glob2rx("bg10e4Bgbias_*.csv") ,full.names=T)
  
  
  for (m in 9:11){
    
    save_dir = paste0("./output/Maxent_5CV/Maxent_result_detail/grid_target/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias")
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    Maxent_output = 
      Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_target_filename[m-3],row.names = 1,stringsAsFactors = F),
                             Env_mth_data=Env[[m-8]],
                             bg_env_data=read.csv(bg_env_target_filename[m-8],row.names = 1,stringsAsFactors = F),
                             save_dir=save_dir)
    
    Maxent_mod = Maxent_output[[2]]
    write.csv(Maxent_output[[1]],paste0("./output/Maxent_5CV/Maxent_result_summary/grid_target/Maxent_Summary_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".csv"))
    save(Maxent_mod,file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_target/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
    write.csv(Maxent_output[[3]],paste0("./output/Maxent_5CV/Maxent_pred/grid_target/pred_bg10e4Bgbias_",formatC(m,width = 2,flag = 0),".csv"))
  }
}


# prediction for whole map -------------------------------------------------------
## both
## --original
all_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_both/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  
  all_map_pred_mth = NULL
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    all_map_pred_fold = dismo::predict(maxent_mod_indiv,Env[[m-8]])
    all_map_pred_mth = c(all_map_pred_mth,all_map_pred_fold)
  }
  all_map_pred = cbind(all_map_pred,all_map_pred_mth)
}
write.csv(all_map_pred,"./output/Maxent_5CV/Maxent_pred_Wholemap/grid_both/pred_bg10e4Bgbias_both_allenv.csv")


## byactch
## --original
all_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_bycatch/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  
  all_map_pred_mth = NULL
    for(f in 1:5){
      maxent_mod_indiv = Maxent_mod[[f]]
      all_map_pred_fold = dismo::predict(maxent_mod_indiv,Env[[m-8]])
      all_map_pred_mth = c(all_map_pred_mth,all_map_pred_fold)
    }
    all_map_pred = cbind(all_map_pred,all_map_pred_mth)
}
write.csv(all_map_pred,"./output/Maxent_5CV/Maxent_pred_Wholemap/grid_bycatch/pred_bg10e4Bgbias_bycatch_allenv.csv")


## target
## --original
all_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_target/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  
  all_map_pred_mth = NULL
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    all_map_pred_fold = dismo::predict(maxent_mod_indiv,Env[[m-8]])
    all_map_pred_mth = c(all_map_pred_mth,all_map_pred_fold)
  }
  all_map_pred = cbind(all_map_pred,all_map_pred_mth)
}
write.csv(all_map_pred,"./output/Maxent_5CV/Maxent_pred_Wholemap/grid_target/pred_bg10e4Bgbias_target_allenv.csv")



# prediction for fishing ground -------------------------------------------------------

# import true occ data name
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)

## both
fg_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_both/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Env_mth_data = Env[[m-8]]
  obs_occ$Month = m
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    obs_occ$pred_fg = dismo::predict(maxent_mod_indiv,occ_env)
    obs_occ$mod_fold = f
    fg_map_pred = rbind(fg_map_pred,obs_occ)
  }
}
write.csv(fg_map_pred,"./output/Maxent_5CV/Maxent_pred_fg_for_AUC/grid_both/pred_bg10e4Bgbias_both_fg.csv")



## byactch
## --original 
fg_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_bycatch/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Env_mth_data = Env[[m-8]]
  obs_occ$Month = m
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }

  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    obs_occ$pred_fg = dismo::predict(maxent_mod_indiv,occ_env)
    obs_occ$mod_fold = f
    fg_map_pred = rbind(fg_map_pred,obs_occ)
  }
}
write.csv(fg_map_pred,"./output/Maxent_5CV/Maxent_pred_fg/grid_bycatch/pred_bg10e4Bgbias_bycatch_fg.csv")


## target 
## --original
fg_map_pred = NULL
for(m in 9:11){
  load(file=paste0("./output/Maxent_5CV/Maxent_result_summary/grid_target/Maxent_mod_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".RData"))
  obs_occ = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  Env_mth_data = Env[[m-8]]
  obs_occ$Month = m
  
  # arrange
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = obs_occ %>%dplyr::filter(Year==y)
    Env_mth_yr = Env_mth_data %>%dplyr::filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  for(f in 1:5){
    maxent_mod_indiv = Maxent_mod[[f]]
    obs_occ$pred_fg = dismo::predict(maxent_mod_indiv,occ_env)
    obs_occ$mod_fold = f
    fg_map_pred = rbind(fg_map_pred,obs_occ)
  }
}
write.csv(fg_map_pred,"./output/Maxent_5CV/Maxent_pred_fg/grid_target/pred_bg10e4Bgbias_target_fg.csv")




