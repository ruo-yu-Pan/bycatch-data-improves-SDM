
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
  Env_mth = Env_data%>%dplyr::filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}

# function for MaxEnt model training and prediction for 5 fold CV
Maxent_train_pred_func = function(occ_data=Cutlassfish_occ,
                                  Env_mth_data=Env_mth,
                                  bg_env_data=bg_env_mth,
                                  save_dir){
  # arrange for model training
  occ_env = NULL
  for(y in 2009:2019){
    occ_yr = occ_data %>%dplyr::filter(occ==1&Year==y)
    Env_mth_yr = Env_mth_data %>%filter(Year==y)
    
    occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
    occ_env=rbind(occ_env,occ_env_yr)
  }
  
  # MaxEnt 

    mod_occ = c(occ_env$occ,rep(0,nrow(bg_env_data)))
    mod_env = rbind(occ_env[,c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")],
                          bg_env_data[,c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")])
  
    # maxent
    
    ifelse(!dir.exists(save_dir), 
           dir.create(file.path(save_dir)),FALSE)
    
    maxent_mod = maxent(x=mod_env,
                        p=mod_occ,
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
                        path=save_dir)
    
    mod_result = maxent_mod@results
    Maxent_model_summary = mod_result
    Maxent_model_list = list(maxent_mod)
    
    # predict for fishing ground -----------------------------------
    
    ## arrange for fishing ground data
    occ_fg_env = NULL
    for(y in 2009:2019){
      occ_yr = occ_data %>%filter(Year==y)
      Env_mth_yr = Env_mth_data %>%filter(Year==y)
      occ_env_yr = merge(x=occ_yr,y=Env_mth_yr, by="grid_lab",all.x=T)
      occ_fg_env=rbind(occ_fg_env,occ_env_yr)
    }

    mod_env_fg = occ_fg_env[,c("Bathy","MLD","SSH","SSS","SST","EKE","GM","chla")]
    
    pred_occ = data.frame(pred = predict(maxent_mod,mod_env_fg))
    pred_occ = cbind(occ_fg_env,pred_occ)

  output_list = list(Maxent_model_summary,Maxent_model_list,pred_occ)
}

# Train Maxent model and predict ------------------------------------------------
foldname=c("grid_both","grid_bycatch","grid_target")
gear = c("both","bycatch","target")


for(fod in 1:3){
  Cutlassfish_filename=
    list.files(paste0("./compiled_data/occ_data_in_ground/",foldname[fod]),full.names=T)
  bg_env_filename=
    list.files(paste0("./compiled_data/bg_file/",foldname[fod]),
               pattern =glob2rx("bg10e4Bgbias_*.csv") ,full.names=T)
  
  for (m in 9:11){
    
    save_dir = paste0("./output/Maxent_allData/Maxent_result_detail/",foldname[fod],"/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias")
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    if(fod==3){
      Maxent_output = 
        Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_filename[m-3],row.names = 1,stringsAsFactors = F),
                               Env_mth_data=Env[[m-8]],
                               bg_env_data=read.csv(bg_env_filename[m-8],row.names = 1,stringsAsFactors = F),
                               save_dir=save_dir)
    }else{
      Maxent_output = 
        Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_filename[m],row.names = 1,stringsAsFactors = F),
                               Env_mth_data=Env[[m-8]],
                               bg_env_data=read.csv(bg_env_filename[m-8],row.names = 1,stringsAsFactors = F),
                               save_dir=save_dir)
    }
  
    Maxent_mod = Maxent_output[[2]]
    write.csv(Maxent_output[[1]],paste0("./output/Maxent_allData/Maxent_result_summary/",foldname[fod],"/Maxent_Summary_bg10e4Bgbias_",formatC(m,width = 2,flag = 0),".csv"))
    save(Maxent_mod,file=paste0("./output/Maxent_allData/Maxent_result_summary/",foldname[fod],"/Maxent_mod_bg10e4Bgbias_",formatC(m,width = 2,flag = 0),".RData"))
    write.csv(Maxent_output[[3]],paste0("./output/Maxent_allData/Maxent_pred_fg/",foldname[fod],"/pred_bg10e4Bgbias_",formatC(m,width = 2,flag = 0),".csv"))
  }
}


### -resample


for (m in 9:11){
  # list the  resample file
  Cutlassfish_resample_filename = list.files(paste0("./compiled_data/occ_data_in_ground/test_smp_range/",
                                                            formatC(m,width = 2,flag = 0)),
                                                     full.names=T)
  # list the bg file of bycatch resample data 
  bg_env_resample_filename=list.files("./compiled_data/bg_file/test_smp_range/bg10e4Bgbias",
                                              pattern =glob2rx(paste0("*",formatC(m,width = 2,flag = 0),"*.csv")) ,full.names=T)
  
  for(file in 1: length(Cutlassfish_resample_filename)){
    dash_pos = gregexpr("_",Cutlassfish_resample_filename[file])[[1]]
    slash_pos = gregexpr("/",Cutlassfish_resample_filename[file])[[1]]
    gear = substr(Cutlassfish_resample_filename[file],slash_pos[[5]]+1,dash_pos[[7]]-1)
    mm = as.numeric(substr(Cutlassfish_resample_filename[file],slash_pos[[4]]+1,slash_pos[[4]]+2))
    SXX = substr(Cutlassfish_resample_filename[file],dash_pos[[7]]+1,dash_pos[[7]]+3)
    NX_rX = substr(Cutlassfish_resample_filename[file],dash_pos[[8]]+1,dash_pos[[9]]+3)
    Nlab = substr(Cutlassfish_resample_filename[file],dash_pos[[8]]+2,dash_pos[[8]]+2)
    r = as.numeric(substr(Cutlassfish_resample_filename[file],dash_pos[[9]]+2,dash_pos[[9]]+3))

    
    save_dir = paste0("./output/Maxent_allData/Maxent_result_detail/test_smp_range/",formatC(m,width = 2,flag = 0),"_bg10e4Bgbias/",gear,"_",SXX,"_",NX_rX)
    # create folder
    ifelse(!dir.exists(save_dir), dir.create(file.path(save_dir)),FALSE)
    
    # Maxent Model
    Maxent_output = 
      Maxent_train_pred_func(occ_data=read.csv(Cutlassfish_resample_filename[file],row.names = 1,stringsAsFactors = F),
                             Env_mth_data=Env[[m-8]],
                             bg_env_data=read.csv(bg_env_resample_filename[file],row.names = 1,stringsAsFactors = F),
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
      
      write.csv(Maxent_model_summary_resample,paste0("./output/Maxent_allData/Maxent_result_summary/test_smp_range/Maxent_Summarybg10e4Bgbias_",formatC(m,width = 2,flag = 0),SXX,"_N",Nlab,".csv"))
      save(Maxent_model_resample,file=paste0("./output/Maxent_allData/Maxent_result_summary/test_smp_range/Maxent_modbg10e4Bgbias_",formatC(m,width = 2,flag = 0),SXX,"_N",Nlab,".RData"))
      write.csv(Maxent_pred_resample,paste0("./output/Maxent_allData/Maxent_pred_fg/test_smp_range/predbg10e4Bgbias_",formatC(m,width = 2,flag = 0),"_",SXX,"_N",Nlab,".csv"))
    }else{
      Maxent_mod = Maxent_output[[2]]
      Maxent_model_summary_resample = cbind(Maxent_model_summary_resample,Maxent_output[[1]])
      Maxent_pred_resample = rbind(Maxent_pred_resample,Maxent_output[[3]])
      Maxent_model_resample = c(Maxent_model_resample,Maxent_mod)
    }
  }
}


