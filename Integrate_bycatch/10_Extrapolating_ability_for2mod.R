
######################################################################
# extrapolating ability
######################################################################
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library(dplyr)


# import true pa data
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)
Cutlassfish_bycatch_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_bycatch",full.names=T)
Cutlassfish_target_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_target",full.names=T)

# ----------------------------------
Thresh_P10_func = function(occPredVals){
  if(length(occPredVals) < 10){
    p10 <- floor(length(occPredVals) * 0.9)
  } else {
    p10 <- ceiling(length(occPredVals) * 0.9)
  }
  thresh <- rev(sort(occPredVals))[p10]
  return(thresh)
}
# ----------------------------------


bycatch_mean_HS = NULL
target_mean_HS = NULL

for(m in 9:11){
  
  # import allData prediction
  fg_pred_both = read.csv(
    paste0("./output/Maxent_allData/Maxent_pred_fg/grid_both/pred_bg10e4Bgbias_",formatC(m,width = 2,flag=0),".csv"),row.names = 1)
  fg_pred_bycatch = read.csv(
    paste0("./output/Maxent_allData/Maxent_pred_fg/grid_bycatch/pred_bg10e4Bgbias_",formatC(m,width = 2,flag=0),".csv"),row.names = 1)
  fg_pred_target = read.csv(
    paste0("./output/Maxent_allData/Maxent_pred_fg/grid_target/pred_bg10e4Bgbias_",formatC(m,width = 2,flag=0),".csv"),row.names = 1)
  
  # import true data
  raw_data = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  raw_data_occ = raw_data%>%filter(occ==1)
  bycatch_data = read.csv(Cutlassfish_bycatch_occ_filename[m],row.names = 1)
  bycatch_data_occ = bycatch_data%>%dplyr::filter(occ==1)
  target_data = read.csv(Cutlassfish_target_occ_filename[m-3],row.names = 1)
  target_data_occ = target_data%>%dplyr::filter(occ==1)
  
  # Threshold
  thresh10P_both = Thresh_P10_func(
    fg_pred_both$pred[which(fg_pred_both$occ==1)])
  thresh10P_bycatch = Thresh_P10_func(
    fg_pred_bycatch$pred[which(fg_pred_bycatch$occ==1)])
  thresh10P_target = Thresh_P10_func(
    fg_pred_target$pred[which(fg_pred_target$occ==1)])
  
  
  
  bycatch_only_pred = NULL
  target_only_pred = NULL
  for(y in 2009:2019){
    fg_pred_bycatch_yr = fg_pred_bycatch %>% dplyr::filter(Year.x==y) 
    fg_pred_target_yr = fg_pred_target %>% dplyr::filter(Year.x==y) 
    
    raw_data_occ_yr=raw_data_occ %>% dplyr::filter(Year==y)
    bycatch_data_occ_yr=bycatch_data_occ %>% dplyr::filter(Year==y)
    target_data_occ_yr=target_data_occ %>% dplyr::filter(Year==y)
    
    fg_pred_bycatch_yr = fg_pred_bycatch_yr[
      which(fg_pred_bycatch_yr$grid_lab %in% unique(raw_data_occ_yr$grid_lab)),]
    fg_pred_target_yr = fg_pred_target_yr[
      which(fg_pred_target_yr$grid_lab %in% unique(raw_data_occ_yr$grid_lab)),]
    
    bycatch_pos_bMod = which(fg_pred_bycatch_yr$grid_lab %in% unique(bycatch_data_occ_yr$grid_lab))
    target_pos_bMod = which(fg_pred_bycatch_yr$grid_lab %in% unique(target_data_occ_yr$grid_lab))
    
    bycatch_pos_tMod = which(fg_pred_target_yr$grid_lab %in% unique(bycatch_data_occ_yr$grid_lab))
    target_pos_tMod = which(fg_pred_target_yr$grid_lab %in% unique(target_data_occ_yr$grid_lab))
    
    bycatch_only_tMod = fg_pred_target_yr[!(bycatch_pos_tMod %in% target_pos_tMod),]
    target_only_bMod = fg_pred_bycatch_yr[!(target_pos_tMod %in% bycatch_pos_tMod),]
    
    if(nrow(bycatch_only_tMod)>0){
      bycatch_only_pred = rbind(bycatch_only_pred,bycatch_only_tMod)
    }
    if(nrow(target_only_bMod)>0){
      target_only_pred = rbind(target_only_pred,target_only_bMod)
    }
  }
  
  # ----------------------------------
  bycatch_only_pred$pred_occ = ifelse(bycatch_only_pred$pred>=thresh10P_target,1,0)
  target_only_pred$pred_occ = ifelse(target_only_pred$pred>=thresh10P_bycatch,1,0)
  
  
  bycatch_mean_HS = rbind(bycatch_mean_HS,
                          c(median(bycatch_only_pred$pred),
                            mean(bycatch_only_pred$pred),
                            sd(bycatch_only_pred$pred),
                            (sum(bycatch_only_pred$pred_occ)/nrow(bycatch_only_pred))))
  target_mean_HS = rbind(target_mean_HS,
                         c(median(target_only_pred$pred),
                           mean(target_only_pred$pred),
                           sd(target_only_pred$pred),
                           (sum(target_only_pred$pred_occ)/nrow(target_only_pred))))
}

bycatch_mean_HS_tmod = as.data.frame(bycatch_mean_HS)
target_mean_HS_bmod = as.data.frame(target_mean_HS)

colnames(bycatch_mean_HS_tmod) = c("med_HS","mu_HS","sd_HS","Sensitivity")
colnames(target_mean_HS_bmod) = c("med_HS","mu_HS","sd_HS","Sensitivity")

write.csv(bycatch_mean_HS_tmod,"./output/Model_performance/Extrapolation_ability/Extrapolate_target_mod.csv")
write.csv(target_mean_HS_bmod,"./output/Model_performance/Extrapolation_ability/Extrapolate_bycatch_mod.csv")

