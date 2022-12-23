
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library(sp)
library("ape")
library(dplyr)
library("spdep")
library("gstat")

fold_name = c("grid_both","grid_bycatch","grid_target")

SAC_val_ls = list()
SAC_pval_ls = list()
SAC_val_mat = NULL
SAC_pval_mat = NULL

for(fod in 1:3){
  pred_fg_filename = list.files(paste0("./output/Maxent_allData/Maxent_pred_fg/",
                                       fold_name[fod]),full.names = T)
  
  SAC_value = NULL
  SAC_pval = NULL
  
  for(m in 9:11){
    pred_fg = read.csv(pred_fg_filename[m-8],row.names = 1)
    pred_fg_mth = pred_fg
    
    
    pred_fg_train = pred_fg_mth %>% dplyr::filter(occ==1)
    resi = data.frame(
      x = pred_fg_train$Lon.x,
      y = pred_fg_train$Lat.x,
      resi = pred_fg_train$occ-pred_fg_train$pred)
    
    # Moran's I
    resi_dists = as.matrix(dist(cbind(resi$x, resi$y)))
    
    resi_dists_inv = 1/resi_dists
    diag(resi_dists_inv) = 0
    resi_dists_inv[is.infinite(resi_dists_inv)] = 0
    
    resi_dists_inv_lsw = mat2listw(resi_dists_inv, 
                                   row.names = NULL, style="M")
    SAC_M = moran.test(resi$resi, resi_dists_inv_lsw)
    main_label = paste0("Moran_I = ",SAC_M$estimate[1])
    png(paste0("./output/Spatial_Auto_COR/",fold_name[fod],"_",month.abb[m],".png"))
    moran.plot(resi$resi, 
               listw = resi_dists_inv_lsw,
               main=main_label)
    dev.off()
    
    
    SAC_value = c(SAC_value,SAC_M$estimate[1])
    SAC_pval=c(SAC_pval,SAC_M$p.value)
    
  }
  
  SAC_val_mat=rbind(SAC_val_mat,SAC_value)
  SAC_pval_mat=rbind(SAC_pval_mat,SAC_pval)
}




# > SAC_val_mat
# Moran I statistic Moran I statistic Moran I statistic
# SAC_value        0.04987565       0.045174780       0.015912481
# SAC_value        0.04831212       0.008699534       0.038979895
# SAC_value        0.03111273       0.049394455       0.006230328

# > SAC_pval_mat
# [,1]         [,2]         [,3]
# SAC_pval 1.733625e-92 2.092410e-31 5.552454e-06
# SAC_pval 2.783170e-78 1.048469e-02 2.868698e-04
# SAC_pval 3.507700e-03 4.437884e-12 1.443724e-02
