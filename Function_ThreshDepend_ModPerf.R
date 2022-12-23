
library(PresenceAbsence)

### function for OR for original data
Thresh_depend_ModPerf_func = function(mth_pred_data=fg_pred_mth,
                                      thresh=Thresh[which(Thresh$gear==gear[fod]&Thresh$m==m),1]){
  
  fg_pred_mth = as.data.frame(mth_pred_data)
  
  ModPerf_info = NULL
  for(f in 1:5){
    fg_pred_mth_fold = fg_pred_mth %>% dplyr::filter(nth_fold==f & mod_fold==f)
    
    ## - confusion matrix
    cm_dat = data.frame(plotID = c(1:nrow(fg_pred_mth_fold)),
                        Observed = fg_pred_mth_fold$occ,
                        Predicted = fg_pred_mth_fold$pred_fg)
    cm_dat = na.omit(cm_dat)
    cm_tab = cmx(cm_dat, threshold = thresh[f])
    
    tss = sensitivity(cm_tab)[,1]+specificity(cm_tab)[,1]-1
    
    ## F measure
    precision = diag(cm_tab)[1] / apply(cm_tab,1,sum)[1]
    recall = diag(cm_tab)[1] / apply(cm_tab,2,sum)[1]
    F2 = (1+4)* precision * recall / (4*precision + recall) 
    
    ModPerf_info = rbind(ModPerf_info,
                               c(tss,F2,1-sensitivity(cm_tab)[,1],
                                 sensitivity(cm_tab)[,1],
                                 specificity(cm_tab)[,1],precision,
                                 cm_tab[1,1],cm_tab[1,2],cm_tab[2,1],cm_tab[2,2]))
  }
  colnames(ModPerf_info) = c("TSS","F2","OR","Sensitive","Specificity","Precision",
                             "N00","N01","N10","N11")
  return(ModPerf_info)
}

### function for OR for resample data
Thresh_depend_ModPerf_resample_func = function(mth_pred_data=fg_pred_resample_onesample,thresh=thresh_r[,2]){
  
  fg_pred_mth = as.data.frame(mth_pred_data)
  
  ModPerf_info = NULL
  for(f in 1:5){
    fg_pred_mth_fold = fg_pred_mth %>% dplyr::filter(resample_nthfold==f & mod_fold==f)
    
    ## - confusion matrix
    cm_dat = data.frame(plotID = c(1:nrow(fg_pred_mth_fold)),
                        Observed = fg_pred_mth_fold$occ,
                        Predicted = fg_pred_mth_fold$pred_fg)
    cm_dat = na.omit(cm_dat)
    cm_tab = cmx(cm_dat, threshold = thresh[f])
    
    tss = sensitivity(cm_tab)[,1]+specificity(cm_tab)[,1]-1
    
    ## F measure
    precision = diag(cm_tab)[1] / apply(cm_tab,1,sum)[1]
    recall = diag(cm_tab)[1] / apply(cm_tab,2,sum)[1]
    F2 = (1+4)* precision * recall / (4*precision + recall) 
    
    ModPerf_info = rbind(ModPerf_info,
                         c(tss,F2,1-sensitivity(cm_tab)[,1],
                           sensitivity(cm_tab)[,1],
                           specificity(cm_tab)[,1],precision,
                           cm_tab[1,1],cm_tab[1,2],cm_tab[2,1],cm_tab[2,2]))
  }
  colnames(ModPerf_info) = c("TSS","F2","OR","Sensitive","Specificity","Precision",
                             "N00","N01","N10","N11")
  return(ModPerf_info)
}