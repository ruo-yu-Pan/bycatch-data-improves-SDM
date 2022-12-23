
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)
Root_dir = "D:/Ruo_data/2020_Assis_FishTactic/"

library("ggplot2")
library("gridExtra")
library("pROC")
library("dplyr")
library("reshape2")
library("sp")

# Calculate AUC based on randomBg ----------------------------------------------------------------------------------
pred_fg_foldname = list.files("./output/Maxent_5CV/Maxent_pred_fg_for_AUC",full.names=T)
foldname = c("grid_both","grid_bycatch","grid_target")
gear = c("both","bycatch","target")

for(fod in 1:3){
  
  #if(fod==1){gear="bycatch"}else{gear="target"}
  
  pred_fg_filename = list.files(pred_fg_foldname[fod],full.names = T)
  
  for(fil in 1:length(pred_fg_filename)){
    for(m in 9:11){
      AUC_value = NULL
      
      pred_fg = read.csv(pred_fg_filename[fil])
      pred_fg_mth = pred_fg %>% dplyr::filter(Month.x==m)
      
      for(mod_f in 1:5){
        pred_fg_train = pred_fg_mth %>% dplyr::filter(mod_fold==mod_f&nth_fold!=mod_f)
        pred_fg_test = pred_fg_mth %>% dplyr::filter(mod_fold==mod_f&nth_fold==mod_f)
        
        # Train AUC
        train_auc_result = roc(response = pred_fg_train$occ, 
                               predictor = pred_fg_train$pred_fg,
                               levels = c(0,1), direction='<')$auc
        
        # Test AUC
        test_auc_result = roc(response = pred_fg_test$occ, 
                              predictor = pred_fg_test$pred_fg,
                              levels = c(0,1), direction='<')$auc
        auc_result = c(train_auc_result,test_auc_result,mod_f)
        
        # save
        AUC_value = rbind(AUC_value,auc_result)
      }
      
      AUC_value = as.data.frame(AUC_value)
      colnames(AUC_value)=c("train_AUC","test_AUC","mod_fold")
      
      write.csv(AUC_value,paste0("./output/Model_performance/AUC_fg_randomBg/",foldname[fod],"/AUC_bg10e4Bgbias_",gear[fod],"_",formatC(m,width = 2,flag = 0),".csv"))
      
    }
  }
}

# Organize AUC and threshold -----------------------------------------------------------
foldname = c("grid_both","grid_bycatch","grid_target")
gear = c("both","bycatch","target")
## 5foldCV allmap training AUC, testing AUC and CBI, 
## and self-mod Threshold data #######

training_AUC = NULL
testing_AUC = NULL
Thresh = NULL

for(fod in 1:3){
  for(m in 9:11){
    AUC_fg = read.csv(paste0("./output/Model_performance/AUC_fg_randomBg/",
                             foldname[fod],
                             "/AUC_bg10e4Bgbias_",gear[fod],"_",
                             formatC(m,width = 2,flag = 0),".csv"))
    training_AUC = c(training_AUC,AUC_fg$train_AUC)
    testing_AUC = c(testing_AUC,AUC_fg$test_AUC)
    
    Cutlassfish_bycatch_5CV_summary=read.csv(
      paste0("./output/Maxent_5CV/Maxent_result_summary/",foldname[fod],
             "/Maxent_Summary_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".csv"),
      row.names = 1)
    
    thresh_10P = unlist(c(Cutlassfish_bycatch_5CV_summary["X10.percentile.training.presence.Cloglog.threshold",]))
    Thresh = rbind(Thresh,cbind(thresh_10P,gear[fod],m))
  }
}
Thresh = as.data.frame(Thresh)
Thresh[,c(1,3)] = apply(Thresh[,c(1,3)],2,function(x){as.numeric(as.character(x))})
Thresh[,2] = as.character(Thresh[,2])
colnames(Thresh)[c(2:3)] = c("gear","month")



## 5foldCV testing Omission rate ########
source(paste0(Root_dir,"Fishtactics_SDM_code_check/Function_ThreshDepend_ModPerf.R"))

# import environment data
Env_data = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv",stringsAsFactors = F)
Env = list()
for(m in 1:3){
  Env_mth = Env_data%>%filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}

# import true occ data
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)

# calculate OR
Thresh_depend_ModPerf = NULL

for(fod in 1:3){
  fg_pred = read.csv(
    paste0("./output/Maxent_5CV/Maxent_pred_fg_for_AUC/",
           foldname[fod],"/pred_bg10e4Bgbias_",gear[fod],"_fg.csv"),
    row.names = 1)
  
  OR_GEAR_mth = matrix(nrow = 15,ncol = 3)
  for(m in 9:11){
    fg_pred_mth = fg_pred %>% dplyr::filter(Month.x==m)
    
    ## Mod Perf
    ModPerf_gear = 
      Thresh_depend_ModPerf_func(mth_pred_data = fg_pred_mth,
                                 thresh = Thresh[which(Thresh$gear==gear[fod]&Thresh$m==m),1])
    Thresh_depend_ModPerf = rbind(Thresh_depend_ModPerf,ModPerf_gear)
  }
}

Thresh_depend_ModPerf = as.data.frame(Thresh_depend_ModPerf)


# ## Variable contribution #####################################

Var_contrib = NULL
for(fod in 1:3){
  for(m in 9:11){
    Cutlassfish_5CV_summary=read.csv(
      paste0("./output/Maxent_5CV/Maxent_result_summary/",
             foldname[fod],
             "/Maxent_Summary_bg10e4Bgbias_5CV_",formatC(m,width = 2,flag = 0),".csv"),
      row.names = 1)
    Var_contrib = rbind(Var_contrib,
                        t(Cutlassfish_5CV_summary[15:22,]))
  }
}

colnames(Var_contrib) = c("Bathy","EKE","TGM","MLD","SSH","SSS","SST","Chla")


## merge all model performance #######################
Mod_perform = cbind(training_AUC,testing_AUC,
                    Thresh_depend_ModPerf,Var_contrib)

Mod_perform$month = c(rep(rep(9:11,each=5),3))
Mod_perform$fold = rep(c(1:5),9)
Mod_perform$gear = c(rep(c("both","bycatch","target"),each=15))

write.csv(Mod_perform,"./output/Model_performance/Mod_performance_Orig3mod_5CV.csv")

# ----------------------------------------------------------------------------------
## Merge data property
data_property_filename = list.files("./compiled_data/data_property/Original_data",pattern = glob2rx("*occ*.csv"),full.names = T)
data_property_colnames = c("Month","effort_smpsize","Extent_smp_range",
                           "kd_CV","kd_rad_auc","No_site")

data_property = NULL
for(i in 1:3){
  data_property_i = read.csv(data_property_filename[i],row.names = 1,stringsAsFactors = F)
  colnames(data_property_i) = data_property_colnames
  data_property_i$gear = gear[i]
  
  for(m in 1:3){
    data_property_m = rbind(data_property_i[m,],data_property_i[m,],
                            data_property_i[m,],data_property_i[m,],
                            data_property_i[m,])
    data_property = rbind(data_property,data_property_m)
  }
}

Mod_perform_data_property = cbind(Mod_perform,data_property)
write.csv(Mod_perform_data_property, "./output/Model_performance/ModPerf_DataPerf_Orig3mod_5CV.csv")

##############################################################################
# COMPARE MODEL PEROFRMANCE ####
Mod_perform = read.csv("./output/Model_performance/Mod_performance_Orig3mod_5CV.csv",row.names = 1,stringsAsFactors = F)
Mod_perform_data_property = read.csv("./output/Model_performance/ModPerf_DataPerf_Orig3mod_5CV.csv",row.names = 1,stringsAsFactors = F)

# change name
Mod_perform$gear[which(Mod_perform$gear=="both")] = "integrated"
Mod_perform_data_property$gear[which(Mod_perform_data_property$gear=="both")] = "integrated"

# wilcox test
library(dplyr)

wicox_result_all = NULL

pair = list(c("integrated","bycatch"),
            c("integrated","target"),
            c("bycatch","target"))
colnames(Mod_perform)[c(2,6)] = c("AUC","Sensitivity")

for(m in 9:11){
  for(mp in c(2:4,6:8)){
    for(p in 1:3){
      ModPerf_DatProp_origTwoG_one = Mod_perform %>% dplyr::filter(month==m & gear %in% pair[[p]])
      ModPerf_DatProp_origTwoG_one = ModPerf_DatProp_origTwoG_one[,c(mp,23)]
      colnames(ModPerf_DatProp_origTwoG_one)=c("ModPerf","gear")
      wicox_result = wilcox.test(ModPerf~gear,ModPerf_DatProp_origTwoG_one,exact=F)
      wicox_result_all = rbind(wicox_result_all, 
                               c(wicox_result$p.value,m,p,colnames(Mod_perform)[mp]))
    }
  }
}

colnames(wicox_result_all) = c("Pval","month","pair","ModPerf")

# Plot
library(ggplot2)
library(ggpubr)

Mod_perform_forAnalyze = Mod_perform[,c(2:4,6:8,21:23)]
Mod_perform_long = reshape2::melt(Mod_perform_forAnalyze,id=c("month","fold","gear"))

Mod_perform_long$gear = factor(Mod_perform_long$gear, 
                               levels = c("integrated",
                                          "bycatch",
                                          "target"))

gear_compar = list(c("integrated","bycatch"),
                   c("integrated","target"),
                   c("bycatch","target"))

mon.labs <- month.name[9:11]
names(mon.labs) <- c("9", "10", "11")

jpeg("./output/ModelPerform_Comp_originData/Comp_Orig_PerfMetric_3mod.jpeg",
     width = 8,height = 6,res=300,units = "in")
ggplot(Mod_perform_long,aes(x=gear,y=value))+
  geom_boxplot()+
  stat_compare_means(comparisons = gear_compar,
                     aes(group=gear),label.y = c(0.95, 1.15, 1.35),
                     label = "p.signif", method = "wilcox.test")+
  facet_grid(month~variable,
             labeller = labeller(month = mon.labs))+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank())+
  ylim(0,1.5)

dev.off()


