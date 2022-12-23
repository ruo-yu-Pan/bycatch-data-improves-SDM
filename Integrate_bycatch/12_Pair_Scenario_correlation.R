wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library("ggplot2")

pred_folder_name = list.files("./output/Maxent_allData/Maxent_pred_fg",full.names=T)

# Organize the pred allData ####################################################
pred_allData_all = list()
for(m in 9:11){
  pred_allData_mth = NULL
  for(dir in c(1,2,3)){
    pred_filename = list.files(pred_folder_name[dir],full.names = T)
    pred_allData = read.csv(pred_filename[m-8],row.names = 1)
    pred_allData_mth = cbind(pred_allData_mth,pred_allData$pred)
  }
  pred_allData_all[[m-8]] = pred_allData_mth
}

# Calculate correlation and arrange it to df ####################################################
Cor_mat_ls = list()
Cor_pair_scen = list()

for(i in 1:3){
  Cor_mat_ls[[i]] = round(cor(pred_allData_all[[i]]),2)
  Cor_mat = Cor_mat_ls[[i]]
  diag(Cor_mat)=NA

  Cor_pair_scen[[i]] =Cor_mat
}

Cor_pair_scen[[1]]
save(Cor_pair_scen,file = "./output/ModelPerform_Comp_originData/CorMat_CorPairDf_list.RData")
