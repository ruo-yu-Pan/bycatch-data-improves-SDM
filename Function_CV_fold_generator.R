
  
generate_CV_fold = function(Cutlassfish_occ_data,nfold=5){
  
  Cutlassfish_occ = Cutlassfish_occ_data
  Cutlassfish_occ$nth_fold = NA
  
  presence_pos = which(Cutlassfish_occ$occ==1)
  absence_pos = which(Cutlassfish_occ$occ==0)
  
  presence_fold_sample = sample(1:length(presence_pos),length(presence_pos),replace = F)
  absence_fold_sample = sample(1:length(absence_pos),length(absence_pos),replace = F)
  
  presence_pos = cbind(presence_pos,presence_fold_sample)
  absence_pos = cbind(absence_pos,absence_fold_sample)
  
  presence_pos = presence_pos[order(presence_pos[,2]),]
  absence_pos = absence_pos[order(absence_pos[,2]),]
  
  pres_numb_in_one_fold = floor(length(presence_pos[,1])/nfold)
  abse_numb_in_one_fold = floor(length(absence_pos[,1])/nfold)
  
  for(i in 1:(nfold-1)){
    Cutlassfish_occ$nth_fold[presence_pos[(pres_numb_in_one_fold*(i-1)+1):(pres_numb_in_one_fold*i),1]] = i
    Cutlassfish_occ$nth_fold[absence_pos[(abse_numb_in_one_fold*(i-1)+1):(abse_numb_in_one_fold*i),1]] = i
  }
  Cutlassfish_occ$nth_fold[which(is.na(Cutlassfish_occ$nth_fold))] = nfold
  
  return(Cutlassfish_occ)
}
  


