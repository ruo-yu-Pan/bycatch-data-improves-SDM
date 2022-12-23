
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)


library("ggplot2")
library("viridis")
library("colorRamps")
library("dplyr")

library("ggplot2")
library("maps")
mp <- map_data("world")

# import environment data
Env_data = read.csv("D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data/raw_data/Environmental data/FillNA/EV_Knn_0.1.csv",stringsAsFactors = F)
Env = list()
for(m in 1:3){
  Env_mth = Env_data%>%dplyr::filter(Month == m+8)
  Env_mth$grid_lab = as.character(Env_mth$grid_lab)
  Env[[m]] = Env_mth
}
# import true pa data
Cutlassfish_both_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_both",full.names=T)
Cutlassfish_bycatch_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_bycatch",full.names=T)
Cutlassfish_target_occ_filename=list.files("./compiled_data/occ_data_in_ground/grid_target",full.names=T)

# original data ################################################################
# import predicted data in fishing ground
path="./output/fig/predict_prob/original/"
n = 16
bycatch_mean_HS = NULL
target_mean_HS = NULL

for(m in 9:11){
  fg_pred_both = read.csv(
    paste0("./output/Maxent_allData/Maxent_pred_fg/grid_both/pred_bg10e4Bgbias_",formatC(m,width = 2,flag=0),".csv"),row.names = 1)
  fg_pred_bycatch = read.csv(
    paste0("./output/Maxent_allData/Maxent_pred_fg/grid_bycatch/pred_bg10e4Bgbias_",formatC(m,width = 2,flag=0),".csv"),row.names = 1)
  fg_pred_target = read.csv(
    paste0("./output/Maxent_allData/Maxent_pred_fg/grid_target/pred_bg10e4Bgbias_",formatC(m,width = 2,flag=0),".csv"),row.names = 1)

  fg_pred_mth = rbind(fg_pred_both,fg_pred_bycatch,fg_pred_target)
  fg_pred_mth$gear = rep(c("integrated","bycatch","target"),each=nrow(fg_pred_bycatch))
  fg_pred_mth$gear = factor(fg_pred_mth$gear, 
                            levels = c("integrated",
                                       "bycatch",
                                       "target"))
  
  raw_data = read.csv(Cutlassfish_both_occ_filename[m],row.names = 1)
  raw_data_occ = raw_data%>%dplyr::filter(occ==1)
  
  colnames(fg_pred_mth)[24]="HS"

  p = ggplot()+ # open a blank plot
    coord_sf(xlim = c(118,126.2), ylim = c(23,30), expand = F)+
    geom_tile(data=fg_pred_mth,aes(x=Lon.x,y=Lat.y,fill=HS))+
    scale_fill_viridis(discrete=FALSE, option="viridis",limits = c(0,1))+
    geom_polygon(data=mp , aes(x = long, y = lat , group=group), fill="lightgray", colour = NA) +
    #geom_point(data=raw_data_occ,aes(x=Lon,y=Lat),fill="white",colour="black",shape=21,size=0.6,stroke=0.1)+
    #scale_fill_manual(values=c("white")) + 
    facet_grid(gear~Year.x)+
    ggtitle(month.name[m])+
    theme_bw()+
    ylab("Latitude")+
    xlab("Longitude")+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
    

  ggsave(p, file = paste0(path,"predProb_and_occ_",formatC(m,width = 2,flag = 0),".jpeg") ,width =16, height =8, dpi=300, units = "in",limitsize = FALSE) 
}

