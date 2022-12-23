
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

library(dplyr)

# Environmental data in ground ###################################

Env_data = read.csv("./raw_data/Environmental data/FillNA/EV_Knn_0.1.csv")
grid_in_ground = read.csv("./compiled_data/map/grid_in_ground.csv",row.names = 1)


## grid information ###########################################################
# env grid

for(m in 9:11){
  grid_in_ground_env_yr = NULL
  for(y in 2009:2019){
    Env_data_yr = Env_data %>% filter(Year == y,Month == m)
    grid_env_yr = merge(x=grid_in_ground,y=Env_data_yr,by="grid_lab",all.x=T)
    grid_in_ground_env_yr = rbind(grid_in_ground_env_yr,grid_env_yr)
  }
  grid_in_ground_env_yr = grid_in_ground_env_yr[,-c(4,7,8)]
  write.csv(grid_in_ground_env_yr,
            paste0("./compiled_data/env_data_in_ground/env_",ifelse(m %/%10 ==0,paste0("0",m),m),".csv"))
}


