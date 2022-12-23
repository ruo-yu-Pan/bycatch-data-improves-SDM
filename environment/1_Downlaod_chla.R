
wd = "D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/"
setwd(wd)

#-------------------------------------------
# Downlaod Chla data

month_days_365 = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
month_days_366 = c(0,31,60,91,121,152,182,213,244,274,305,335,366)

for(y in 2012:2019){
  for(m in 1:12) { 
    if(y!=2012&y!=2016){
      appkey = "c5d4c690ddd9fcc73f083f5c7c95aeeeb2febd7a"
      filename = paste0("A",y,formatC(month_days_365[m]+1, width = 3, flag = 0),
                        y,formatC(month_days_365[m+1], width = 3, flag = 0),
                        ".L3m_MO_CHL_chlor_a_4km.nc")
      url <- paste0("https://oceandata.sci.gsfc.nasa.gov/ob/getfile/",filename,"?appkey=",appkey)
      destfile <- paste0("D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/raw_data/Environmental Data/Download/Chla_Monthly_4km/",filename)
      download.file(url, destfile, mode="wb")
    }
    
    if(y==2012|y==2016){
      filename = paste0("A",y,formatC(month_days_366[m]+1, width = 3, flag = 0),
                        y,formatC(month_days_366[m+1], width = 3, flag = 0),
                        ".L3m_MO_CHL_chlor_a_4km.nc")
      url <- paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/",filename,"?appkey=",appkey)
      destfile <- paste0("D:/Ruo_data/2020_Assis_FishTactic/Fishtactics_SDM_data_check/raw_data/Environmental Data/Download/Chla_Monthly_4km/",filename)
      download.file(url, destfile, mode="wb")
    }
  }
}
