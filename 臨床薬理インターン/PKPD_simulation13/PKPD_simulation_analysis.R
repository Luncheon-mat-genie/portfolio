# シミュレーション結果の解析をする。

library(tidyverse)
library(data.table)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
datapath = file.path(basepath,"output_ke_eta1_exp","csv")
outpath = file.path(basepath,"output_ke_eta1_exp","csv")


# メモリが足りないので一つ一つのデータを読み込んで処理していく

# 24hに一回投与する群
# timing 22hr
# days 84
# dose 1,10,20,40,60,80,100,120,140,160,200,240,280,320,360

dose_list = c(1,10,20,40,60,80,100,120,140,160,200,240,280,320,360)
days = 84
timing = c(22)

# Cmaxを算出する関数
calc_Cmax = function(data){
  Cmax = max(data$C)
  Cmax
}

# 24hごとのAUCを計算する関数う
# 24hごとのAUCを計算する。
calc_24AUC2 = function(time,data){
  time_range = c(time:(time+23))
  C = data$C[time_range]
  AUC = sum(C)
  AUC
}

calc_24AUC = function(data){
  # 投与開始から2016hまで24hごとにAUCを計算する。
  time_list = seq(min(data$t),2016,24)
  calc_time_list = time_list[-(length(time_list))]-(min(data$t)-1)
  AUC = as.numeric(unlist(mapply(calc_24AUC2,
                                 calc_time_list,list(data))))
  # AUCの最大値を出力
  AUC_max = max(AUC)
  AUC_max
}

# 各投与ごとのtmaxを計算する。
# tmaxを計算する。
calc_tmax2 = function(time,data){
  time_range = c(time:(time+23))
  start_time = data$t[time]
  tibble_f1 = data[time_range,]
  tmax = tibble_f1$t[which.max(tibble_f1$C)]
  tmax = tmax-start_time
}

calc_tmax = function(data,time_gap){
  # 投与開始から2016hまでtime_gapごとにtmaxを計算する。
  time_list = seq(min(data$t),2016,24)
  calc_time_list = time_list[-(length(time_list))]-(min(data$t)-1)
  tmax = as.numeric(unlist(mapply(calc_tmax2,
                                 calc_time_list,list(data))))
  list(tmax)
}

# 5%, 50%, 95%, aveを取り出す
pull_rv = function(data,category){ # categoryは取り出したい列名
  n = nrow(data)
  s = paste("value = data$",category,sep="")
  eval(parse(text =s))
  value=value %>% sort()
  value_rv = value[c(n*0.05,n*0.5,n*0.95)]
  ave = mean(value)
  tibble_f1 = tibble(
    value = c(value_rv,ave),
    percent = c("5%","50%","95%","ave")
  )
  list(tibble_f1)
}

tibble_Cmax_auc = list()
tibble_tmax = list()
for(dose in dose_list){
  print(dose)
  
  # データを読み込む
  csvname = paste("result_",timing,"h_",dose,"mg_",days,"d.csv",sep="")
  setwd(datapath)
  tibble_p1 = fread(csvname,data.table=FALSE)
  tibble_p1_ID = tibble_p1 %>% 
    group_by(ID) %>% 
    nest() %>% 
    ungroup()
  tibble_p2 = tibble(
    dose=dose,
    ID = tibble_p1_ID$ID,
    Cmax = as.numeric(unlist(mapply(calc_Cmax,tibble_p1_ID$data))),
    AUC_max = as.numeric(unlist(mapply(calc_24AUC,tibble_p1_ID$data)))
  )
  
  tibble_p1_t = tibble_p1 %>% 
    group_by(t) %>% 
    nest() %>% 
    ungroup()
  tibble_p1_rv = tibble(
    t = tibble_p1_t$t,
    mapply(pull_rv,tibble_p1_t$data,"C")
    )
  tibble_p1_rv = tibble_p1_rv %>% 
    group_by(t) %>% 
    unnest() %>% 
    ungroup()
  
  ggplot(tibble_p1_rv,aes(t,value,color=percent))+
    geom_line()
  
  tibble_Cmax_rv = mapply(pull_rv,
                          tibble_p2,"Cmax")
  
  # ID100抜き出してプロット
  tibble_p3 = tibble_p1_ID[c(1:2),] %>% 
    group_by(ID) %>% 
    unnest() %>% 
    ungroup() %>% 
    filter(t%in%c(0:5000))
  
  tibble_test = tibble_p2[c(1:1000),]
  
  ggplot(tibble_p3,aes(t,C,color=factor(ID),group=ID))+
    geom_line(linewidth=1)
  setwd(outpath)
  png("C_time.png",res=300,width=2500,height = 1500)
  ggplot(tibble_p3,aes(t,C,color=factor(ID),group=ID))+
    geom_line(linewidth=1)
  dev.off()
  
  
  # tmaxを計算
  if(length(timing)==1){
    time_gap=24
  } else if(lenth(timing)==2){
    time_gap=12
  }
  
  tibble_p_tmax =  tibble(
    dose=dose,
    tmax = unlist(mapply(calc_tmax,tibble_p1_ID$data,time_gap))
  )
  
  tibble_Cmax_auc[[length(tibble_Cmax_auc)+1]] = tibble_p2
  tibble_tmax[[length(tibble_tmax)+1]] = tibble_p_tmax
}

tibble_Cmax_auc = tibble_Cmax_auc %>% 
  bind_rows()
tibble_tmax = tibble_tmax %>% 
  bind_rows()

# doseごとにCmaxをプロット

ggplot(tibble_Cmax_auc,aes(dose,Cmax,group=dose))+
  geom_boxplot()+
  scale_y_continuous(breaks=seq(0,1000,100))

ggplot(tibble_Cmax_auc,aes(dose,AUC_max,group=dose))+
  geom_boxplot()




dose_list = c(1,10,20,40,60,80,100,120,140,160,200,240,280,320,360)
days = 84
timing = c(9,19)

csv_list_2t=c()
for(dose in dose_list){
  t = paste(timing,sep="",collapse = "_")
  do = paste(dose,"_",dose,sep="")
  csvname = paste("result_",t,"h_",do,"mg_",days,"d.csv",sep="")
  csv_list_2t=c(csv_list_2t,csvname)
}