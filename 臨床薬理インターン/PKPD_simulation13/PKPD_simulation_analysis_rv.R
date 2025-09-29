# シミュレーション結果の解析をする。

library(tidyverse)
library(data.table)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
datapath = file.path(basepath,"output_ke_eta1_exp","csv")
outpath = file.path(basepath,"output_ke_eta1_exp","csv")

Cmax_criteria = 401.439399
AUC_criteria = 9137.985073

# 24hに一回投与する群
# timing 22hr
# days 84
# dose 1,10,20,40,60,80,100,120,140,160,200,240,280,320,360

dose_list = c(1,10,20,40,60,80,100,120,140,160,200,240,280,320,360)
days = 84
timing = c(22)

# データ読み込む関数
pulldata = function(csv,dose){
  setwd(datapath)
  tibble_f1 = read_csv(csv,col_types=cols()) %>% 
    mutate(dose=dose)
  tibble_f1
}

tibble1 = list()
for(dose in dose_list){
  csvname = paste("rv_",timing,"h_",dose,"mg_",days,"d.csv",sep="")
  tibble_p = pulldata(csvname,dose)
  tibble1[[length(tibble1)+1]] = tibble_p
}

# それぞれのpercentについてCmaxを計算する。
calc_Cmax = function(data){
  tibble_f1 = data %>% 
    group_by(percent) %>% 
    nest() %>% 
    ungroup()
  tibble_f2 = NULL
  for(i in 1:nrow(tibble_f1)){
    percent = tibble_f1$percent[i]
    tibble_p1 = tibble_f1$data[[i]]
    Cmax = max(tibble_p1$C)
    Cmax_ratio = Cmax/Cmax_criteria
    s = paste("tibble_p2 = tibble('",percent,"'=Cmax,'",percent,"ratio'=Cmax_ratio)",sep="")
    eval(parse(text=s))
    tibble_f2 = bind_cols(tibble_f2,tibble_p2)
  }
  list(tibble_f2)
}

tibble1_dose = tibble1 %>% 
  bind_rows() %>% 
  group_by(dose) %>% 
  nest() %>% 
  ungroup()
tibble2 = tibble1_dose %>% 
  mutate(data = mapply(calc_Cmax,tibble1_dose$data))
tibble2 = tibble2 %>% 
  group_by(dose) %>% 
  unnest() %>% 
  ungroup()

colnames(tibble2) = c("dose","five","five_ratio","fifty","fifty_ratio","nf","nf_ratio","ave","ave_ratio")

# プロットしてみる。
ggplot()+
  geom_line(data=tibble2,aes(dose,fifty),linewidth=1)+
  geom_ribbon(data=tibble2,aes(x=dose,ymax=nf,ymin=five),alpha=0.2)+
  labs(x="dose(mg)",
       y="Concentration(ng/mL)")+
  scale_x_continuous(breaks=seq(0,400,20))+
  theme_gray()+
  geom_abline(intercept=Cmax_criteria,slope=0,color="red",linewidth=1)

# 割合でプロット
plot1 = ggplot()+
  geom_line(data=tibble2,aes(dose,fifty_ratio),linewidth=1)+
  geom_ribbon(data=tibble2,aes(x=dose,ymax=nf_ratio,ymin=five_ratio),alpha=0.2)+
  labs(x="dose(mg)",
       y="Cmax",
       title="")+
  scale_x_continuous(limits=c(0,240),breaks=seq(0,240,40))+
  theme_gray()+
  geom_abline(intercept=1,slope=0,color="orange",linewidth=1)+
  geom_abline(intercept=1.5,slope=0,color="red",linewidth=1)+
  scale_y_continuous(limits=c(0,2),breaks=seq(0,2,0.25))+
  geom_vline(xintercept=160,color="blue")
plot1

# 出力
setwd(outpath)
png("singledose_Cmax.png",res=300,width=1800,height = 1400)
plot1
dev.off()

# AUC_maxを計算
# 24hごとのAUCを計算
calc_AUC_max = function(data){
  tibble_f1 = data %>% 
    group_by(percent) %>% 
    nest() %>% 
    ungroup()
  tibble_f2 = NULL
  for(i in 1:nrow(tibble_f1)){
    percent=tibble_f1$percent[i]
    tibble_p1 = tibble_f1$data[[i]]
    time_list = seq(22,84*24,24)
    time_list=time_list[-length(time_list)]
    # それぞれ24hでAUCを計算する。
    AUC_list=c()
    for(k in time_list){
      time_range = c(k:(k+23))
      tibble_p2 = tibble_p1 %>% 
        filter(t%in%time_range)
      AUC = sum(tibble_p2$C)
      AUC_list=c(AUC_list,AUC)
    }
    # 後半の30個を平均する。
    AUC = mean(AUC_list[c((length(AUC_list)-30):length(AUC_list))])
    AUC_ratio = AUC/AUC_criteria
    s = paste("tibble_p3 = tibble('",percent,"'=AUC,'",percent,"ratio'=AUC_ratio)",sep="")
    eval(parse(text=s))
    tibble_f2 = bind_cols(tibble_f2,tibble_p3)
  }
  list(tibble_f2)
}

tibble_AUC = tibble1_dose %>% 
  mutate(data=mapply(calc_AUC_max,tibble1_dose$data))
tibble_AUC = tibble_AUC %>% 
  group_by(dose) %>% 
  unnest() %>% 
  ungroup()

colnames(tibble_AUC) = c("dose","five","five_ratio","fifty","fifty_ratio","nf","nf_ratio","ave","ave_ratio")

# プロットしてみる。
ggplot()+
  geom_line(data=tibble_AUC,aes(dose,fifty),linewidth=1)+
  geom_ribbon(data=tibble_AUC,aes(x=dose,ymax=nf,ymin=five),alpha=0.2)+
  labs(x="dose(mg)",
       y="Concentration(ng/mL)")+
  scale_x_continuous(breaks=seq(0,400,20))+
  theme_gray()+
  geom_abline(intercept=AUC_criteria,slope=0,color="red",linewidth=1)

# 割合でプロット
plot2=ggplot()+
  geom_line(data=tibble_AUC,aes(dose,fifty_ratio),linewidth=1)+
  geom_ribbon(data=tibble_AUC,aes(x=dose,ymax=nf_ratio,ymin=five_ratio),alpha=0.2)+
  labs(x="dose(mg)",
       y="AUC0-24h",
       title="")+
  scale_x_continuous(limits=c(0,240),breaks=seq(0,240,40))+
  theme_gray()+
  geom_abline(intercept=1,slope=0,color="orange",linewidth=1)+
  geom_abline(intercept=1.5,slope=0,color="red",linewidth=1)+
  scale_y_continuous(limits=c(0,2),breaks=seq(0,2,0.5))+
  geom_vline(xintercept=160,color="blue")
plot2

# 出力
setwd(outpath)
png("singledose_AUC.png",res=300,width=1100,height = 850)
plot2
dev.off()





# 12hに一回投与する群
# timing 9,19hr
# days 84
# dose 1,10,20,40,60,80,100,120,140,160,200,240,280,320,360

dose_list = c(1,10,20,40,60,80,100,120,140,160,200,240,280,320,360)
days = 84
timing = c(9,19)

# データ読み込む関数
pulldata = function(csv,dose){
  setwd(datapath)
  tibble_f1 = read_csv(csv,col_types=cols()) %>% 
    mutate(dose=dose)
  tibble_f1
}

tibble1 = list()
for(dose in dose_list){
  timing_p=paste(timing,sep="",collapse="_")
  dose_p = paste(dose,"_",dose,sep="")
  csvname = paste("rv_",timing_p,"h_",dose_p,"mg_",days,"d.csv",sep="")
  tibble_p = pulldata(csvname,dose)
  tibble1[[length(tibble1)+1]] = tibble_p
}


tibble1_dose = tibble1 %>% 
  bind_rows() %>% 
  group_by(dose) %>% 
  nest() %>% 
  ungroup()
tibble2 = tibble1_dose %>% 
  mutate(data = mapply(calc_Cmax,tibble1_dose$data))
tibble2 = tibble2 %>% 
  group_by(dose) %>% 
  unnest() %>% 
  ungroup()

colnames(tibble2) = c("dose","five","five_ratio","fifty","fifty_ratio","nf","nf_ratio","ave","ave_ratio")

# プロットしてみる。
ggplot()+
  geom_line(data=tibble2,aes(dose,fifty),linewidth=1)+
  geom_ribbon(data=tibble2,aes(x=dose,ymax=nf,ymin=five),alpha=0.2)+
  labs(x="dose(mg)",
       y="Concentration(ng/mL)")+
  scale_x_continuous(breaks=seq(0,400,20))+
  theme_gray()+
  geom_abline(intercept=Cmax_criteria,slope=0,color="orange",linewidth=1)+
  geom_abline(intercept=(Cmax_criteria*1.5),slope=0,color="red",linewidth=1)
  

# 割合でプロット
plot3=ggplot()+
  geom_line(data=tibble2,aes(dose,fifty_ratio),linewidth=1)+
  geom_ribbon(data=tibble2,aes(x=dose,ymax=nf_ratio,ymin=five_ratio),alpha=0.2)+
  labs(x="dose(mg)",
       y="Cmax",
       title="multidose")+
  scale_x_continuous(breaks=seq(0,400,20))+
  scale_x_continuous(limits=c(0,240),breaks=seq(0,240,40))+
  theme_gray()+
  geom_abline(intercept=1,slope=0,color="orange",linewidth=1)+
  geom_abline(intercept=1.5,slope=0,color="red",linewidth=1)+
  scale_y_continuous(limits=c(0,2),breaks=seq(0,2,0.25))+
  geom_vline(xintercept=100,color="blue")
plot3

# 出力
setwd(outpath)
png("multidose_Cmax.png",res=300,width=1800,height = 1400)
plot3
dev.off()

# AUC_maxを計算
# 24hごとのAUCを計

tibble_AUC = tibble1_dose %>% 
  mutate(data=mapply(calc_AUC_max,tibble1_dose$data))
tibble_AUC = tibble_AUC %>% 
  group_by(dose) %>% 
  unnest() %>% 
  ungroup()


colnames(tibble_AUC) = c("dose","five","five_ratio","fifty","fifty_ratio","nf","nf_ratio","ave","ave_ratio")

# プロットしてみる。
ggplot()+
  geom_line(data=tibble_AUC,aes(dose,fifty),linewidth=1)+
  geom_ribbon(data=tibble_AUC,aes(x=dose,ymax=nf,ymin=five),alpha=0.2)+
  labs(x="dose(mg)",
       y="Concentration(ng/mL)")+
  scale_x_continuous(breaks=seq(0,400,20))+
  theme_gray()+
  geom_abline(intercept=AUC_criteria,slope=0,color="orange",linewidth=1)+
  geom_abline(intercept=(AUC_criteria*1.5),slope=0,color="red",linewidth=1)

# 割合でプロット
plot4=ggplot()+
  geom_line(data=tibble_AUC,aes(dose,fifty_ratio),linewidth=1)+
  geom_ribbon(data=tibble_AUC,aes(x=dose,ymax=nf_ratio,ymin=five_ratio),alpha=0.2)+
  labs(x="dose(mg)",
       y="",
       title="")+
  scale_x_continuous(breaks=seq(0,400,20))+
  theme_gray()+
  scale_x_continuous(breaks=seq(0,400,20))+
  scale_x_continuous(limits=c(0,240),breaks=seq(0,240,40))+
  theme_gray()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_abline(intercept=1,slope=0,color="orange",linewidth=1)+
  geom_abline(intercept=1.5,slope=0,color="red",linewidth=1)+
  scale_y_continuous(limits=c(0,2),breaks=seq(0,2,0.5))+
  geom_vline(xintercept=100,color="blue")
plot4

# 出力
setwd(outpath)
png("multidose_AUC.png",res=300,width=1100,height = 850)
plot4
dev.off()
