# 血中濃度の複数回投与から阻害率の複数回答夜を再現する。

library(tidyverse)
library(drc)
library(gridExtra)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath,"output")
datapath = file.path(basepath,"multi_dose_conc")

# 阻害率計算のためのフィッティングパラメーター読み込み
setwd(file.path(basepath,"..","conc_vs_inhibition"))
tibble_coef = read_csv("PKPD_coef.csv",col_types=cols())
b = tibble_coef$slope[1]
c = tibble_coef$lower_limit[1]
d = tibble_coef$Emax[1]
e = tibble_coef$EC50[1]

conc_to_inh = function(conc){
  y = c+(d-c)/(1+exp(b*(log(conc)-log(e))))
  return(y)
}

# データの読み込み

setwd(datapath)
csv_list = list.files(pattern="csv")

pulldata = function(csv){
  tibble_f1 = read_csv(csv,col_types=cols())
  tibble_f1 = tibble_f1 %>% 
    mutate(inh = as.numeric(unlist(mapply(conc_to_inh,tibble_f1$CP))))
  colnames(tibble_f1) = c("TIME_min","ID","TIME","CP","DOSE","inh")
  list(tibble_f1)
}

tibble1 = mapply(pulldata,csv_list) %>% 
  bind_rows() 

# 出力
setwd(outpath)
write.csv(tibble1,"multi_conc_inh.csv",row.names=FALSE)


# 血中濃度と阻害率の複数回投与のデータを読み込む

setwd(file.path(basepath,".."))
tibble_multi_conc = read_csv("Individual_concentration_multiple_dose.csv",col_types=cols())
tibble_multi_inh = read_csv("Individual_IL31inhibition.csv",col_types=cols())

# VAS50の解析結果を読み込む
setwd(file.path(basepath,"..","inhibition_vs_endpoint"))
tibble_multi_VAS50 = read_csv("endpoint.csv",col_types=cols())%>% 
  arrange(ID,DOSE,TIME)
  

# プロットを作成する。まずは個人で
plot_1 = ggplot()+
  geom_line(data=tibble1,aes(TIME,CP),color="blue")+
  geom_point(data=tibble_multi_conc,aes(TIME,CP),color="blue",size=4)+
  theme_classic()+
  facet_grid(ID~DOSE)
plot_1
    
plot_2 = ggplot()+
  geom_line(data=tibble1,aes(TIME,inh),color="red")+
  geom_point(data=tibble_multi_inh,aes(TIME,IL31inh),color="red",size=4)+
  theme_classic()+
  facet_grid(ID~DOSE)
plot_2

plot_3 = ggplot(tibble_multi_VAS50,aes(TIME,ratio))+
  stat_summary(fun="mean",geom="bar")+
  theme_light()+
  facet_grid(ID~DOSE)+
  scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=seq(0,1,0.1))
plot_3

# 阻害率と改善率を経時的に見てみる。
# plotの都合上ここだけVAS50の時間を-6時間する。
tibble_multi_VAS50_2 = tibble_multi_VAS50 %>% 
  mutate(TIME=tibble_multi_VAS50$TIME-6)
plot_4 =  ggplot()+
  geom_line(data=tibble1,aes(TIME,inh),color="red")+
  geom_point(data=tibble_multi_inh,aes(TIME,IL31inh),color="red",size=4)+
  theme_classic()+
  facet_grid(ID~DOSE)+
  stat_summary(data=tibble_multi_VAS50_2,aes(TIME,ratio),fun="mean",geom="bar")
plot_4

setwd(outpath)
png("inh_imp_ID.png",res=300,width=3000,height = 5000)
plot_4
dev.off()

  
  

# inhibition rateの12hごとの平均を計算する。
tibble2 = tibble1 %>% 
  group_by(ID,DOSE) %>% 
  nest() %>% 
  ungroup()

inh_rate_ave = function(data){
  
  tibble_f3 = NULL
  for(i in 1:14){
    start = (i-1)*720+1
    end = i*720
    range = c(start:end)
    tibble_f1 = data[range,] 
    inh_ave = mean(tibble_f1$inh)
    tibble_f2 = tibble(
      TIME = (end/60),
      inh_ave=inh_ave
    )
    tibble_f3 = bind_rows(tibble_f3,tibble_f2)
  }
  list(tibble_f3)
}

tibble2 = tibble2 %>% 
  mutate(ave_inh = mapply(inh_rate_ave,tibble2$data))
tibble2 = tibble2 %>% 
  dplyr::select(-data) %>% 
  group_by(ID,DOSE) %>% 
  unnest() %>% 
  ungroup() %>% 
  arrange(ID,DOSE,TIME)

# 阻害率の連続プロットと12hごとの平均をプロット
# 阻害率と改善率を経時的に見てみる。
plot_5 =  ggplot()+
  geom_line(data=tibble1,aes(TIME,inh),color="red")+
  geom_point(data=tibble_multi_inh,aes(TIME,IL31inh),color="red",size=4)+
  theme_classic()+
  facet_grid(ID~DOSE)+
  stat_summary(data=tibble2,aes(TIME,inh_ave),fun="mean",geom="bar",alpha=0.3)
plot_5


# tibble2にVAS50追加

tibble2 = tibble2 %>% 
  mutate(VAS50_rate = tibble_multi_VAS50$ratio)

# 出力
setwd(outpath)
write.csv(tibble2,"inh_VAS50.csv",row.names=FALSE)

# inh vs vas50でプロット

plot_4 = ggplot()+
  geom_point(data=tibble2,aes(inh_ave,VAS50_rate),size=2)+
  theme_classic()+
  facet_wrap(~ID)
plot_4


setwd(outpath)
png("inh_imp_12have.png",res=300,width=3000,height = 2000)
plot_4
dev.off()



# よくわからんので合わせてプロット

plot_5 = ggplot(data=tibble2,aes(inh_ave,VAS50_rate))+
  stat_summary(fun="mean",geom="point")+
  stat_summary(fun.data="mean_se",geom="errorbar")+
  theme_classic()+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.1))
plot_5

setwd(outpath)
png("inh_imp_12have_all.png",res=300,width=1200,height = 1000)
plot_5
dev.off()

# 全体でフィッティング

fitting_drc = function(data,dose){
  if(dose=="all"){
    tibble_f1 = data
  } else{  tibble_f1 = data %>% 
    filter(DOSE==dose)
  }
  
  model = drm(VAS50_rate~inh_ave,data = tibble_f1,
              fct = LL.4(names = c("slope","lower_limit","Emax","EC50"),
                         fixed = c(NA,NA,NA,NA)))
  coef = model$coefficients
  result = data.frame(
    slope = coef[1],
    lower_limit = coef[2],
    Emax = coef[3],
    EC50 = coef[4]
  )
  
  # 残差のdataとプロットも返すようにする。
  tibble_f1 = data %>% 
    mutate(residuals=residuals(model))
  
  plot = ggplot(tibble_f1,aes(inh_ave,residuals))+
    geom_point()+
    theme_minimal()
  plot
  
  output_list = list(
    result,tibble_f1,plot
  )
  
  list(output_list)
}

# 直接mutateにmapplyを置くとエラーが発生したのでtmpで仮置きしている。
tmp = mapply(fitting_drc,list(tibble2),"all")

tibble3 = tibble(
  dose = "all",
  model_result = tmp
)

# それぞれのdoseでグラフを作成する。

make_fitting_plot = function(data,dose,coef){
  if(dose=="all"){
    tibble_f1 = data
  } else{  tibble_f1 = data %>% 
    filter(DOSE==dose)
  }
  
  b = coef$slope
  c = coef$lower_limit
  d = coef$Emax
  e = coef$EC50
  
  # inhibition rate 1 まで表示する。
  plot_f2 = ggplot(data=data.frame(X=c(0,1)),aes(x=X))+
    theme_light()+
    geom_point(data=tibble_f1,aes(x=inh_ave,y=VAS50_rate),size=2)+
    stat_function(fun = function(x) (d-c)/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
    labs(x = "inhibition rate",
         y = "improvement rate",
         title = paste("dose ",dose," mg",sep=""),
         subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=",format(d,nsmall=1),", EC50=",format(e,nsmall=1),sep=""))+
    scale_y_continuous(limits=c(0,1))
  plot_f2
  
  # listにして出力
  list(plot_f2)
}

# 実行

tibble3 = tibble3 %>% 
  mutate(plot = mapply(make_fitting_plot,list(tibble2),"all",list(tibble3$model_result[[1]][[1]])))


# Emaxを指定しなかったらEmaxがぶっとんだのでEmaxを1に指定してもう一度フィッティングを行う。
fitting_drc_fixEmax = function(data,dose){
  if(dose=="all"){
    tibble_f1 = data
  } else{  tibble_f1 = data %>% 
    filter(DOSE==dose)
  }
  
  model = drm(VAS50_rate~inh_ave,data = tibble_f1,
              fct = LL.4(names = c("slope","lower_limit","Emax","EC50"),
                         fixed = c(NA,0,1,NA)))
  summary = summary(model)
  
  coef = model$coefficients
  result = data.frame(
    slope = coef[1],
    lower_limit = 0,
    Emax = 1,
    EC50 = coef[2],
    slope_25 = confint(model)[1],
    slope_975 = confint(model)[4],
    lower_limit_25 = confint(model)[2],
    lower_limit_975 = confint(model)[5],
    EC50_25 = confint(model)[3],
    EC50_975 = confint(model)[6]
    
  )
  list(result)
  
  # 残差のdataとプロットも返すようにする。
  tibble_f1 = data %>% 
    mutate(residuals=residuals(model))
  
  plot = ggplot(tibble_f1,aes(inh_ave,residuals))+
    geom_point()+
    theme_minimal()
  plot
  
  output_list = list(
    result,tibble_f1,plot
  )
  
  list(output_list)
  
}

# 直接mutateにmapplyを置くとエラーが発生したのでtmpで仮置きしている。
tmp = mapply(fitting_drc_fixEmax,list(tibble2),"all")

tibble3 = tibble3 %>% 
  mutate(model_result_fixEmax = tmp)

tibble3 = tibble3 %>% 
  mutate(plot_fix = mapply(make_fitting_plot,list(tibble2),"all",list(tibble3$model_result_fixEmax[[1]][[1]])))

# grid.arrangeでまとめる。

gridExtra::grid.arrange(tibble3$plot[[1]],tibble3$plot_fix[[1]],
                        tibble3$model_result[[1]][[3]],tibble3$model_result_fixEmax[[1]][[3]],
                        nrow=2)



# 出力
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive=TRUE)
}
setwd(outpath)

png("improvement_fitting.png",width=2800,height=2600,res=300)
gridExtra::grid.arrange(tibble3$plot[[1]],tibble3$plot_fix[[1]],
                        tibble3$model_result[[1]][[3]],tibble3$model_result_fixEmax[[1]][[3]],
                        nrow=2)
dev.off()


# 非線形混合効果モデルにより個人間変動を考慮する。

library(nlme)

# b=slope,c=0,d=1,e=EC50
# fun=1/(1+exp(b*(log(x)-log(e))))


model=nlme(VAS50_rate~1/(1+exp(b*(log(inh_ave)-log(e)))),
           fixed=list(b|e~1),
           random=list(ID=pdDiag(list(b|e~1))),
           start=c(b=-6,e=0.5),
           data=tibble2)
residuals=residuals(model)
tibble2=tibble2 %>% 
  mutate(residuals=residuals)

# 固定効果部分のプロットとランダム効果を含めたプロットをなんとかしたい。
color_list=c("red","black","blue","yellow","purple"
             ,"green4","grey50","skyblue3","salmon4","orchid3")

# 固定効果のプロット
b = model$coefficients$fixed[1]
e = model$coefficients$fixed[2]
plot_fixed = ggplot(tibble2,aes(inh_ave,VAS50_rate))+
  geom_point()+
  stat_function(fun = function(x) 1/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
  theme_minimal()+
  labs(x = "inhibition rate",
       y = "improvement rate",
       title = "fixed model",
       subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=1, EC50=",format(e,nsmall=1),sep=""))+
  scale_y_continuous(limits=c(0,1))
plot_fixed

# ランダム効果のプロット
b_list = model$coefficients$random$ID[,1]+b
e_list = model$coefficients$random$ID[,2]+e
s = paste("stat_function(fun = function(x) 1/(1+exp(",b_list,"*(log(x)-log(",e_list,")))),color='",color_list,"',linewidth=1)",
          sep="",collapse="+")
tibble2_2 = tibble2 %>% 
  mutate(ID = as.character(tibble2$ID))
plot_randomed = ggplot(tibble2_2,aes(inh_ave,VAS50_rate,color=ID))+
  geom_point(size=3)+
  theme_minimal()+
  labs(x = "inhibition rate",
       y = "improvement rate",
       title = "非線形混合効果モデル",
       subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=1, EC50=",format(e,nsmall=1),sep=""))+
  scale_y_continuous(limits=c(0,1))+
  scale_color_manual(values=color_list)
plot_randomed
ss = paste("plot_randomed = plot_randomed+",s,sep="")
eval(parse(text = ss))
plot_randomed

# 残差のプロット
plot_residuals = ggplot(tibble2,aes(inh_ave,residuals))+
  geom_point()+
  theme_minimal()+
  labs(title="residuals")
plot_residuals

# 残差のsdをinh_ave 0.1ずつ計算

inh_criteria=seq(0.1,1,0.1)
sd_inh = function(data,criteria){
  start=criteria-0.1
  end=criteria
  tibble_f1 = data %>% 
    filter(inh_ave>start) %>% 
    filter(inh_ave<end)
  if(nrow(tibble_f1)>0){
    sd=sd(tibble_f1$residuals)
  } else{
    sd=NA
  }
  return(sd)
}
tibble_residual_sd = tibble(
  criteria=inh_criteria,
  data=as.numeric(unlist(mapply(sd_inh,list(tibble2),inh_criteria)))
)
plot_residual_sd = ggplot(tibble_residual_sd,aes(criteria,data))+
  geom_point(size=4)+
  theme_minimal()+
  labs(x="inhibition rate",
       y = "residuals_sd",
       title = "residuals SD")
plot_residual_sd

gridExtra::grid.arrange(plot_fixed,plot_randomed,
                        plot_residuals,plot_residual_sd,
                        nrow=2)



setwd(outpath)

png("NMEM_fitting.png",width=2800,height=2600,res=300)
gridExtra::grid.arrange(plot_fixed,plot_randomed,
                        plot_residuals,plot_residual_sd,
                        nrow=2)
dev.off()




# 残差のSDを線形回帰して一般化する。
# inhが0.2以上のデータを利用して線形回帰する。
tibble_residual_lm = tibble_residual_sd %>% 
  filter(criteria>=0.2)
lm_result = lm(data~criteria,tibble_residual_lm)
summary = summary(lm_result)

# 阻害率に応じて残差のSDを計算する関数
b = lm_result$coefficients[1]
a = lm_result$coefficients[2]
residuals_SD = function(x){ # xは阻害率
  y = a*x+b
  if(y<0){
    y=0
  }
  return(y)
}

# 残差のSDについてプロットを作成する。

plot_residual_sd_2 = ggplot(tibble_residual_sd,aes(criteria,data))+
  geom_point(size=4)+
  theme_minimal()+
  labs(x="inhibition rate",
       y = "residuals_sd",
       title = "residuals SD",
       subtitle = paste("y = ",format(a,digits=3),"x+",format(b,digits=3)," R^2=",format(summary$r.squared,digits=2),sep=""))+
  stat_function(fun=function(x) a*x+b,color="blue", linewidth=1)+
  geom_abline(intercept=0,slope=0,color="red",linewidth=1)
plot_residual_sd_2


setwd(outpath)
png("residuals_SD_lm.png",width=1000,height=800,res=200)
plot_residual_sd_2
dev.off()


