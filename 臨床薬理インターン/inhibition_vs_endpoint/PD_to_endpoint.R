# inhibition vs endpoint(VAS score)の解析

library(tidyverse)
library(drc)
library(gridExtra)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath)
datapath = file.path(basepath,"..")
time_list = c(0,2,4,6,8,12,24,30,36,48,72,96,120,144,168)

# Individual_IL31inhibition.csvの読み込み
setwd(datapath)
tibble1 = read_csv("Individual_IL31inhibition.csv",col_types=cols())

# Individual_Endpoint.csvの読み込み
setwd(datapath)
tibble2 = read_csv("Individual_Endpoint.csv",col_types=cols())

# VAS scoreのFALSE=0,TRUE=1とする。

make_flag = function(VAS50){
  if(VAS50=="FALSE"){
    flag=0
  } else{
    flag=1
  }
  return(flag)
}

tibble2 = tibble2 %>% 
  mutate(flag = as.numeric(unlist(mapply(make_flag,tibble2$VAS50))))

# IDごとにグループ化して12hごとのflagの割合を計算する。

tibble3 = tibble2 %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup()

flag_percent = function(ID,data){
  tibble_f4 = NULL
  for(dose in c(80,160)){
    tibble_f1 = data %>% 
      filter(DOSE==dose) %>% 
      filter(!(TIME==0))
    tibble_f3 = NULL
    for(i in 1:(168/12)){
      tibble_f2 = tibble_f1[c((i*12-11):(i*12)),]
      sum = sum(tibble_f2$flag)
      ratio = sum/12
      tibble_p = tibble(
        TIME = i*12,
        DOSE = dose,
        ratio = ratio
      )
      tibble_f3 = bind_rows(tibble_f3,tibble_p)
    }
    tibble_f4 = bind_rows(tibble_f4,tibble_f3)
  }
  tibble_f4 = tibble_f4 %>% 
    mutate(ID = ID)
  list(tibble_f4)
}


tibble3 = tibble3 %>% 
  mutate(data_ratio = mapply(flag_percent,tibble3$ID,tibble3$data))

# 統合する

tibble4 = bind_rows(tibble3$data_ratio)

# tibble4をtimeごとにプロットしてみる。

plot1 = ggplot(tibble4,aes(TIME,ratio))+
  stat_summary(fun="mean",geom="bar")+
  stat_summary(fun.data="mean_se",geom="errorbar",width=5)+
  geom_point(size=2)+
  theme_light()+
  facet_wrap(~DOSE)+
  scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=seq(0,1,0.1))
plot1

# tibble4を出力する。
setwd(basepath)
write.csv(tibble4,"endpoint.csv",row.names=FALSE)

# tibble1のID,TIME,DOSEの情報から、該当する列をtibble4から抜きだし、ratioを追加する。

add_ratio = function(id,time,dose){
  ratio = tibble4 %>% 
    filter(ID==id&TIME==time&DOSE==dose)
  if(nrow(ratio)==0){
    ratio = NA
  } else{
    ratio = ratio$ratio
  }
  return(ratio)
}

tibble1 = tibble1 %>% 
  mutate(ratio = as.numeric(unlist(mapply(add_ratio,
                                            tibble1$ID,tibble1$TIME,tibble1$DOSE))))

tibble5 = tibble1 %>% 
  filter(!(is.na(ratio)))


# inhibition vs improvementでプロットしてみる。

plot2 = ggplot(tibble5,aes(IL31inh,ratio))+
  geom_point(size=3)+
  theme_classic()
plot2


# fittingする。
tibble6 = tibble(
  dose = c(80,160,"all")
)

fitting_drc = function(dose){
  if(dose=="all"){
    tibble_f1 = tibble5
  } else{  tibble_f1 = tibble5 %>% 
    filter(DOSE==dose)
  }
  
  model = drm(ratio~IL31inh,data = tibble_f1,
              fct = LL.4(names = c("slope","lower_limit","Emax","EC50"),
                         fixed = c(NA,NA,NA,NA)))
  coef = model$coefficients
  result = data.frame(
    slope = coef[1],
    lower_limit = coef[2],
    Emax = coef[3],
    EC50 = coef[4]
  )
  list(result)
}

# 直接mutateにmapplyを置くとエラーが発生したのでtmpで仮置きしている。
tmp = mapply(fitting_drc,tibble6$dose)

tibble6 = tibble6 %>% 
  mutate(model_coef = tmp)

# それぞれのdoseでグラフを作成する。

make_fitting_plot = function(dose,coef){
  if(dose=="all"){
    tibble_f1 = tibble5
  } else{  tibble_f1 = tibble5 %>% 
    filter(DOSE==dose)
  }
  
  b = coef$slope
  c = coef$lower_limit
  d = coef$Emax
  e = coef$EC50
  
  plot_f1 = ggplot(tibble_f1,aes(x=IL31inh,y=ratio))+
    stat_summary(fun="mean",geom="point",size=3)+
    stat_summary(fun.data="mean_se",geom="errorbar",width=.1)+
    theme_classic()+
    stat_function(fun = function(x) (d-c)/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
    labs(x = "inhibition rate",
         y = "improvement rate",
         title = paste("dose ",dose," mg",sep=""),
         subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=",format(d,nsmall=1),", EC50=",format(e,nsmall=1),sep=""))+
    scale_y_continuous(limits=c(0,1))
  plot_f1
  
  # inhibition rate 1 まで表示する。
  plot_f2 = ggplot(data=data.frame(X=c(0,1)),aes(x=X))+
    theme_classic()+
    geom_point(data=tibble_f1,aes(x=IL31inh,y=ratio))+
    stat_function(fun = function(x) (d-c)/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
    labs(x = "inhibition rate",
         y = "improvement rate",
         title = paste("dose ",dose," mg",sep=""),
         subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=",format(d,nsmall=1),", EC50=",format(e,nsmall=1),sep=""))+
    scale_y_continuous(limits=c(0,1))
  plot_f2
  
  # listにして出力
  plot_list = list(
    plot_f1,plot_f2
  )
  list(plot_list)
}

# 実行

tibble6 = tibble6 %>% 
  mutate(plot = mapply(make_fitting_plot,tibble6$dose,tibble6$model_coef))

# grid.arrangeでまとめる。

gridExtra::grid.arrange(tibble6$plot[[1]][[1]],tibble6$plot[[1]][[2]],
                        tibble6$plot[[2]][[1]],tibble6$plot[[2]][[2]],
                        tibble6$plot[[3]][[1]],tibble6$plot[[3]][[2]],
                        nrow=3)

# グラフを出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
setwd(outpath)

# 棒グラフ
png("improvement.png",width=2000,height = 1000,res=300)
plot1
dev.off()

# fittingしたプロット
png("conc_vs_inhibition.png",width=2400,height = 3000,res=300)
gridExtra::grid.arrange(tibble6$plot[[1]][[1]],tibble6$plot[[1]][[2]],
                        tibble6$plot[[2]][[1]],tibble6$plot[[2]][[2]],
                        tibble6$plot[[3]][[1]],tibble6$plot[[3]][[2]],
                        nrow=3)
dev.off()