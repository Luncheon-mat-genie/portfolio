# inhibition vs endpoint(VAS score)の解析
# VASの解析方法として、ある時間の10人いる患者のTRUE割合を出したら1時間おきのデータがとれるのではないかと考えた。


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

# DOSE,TIMEごとにグループ化して時間ごとのflagの割合(10人の患者の中での）を計算する。

tibble3 = tibble2 %>% 
  group_by(DOSE,TIME) %>% 
  nest() %>% 
  ungroup()

flag_percent = function(data){
  sum = sum(data$flag)
  ratio = sum/10
  return(ratio)
}


tibble3 = tibble3 %>% 
  mutate(ratio = as.numeric(unlist(mapply(flag_percent,tibble3$data))))

# 統合する

tibble4 = tibble3 %>% 
  arrange(TIME,DOSE)

# tibble4をtimeごとにプロットしてみる。

plot1 = ggplot(tibble4,aes(TIME,ratio))+
  geom_point()+
  geom_line()+
  geom_point(size=2)+
  theme_light()+
  facet_wrap(~DOSE)+
  scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=seq(0,1,0.1))
plot1

# 複数回投与時の阻害率の予測値をもってくる。
datapath_multi = file.path(basepath,"..","conc_to_inhibiton","output")
setwd(datapath_multi)
tibble5 = read_csv("multi_conc_inh.csv",col_types=cols())

# 阻害率と改善率を合わせてプロット

plot2 = ggplot()+
  geom_point(data=tibble4,aes(TIME,ratio))+
  geom_line(data=tibble4,aes(TIME,ratio))+
  stat_summary(data=tibble5,aes(TIME,inh),fun="mean",geom="line",color="red")+
  stat_summary(data=tibble5,aes(TIME,inh),fun.data="mean_se",geom="errorbar",color="red")+
  facet_wrap(~DOSE)+
  theme_light()
plot2

# 1h時点での阻害率を抜きだすパターンと0~1hの阻害率の平均を計算するパターンを試す。

# 1h時点の阻害率を抜き出す
tibble6 = tibble5 %>% 
  filter(TIME_min%in%seq(0,10079,60))
tibble7 = tibble6 %>% 
  group_by(TIME,DOSE) %>% 
  nest() %>% 
  ungroup() %>% 
  arrange(TIME,DOSE)
tibble8 = tibble4 %>% 
  dplyr::select(-data) %>% 
  filter(!(TIME==168))%>% 
  arrange(TIME,DOSE)
tibble8 = tibble8 %>% 
  dplyr::select(ratio)
tibble8 = bind_cols(tibble7,tibble8) %>% 
  group_by(TIME,DOSE,ratio) %>% 
  unnest() %>% 
  ungroup()

plot3 = ggplot(tibble8,aes(inh,ratio))+
  geom_point()+
  theme_light()
plot3


# 全体でフィッティング

fitting_drc = function(data,dose){
  if(dose=="all"){
    tibble_f1 = data
  } else{  tibble_f1 = data %>% 
    filter(DOSE==dose)
  }
  
  model = drm(ratio~inh,data = tibble_f1,
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
tmp = mapply(fitting_drc,list(tibble8),"all")

tibble9 = tibble(
  dose = "all",
  model_coef = tmp
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
    geom_point(data=tibble_f1,aes(x=inh,y=ratio),size=2)+
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

tibble9 = tibble9 %>% 
  mutate(plot = mapply(make_fitting_plot,list(tibble8),"all",tibble3$model_coef))

# ばらつきが大きすぎるし意味合い的にも歪になっている。

# 1hごとに平均して考えてみる。

tibble10 = tibble5 %>% 
  group_by(DOSE) %>% 
  nest() %>% 
  ungroup()
inh_ave_1h = function(data){
  time_list = seq(60,10080,60)
  tibble_f3 = NULL
  for(i in time_list){
    start = i-60
    end=i
    tibble_f1 = data %>% 
      filter(TIME_min>=start) %>% 
      filter(TIME_min<=end)
    mean = mean(tibble_f1$inh)
    tibble_f2 = tibble(
      TIME = i/60,
      inh_ave = mean
    )
    tibble_f3 = bind_rows(tibble_f3,tibble_f2)
  }
  list(tibble_f3)
}
tibble11 = tibble10 %>% 
  mutate(data = mapply(inh_ave_1h,tibble10$data))
tibble11 = tibble11 %>% 
  group_by(DOSE) %>% 
  unnest() %>% 
  ungroup() %>% 
  arrange(TIME,DOSE)
tibble12= tibble4 %>% 
  filter(!(TIME==0)) 
tibble11 = tibble11 %>% 
  mutate(imp_ratio = tibble12$ratio)

plot4 = ggplot(tibble11,aes(inh_ave,imp_ratio))+
  geom_point()+
  theme_light()
plot4
