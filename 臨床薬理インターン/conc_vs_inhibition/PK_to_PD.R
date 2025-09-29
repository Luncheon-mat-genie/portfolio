# inhibition vs concの解析

library(tidyverse)
library(drc)
library(gridExtra)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath)
datapath = file.path(basepath,"..")
time_list = c(0,2,4,6,8,12,24,30,36,48,72,96,120,144,168)

# Individual_IL31inhibition.csvから解析
# fileの読み込み
setwd(datapath)
tibble1 = read_csv("Individual_IL31inhibition.csv",col_types=cols())

# doseごとにconc vs inhibitionでプロット

plot1 = ggplot(tibble1,aes(x=CP,y=IL31inh))+
  stat_summary(fun="mean",geom="point",size=3)+
  stat_summary(fun.data="mean_se",geom="errorbar",width=.1)+
  theme_classic()+
  facet_wrap(~DOSE)
plot1

# シグモイドカーブでフィッティングできそう
tibble2 = tibble(
  dose = c(80,160,"all")
)


fitting_drc = function(dose){
  if(dose=="all"){
    tibble_f1 = tibble1
  } else{  tibble_f1 = tibble1 %>% 
    filter(DOSE==dose)
  }
  
  model = drm(IL31inh~CP,data = tibble_f1,
              fct = LL.4(names = c("slope","lower_limit","Emax","EC50"),
                         fixed = c(NA,0,NA,NA)))
  list(model)
}

pull_coef = function(model){
  coef = model$coefficients
  result = data.frame(
    slope = coef[1],
    lower_limit = 0,
    Emax = coef[2],
    EC50 = coef[3]
  )
  list(result)
}

# 直接mutateにmapplyを置くとエラーが発生したのでtmpで仮置きしている。
tmp = mapply(fitting_drc,tibble2$dose)

tibble2 = tibble2 %>% 
  mutate(model = tmp)

tmp = mapply(pull_coef,tibble2$model)
tibble2 = tibble2 %>% 
  mutate(model_coef = tmp)

# それぞれのdoseでグラフを作成する。

make_fitting_plot = function(dose,coef){
  if(dose=="all"){
    tibble_f1 = tibble1
  } else{  tibble_f1 = tibble1 %>% 
    filter(DOSE==dose)
  }
  
  b = coef$slope
  c = coef$lower_limit
  d = coef$Emax
  e = coef$EC50
  
  plot_f1 = ggplot(tibble_f1,aes(x=CP,y=IL31inh))+
    stat_summary(fun="mean",geom="point",size=3)+
    stat_summary(fun.data="mean_se",geom="errorbar",width=.1)+
    theme_classic()+
    stat_function(fun = function(x) (d-c)/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
    labs(x = "concentration(ng/mL)",
         y = "inhibition rate",
         title = paste("dose ",dose," mg",sep=""),
         subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=",format(d,nsmall=1),", EC50=",format(e,nsmall=1),sep=""))
  plot_f1
  
  # conc 1000まで表示する。
  plot_f2 = ggplot(data=data.frame(X=c(0,750)),aes(x=X))+
    theme_classic()+
    stat_function(fun = function(x) (d-c)/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
    geom_point(data=tibble_f1,aes(CP,IL31inh))+
    labs(x = "concentration(ng/mL)",
         y = "inhibition rate",
         title = paste("dose ",dose," mg",sep=""),
         subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=",format(d,nsmall=1),", EC50=",format(e,nsmall=1),sep=""))
  plot_f2
  
  # listにして出力
  plot_list = list(
    plot_f1,plot_f2
  )
  list(plot_list)
}

# 実行

tibble2 = tibble2 %>% 
  mutate(plot = mapply(make_fitting_plot,tibble2$dose,tibble2$model_coef))

# grid.arrangeでまとめる。

gridExtra::grid.arrange(tibble2$plot[[1]][[1]],tibble2$plot[[1]][[2]],
                        tibble2$plot[[2]][[1]],tibble2$plot[[2]][[2]],
                        tibble2$plot[[3]][[1]],tibble2$plot[[3]][[2]],
                        nrow=3)

# グラフを出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
setwd(outpath)
png("conc_vs_inhibition.png",width=2400,height = 3000,res=300)
gridExtra::grid.arrange(tibble2$plot[[1]][[1]],tibble2$plot[[1]][[2]],
                        tibble2$plot[[2]][[1]],tibble2$plot[[2]][[2]],
                        tibble2$plot[[3]][[1]],tibble2$plot[[3]][[2]],
                        nrow=3)
dev.off()

# fittingで算出したパラメーターを出力する。
setwd(outpath)
write.csv(tibble2$model_coef[[3]],"PKPD_coef.csv",row.names = FALSE)


output_plot = tibble2$plot[[3]][[2]]
output_plot

setwd(outpath)
png("conc_to_inh.png",res=300,width=1000,height=800)
print(output_plot)
dev.off()



# 残差を計算してtibble1に加える。
calc_residual = function(model,data){
  residual = residuals(model)
  list(residual)
}

tibble2 = tibble2 %>% 
  mutate(residuals = mapply(calc_residual,tibble2$model,list(tibble1)))

# 自由度修正済み決定係数を計算する。
r2_coef<- function(dose,data,residuals){
  if(dose=="all"){
    tibble_f1 = data
  } else{
    tibble_f1 = data %>% 
      filter(DOSE==dose)
  }
  
 
  RSS = sum(residuals^2)
  ave = mean(tibble_f1$IL31inh)
  hensa = tibble_f1$IL31inh-ave
  hensa_2 = sum(hensa^2)
  
  R_2 = 1-RSS/hensa_2
  
  N = nrow(tibble_f1)
  
  adjusted_r2 = 1-(RSS/(N-3-1))/(hensa_2/N-1)
  
  
}

# 残差のヒストグラムを作成する。
residual_histogram = function(residuals){
  tibble_f1 = tibble(
    residuals=residuals
  )
 plot = ggplot(tibble_f1,aes(y=residuals))+
   geom_boxplot()+
   theme_classic()+
   theme(
     axis.text.x=element_blank()
   )+
   scale_y_continuous(limits=c(-8e-06,1e-06))+
   labs(title="残差の分布",
        y="残差(実測値-モデルからの予測値)")
 list(plot)
}

tibble2 = tibble2 %>% 
  mutate(residual_plot = mapply(residual_histogram,tibble2$residuals))

setwd(outpath)
png("residuals.png",res=300,width=500,height = 800)
print(tibble2$residual_plot[[3]])
dev.off()
