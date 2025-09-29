# 斎藤さんが算出してくれたPKの代表値を解析する。

library(tidyverse)
library(drc)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath,"figures")

# すべてのファイルの読み込み
setwd(basepath)
csv_list = list.files(pattern="csv")

tibble1 = tibble(
  csv = csv_list
)

readdata = function(file){
  tibble_f1 = read_csv(file,col_types=cols())
  list(tibble_f1)
}

tibble1 = tibble1 %>% 
  mutate(data = mapply(readdata,tibble1$csv))

# データを処理しやすいように変形する。

data_process = function(data){
  tibble_f1 = NULL
  for(i in 2:9){
    colname = colnames(data)[i]
    tibble_f2 = data %>% 
      dplyr::select(c(ID,colname))
    colnames(tibble_f2) = c("ID","value")
    tibble_f2 = tibble_f2 %>% 
      mutate(dose = as.numeric(colname))
    tibble_f1 = bind_rows(tibble_f1,tibble_f2)
  }
  list(tibble_f1)
}

tibble1 = tibble1 %>% 
  mutate(data_process = mapply(data_process,tibble1$data))

# すべてのパラメーターについて、dose vs valueでプロットしてみる。

dose_value_plot = function(file,data){
  plot = ggplot(data,aes(dose,value))+
    geom_point()+
    stat_summary(fun="mean",geom="point",size=4,color="red")+
    theme_classic()+
    labs(title=file)
  plot2 = ggplot(data,aes(dose,value))+
    geom_point()+
    stat_summary(fun="mean",geom="point",size=4,color="red")+
    theme_classic()+
    labs(title=file)+
    scale_x_log10()
  plot_list = list(plot,plot2)
  list(plot_list)
}

tibble1 = tibble1 %>% 
  mutate(dose_vs_value_plot = mapply(dose_value_plot,
                                     tibble1$csv,tibble1$data_process))

gridExtra::grid.arrange(tibble1$dose_vs_value_plot[[1]][[1]],
                        tibble1$dose_vs_value_plot[[2]][[1]],
                        tibble1$dose_vs_value_plot[[3]][[1]],
                        tibble1$dose_vs_value_plot[[4]][[1]],
                        tibble1$dose_vs_value_plot[[5]][[1]],
                        tibble1$dose_vs_value_plot[[6]][[1]],
                        ncol=3)

gridExtra::grid.arrange(tibble1$dose_vs_value_plot[[1]][[2]],
                        tibble1$dose_vs_value_plot[[2]][[2]],
                        tibble1$dose_vs_value_plot[[3]][[2]],
                        tibble1$dose_vs_value_plot[[4]][[2]],
                        tibble1$dose_vs_value_plot[[5]][[2]],
                        tibble1$dose_vs_value_plot[[6]][[2]],
                        ncol=3)

# AとAUCとCmaxをフィッティングしてみる。

fitting_list = c("A_Phase1a.csv","AUC_Phase1a.csv","Cmax_Phase1a.csv")

tibble2 = tibble1 %>% 
  filter(csv %in% fitting_list)

fitting_drc = function(data){
  tibble_f1 = data
  
  model = drm(value~dose,data = tibble_f1,
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

tmp = mapply(fitting_drc,tibble2$data_process)

tibble2 = tibble2 %>% 
  mutate(model_coef = tmp)


# グラフを作成する。

make_fitting_plot = function(csv,data,coef){
  tibble_f1 = data
  
  b = coef$slope
  c = coef$lower_limit
  d = coef$Emax
  e = coef$EC50
  
  plot_f1 = ggplot(tibble_f1,aes(x=dose,y=value))+
    stat_summary(fun="mean",geom="point",size=3)+
    geom_point(size=1.5,color="blue")+
    stat_summary(fun.data="mean_se",geom="errorbar",width=.1)+
    theme_classic()+
    scale_x_log10()+
    stat_function(fun = function(x) c+(d-c)/(1+exp(b*(log(x)-log(e)))),color="red",linewidth=1)+
    labs(x = "log10(dose)",
         y = "value",
         title = csv,
         subtitle = paste("slope=",format(b,nsmall=1),",lower_limit=0,
Emax=",format(d,nsmall=1),", EC50=",format(e,nsmall=1),sep=""))
  plot_f1
  
  list(plot_f1)
}

# 実行

tibble2 = tibble2 %>% 
  mutate(plot = mapply(make_fitting_plot,
                       tibble2$csv,tibble2$data_process,tibble2$model_coef))

# grid.arrangeでまとめる。

gridExtra::grid.arrange(tibble2$plot[[1]],tibble2$plot[[2]],
                        tibble2$plot[[3]],
                        ncol=3)


# グラフを出力
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
setwd(outpath)

# dose vs value 線形scale
png("dose_vs_value.png",width=2800,height=1600,res=300)
gridExtra::grid.arrange(tibble1$dose_vs_value_plot[[1]][[1]],
                        tibble1$dose_vs_value_plot[[2]][[1]],
                        tibble1$dose_vs_value_plot[[3]][[1]],
                        tibble1$dose_vs_value_plot[[4]][[1]],
                        tibble1$dose_vs_value_plot[[5]][[1]],
                        tibble1$dose_vs_value_plot[[6]][[1]],
                        ncol=3)
dev.off()

# dose vs value log
png("dose_vs_value_log.png",width=2800,height=1600,res=300)
gridExtra::grid.arrange(tibble1$dose_vs_value_plot[[1]][[2]],
                        tibble1$dose_vs_value_plot[[2]][[2]],
                        tibble1$dose_vs_value_plot[[3]][[2]],
                        tibble1$dose_vs_value_plot[[4]][[2]],
                        tibble1$dose_vs_value_plot[[5]][[2]],
                        tibble1$dose_vs_value_plot[[6]][[2]],
                        ncol=3)
dev.off()

# fittingしたプロット
png("dose_vs_value_fitting.png",width=3200,height=900,res=300)
gridExtra::grid.arrange(tibble2$plot[[1]],tibble2$plot[[2]],
                        tibble2$plot[[3]],
                        ncol=3)
dev.off()


# ID, value, doseでtibbleを整理

data_process2 = function(csv,data_process){
  name = str_sub(csv,end=-5)
  data_process = data_process %>% 
    mutate(category=name)
  list(data_process)
}

tibble3 = mapply(data_process2,tibble1$csv,tibble1$data_process) %>% 
  bind_rows()

# 変形
tibble4 = NULL
for(id in 1:6){
  tibble_f1 = tibble3 %>% 
    filter(ID==id)
  for(dose_f in c(0.1,1,10,20,40,80,160,320)){
    tibble_f2 = tibble_f1 %>% 
      filter(dose==dose_f)
    tibble_f3 = tibble(
      ID = id,
      dose=dose_f,
      A = filter(tibble_f2,category=="A_Phase1a")$value,
      AUC = filter(tibble_f2,category=="AUC_Phase1a")$value,
      Cmax = filter(tibble_f2,category=="Cmax_Phase1a")$value,
      ka = filter(tibble_f2,category=="ka_Phase1a")$value,
      ke = filter(tibble_f2,category=="ke_Phase1a")$value,
      tmax = filter(tibble_f2,category=="tmax_Phase1a")$value,
    )
    tibble4 = bind_rows(tibble4,tibble_f3)
  }
}

# csvで出力する。

setwd(basepath)
write.csv(tibble4,"PK_parameters.csv",row.names=FALSE)
