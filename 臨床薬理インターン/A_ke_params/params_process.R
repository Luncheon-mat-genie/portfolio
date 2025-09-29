# 斎藤さんが算出してくれたPKの代表値を解析する。

library(tidyverse)

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
  
  if(file%in%c("A_Phase1b.csv","ke_Phase1b.csv")){
    tibble_f1  =tibble_f1 %>% 
      mutate(ID = (tibble_f1$ID+6))
  }
  
  list(tibble_f1)
}

tibble1 = tibble1 %>% 
  mutate(data = mapply(readdata,tibble1$csv))

# データを処理しやすいように変形する。

data_process = function(data){
  tibble_f1 = NULL
  for(i in 2:ncol(data)){
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
                        ncol=3)

gridExtra::grid.arrange(tibble1$dose_vs_value_plot[[1]][[2]],
                        tibble1$dose_vs_value_plot[[2]][[2]],
                        tibble1$dose_vs_value_plot[[3]][[2]],
                        tibble1$dose_vs_value_plot[[4]][[2]],
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
                        ncol=3)
dev.off()

# dose vs value log
png("dose_vs_value_log.png",width=2800,height=1600,res=300)
gridExtra::grid.arrange(tibble1$dose_vs_value_plot[[1]][[2]],
                        tibble1$dose_vs_value_plot[[2]][[2]],
                        tibble1$dose_vs_value_plot[[3]][[2]],
                        tibble1$dose_vs_value_plot[[4]][[2]],
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
for(id in 1:16){
  tibble_f1 = tibble3 %>% 
    filter(ID==id)
  for(dose_f in c(0.1,1,10,20,40,80,160,320)){
    tibble_f2 = tibble_f1 %>% 
      filter(dose==dose_f)
    if(nrow(tibble_f2)!=0){
      tibble_f3 = tibble(
        ID = id,
        dose=dose_f,
        A = filter(tibble_f2,category%in%c("A_Phase1a","A_Phase1b"))$value,
        ke = filter(tibble_f2,category%in%c("ke_Phase1a","ke_Phase1b"))$value
      )
      tibble4 = bind_rows(tibble4,tibble_f3)
      print(tibble4)
    }
  }
}

# csvで出力する。

setwd(basepath)
write.csv(tibble4,"PK_parameters.csv",row.names=FALSE)
