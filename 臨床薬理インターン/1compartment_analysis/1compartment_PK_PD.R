# 単回投与のモデルから複数回投与の血中濃度の再現
# そこから、複数回答夜の阻害率の再現

library(tidyverse)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath,"figures")
datapath = file.path(basepath,"..","conc_vs_inhibition") # conc_inhibitonのcoef.csvがあるpath
t.obs = seq(0,168,1)

# PK_parameters.csvを読み込む

setwd(basepath)
tibble1 = read_csv("PK_parameters.csv",col_types=cols()) %>% 
  filter(dose %in% c(80,160))
tibble2 = tibble1 %>% 
  group_by(dose) %>% 
  summarise(A=mean(A),
            ka=mean(ka),
            ke=mean(ke)) %>% 
  ungroup()

# 関数を結合して複数回投与の式を出す

multiple_model = function(A,ka,ke){
  # 関数を結合して複数回投与の式を出す
  for(i in 1:7){
    lag = (i-1)*24
    
    s = paste("y_p",i," = function(t){
              if((t-",lag,")<0){
              result=0
              return(result)
              } else{
              result=",A,"*(exp(-",ke,"*(t-",lag,"))-exp(-",ka,"*(t-",lag,")))
              return(result)
              }
              } ",sep="")
    eval(parse(text=s))
  }
  
  y = function(t){
    result = y_p1(t)+y_p2(t)+y_p3(t)+y_p4(t)+y_p5(t)+y_p6(t)+y_p7(t)
    return(result)
  }
  list(y)
}

tibble2 = tibble2 %>% 
  mutate(multiple_model = mapply(multiple_model,tibble2$A,tibble2$ka,tibble2$ke))

# 数値データを計算する。

multiple_value = function(multiple_model){
  tibble_f1 = tibble(
    t = t.obs,
    C = as.numeric(unlist(mapply(multiple_model,t.obs)))
  )
  list(tibble_f1)
}

tibble2 = tibble2 %>% 
  mutate(multiple_value = mapply(multiple_value,tibble2$multiple_model))

# 複数回投与のプロットを作成する。
PK_multiple_plot = function(dose,multiple_value){
  
  plot = ggplot(multiple_value,aes(x=t,y=C))+
    geom_line()+
    theme_classic()+
    labs(title = paste("DOSE",dose,"mg"))
  list(plot)
}

tibble2 = tibble2 %>% 
  mutate(multiple_plot = mapply(PK_multiple_plot,tibble2$dose,tibble2$multiple_value))

# 複数回投与のデータを読み込む

setwd(file.path(basepath,".."))
tibble3 = read_csv("Individual_concentration_multiple_dose.csv",col_types=cols()) %>% 
  group_by(DOSE) %>% 
  nest()

tibble2 = tibble2 %>% 
  mutate(raw_data = tibble3$data)

# plot

PK_multiple_plot  = function(dose,multiple_value,raw_data){
  
  plot = ggplot(NULL)+
    geom_point(data=raw_data,aes(TIME,CP))+
    stat_summary(data=raw_data,aes(TIME,CP,color=ID),fun="mean",geom="point",size=4)+
    geom_line(data=multiple_value,aes(t,C))+
    theme_classic()
  plot
}

# あわんやないかい。斎藤さんのデータでやろう。
