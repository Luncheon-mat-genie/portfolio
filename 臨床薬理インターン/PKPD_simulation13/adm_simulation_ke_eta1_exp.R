# 薬物を投与した時の血漿中薬物濃度、阻害率、改善確率を計算する。

library(tidyverse)
library(gridExtra)
library(data.table)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
datapath = file.path(basepath,"output_ke_eta1_exp")
outpath = file.path(basepath,"output_ke_eta1_exp")


# 各種パラメータの取り出し
setwd(datapath)
tibble_params = read_csv("estimate_params.csv",col_types=cols())
param_list = colnames(tibble_params)
for(i in param_list){
  s = paste(i," = tibble_params$",i,sep="")
  eval(parse(text = s))
}

ka=0.5
b0=-6.834934188610949
b1=9.798149516668506
eta0_sd=0.170910378387596


n_simulations_inner=1 # 個人内変動を生成する回数 
n_simulations_inter=5000 # 個人間変動を生成する回数　基本的にこっちの回数を増やすことでサンプル数を増やす。
monte_gap = 1 # モンテカルロシミュレーションで予測区間を推定するときにどういう時間の刻みでデータを生成していくか



# 投与終了までにどのタイミング・投与量で投与するかをtibbleにする関数
admini_timing = function(timing,dose,days){
  days_h = days*24
  timing_list = c(timing)
  while(timing<days_h){
    timing = timing + 24
    timing_list = c(timing_list,timing)
  }
  timing_list=timing_list[-length(timing_list)]
  tibble_f1 = tibble(
    timing = timing_list,
    dose = dose
  )
  list(tibble_f1)
}


# 個人内変動のみを考慮した一人分のデータを再現する関数
make_C_all_id_last = function(timing_list,days,n){
  
  print(n_simulations_inter+1-n)
  
  # 個人間変動を1つ乱数生成する。
  ke_eta1 = rnorm(1,0,ke_eta1_sd)
  A_eta2 = rnorm(1,0,A_eta2_sd)
  A_eta4 = rnorm(1,0,A_eta4_sd)
  # A_eta4_eta1 = rnorm(1,0,A_eta4_eta1_sd)
  # A_eta4 = A_eta4_a*ke_eta1 + A_eta4_b + A_eta4_eta1
  
  # calc_C_all_id2をforループにする。
  tibble_f3 = list()
  for(i in 1:nrow(timing_list)){
    timing = timing_list$timing[i]
    dose = timing_list$dose[i]
    time_range = seq(timing,(days*24),monte_gap)
    time_range_calc = time_range-timing
    
    # calc_C_all_idをforループにする。
    tibble_f2 = list()
    for(k in 1:n_simulations_inner){
      # ある服用タイミング、投与量において個人内変動のみを考慮した血中濃度をシミュレーションする
        
      # keを乱数生成する。
      eta2 = rnorm(1,0,ke_eta2_sd)
      ke = ke_fixed * exp(ke_eta1) * exp(eta2)
      
      # Aを乱数生成する。
      A_eta3 = rnorm(1,0,A_eta3_sd)
      logA_rv = A_a*log10(dose)+ (A_b+A_eta2)
      c_rv = A_c + A_eta4
      if(logA_rv>c_rv){
        logA_rv = c_rv
      }
      logA = logA_rv + A_eta3
      A = exp(logA)
      
      # 血中濃度を計算する。
      C = A*(exp(-ke*time_range_calc)-exp(-ka*time_range_calc))
      tibble_f1 = tibble(
        t = time_range,
        C = C,
        number = c(1:n_simulations_inner)
      )
      tibble_f2[[length(tibble_f2)+1]] = tibble_f1
    }
    tibble_f2 = tibble_f2 %>% 
      bind_rows()
    tibble_f3[[length(tibble_f3)+1]] = tibble_f2
  }
  
  tibble_f3 = tibble_f3 %>% 
    bind_rows() %>% 
    group_by(t) %>% 
    summarise(C = sum(C))
  
  # 阻害率もけいさんする。IDもつけておく
  tibble_f3 = tibble_f3 %>% 
    mutate(inh = 1/(1+exp(-1.999967*(log(tibble_f3$C)-log(200.006)))),
           ID = n)
  list(tibble_f3)
}

# それぞれの患者さんについてひとつlogitを生成し、改善率を計算するスクリプト
calc_imp_p_id = function(data){
  # 阻害率から改善率を計算する。
  # 改善率は個人間変動のみを考える。
  # 斎藤さんが推定したパラメータを用いる。
  # logit = (b0 + eta0) + b1*inh
  # eta0は正規分布に従うと仮定する。
  # b0=-6.834934188610949,b1=9.798149516668506,eta0_sd=0.170910378387596
  
  # eta0を一つ乱数生成する。
  eta0 = rnorm(1,0,eta0_sd)
  logit = (b0 + eta0) + b1 * data$inh
  
  # 阻害率から改善確率をけいさんする
  tibble_f1 = data %>% 
    mutate(imp_p = 1/(1+exp(-logit)))
  list(tibble_f1)
}

# 5%,50%,95%,averageをそれぞれの時間で計算する
calc_percent = function(data,category){ # categoryはC or inh or imp_p
  n_all = n_simulations_inner*n_simulations_inter
  s = paste("C = data$",category," %>%  sort()",sep="")
  eval(parse(text=s))
  C_percent = C[c(n_all*0.05,n_all*0.5,n_all*0.95)]
  C_ave = mean(C)
  tibble_f1 = tibble(
    percent = c("5%","50%","95%","average"),
    value = c(C_percent,C_ave)
  )
  list(tibble_f1)
}

# 個人のデータをn個作成し、最終的な信頼区間を推定する
PKPD_simulation = function(timing,dose,days){
  
  timing_list = mapply(admini_timing,timing,dose,days) %>% 
    bind_rows() %>% 
    arrange(timing)
  
  print(timing_list)
  
  n_range = c(1:n_simulations_inter)
  tibble_Clast = mapply(make_C_all_id_last,list(timing_list),days,n_range)
  
  tibble_Clast2 = tibble_Clast %>% 
    bind_rows()
  tibble_Clast2 = tibble_Clast2 %>% 
    group_by(t) %>% 
    unnest() %>% 
    nest() %>% 
    ungroup()
  
  
  tibble_C = tibble_Clast2 %>% 
    mutate(data = mapply(calc_percent,tibble_Clast2$data,"C"))
  
  # 血中濃度のシミュレーション結果
  tibble_C = tibble_C %>% 
    group_by(t) %>% 
    unnest() %>% 
    ungroup()
  
  tibble_inh = tibble_Clast2 %>% 
    mutate(data = mapply(calc_percent,tibble_Clast2$data,"inh"))
  
  # 阻害率のシミュレーション結果
  tibble_inh = tibble_inh %>% 
    group_by(t) %>% 
    unnest() %>% 
    ungroup()
  
  # 阻害率から改善率を計算する。
  tibble_imp = mapply(calc_imp_p_id,tibble_Clast) %>% 
    bind_rows()
  
  tibble_imp = tibble_imp %>% 
    group_by(t) %>% 
    nest() %>% 
    ungroup()
  
  tibble_imp2 = tibble_imp %>% 
    mutate(data = mapply(calc_percent,tibble_imp$data,"imp_p"))
  
  # 改善率のシミュレーションデータ生成完了
  tibble_imp2 = tibble_imp2 %>% 
    group_by(t) %>% 
    unnest() %>% 
    ungroup()
  
  # 結果の出力
  outpath_csv = file.path(outpath,"csv")
  if(!(dir.exists(outpath_csv))){
    dir.create(outpath_csv)
  }
  setwd(outpath_csv)
  result = tibble_imp %>% 
    bind_rows() %>% 
    group_by(t) %>% 
    unnest()
  csvname_timing = paste(timing,sep="",collapse="_")
  csvname_dose = paste(dose,sep="",collapse="_")
  csvname = paste("result_",csvname_timing,"h_",csvname_dose,"mg_",days,"d.csv",sep="")
  fwrite(result,csvname)
  
  # 代表値の出力
  result_rv = tibble_C %>% 
    dplyr::rename("C"=value) %>% 
    mutate(inh = tibble_inh$value,
           imp = tibble_imp2$value)
  csvname_timing = paste(timing,sep="",collapse="_")
  csvname_dose = paste(dose,sep="",collapse="_")
  csvname = paste("rv_",csvname_timing,"h_",csvname_dose,"mg_",days,"d.csv",sep="")
  fwrite(result_rv,csvname)
  

  list(
    C = tibble_C,
    inh = tibble_inh,
    imp_p = tibble_imp2
  )
}


# PK_simulation関数を実行するとシミュレーションが実行されます。

# 22時に24時間ごとに1,10,20,40,60,80,100,120,140,160,200,240,280,320,360 mgをシミュレーションする。12*7日間
for(dose in c(160)){
  PKPD_simulation(22,dose,(12*7))
}

# 9時と19時に1,10,20,40,60,80,100,120,140,160,200,240,280,320,360 mgをシミュレーションする。 12*7日間
for(dose in c(70,110)){
  dose_list = c(dose,dose)
  PKPD_simulation(timing=c(7,19),dose_list,(12*7))
}


tibble_80 = PKPD_simulation(0,80,7)
tibble_160 = PKPD_simulation(0,160,7)


# PKPD_simulationで出力された結果をプロットにする関数
make_plot = function(data,category,dose){ # categoryはC or inh or imp_p 
  s = paste("tibble_f1 = data$",category,sep="")
  eval(parse(text = s))
  
  plot_f1 = ggplot()+
    geom_line(data=tibble_f1,aes(t,value,color=percent),linewidth=1)+
    labs(y=category)+
    labs(x = paste("time dose:",dose,"mg",sep=""))
  if(category%in%c("inh","imp_p")){
    plot_f1 = plot_f1 +
      scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))
  }
  plot_f1
}

# 3つのプロットをまとめて出力する関数
arrange_plot = function(data,dose){
  plot_C = make_plot(data,"C",dose)
  plot_inh = make_plot(data,"inh",dose)
  plot_imp_p = make_plot(data,"imp_p",dose)
  grid.arrange(plot_C,plot_inh,plot_imp_p)
}

plot_80 = arrange_plot(tibble_80,80)
plot_160 = arrange_plot(tibble_160,160)

# 出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
setwd(outpath)
png("simulation_80mg.png",res=300,width=2400,height = 2400)
arrange_plot(tibble_80,80)
dev.off()

png("simulation_160mg.png",res=300,width=2400,height = 2400)
arrange_plot(tibble_160,160)
dev.off()




# 血中濃度、阻害率のデータを読み込む
setwd(basepath)
tibble_c_data = read_csv("Individual_concentration_multiple_dose.csv",col_types=cols())
tibble_inh_data = read_csv("Individual_IL31inhibition.csv",col_types=cols())

# PKPD_simulationで出力された80mg, 160mgのプロットを作成する。
make_plot_80160 = function(category){ # categoryはC or inh or imp_p 
  s = paste("tibble_f1_80 = tibble_80$",category,sep="")
  eval(parse(text = s))
  s = paste("tibble_f1_160 = tibble_160$",category,sep="")
  eval(parse(text = s))
  tibble_f2 = bind_rows(tibble_f1_80 %>% mutate(DOSE=80),tibble_f1_160 %>% mutate(DOSE=160))
  
  plot_f1 = ggplot()+
    geom_line(data=tibble_f2,aes(t,value,color=percent),linewidth=1)+
    labs(y=category)+
    theme_grey()+
    labs(x = "time(hr)")+
    facet_wrap(~DOSE)
    scale_x_continuous(breaks=seq(0,168,24))
  if(category=="C"){
    plot_f1 = plot_f1+
      geom_point(data=tibble_c_data,aes(TIME,CP),alpha=0.5)+
      stat_summary(data=tibble_c_data,aes(TIME,CP),fun="mean",geom="point",size=4,alpha=0.5)
  } else if(category%in%c("inh")){
    plot_f1 = plot_f1 +
      geom_point(data=tibble_inh_data,aes(TIME,IL31inh),alpha=0.5)+
      stat_summary(data=tibble_inh_data,aes(TIME,IL31inh),fun="mean",geom="point",size=4,alpha=0.5)+
      scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25))
  } else{
    plot_f1 = plot_f1 +
      scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25))
  }
  plot_f1
}

# 3つのプロットをまとめて出力する関数
arrange_plot_80160 = function(){
  plot_C = make_plot_80160("C")
  plot_inh = make_plot_80160("inh")
  plot_imp_p = make_plot_80160("imp_p")
  grid.arrange(plot_C,plot_inh,plot_imp_p)
}

arrange_plot_80160()

# 出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
setwd(outpath)
png("simulation_80_160.png",res=300,width=2700,height =3000)
arrange_plot_80160()
dev.off()






