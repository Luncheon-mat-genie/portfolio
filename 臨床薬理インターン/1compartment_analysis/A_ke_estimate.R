# keとAのパラメータ推定をする。

library(tidyverse)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath,"figures")
n_simulations = 3000

# データを読み込む

setwd(basepath)
tibble1 = read_csv("PK_parameters.csv",col_types=cols())


# keのパラメータ推定を行う

tibble_ke = tibble1 %>% dplyr::select(c(ID,dose,ke))

ggplot(tibble_ke,aes(ID,ke))+geom_point()+stat_summary(fun="mean",geom="point",size=3)
ke_ave = mean(tibble_ke$ke)

# ke = ave*exp(eta)と考え、etaの分布を推定する。
# それぞれの患者さんの平均=全体平均*exp(eta)と考えてそれぞれの患者さんのetaを算出する。

tibble_ke2 = tibble_ke %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup()

# それぞれの平均値の算出
calc_ke_idave = function(data){
  id_ave = mean(data$ke)
  return(id_ave)
}

tibble_ke2 = tibble_ke2 %>% 
  mutate(id_ave = as.numeric(unlist(mapply(calc_ke_idave,tibble_ke2$data))))

# eta1の算出
calc_ke_eta1 = function(data){
  id_ave = mean(data$ke)
  eta = log(id_ave/ke_ave)
  return(eta)
}

tibble_ke2 = tibble_ke2 %>% 
  mutate(eta1 = mapply(calc_ke_eta1,tibble_ke2$data))

# etaの普遍標準偏差を計算
eta1_sd = sd(tibble_ke2$eta1)

# 次にkeの個人内変動eta2を計算する。そのために残差を計算する。
calc_ke_eta2 = function(data){
  id_ave = mean(data$ke)
  residuals = data$ke-id_ave
  list(residuals)
}

residuals = unlist(mapply(calc_ke_eta2,tibble_ke2$data))
eta2_sd = sd(residuals)

# 確認のため分布を再現してみる。
# ke = ke_ave * exp(eta1) + eta2
# keを乱数生成する関数
make_ke = function(n_simulations){
  eta1 = rnorm(n_simulations,0,eta1_sd)
  eta2 = rnorm(n_simulations,0,eta2_sd)
  
  ke = ke_ave * exp(eta1) + eta2
  ke
}

tibble_ke_simulate = tibble(
  category = "simulate",
  ke = make_ke(3000)
)

tibble_ke = tibble_ke %>% 
  mutate(category="data")

ggplot()+geom_point(data=tibble_ke,aes(category,ke))+geom_point(data=tibble_ke_simulate,aes(category,ke))

# 本来のデータよりも結構ばらつきが大きい
# keの個人間変動のみを考慮するとどうなる
make_ke_2 = function(n_simulations){
  eta1 = rnorm(n_simulations,0,eta1_sd)
  
  ke = ke_ave * exp(eta1)
  ke
}

ke_simulate = make_ke_2(n_simulations) %>% 
  sort()

tibble_ke_p = tibble(
  category="simulate",
  percent = c("10%","50%","90%"),
  ke = ke_simulate[c(n_simulations*0.1,n_simulations*0.5,n_simulations*0.9)]
)

tibble_ke_simulate = tibble(
  category = "simulate",
  ke = make_ke_2(n_simulations)
)

tibble_ke2 = tibble_ke2 %>% 
  mutate(category="data")

ggplot()+geom_point(data=tibble_ke2,aes(category,id_ave))+geom_point(data=tibble_ke_simulate,aes(category,ke))+
  geom_violin(data=tibble_ke_simulate,aes(category,ke))+
  geom_point(data=tibble_ke_p,aes(category,ke,color=percent))

# eta2が加算誤差だとkeが負の値をとる可能性があり不適切。指数誤差で考える。
# 次にkeの個人内変動eta2を計算する。それぞれのデータの指数誤差を計算する。
calc_ke_eta2 = function(data){
  id_ave = mean(data$ke)
  residuals = log(data$ke)-log(id_ave)
  list(residuals)
}
residuals = unlist(mapply(calc_ke_eta2,tibble_ke2$data))
eta2_sd = sd(residuals)

# 確認のため分布を再現してみる。
# ke = ke_ave * exp(eta1) + eta2

# keを乱数生成する関数
make_ke = function(n_simulations){
  eta1 = rnorm(n_simulations,0,eta1_sd)
  eta2 = rnorm(n_simulations,0,eta2_sd)
  
  ke = ke_ave * exp(eta1) * exp(eta2)
  ke
}


tibble_ke_simulate = tibble(
  category = "simulate",
  ke = make_ke(n_simulations)
)


ke_simulate = make_ke(n_simulations) %>% 
  sort()

tibble_ke_p = tibble(
  category="simulate",
  percent = c("10%","50%","90%"),
  ke = ke_simulate[c(n_simulations*0.1,n_simulations*0.5,n_simulations*0.9)]
)


ggplot()+geom_point(data=tibble_ke,aes(category,ke))+geom_point(data=tibble_ke_simulate,aes(category,ke))+
  geom_violin(data=tibble_ke_simulate,aes(category,ke))+
  geom_point(data=tibble_ke_p,aes(category,ke,color=percent))

# なんかぶっ飛ぶ値がでてきたけど、いい具合じゃない？


# Aのパラメータを推定
# 個人間変動を算出する。
# 全体のデータからslope,emax,ec50の値を求める。
# それぞれの患者さんからslope,ec50を固定してemaxを算出する。
# emaxの誤差を対数正規分布として分布を求める。

library(drc)

tibble_A = tibble1 %>% dplyr::select(c(ID,dose,A))

# drcのフィッティングをする関数
fitting_drc = function(data){
  
  model = drm(A~dose,data = data,
              fct = LL.4(names = c("slope","lower_limit","Emax","EC50"),
                         fixed = c(NA,0,NA,NA)))
  list(model)
}

# 全体でフィッティング
model_A = unlist(fitting_drc(tibble_A))
slope = model_A$`coefficients.slope:(Intercept)`
Emax = model_A$`coefficients.Emax:(Intercept)`
EC50 = model_A$`coefficients.EC50:(Intercept)`

# slopeとec50も固定した状態でフィッティングする関数

# drcのフィッティングをする関数
fitting_drc_fixed = function(data){
  
  model = drm(A~dose,data = data,
              fct = LL.4(names = c("slope","lower_limit","Emax","EC50"),
                         fixed = c(slope,0,NA,EC50)))
  list(model)
}


# それぞれの患者さんでfitting
tibble_A_id = tibble_A %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup()

tibble_A_id = tibble_A_id %>% 
  mutate(model_A = mapply(fitting_drc_fixed,tibble_A_id$data))

# Emaxを算出
pull_Emax = function(model){
  Emax = model$coefficients[1]
  return(Emax)
}

tibble_A_id = tibble_A_id %>% 
  mutate(Emax = as.numeric(unlist(mapply(pull_Emax,tibble_A_id$model_A))))

# それぞれの患者さんについてプロット
plot_id = function(ID,data,Emax){
  plot = ggplot(data,aes(dose,A))+
    geom_point()+
    stat_function(fun = function(x) Emax/(1+exp(slope*(log(x)-log(EC50)))))+
    scale_x_log10()+
    labs(ttle="ID")
  plot
  list(plot)
}

tibble_A_id = tibble_A_id %>% 
  mutate(plot = mapply(plot_id,tibble_A_id$ID,tibble_A_id$data,tibble_A_id$Emax))

gridExtra::grid.arrange(tibble_A_id$plot[[1]],tibble_A_id$plot[[2]],tibble_A_id$plot[[3]],
                        tibble_A_id$plot[[4]],tibble_A_id$plot[[5]],tibble_A_id$plot[[6]],
                        ncol=3)

# Emaxの誤差を求める。
# Emax_id = Emax_all * exp(eta1)

calc_eta1_Emax = function(Emax_id){
  Emax_all = Emax
  eta1 = log(Emax_id/Emax_all)
  return(eta1)
}

tibble_A_id = tibble_A_id %>% 
  mutate(eta1 = as.numeric(unlist(mapply(calc_eta1_Emax,tibble_A_id$Emax))))

eta1_sd = sd(tibble_A_id$eta1)

# 乱数生成したEmaxとEmaxの分布を比べる。
# Emaxを乱数生成
make_Emax = function(n_simulations){
  eta1 = rnorm(n_simulations,0,eta1_sd)
  Emax_make = Emax * exp(eta1)
  Emax_make
}


Emax_simulate = make_Emax(n_simulations) %>% 
  sort()

tibble_Emax_simulation = tibble(
  category="simulation",
  Emax = Emax_simulate
)

tibble_Emax_p = tibble(
  category="simulation",
  percent = c("10%","50%","90%"),
  Emax = Emax_simulate[c(n_simulations*0.1,n_simulations*0.5,n_simulations*0.9)]
)


tibble_A_id = tibble_A_id %>% 
  mutate(category="data")

ggplot()+
  geom_point(data=tibble_A_id,aes(category,Emax))+
  geom_point(data=tibble_Emax_simulation,aes(category,Emax))+
  geom_violin(data=tibble_Emax_simulation,aes(category,Emax))+
  geom_point(data=tibble_Emax_p,aes(category,Emax,color=percent))
  

# 乱数生成したEmaxでの曲線と実データを比べてみる

Emax_simulate_plot = make_Emax(200) %>% 
  sort()

s = paste("stat_function(fun = function(x) ",Emax_simulate_plot,"/(1+exp(slope*(log(x)-log(EC50)))),color='grey50')",
          sep="",collapse = "+")
plot_Emax_simulate = ggplot(tibble_A,aes(dose,A))
ss = paste("plot_Emax_simulate = plot_Emax_simulate +",s,sep="")
eval(parse(text = ss))
plot_Emax_simulate = plot_Emax_simulate+
  geom_point()+
  stat_summary(fun="mean",geom="point",size=4)+
  scale_x_log10()

plot_Emax_simulate

# さらに、Aの個人内変動eta2を推定する。そのために各患者でのeta2を計算する。

calc_eta2_Emax = function(data,Emax_id){
  tibble_f2  = NULL
  for(i in 1:nrow(data)){
    dose = data$dose[i]
    A = data$A[i]
        eta2 = log((A*(1+exp(slope*log(dose)-log(EC50))))/Emax_id)
    tibble_f1 = tibble(
      dose=dose,
      A=A,
      eta2 = eta2
    )
    tibble_f2 = bind_rows(tibble_f2,tibble_f1)
  }
  list(tibble_f2)
}

tibble_eta2 = mapply(calc_eta2_Emax,tibble_A_id$data,tibble_A_id$Emax) %>% 
  bind_rows()

ggplot()+
  geom_point(data=tibble_eta2,aes(dose,eta2))+
  scale_x_log10()

# eta2は一定の残差を持ち、160mgから頭打ちになる
# eta2はdoseを常用対数変換した時線形回帰することができる。
# 0.1 ~ 160mgのデータを用いて線形回帰、160mg ~ 320mgの平均で傾き0の直線
# 二つの線が交わるラインで使用する回帰直線を変更する。
# 残差はdoseにかかわらず一定の正規分布に従うものとする。

tibble_eta2 = tibble_eta2 %>% 
  mutate(dose_log10 = log10(tibble_eta2$dose))
tibble_01_160 = tibble_eta2 %>% 
  filter(dose<=160)
tibble_160_320 = tibble_eta2 %>% 
  filter(dose>=320)
lm1 = lm(eta2~dose_log10,tibble_01_160)
a1 = lm1$coefficients[2]
b1 = lm1$coefficients[1]
top_intercept = mean(tibble_160_320$eta2)

ggplot(data=tibble_eta2,aes(dose_log10,eta2))+
  geom_point()+
  stat_function(fun=function(x) lm1$coefficients[2]*x+lm1$coefficients[1],color="red")+
  geom_abline(intercept=top_intercept,slope=0,color="blue")

# doseからeta2の平均値を計算する関数
calc_eta2_ave = function(dose){ # doseは対数変換していないもの
  if(dose>0){
    dose_log10 = log10(dose)
    x_switch = (top_intercept-b1)/a1
    if(dose_log10<x_switch){
      y = a1*dose_log10+b1
    } else{
      y = top_intercept
    }
    return(y)
  } else{
    eta2 = -100
  }
}

# それぞれの直線の残差を計算する。
tibble_01_160 = tibble_01_160 %>% 
  mutate(residuals = lm1$residuals)
tibble_160_320 = tibble_160_320 %>% 
  mutate(residuals = tibble_160_320$eta2-top_intercept)
tibble_eta2_all = bind_rows(tibble_01_160,tibble_160_320) %>% 
  mutate(label="all")

ggplot(data=tibble_eta2_all,aes(label,residuals))+
  geom_point()

eta2_residual_sd = sd(tibble_eta2_all$residuals)

# Aの分布をシミュレーションしてみる。
# 100個のモデルを作成し、dose 1ごとに50個プロットを作成する。

# デバッグ用
model_number=100
plot_number=50
gap=1
make_A = function(model_number,plot_number,gap){
  dose_list = seq(0,320,gap)
  tibble_f2 = NULL
  for(dose in dose_list){
    # model_numberの数だけモデルを作成する。
    eta1 = rnorm(model_number,0,eta1_sd)
    Emax_list = Emax * exp(eta1)
    A_list = c()
    for(Emax_id in Emax_list){
      eta2_residual = rnorm(plot_number,0,eta2_residual_sd)
      eta2_ave = calc_eta2_ave(dose)
      eta2 = eta2_ave+eta2_residual
      A = (Emax_id/(1+exp(slope*log(dose)-log(EC50)))) * exp(eta2)
      A_list = c(A_list,A)
    }
    tibble_f1 = tibble(
      dose = dose,
      A = A_list
    )
    tibble_f2 = bind_rows(tibble_f2,tibble_f1)
  }
  tibble_f2
}

tibble_A_simulation = make_A(50,50,10)

tibble_A_sim_id = tibble_A_simulation %>% 
  group_by(dose) %>% 
  nest() %>% 
  ungroup()
# 各doseで90%信頼区間に該当する点を選択
percent_90 = function(dose,data){
  num1 = data$A %>% sort()
  num1 = num1[c(2500*0.1,2500*0.9)]
  tibble_f1 = tibble(
    dose = dose,
    percent = c("10%","90%"),
    A = num1
  )
  list(tibble_f1)
}

tibble_A_percent = mapply(percent_90,tibble_A_sim_id$dose,tibble_A_sim_id$data) %>% 
  bind_rows()

ggplot()+
  geom_point(data=tibble_A_simulation,aes(dose,A))+
  geom_point(data=tibble_A,aes(dose,A),color="blue",size=3)+
  stat_summary(data=tibble_A,aes(dose,A),fun="mean",geom="point",size=5,color="blue")+
  geom_point(data=tibble_A_percent,aes(dose,A,color=percent),size=4)+
  scale_x_log10()




# 実際に血中濃度→阻害率→改善率をプロットしてみる。

# 服用タイミング、投与量の指定、投与期間
timing = c(9,19)
dose = c(60,60)
days = 10

# 各種パラメータ
ka = 0.5
ke_ave = 0.02999667
ke_eta1_sd = 0.08947654
ke_eta2_sd = 0.3202013
A_slope = -1.72909
A_EC50 = 63.85906
Emax_ave = 203.2953
Emax_eta1_sd = 0.09366983
a1 = 2.07656 
b1 = -4.657897 
top_intercept = -0.1177216
A_eta2_residual_sd = 0.2933491

n_simulations_eta1=150
n_simulations_eta2=150

# 血中濃度
# あるke, Aから血中濃度を計算する関数
calc_C = function(t,ke,A){
  C = A*(exp(-ke*t)-exp(-ka*t))
  return(C)
}

# keを乱数生成する関数
make_ke = function(n_simulations){
  eta1 = rnorm(n_simulations,0,ke_eta1_sd)
  eta2 = rnorm(n_simulations,0,ke_eta2_sd)
  
  ke = ke_ave * exp(eta1) * exp(eta2)
  ke
}

# Emaxを乱数生成する関数
make_Emax = function(n_simulations){
  eta1 = rnorm(n_simulations,0,Emax_eta1_sd)
  Emax_make = Emax_ave * exp(eta1)
  Emax_make
}

# dose(投与量)からA_eta2_aveを計算する関数
calc_eta2_ave = function(dose){ # doseは対数変換していないもの
  if(dose>0){
    dose_log10 = log10(dose)
    x_switch = (top_intercept-b1)/a1
    if(dose_log10<x_switch){
      y = a1*dose_log10+b1
    } else{
      y = top_intercept
    }
  } else{
    y = -100
  }
  return(y)
}

# Aを乱数生成する関数
make_A = function(dose,n_simulations){
  
  Emax = make_Emax(n_simulations)
  eta2_residual = rnorm(plot_number,0,A_eta2_residual_sd)
  eta2_ave = calc_eta2_ave(dose)
  eta2 = eta2_ave+eta2_residual
  A = (Emax/(1+exp(slope*log(dose)-log(EC50)))) * exp(eta2)
  return(A)
}


# 投与終了までにどのタイミング・投与量で投与するかをtibbleにする。
admini_timing = function(timing,dose){
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

timing_list = mapply(admini_timing,timing,dose) %>% 
  bind_rows() %>% 
  arrange(timing)


# ある投与量の時に、n_simulationsの数だけ1時間ごとに1000個乱数を生成し血中濃度を取得する関数
# days*24-timingの期間1hずつ血中濃度を計算する。
calc_C24 = function(time,dose){
  time_range = c(time,days*24)
  time_range_calc = time_range-time
  # keを乱数生成する
  ke = make_ke(n_simulations)
  # Aを乱数生成する
  A = make_A(dose,n_simulations)
}

# 個人内変動のみを考慮して1人の通算のデータを計算する。

make_C_all_id_last = function(timing_list){
  calc_C_all_id = function(time,dose,n_simulations_eta2){
    # keを乱数生成する関数
    make_ke_id = function(n_simulations){
      eta1 = rnorm(1,0,ke_eta1_sd)
      eta2 = rnorm(n_simulations,0,ke_eta2_sd)
      
      ke = ke_ave * exp(eta1) * exp(eta2)
      ke
    }
    # keを乱数生成する。
    ke = make_ke_id(n_simulations_eta2)
    # Aを乱数生成する関数
    make_A_id = function(dose,n_simulations){
      # Emaxを1つ生成する関数
      make_Emax_id = function(){
        eta1 = rnorm(1,0,Emax_eta1_sd)
        Emax_make = Emax_ave * exp(eta1)
        Emax_make
      }
      Emax = make_Emax_id()
      eta2_residual = rnorm(n_simulations,0,A_eta2_residual_sd)
      eta2_ave = calc_eta2_ave(dose)
      eta2 = eta2_ave+eta2_residual
      A = (Emax/(1+exp(slope*log(dose)-log(EC50)))) * exp(eta2)
      return(A)
    }
    # Aを乱数生成する。
    A = make_A_id(dose,n_simulations_eta2)
    
    C = calc_C(time,ke,A)
    tibble_f1 = tibble(
      t = time,
      C = C,
      number = c(1:n_simulations_eta2)
    )
    list(tibble_f1)
  }
  
  # timing_listにそってまずは個人内変動のみを考慮した血中濃度の再現をする。
  calc_C_all_id2 = function(time,dose){
    time_range = c(time:(days*24))
    time_range_calc = time_range-time
    tibble_f1 = tibble(
      real_time = time_range,
      C = mapply(calc_C_all_id,time_range_calc,dose,n_simulations_eta2)
    )
    list(tibble_f1)
    
  }
  
  tibble_sim1 = timing_list %>% 
    mutate(data = mapply(calc_C_all_id2,timing_list$timing,timing_list$dose))
  tibble_sim1 = tibble_sim1 %>% 
    group_by(timing,dose) %>% 
    unnest() %>% 
    ungroup()
  
  tibble_sim1 = tibble_sim1 %>% 
    group_by(real_time) %>% 
    nest() %>% 
    ungroup()
  
  # 各投与の血中濃度を線形結合する。
  bind_C = function(data){
    # 線形結合する。
    tibble_f1 = data$C %>% bind_rows() %>% 
      group_by(number) %>% 
      summarise(C = sum(C)) %>% 
      ungroup()
    list(tibble_f1)
  }
  
  tibble_sim1 = tibble_sim1 %>% 
    mutate(data = mapply(bind_C,tibble_sim1$data))
  tibble_sim1
}

# これで個人内変動のみを考慮した一人分のデータを再現する関数ができた
# 個人のデータを100個作成し、最終的な信頼区間を推定する


tibble_Clast = NULL
for(i in 1:n_simulations_eta1){
  tibble_p = make_C_all_id_last(timing_list)
  tibble_Clast[[length(tibble_Clast)+1]] = tibble_p
  print(i)
}

tibble_Clast2 = tibble_Clast %>% 
  bind_rows()
tibble_Clast2 = tibble_Clast2 %>% 
  group_by(real_time) %>% 
  unnest() %>% 
  nest() %>% 
  ungroup()


# 10%,50%,90%,averageをそれぞれの時間で計算する
calc_percent = function(data){
  n_all = n_simulations_eta1*n_simulations_eta2
  C = data$C %>% sort()
  C_percent = C[c(n_all*0.1,n_all*0.5,n_all*0.9)]
  C_ave = mean(C)
  tibble_f1 = tibble(
    percent = c("10%","50%","90%","average"),
    C = c(C_percent,C_ave)
  )
  list(tibble_f1)
}

tibble_Clast2 = tibble_Clast2 %>% 
  mutate(data = mapply(calc_percent,tibble_Clast2$data))

# シミュレーションデータ生成完了
tibble_Clast2 = tibble_Clast2 %>% 
  group_by(real_time) %>% 
  unnest() %>% 
  ungroup()

# プロットする。
ggplot()+
  geom_line(data=tibble_Clast2,aes(real_time,C,color=percent))+
  scale_y_continuous(limits=c(0,400))




# ここ手動で振り分けないといけないです。
tibble_80 = tibble_Clast2
tibble_160 = tibble_Clast2
tibble_60_60 = tibble_Clast2


tibble_80 = tibble_80 %>% mutate(DOSE=80)
tibble_160 = tibble_160 %>% mutate(DOSE=160)
tibble_sim_80_160 = bind_rows(tibble_80,tibble_160)


# Individual_concentration_multiple_dose.csvを読み込む

setwd(basepath)
tibble_data = read_csv("Individual_concentration_multiple_dose.csv",col_types=cols())

# プロットする。
plot_simulate = ggplot()+
  geom_line(data=tibble_sim_80_160,aes(real_time,C,color=percent),linewidth=1)+
  geom_point(data=tibble_data,aes(TIME,CP),alpha=0.3)+
  stat_summary(data=tibble_data,aes(TIME,CP),fun="mean",geom="point",size=4,alpha=0.3)+
  scale_y_continuous(limits=c(0,400))+
  facet_wrap(~DOSE)


# 出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
setwd(outpath)
png("simulation_result.png",res=300,width=2800,height = 1200)
plot_simulate
dev.off()


