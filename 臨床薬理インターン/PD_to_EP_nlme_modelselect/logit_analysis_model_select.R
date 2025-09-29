# logitを用いて二値データをそのまま解析する。
# どのモデルがいいのか全パターン非線形混合効果モデルを作成してAICを算出する。

library(tidyverse)
library(gt)
library(nlme)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
datapath = file.path(basepath,"P1b_Inh_EP")
outpath = file.path(basepath,"figures")


ID = c(1:10)
Doses = c(80,160)

# データの読み込み
setwd(datapath)
tibble1 = NULL
for(id in ID){
  for(dose in Doses){
    csv = paste("ID=",id,"_DOSE=",dose,"_Inh_EP.csv",sep="")
    tibble_p1 = read_csv(csv,col_types=cols()) %>% 
      mutate(ID=id,
             DOSE=dose)
    tibble1 = bind_rows(tibble1,tibble_p1)
  }
}

# TRUE=1,FALSE=0に変換
makeflag = function(VAS50){
  if(VAS50=="TRUE"){
    flag=1
  } else{
    flag=0
  }
  return(flag)
}

tibble1 = tibble1 %>% 
  mutate(flag = as.numeric(unlist(mapply(makeflag,tibble1$VAS50))))



# 非線形混合効果モデルを適用する。
# modelのパターンを指定する。
model_list = list(
  linear = flag~1/(1+exp(-(b0+b1*Inh))),
  linear = flag~1/(1+exp(-(b0+b1*Inh))),
  linear = flag~1/(1+exp(-(b0+b1*Inh))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh)))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh)))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh)))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh)))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh)))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh)))),
  Emax = flag~1/(1+exp(-(b0+(b1*Inh)/(b2*Inh))))
  
)
# modelの固定効果
fixed_list = list(
  linear = list(b0|b1~1), # b0,b1を固定
  linear = list(b0|b1~1), 
  linear = list(b0|b1~1), 
  Emax = list(b0|b1|b2~1),
  Emax = list(b0|b1|b2~1),
  Emax = list(b0|b1|b2~1),
  Emax = list(b0|b1|b2~1),
  Emax = list(b0|b1|b2~1),
  Emax = list(b0|b1|b2~1),
  Emax = list(b0|b1|b2~1)

)

# modelの変量効果
randomed_list = list(
  l1 = b0~1|ID,
  l2 = b1~1|ID,
  l3 = b0+b1~1|ID,
  E1 = b0~1|ID,
  E2 = b1~1|ID,
  E3 = b2~1|ID,
  E4 = b0+b1~1|ID,
  E5 = b0+b2~1|ID,
  E6 = b1+b2~1|ID,
  E7 = b0+b1+b2~1|ID
)

# 初期値のリスト
start_list = list(
  linear = c(b0=0,b1=1),
  linear = c(b0=0,b1=1),
  linear = c(b0=0,b1=1),
  Emax = c(b0=0,b1=0,b2=1),
  Emax = c(b0=0,b1=0,b2=1),
  Emax = c(b0=0,b1=0,b2=1),
  Emax = c(b0=0,b1=0,b2=1),
  Emax = c(b0=0,b1=0,b2=1),
  Emax = c(b0=0,b1=0,b2=1),
  Emax = c(b0=0,b1=0,b2=1)
)

# modelの名前
name_list = list(
  l1 = "l1",
  l2 = "l2",
  l3 = "l3",
  E1 = "E1",
  E2 = "E2",
  E3 = "E3",
  E4 = "E4",
  E5 = "E5",
  E6 = "E6",
  E7 = "E7"
  
)


# 全てのパターンのモデルをフィッティングする。
tibble2 = NULL
for(i in 1:length(model_list)){
  model = try(nlme(model_list[[i]],
               fixed=fixed_list[[i]],
               random=randomed_list[[i]],
               start=start_list[[i]],
               data=tibble1))
  # if(class(model)=="try-error"){
  #   AIC <- AIC(model)
  #   # パラメータ数の取得
  #   p <- length(fixef(model)) + length(ranef(model))
  #   
  #   # cAIC を計算
  #   cAIC <- AIC + 2*p*(p+1)/(nrow(tibble1)-p-1)
  #   # BICを計算
  #   BIC = BIC(model)
  # } else{
  #   AIC=NA
  #   cAIC=NA
  #   BIC=NA
  # }
  
  tibble_p1 = tibble(
    name = name_list[[i]],
    model = list(model)
  )
  tibble2 = bind_rows(tibble2,tibble_p1)
                    
}

# tibble2からNAを取り除く
tibble2 = tibble2[c(1:3),]

AIC_model = function(model){
  return(AIC(model))
}

BIC_model = function(model){
  return(BIC(model))
}

tibble2 = tibble2 %>% 
  mutate(AIC = as.numeric(unlist(mapply(AIC_model,tibble2$model))),
         BIC = as.numeric(unlist(mapply(BIC_model,tibble2$model))))

# 予測区間を推定する。
inh_list = seq(0,1,0.001)
n_simulations=100

l1_plot = function(Inh){
  model = tibble2$model[[1]]
  b0 = model$coefficients$fixed["b0"]
  b1 = model$coefficients$fixed["b1"]
  b0_eta_sd = as.numeric(VarCorr(model)["b0","StdDev"])
  # 乱数生成
  b0_eta = rnorm(n_simulations,0,b0_eta_sd)
  logit = b0+b0_eta+b1*Inh
  p = 1/(1+exp(-logit))
  tibble_f1 = tibble(
    inh = Inh,
    p  = p
  )
  list(tibble_f1)
}

l2_plot = function(Inh){
  model = tibble2$model[[2]]
  b0 = model$coefficients$fixed["b0"]
  b1 = model$coefficients$fixed["b1"]
  b1_eta_sd = as.numeric(VarCorr(model)["b1","StdDev"])
  # 乱数生成
  b1_eta = rnorm(n_simulations,0,b1_eta_sd)
  logit = b0+(b1_eta+b1)*Inh
  p = 1/(1+exp(-logit))
  tibble_f1 = tibble(
    inh = Inh,
    p  = p
  )
  list(tibble_f1)
}

l3_plot = function(Inh){
  model = tibble2$model[[3]]
  b0 = model$coefficients$fixed["b0"]
  b1 = model$coefficients$fixed["b1"]
  b0_eta_sd = as.numeric(VarCorr(model)["b0","StdDev"])
  b1_eta_sd = as.numeric(VarCorr(model)["b1","StdDev"])
  # 乱数生成
  b0_eta = rnorm(n_simulations,0,b0_eta_sd)
  b1_eta = rnorm(n_simulations,0,b1_eta_sd)
  logit = b0+b0_eta+(b1+b1_eta)*Inh
  p = 1/(1+exp(-logit))
  tibble_f1 = tibble(
    inh = Inh,
    p  = p
  )
  list(tibble_f1)
}

s1_plot = function(Inh){
  b0 = -6.8349
  b1 = 9.798149516668506
  b0_eta_sd = 0.170910378387596
  # 乱数生成
  b0_eta = rnorm(n_simulations,0,b0_eta_sd)
  logit = b0+b0_eta+b1*Inh
  p = 1/(1+exp(-logit))
  tibble_f1 = tibble(
    inh = Inh,
    p  = p
  )
  list(tibble_f1)
}

# plotを作成する。
make_plot  = function(fun){
  tibble_f1 = mapply(fun,inh_list) %>% 
    bind_rows()
  plot_f1 = ggplot(tibble_f1,aes(inh,p))+
    geom_point()+
    theme_classic()
  plot_f1
}

make_plot(l1_plot)
make_plot(l2_plot)
make_plot(l3_plot)
make_plot(s1_plot)


# 90%信頼区間を推定する。
n_simulations = 1000
estimate_90 = function(fun){
  tibble_f1 = mapply(fun,inh_list)
  # それぞれの阻害率で上位5%, 下位5%を削る。
  n = n_simulations*0.05
  estimate_90_2 = function(data){
    tibble_f1 = data %>% 
      arrange(p)
    tibble_f1 = tibble_f1[c(n,(n_simulations-n)),]
    list(tibble_f1)
  }
  tibble_f2 = mapply(estimate_90_2,tibble_f1) %>% 
    bind_rows()
  plot_f1 = ggplot(tibble_f2,aes(inh,p))+
    geom_point(size=1)+
    theme_classic()
  plot_f1
}

estimate_90(l1_plot)
estimate_90(l2_plot)
estimate_90(l3_plot)
estimate_90(s1_plot)

# 出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath)
}
setwd(outpath)
png("model_l1.png",res=300,width=1200,height = 800)
estimate_90(l1_plot)
dev.off()

png("model_l2.png",res=300,width=1200,height = 800)
estimate_90(l2_plot)
dev.off()

png("model_l3.png",res=300,width=1200,height = 800)
estimate_90(l3_plot)
dev.off()

png("model_s1.png",res=300,width=1200,height = 800)
estimate_90(s1_plot)
dev.off()

library(gt)
tibble_table = tibble2 %>% dplyr::select(-model)
gtsave(gt(tibble_table),"table.png")
