# logitを用いて二値データをそのまま解析する。
# 斎藤さんの行った手法と非線形混合効果のnlme関数を用いた結果を比較する。

library(tidyverse)
library(gt)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
datapath = file.path(basepath,"P1b_Inh_EP")
outpath = file.path(basepath,"figures")

# とりあえず、斎藤さんの行った方法をたどってみる。

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

# 全てのデータを合わせてモデルにフィッティング
# linear modelの対数尤度関数を
linear_logit = function(params,x,y){
  beta0 = params[1]
  beta1 = params[2]
  logit = beta0+beta1*x
  p = 1/(1+exp(-logit))
  # 負の対数尤度を計算
  likelihood = -sum(y*log(p)+(1-y)*log(1-p))
  return(likelihood)
}

Emax_logit = function(params,x,y){
  beta0 = params[1]
  beta1 = params[2]
  beta2 = params[3]
  logit = beta0+(beta1*x/(beta2+x))
  p = 1/(1+exp(-logit))
  # 負の対数尤度を計算
  likelihood = -sum(y*log(p)+(1-y)*log(1-p))
  return(likelihood)
}

x = tibble1$Inh
y = tibble1$flag

# 初期パラメータの設定
initial_params_lin = c(0,1) # beta0,beta1
initial_params_Emax = c(0,0,1) # beta0,beta1

# 最小化問題を解く
result_lin = optim(par=initial_params_lin,fn=linear_logit,x=x,y=y)
result_Emax = optim(par=initial_params_Emax,fn=Emax_logit,x=x,y=y)

# 推定されたパラメータを取得
beta0 = result_lin$par[1]
beta1 = result_lin$par[2]
gammma0 = result_Emax$par[1]
gammma1 = result_Emax$par[2]
gammma2 = result_Emax$par[3]


# 固定効果のみでプロットしてみる。

linear_plot = function(x){
  logit = beta0+beta1*x
  p = 1/(1+exp(-logit))
  return(p)
}

plot1 = ggplot(tibble1,aes(Inh,flag))+
  geom_point()+
  stat_function(fun=linear_plot,linewidth=1)+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  labs(x = "Inhibition rate",
       y = "probability")+
  theme_light()


# 非線形混合効果モデルを適用する。

library(nlme)

# linear modelでlogit=beta0+beta1*x+eta2 eta2~N(0,s^2)
model_l1=nlme(flag~1/(1+exp(-(b0+b1*x+eta2))),
           fixed=list(b0|b1~1), # b0,b1を固定
           random=eta2~1|ID, # etaをランダムパラメータとして推定
           start=c(b0=0,b1=1),
           data=tibble1)



# linear modelでlogit=beta0+(beta1+eta1)*x eta1~N(0,s^2)
model_l2=nlme(flag~1/(1+exp(-(b0+(b1+eta1)*x))),
           fixed=list(b0|b1~1),
           random=list(ID=pdDiag(list(eta1~1))),
           start=c(b0=0,b1=1),
           data=tibble1)

# linear modelでlogit=beta0+(beta1+eta1)*x+eta2 eta1,eta2~N(0,s^2)
model_l3=nlme(flag~1/(1+exp(-(b0+(b1+eta1)*x+eta2))),
              fixed=list(b0|b1~1),
              random=list(ID=pdDiag(list(eta1|eta2~1))),
              start=c(b0=0,b1=1),
              data=tibble1)


# link関数にロジットを指定した線形混合効果モデルも試してみる。
library(lme4)
model_gl1 = glmer(flag~Inh + (1|ID),
                 data=tibble1,family=binomial,
                 control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))
# boundary (singular) fit: see help('isSingular')
# ランダム効果が小さすぎて異常なフィット

model_gl = glmer(flag~Inh + (Inh|ID),
                   data=tibble1,family=binomial,
                 control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))



# linear modelでAIC,BICを比較
# 必要な情報を取り出す。
pull_info = function(name,model){
  AIC = AIC(model)
  BIC = BIC(model)
  sd_list = VarCorr(model)
  if(name%in%c("model_l3")){
    eta1_sd = sd_list["eta1","StdDev"]
    eta2_sd = sd_list["eta2","StdDev"]
  } else if(name=="model_l2"){
    eta1_sd = sd_list["eta1","StdDev"]
    eta2_sd=NA
  } else if(name=="model_l1"){
    eta1_sd=NA
    eta2_sd = sd_list["eta2","StdDev"]
  } else{
    eta1_sd = "1.323e-07"
    eta2_sd = "0"
  }

  tibble = tibble(
    model = name,
    AIC = AIC,
    BIC = BIC,
    eta1_sd = eta1_sd,
    eta2_sd = eta2_sd
  )
  list(tibble)
}

name_list = c("model_l1","model_l2","model_l3","model_gl")
model_list = list(model_l1,model_l2,model_l3,model_gl)
tibble_model = mapply(pull_info,name_list,model_list) %>% 
  bind_rows()

tibble_model # model1,2がAIC, BIC同じぐらいでmodel3より小さい


# Emaxモデル
# logit=beta0+(beta1*x)/(beta2*x)+eta3 eta3~N(0,s^2)

# model_l1=nlme(flag~1/(1+exp(-(g0+(g1*x)/(g2*x)+eta3))),
#               fixed=list(g0|g1|g2~1), # b0,b1を固定
#               random=eta3~1|ID, # etaをランダムパラメータとして推定
#               start=c(g0=0,g1=0,g2=1),
#               data=tibble1)

# Error in MEEM(object, conLin, control$niterEM) : Singularity in backsolve at level 0, block 1
# 上のエラーが発生する。すぐに解決はできなさそう。



# 阻害率と改善確率のプロット
# 全体で阻害率0.1ごとに改善確率を計算する。
inh_interval = seq(0.1,1,0.1)
tibble3 = NULL
for(i in inh_interval){
  start = i-0.1
  end=i
  tibble_p1 = tibble1 %>% 
    filter(Inh>start) %>% 
    filter(Inh<end)
  # flagが1の割合を計算
  percent = sum(tibble_p1$flag)/nrow(tibble_p1)
  tibble_p2 = tibble(
    inh = i,
    rate = percent
  )
  tibble3 = bind_rows(tibble3,tibble_p2)
}



# 阻害率0.1ごとに改善率を計算する。
# IDごとに計算する。
inh_interval = seq(0.1,1,0.1)
inh_rate = function(data){
  tibble_f3 = NULL
  for(i in inh_interval){
    start = i-0.1
    end=i
    tibble_f1 = data %>% 
      filter(Inh>start) %>% 
      filter(Inh<end)
    # flagが1の割合を計算
    percent = sum(tibble_f1$flag)/nrow(tibble_f1)
    tibble_f2 = tibble(
      inh = i,
      rate = percent
    )
    tibble_f3 = bind_rows(tibble_f3,tibble_f2)
  }
  list(tibble_f3)
}
tibble4 = tibble1 %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup
tibble4 = tibble4 %>% 
  mutate(data = mapply(inh_rate,tibble4$data)) %>% 
  group_by(ID) %>% 
  unnest() %>% 
  ungroup()

tibble4 = tibble4 %>% 
  mutate(ID = factor(tibble4$ID,levels=c(1:10)))
color_list = c("black","red","blue","green","yellow","purple","grey40","magenta","green4","red4")

plot2 = ggplot(tibble4,aes(inh,rate,color=ID))+
  geom_point(size=2)+
  geom_point(data=tibble3,aes(inh,rate),color="black",size=4)+
  theme_light()+
  labs(x = "阻害率",y="改善確率")+
  scale_color_manual(values=color_list)
plot2

# model_l1
b0_fix = model_l1$coefficients$fixed["b0"]
b1_fix = model_l1$coefficients$fixed["b1"]
eta_list = model_l1$coefficients$random$ID
s = paste("stat_function(fun=function(x) 1/(1+exp(-(",b0_fix,"+",b1_fix,"*x+",eta_list,"))),color='",color_list,"')",sep="",collapse = "+")
ss = paste("plot3 = plot2 + ",s,sep="")
eval(parse(text=ss))
plot3

# model_l2
b0_fix = model_l2$coefficients$fixed["b0"]
b1_fix = model_l2$coefficients$fixed["b1"]
eta_list = model_l2$coefficients$random$ID
s = paste("stat_function(fun=function(x) 1/(1+exp(-(",b0_fix,"+(",b1_fix,"+",eta_list,")*x))),color='",color_list,"')",sep="",collapse = "+")
ss = paste("plot4 = plot2 + ",s,sep="")
eval(parse(text=ss))
plot4

# model_gl
b0_mix = coef(model_gl)$ID[,1]
b1_mix = coef(model_gl)$ID[,2]
s = paste("stat_function(fun=function(x) 1/(1+exp(-(",b0_mix,"+(",b1_mix,")*x))),color='",color_list,"')",sep="",collapse = "+")
ss = paste("plot4 = plot2 + ",s,sep="")
eval(parse(text=ss))
plot4_2


# それぞれのIDについて誤差項を最尤法により求め分布を求めたモデルか非線形混合効果のnlme関数により求められたモデルのどちらがよいか考える。
# 最尤法modelの値を算出する関数
model_s1 = function(inh){
  # パラメータ
  beta0 = -6.8349
  beta1 = 9.79815
  # 誤差項を乱数生成
  eta = rnorm(1,0,0.16214)
  logit = beta0 + beta1*inh+eta
  p = 1/(1+exp(-logit))
  return(p)
}

# 非線形混合効果によるモデル
model_l1_calc =  function(inh){
  # パラメータ
  b0_fix = model_l1$coefficients$fixed["b0"]
  b1_fix = model_l1$coefficients$fixed["b1"]
  eta_sd = as.numeric(VarCorr(model_l1)["eta2","StdDev"])
  # 誤差項を乱数生成
  eta = rnorm(1,0,eta_sd)
  logit = beta0 + beta1*inh+eta
  p = 1/(1+exp(-logit))
  return(p)
}

model_l2_calc =  function(inh){
  # パラメータ
  b0_fix = model_l2$coefficients$fixed["b0"]
  b1_fix = model_l2$coefficients$fixed["b1"]
  eta_sd = as.numeric(VarCorr(model_l2)["eta1","StdDev"])
  # 誤差項を乱数生成
  eta = rnorm(1,0,eta_sd)
  logit = beta0 + (beta1+eta)*inh
  p = 1/(1+exp(-logit))
  return(p)
}

model_gl_calc = function(inh){
  # パラメータ
  b0_fix = -6.899
  b1_fix = 9.894
  eta0_sd = 0
  eta1_sd = 1.323e-07
  # 誤差項を乱数生成
  eta0 = rnorm(1,0,eta0_sd)
  eta1 = rnorm(1,0,eta1_sd)
  logit = (beta0+eta0) + (beta1+eta1)*inh
  p = 1/(1+exp(-logit))
  return(p)
}

n_simulations=100

# inh 0.01ごとに100回ずつ乱数生成をして、どんな分布になるかを確認する。

simulation_mon = function(inh,model){
  inh_list_f = rep(inh,n_simulations)
  simulation = as.numeric(unlist(mapply(model,inh)))
  # 上位2.5%と下位２.５％の数値を除去する。
  upper_limit = 
  tibble_f1 = tibble(
    inh =inh_list_f,
    simulation = as.numeric(unlist(mapply(model,inh)))
  )
  list(tibble_f1)
}

inh_list = seq(0,1,0.01)
tibble_simu1 = mapply(simulation_mon,inh_list,list(model_s1)) %>% 
  bind_rows()

plot5 = ggplot(tibble_simu1,aes(inh,simulation))+
  geom_point(color="blue",alpha=0.5)+
  geom_point(data=tibble3,aes(inh,rate),size=4)+
  geom_point(data=tibble4,aes(inh,rate),size=2)+
  theme_light()+
  labs(x="阻害率",
       y="改善確率")
  
  

tibble_simu2 = mapply(simulation_mon,inh_list,list(model_l1_calc)) %>% 
  bind_rows()

plot6 = ggplot(tibble_simu2,aes(inh,simulation))+
  geom_point(color="blue",alpha=0.5)+
  geom_point(data=tibble3,aes(inh,rate),size=4)+
  geom_point(data=tibble4,aes(inh,rate),size=2)+
  theme_light()+
  labs(x="阻害率",
       y="改善確率")

tibble_simu3 = mapply(simulation_mon,inh_list,list(model_l2_calc)) %>% 
  bind_rows()

plot7 = ggplot(tibble_simu3,aes(inh,simulation))+
  geom_point(color="blue",alpha=0.5)+
  geom_point(data=tibble3,aes(inh,rate),size=4)+
  geom_point(data=tibble4,aes(inh,rate),size=2)+
  theme_light()+
  labs(x="阻害率",
       y="改善確率")

tibble_simu4 = mapply(simulation_mon,inh_list,list(model_gl_calc)) %>% 
  bind_rows()

plot8 = ggplot(tibble_simu4,aes(inh,simulation))+
  geom_point(color="blue",alpha=0.5)+
  geom_point(data=tibble3,aes(inh,rate),size=4)+
  geom_point(data=tibble4,aes(inh,rate),size=2)+
  theme_light()+
  labs(x="阻害率",
       y="改善確率")
plot8

# plot1~plot7を出力する。
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive=TRUE)
}
setwd(outpath)
plotname_list = c("fixed_model.png","probability.png","model_l1.png","model_l2.png",
                  "prediction_model_s1.png","prediction_model_l1.png","prediction_model_l2.png")
plot_list = list(plot1,plot2,plot3,plot4,plot5,plot6,plot7)
for(i in 1:length(plot_list)){
  png(plotname_list[i],res=300,width=1300,height = 900)
  print(plot_list[[i]])
  dev.off()
  print(i)
}

# tibble_modelを表として出力する。
table = gt(tibble_model)
gtsave(table,"model_compare.png")




