# PKのパラメータを推定する。
# 流れ
# doseとlog(A)の線形回帰 log(A) = 

library(tidyverse)
library(drc)
library(gridExtra)
library(lme4)
library(ggpmisc)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
outpath = file.path(basepath,"output_ke_eta1_exp")
datapath = file.path(basepath,"output_ke_eta1_exp")
if(!(dir.exists(outpath))){
  dir.create(outpath,recursive = TRUE)
}
n_simulations = 100000 # 信頼区間を推定する際のモンテカルロシミュレーションの反復回数

# データを読み込む

setwd(basepath)
tibble1 = read_csv("PK_parameters_a.csv",col_types=cols()) 
tibble1 = tibble1 %>% 
  mutate(logA = log(tibble1$A),
         log10dose = log10(tibble1$dose))

setwd(outpath)
png("log10dose_logA.png",res=300,width=1200,height = 1000)
ggplot(tibble1,aes(log10dose,logA))+
  geom_point()
dev.off()

# log(A) = A_a*log10(dose) + (A_b+A_eta2) +A_eta3としてlog(A)の分布を推定する。
# eta1,eta2: 個人間変動　eta3: 個人内変動　全て正規分布に従う
# 線形混合効果モデルで個人間変動を推定する。
# 320mgを除く
tibble_u160 = tibble1 %>% 
  filter(!(dose==320))
fit.lmer = lmer(logA~log10dose + (1|ID),data=tibble_u160)
# 傾きにランダム効果を入れようとするとエラーが発生するので切片のみにランダム効果を考える。
summary.fit.lmer = summary(fit.lmer)
A_a = summary.fit.lmer$coefficients[2]
A_b = summary.fit.lmer$coefficients[1]
A_eta2_sd = 0.05574

# 線形混合効果モデルをプロットしてみる。
ggplot(tibble1,aes(log10dose,logA))+
  geom_point()+
  stat_function(fun=function(x) A_a*x+A_b)+
  stat_function(fun=function(x) A_a*x+A_b+2*A_eta2_sd,color="blue")+
  stat_function(fun=function(x) A_a*x+A_b-2*A_eta2_sd,color="blue")


# それぞれのIDでプロットする。
tibble_u160_ID = tibble_u160 %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(eta2 = ranef(fit.lmer)$ID[,1])

plot_doseA_ID = function(ID,data,A_eta2){
  plot = ggplot(data,aes(log10dose,logA))+
    geom_point()+
    stat_function(fun=function(x) A_a*x+A_b+A_eta2)+
    labs(title=paste("ID",ID,sep=""))
  list(plot)
}

tibble_u160_ID = tibble_u160_ID %>% 
  mutate(plot1 = mapply(plot_doseA_ID,
                        tibble_u160_ID$ID,tibble_u160_ID$data,tibble_u160_ID$eta2))

grid.arrange(tibble_u160_ID$plot1[[1]],tibble_u160_ID$plot1[[2]],tibble_u160_ID$plot1[[3]],
             tibble_u160_ID$plot1[[4]],tibble_u160_ID$plot1[[5]],tibble_u160_ID$plot1[[6]],
             ncol=3)

# 残差を計算する。
calc_residuals = function(data,A_eta2){
  logA_sim = A_a*data$log10dose+A_b+A_eta2
  residuals = data$logA-logA_sim
  tibble_f1 = data %>% 
    mutate(residuals = residuals)
  list(tibble_f1)
}

tibble_u160_ID = tibble_u160_ID %>% 
  mutate(data=mapply(calc_residuals,tibble_u160_ID$data,tibble_u160_ID$eta2))

tmp1 = bind_rows(tibble_u160_ID$data)
ggplot(tmp1,aes(log10dose,residuals))+
  geom_point()

# 頭打ちになるところのパラメータ推定と残差の推定をする。
# 160mgと320mgについてlogA = (c + eta4) + eta3
# eta4: 個人間変動, eta3: 個人内変動
# 全ての160mgと320mgのlogAの平均
tibble_160_320 = tibble1 %>% 
  filter(dose%in%c(160,320))
A_c = mean(tibble_160_320$logA)

# それぞれの患者さんでeta4をけいさんする。
tibble_160_320_ID = tibble_160_320 %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup()

calc_eta4 = function(data){
  c_id = mean(data$logA)
  eta4 = c_id-A_c
  eta4
}

tibble_160_320_ID = tibble_160_320_ID %>% 
  mutate(eta4 = as.numeric(unlist(mapply(calc_eta4,tibble_160_320_ID$data))))
A_eta4_sd = sd(tibble_160_320_ID$eta4)

# 残差を計算する。
calc_residuals_c = function(data,eta4){
  logA_sim = A_c+eta4
  residuals = data$logA-logA_sim
  data = data %>% 
    mutate(residuals=residuals)
  list(data)
}

tibble_160_320_ID = tibble_160_320_ID %>% 
  mutate(data = mapply(calc_residuals_c,tibble_160_320_ID$data,tibble_160_320_ID$eta4))
tmp2 = bind_rows(tibble_160_320_ID$data)
ggplot(tmp2,aes(log10dose,residuals))+
  geom_point()+
  geom_point(data=tmp1,aes(log10dose,residuals))

# 全ての容量で同じ個体内変動に従うと仮定する。
# Aに関しては160mgは完全に飽和していると考えて、線形回帰の方の160mgは除去する。
tmp1 = tmp1 %>% 
  filter(!(dose=160))

tibble_residuals = bind_rows(tmp1,tmp2)
A_eta3_sd = sd(tibble_residuals$residuals)

# 残差込みでプロットしてみる。
plot_doseA_ID_eta3 = function(ID,data,A_eta2){
  plot = ggplot(data,aes(log10dose,logA))+
    geom_point()+
    stat_function(fun=function(x) A_a*x+A_b+A_eta2)+
    stat_function(fun=function(x) A_a*x+A_b+A_eta2+2*A_eta3_sd,color="blue")+
    stat_function(fun=function(x) A_a*x+A_b+A_eta2-2*A_eta3_sd,color="blue")+
    labs(title=paste("ID",ID,sep=""))
  list(plot)
}

tibble_u160_ID = tibble_u160_ID %>% 
  mutate(plot2 = mapply(plot_doseA_ID_eta3,
                        tibble_u160_ID$ID,tibble_u160_ID$data,tibble_u160_ID$eta2))

grid.arrange(tibble_u160_ID$plot2[[1]],tibble_u160_ID$plot2[[2]],tibble_u160_ID$plot2[[3]],
             tibble_u160_ID$plot2[[4]],tibble_u160_ID$plot2[[5]],tibble_u160_ID$plot2[[6]],
             ncol=3)

# 十分に再現できている感じがする。

# Aの分布を再現してみる。
# まず個体間変動をn個乱数生成し、その一つ一つについて個体内変動をつける。
# doseはlog10(dose)について-1~2.5を0.01ずつ
log10dose_list = seq(-1,2.5,0.01)
n_inter = 3000

check_c = function(logA,eta4){
  c_id = A_c+eta4
  if(logA>c_id){
    logA=c_id
  }
  return(logA)
}

make_logA = function(n_inter,log10dose){
  # 個体間変動
  eta2 = rnorm(n_inter,0,A_eta2_sd)
  eta4 = rnorm(n_inter,0,A_eta4_sd)
  
  logA = A_a*log10dose + (A_b+eta2)
  
  # logAがc_idを超えていたらc_idに
  logA = as.numeric(unlist(mapply(check_c,logA,eta4)))
  
  eta3 = rnorm(n_inter,0,A_eta3_sd)
  # 個体内変動
  logA_sim = logA + eta3
  
  tibble_f1 = tibble(
    ID = 1:n_inter,
    dose=10^log10dose,
    log10dose=log10dose,
    A = exp(logA_sim),
    logA = logA_sim
  )
  list(tibble_f1)
}

tibble_A_sim = mapply(make_logA,n_inter,log10dose_list) %>% 
  bind_rows()

# それぞれのdoseで5%, 50%, 95%, aveを抽出する。
tibble_A_sim_dose = tibble_A_sim %>% 
  group_by(dose) %>% 
  nest() %>% 
  ungroup()

pullrv = function(dose,data){
  A = data$A
  A = A %>% sort()
  A_rv = A[c(n_inter*0.05,n_inter*0.5,n_inter*0.95)]
  ave = mean(A)
  tibble_f1 = tibble(
    dose = dose,
    A = c(A_rv,ave),
    logA = log(A),
    category = c("5%","50%","95%","ave")
  )
  list(tibble_f1)
}

tibble_A_sim_rv = mapply(pullrv,tibble_A_sim_dose$dose,tibble_A_sim_dose$data) %>% 
  bind_rows()

ggplot()+
  geom_line(data=tibble_A_sim_rv,aes(dose,A,color=category),linewidth=1)+
  geom_point(data=tibble1,aes(dose,A))+
  stat_summary(data=tibble1,aes(dose,A),fun="mean",geom="point",color="blue",size=4)+
  scale_x_log10()



# keのパラメータ推定を行う
# eta1 個人間変動、eta2 個人内変動
# eta1,eta2自体は正規分布に従うが、指数変換することで対数正規分布に従う誤差を想定している。
# ke = (ke_fixed*exp(eta1)) * exp(eta2)と考えてke_fixed, eta1, eta2の推定を行う

# keのデータを取り出す
tibble_ke = tibble1 %>% dplyr::select(c(ID,dose,ke))

ggplot(tibble_ke,aes(ID,ke))+geom_point()+stat_summary(fun="mean",geom="point",size=3)

# ke_fixedの推定。

tibble_ke2 = tibble_ke %>% 
  group_by(ID) %>% 
  nest() %>% 
  ungroup()

# それぞれの患者さんのkeの平均値の算出
calc_ke_idave = function(data){
  id_ave = mean(data$ke)
  return(id_ave)
}

tibble_ke2 = tibble_ke2 %>% 
  mutate(id_ave = as.numeric(unlist(mapply(calc_ke_idave,tibble_ke2$data))))

ke_fixed = exp(mean(log(tibble_ke2$id_ave)))

# eta1をすいていする。
# それぞれの患者さんの平均=ke_fixed*exp(eta1)と考えてそれぞれの患者さんのeta1を算出する。
# eta1の算出
calc_ke_eta1 = function(data){
  id_ave = mean(data$ke)
  eta = log(id_ave/ke_fixed)
  return(eta)
}

tibble_ke2 = tibble_ke2 %>% 
  mutate(eta1 = mapply(calc_ke_eta1,tibble_ke2$data))

# eta1は正規分布に従っているので普遍標準偏差を計算
ke_eta1_sd = sd(tibble_ke2$eta1)

# eta1の推定がどんなもんか検証する。
# eta1を乱数生成する関数
make_ke_id_ave = function(n_simulations){
  eta1 = rnorm(n_simulations,0,ke_eta1_sd)
  ke_id_ave = ke_fixed*exp(eta1)
  ke_id_ave
}
tibble_ke_id_ave_sim = tibble(
  id_ave = make_ke_id_ave(n_simulations)
)
plot_ke_id_ave = ggplot()+
  geom_violin(data=tibble_ke2,aes("data",id_ave))+
  geom_boxplot(data=tibble_ke2,aes("data",id_ave),width=.5)+
  geom_point(data=tibble_ke2,aes("data",id_ave))+
  geom_violin(data=tibble_ke_id_ave_sim,aes("simulation",id_ave))+
  geom_boxplot(data=tibble_ke_id_ave_sim,aes("simulation",id_ave),width=.5)+
  geom_point(data=tibble_ke_id_ave_sim,aes("simulation",id_ave))
plot_ke_id_ave

# 出力
setwd(outpath)
png("ke_id_ave_check.png",res=300,width=2000,height=1600)
print(plot_ke_id_ave)
dev.off()


# 次にkeの個人内変動eta2を計算する。そのために各データの残差を計算する。
calc_ke_eta2 = function(data,eta1){
  id_ave = mean(data$ke)
  residuals = log(data$ke)-log(id_ave)
  list(residuals)
}
residuals = unlist(mapply(calc_ke_eta2,tibble_ke2$data))
# eta2自体は正規分布に従っているので普遍標準偏差を計算
ke_eta2_sd = sd(residuals)

# 確認のため分布を再現してみる。
# ke = ke_ave * exp(eta1) * exp(eta2)
# keを乱数生成する関数
# まずke_fixed+eta1で個人のke平均値(ke_id_ave)をけいさんする。その後、ke_id_ave*exp(eta2)により残差を付ける
make_ke = make_ke = function(n_simulations){
  eta1 = rnorm(300,0,ke_eta1_sd)
  eta2 = rnorm(300,0,ke_eta2_sd)
  
  ke_id_ave = ke_fixed *exp(eta1)
  ke_list = c()
  for(i in ke_id_ave){
    ke = ke_id_ave * exp(eta2)
    ke_list = c(ke_list,ke)
  }
  ke_list
}


tibble_ke_simulate = tibble(
  ke = make_ke(n_simulations)
)


ke_simulate = make_ke(n_simulations) %>% 
  sort()

# 10%,50%,90%の区間のデータを取り出す
tibble_ke_p = tibble(
  category="simulate",
  percent = c("10%","50%","90%"),
  ke = ke_simulate[c(n_simulations*0.1,n_simulations*0.5,n_simulations*0.9)]
)

# 実際の数値と比較する。
plot_ke = ggplot()+
  geom_violin(data=tibble_ke,aes("data",ke))+
  geom_boxplot(data=tibble_ke,aes("data",ke),width=.5)+
  geom_point(data=tibble_ke,aes("data",ke))+
  geom_violin(data=tibble_ke_simulate,aes("simulate",ke))+
  geom_boxplot(data=tibble_ke_simulate,aes("simulate",ke),width=.5)+
  geom_point(data=tibble_ke_simulate,aes("simulate",ke))
plot_ke

# 出力
setwd(outpath)
png("ke_check.png",res=300,width=2000,height=1600)
print(plot_ke)
dev.off()




# ke_eta1とeta2,eta4の関係

tibble_eta1_2_4 = tibble(
  ID = tibble_ke2$ID,
  ke_eta1 = tibble_ke2$eta1,
  A_eta2 = tibble_160_320_ID$eta4,
  A_eta4 = tibble_u160_ID$eta2
)

ggplot(tibble_eta1_2_4,aes(ke_eta1,A_eta2))+
  geom_point()+
  stat_smooth(method="lm")+
  stat_poly_eq(formula=y~x,
               aes(label=paste(stat(eq.label),
                               stat(adj.rr.label),
                               stat(AIC.label),
                               sep="~~~")),
               parse=TRUE)

ggplot(tibble_eta1_2_4,aes(ke_eta1,A_eta4))+
  geom_point()+
  stat_smooth(method="lm")+
  stat_poly_eq(formula=y~x,
               aes(label=paste(stat(eq.label),
                               stat(adj.rr.label),
                               stat(AIC.label),
                               sep="~~~")),
               parse=TRUE)

# keのeta1からA_eta2, A_eta4を計算するようにしてみる。
# A_eta2 = A_eta2_a * ke_eta1 + A_eta2_b + A_eta2_eta1
# A_eta2 = A_eta4_a * ke_eta1 + A_eta4_b + A_eta4_eta1
lm_ke_eta2 = lm(A_eta2~ke_eta1,tibble_eta1_2_4)
A_eta2_a = lm_ke_eta2$coefficients[2]
A_eta2_b = lm_ke_eta2$coefficients[1]
A_eta2_eta1_sd = sd(lm_ke_eta2$residuals)

lm_ke_eta4 = lm(A_eta4~ke_eta1,tibble_eta1_2_4)
A_eta4_a = lm_ke_eta4$coefficients[2]
A_eta4_b = lm_ke_eta4$coefficients[1]
A_eta4_eta1_sd = sd(lm_ke_eta4$residuals)


# パラメータを出力する。
tibble_out = tibble(
  A_a = A_a,
  A_b = A_b,
  A_c = A_c,
  A_eta2_sd = A_eta2_sd,
  A_eta3_sd = A_eta3_sd,
  A_eta4_sd = A_eta4_sd,
  ke_fixed = ke_fixed*0.5825,
  ke_eta1_sd = ke_eta1_sd,
  ke_eta2_sd = ke_eta2_sd,
  A_eta2_a = A_eta2_a,
  A_eta2_b = A_eta2_b,
  A_eta2_eta1_sd = A_eta2_eta1_sd,
  A_eta4_a = A_eta4_a,
  A_eta4_b = A_eta4_b,
  A_eta4_eta1_sd = A_eta4_eta1_sd
)

setwd(outpath)
write.csv(tibble_out,"estimate_params.csv",row.names=FALSE)


