# 阻害率から改善率のsimulateをするスクリプト
# 流れ
# 1. 血中濃度から阻害率を計算 ここは管轄外。阻害率までは計算できているという体。
# 2. slopeとEC50を乱数生成(slope~N(-6.167064,1.281413), EC50~N(0.6909326,0.01880815)) 実際使うのがSDなので分散の代わりにSDを書いている。
# 3. 阻害率に応じて残差のSDを算出し、残差を乱数生成
# 4. ロジスティックモデルとその残差からデータを生成

library(tidyverse)

basepath = dirname(rstudioapi::getSourceEditorContext()$path)
datapath = file.path(basepath,"output")
outpath = datapath = file.path(basepath,"output")


# 阻害率から残差のSDを計算する関数
residuals_SD = function(x){ # xは阻害率
  a = 0.223
  b = -0.0396
  y = a*x+b
  if(y<0){
    y=0
  }
  return(y)
}

# 残差を計算する関数
residuals = function(residuals_SD){
  return(rnorm(1,0,residuals_SD))
}

# slope, EC50と残差からVAS50の改善率を計算する関数
inh_to_VAS50 = function(inh_ave,slope,EC50,residual){
  y = 1/(1+exp(slope*(log(inh_ave)-log(EC50)))) # モデルから計算
  y = y+residual # 残差をつける
  
  # yが0以下の場合は0に
  if(y<0){
    y=0
  }
  return(y)
}

# 阻害率のデータが入ったデータフレーム
# ID, TIME, inh_aveを列名に。

# デモとしてPhase1 bのデータを用いる。
setwd(datapath)
tibble_inh = read_csv("inh_VAS50.csv",col_types=cols())


# VAS50を各データについて50個計算する。
simulate_VAS50 = function(inh_ave,number){ # numberは生成したい乱数の数
  # slopeとEC50を乱数生成
  slope = rnorm(number,-6.167064,1.281413)
  EC50 = rnorm(number,0.6909326,0.01880815)
  # 阻害率から残差を計算
  residuals_sd = as.numeric(unlist(mapply(residuals_SD,inh_ave))) # 残差のSD
  residuals = rnorm(number,0,residuals_sd)
  tibble_f1 = tibble(
    slope=slope,
    EC50=EC50,
    residuals = residuals
  )
  VAS50 = as.numeric(unlist(mapply(inh_to_VAS50,inh_ave,
                                   tibble_f1$slope,tibble_f1$EC50,tibble_f1$residuals)))
  tibble_f1 = tibble_f1 %>% 
    mutate(improvement_ratio=VAS50)
  list(tibble_f1)
}

tibble_simulate = tibble_inh %>% 
  mutate(simulate = mapply(simulate_VAS50,tibble_inh$inh_ave,50)) %>% 
  dplyr::select(-VAS50_rate)

tibble_simulate = tibble_simulate %>% 
  group_by(ID,TIME,inh_ave) %>% 
  unnest() %>% 
  ungroup()

# プロットしてみる。実際の測定点と比較
plot_simulate = ggplot()+
  geom_point(data=tibble_simulate,aes(inh_ave,improvement_ratio),size=1)+
  geom_point(data=tibble_inh,aes(inh_ave,VAS50_rate),color="red",size=3.5)+
  theme_minimal()+
  labs(title = "simulateしたデータ(黒)と実際の値(赤)")+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.1))
plot_simulate

setwd(outpath)
png("simulate_VAS50.png",res=300,width=2000,height = 1400)
plot_simulate
dev.off()





