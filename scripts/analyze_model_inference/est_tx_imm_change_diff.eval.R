# est changes in tx and imm
state.names = c('S1', 'E1','I1',
                'death1', 'newIobs1','newItot1','beta',
                'Tei','Tir','Trs','Td.mean','Td.sd',
                'p.mob','alpha','ifr')

for(ir in 1:num_runs){
  
  tmp = try(load(paste0(dir_res, loc.t, '_train_r',ir,'.RData')))
  if(class(tmp) == "try-error")
    next
  
  
  
  Rtx_stats = res.train$Rtx_stats %>% setnames('Week.start','date')
  Rtx_ens = res.train$Rtx_ens 
  percS_stats = res.train$Susceptibility_stats %>% setnames('Week.start','date')
  percS_ens = res.train$Susceptibility_ens
  cumI_stats = res.train$cumIperc_stats %>% setnames('Week.start','date')
  cumI_ens = res.train$cumIperc_ens
  xpost_mean = res.train$xpost_mean
  
  best.hyp.t = res.train$hyp.best_all[eval == tag.eval]
  Rtx_stats = Rtx_stats[eval == tag.eval]
  Rtx_ens = Rtx_ens[eval == tag.eval]
  percS_stats = percS_stats[eval == tag.eval]
  percS_ens = percS_ens[eval == tag.eval]
  cumI_stats = cumI_stats[eval == tag.eval]
  cumI_ens = cumI_ens[eval == tag.eval]
  xpost_mean = xpost_mean[eval == tag.eval]
  
  da.t$date = da.t$date %>% as.Date()
  Rtx_stats$date = Rtx_stats$date %>% as.Date()
  percS_stats$date = percS_stats$date %>% as.Date()
  cumI_stats$date =cumI_stats$date %>% as.Date()
  tda = merge(da.t, Rtx_stats, by ='date')
  idxW1 = 1: end1stWave; idxW2 = (end1stWave + 1): nrow(da.t)
  idxW1main = 4: (end1stWave)  # exclude 1st and last 3 weeks
  idxW2main = (main2ndWave + 3) : nrow(da.t)
  idx1 = which.max(tda[idxW1]$mean)
  # idx2 = which.max(tda[idxW2]$mean) + idxW2[1] - 1
  # exclude weeks with low cases that are less accurate, those may be drifting
  tda2 = tda[idxW2]; tda2=tda2[case > 50]
  idx2 = which.max(tda2$mean) 
  idx2 = which(da.t$date == tda2[idx2]$date)
  # also record the dates
  Tmax1.Rtx = da.t$date[idx1]
  Tmax2.Rtx = da.t$date[idx2]
  
  # Rtx1 = Rtx_stats[idx1 + tm_window, c("mean", 'sd', "median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr"), with=F] %>% colMeans()
  # Rtx2 = Rtx_stats[pmin(idx2 + tm_window, nrow(da.t)), c("mean", 'sd',"median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr"), with=F] %>% colMeans()
  
  # 4/15/21
  # based on the mean during the 1st wave and the mean during the main 2nd wave (due to the new variant)
  Rtx1 = Rtx_stats[idxW1main, c("mean", 'sd', "median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr"), with=F] %>% colMeans()
  Rtx2 = Rtx_stats[idxW2main, c("mean", 'sd',"median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr"), with=F] %>% colMeans()
  
  da.t[idx2]
  
  
  
  dRtx = data.table(x.mean = (Rtx2['mean'] - Rtx1['mean'])/Rtx1['mean'], x.median = (Rtx2['median'] - Rtx1['median'])/Rtx1['median'])
  dRtx
  
  
  # change in %S
  # 4/16/21 simply based on the time interval with the highest and lowest S and the cumulative infection rate in between
  idx1 = end1stWave+6 # which.max(percS_stats[idxW2]$mean) + idxW2[1] - 1  # main2ndWave + 1 # 
  idx2 = which.min(percS_stats[idxW2]$mean) + idxW2[1] - 1
  
  S1 = percS_stats[idx1,c("mean", 'median'), with=F] %>% unlist
  S2 = percS_stats[idx2,c("mean", 'median'), with=F] %>% unlist
  dI = (cumI_stats[idx2,c("mean", 'median'), with=F] - cumI_stats[idx1,c("mean", 'median'), with=F]) %>% unlist
  
  # S.w1 = percS_stats[end1stWave,c("mean", 'median'), with=F] %>% unlist
  # dImm.w1 = cumI_stats[end1stWave,c("mean", 'median'), with=F] %>% unlist
  S.w1 = percS_stats[end1stWave+5,c("mean", 'median'), with=F] %>% unlist
  dI.w1 = cumI_stats[end1stWave+5,c("mean", 'median'), with=F] %>% unlist
  # S.w1 = percS_stats[idx1,c("mean", 'median'), with=F] %>% unlist
  # dImm.w1 = cumI_stats[idx1,c("mean", 'median'), with=F] %>% unlist
  dS.w1 = (100 - S.w1)
  p = dS.w1 / dI.w1 
  # the filter could overshot, use the level from wave 1 to scale it
  ImmLoss = (dI * p - 
               (S1 - S2) # /p # should we adjust dS too?
             ) %>% unlist
  dImm.t = data.table(mean = ImmLoss['mean'], median = ImmLoss['median'])
  
  dImm = (dImm.t / (100 - S.w1))  %>% unlist
  dImm
  
  res = rbind(res, 
              data.table(loc = loc.t, irun = ir, hyp.best = best.hyp.t$hyp.test,  state = 'dRtx', tmax = Tmax2.Rtx, dRtx),
              data.table(loc = loc.t, irun = ir, hyp.best = best.hyp.t$hyp.test,  state = 'dImm', tmax = Tmax2.S, dImm),
              fill = T
  ) 
} # end run
