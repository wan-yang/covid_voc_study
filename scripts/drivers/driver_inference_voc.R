# to run model-inference of three VOCs (B.1.1.7, B.1.351, P1) using data from the UK, S.Africa, and Brazil
# 5/5/2021


num_runs = 2  # number of runs
tno = 1:num_runs

num_ens = 500 # number of ensemble members
epi.model = 'SEIRSV' # susceptible-exposed-infectious-recovered-susc

# tag.evals = c('eq','rank', 'obs.more', 'obs.most', 'obs.only') # 'eq','rank', 'obs.only' - not good
tag.evals = c('obs.more', 'obs.most','obs.comb')


if(T){
  dir_data = '../../data/'
  dir_code = '../model_inference/'
  
  dir_res = paste0('../../results/voc/')
  
  
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(xlsx)
  library(MMWRweek)
  library(readr); library(readxl); library(writexl); library(stringi)
  library(lemon)
  library(tgp)
  library(scales) # for transparency
  library(R0)
  library(msm) # for the deltamethod function to compute variance
  
}

source(paste0(dir_code,'setSRparm.R'))
source(paste0(dir_code,'SEIRS.R'))
source(paste0(dir_code,'EAKF.R'))
source(paste0(dir_code,'get_relR0.R'))
source(paste0(dir_code,'set_tm2event.R'))

if(F){
  # if you're using a cluster for computing
  args=commandArgs(TRUE)
  tno=as.integer(args[1])
  
}

if(! file.exists(dir_res))  dir.create(dir_res,recursive = T)


da = read.csv(paste0(dir_data,'da_case.death.mob_uk.sa.br.csv'), stringsAsFactors = F)  %>% data.table()
da = da[!is.na(case.uk)]
da[date < as.Date('2020-02-09') & is.na(da)] = 0
# set population size to 
N = 1e6; # per 1 M
num_gr = num_obs = length(N); # no age structure



stoch = T


seasonality = T; # include seasonality

# parms from flu hk paper, AH/T model, from Yuan et al. 2021
Rwea_parm.bounds = rbind(
  c(2.34, 2.93), # R0max
  c(.86,1.18), # R0diff - this determines the magnitude of seasonality
  c(2.2,4.0)/1000, # qmin
  c(17,20)/1000, # qmax
  c(10.2,11.3)/1000, # qmid
  c(20.2,24), # Tc
  c(.4, 5.1), # Tdiff
  c(.95,1.54) # Texp
)
# use the mean instead to reduce uncertainty?
Rwea_parm.bounds = cbind(rowMeans(Rwea_parm.bounds),rowMeans(Rwea_parm.bounds))
rownames(Rwea_parm.bounds) = c('R0max','R0diff','qmin','qmax','qmid','Tc','Tdiff','Texp')


source(paste0(dir_code,'set_tm2event.R'))


seed = .1

# SR
doSR = T
percSRmajor = .1
SR.perc = .05
SR.perc.local = .03;
SR.perc.full = 0.03
SR.var.local= c('beta','Tei','Tir','Trs','Td.mean','Td.sd','p.mob','alpha','ifr')
SR.var.full= c('beta','Tir','ifr') # ,'alpha','ifr'
SR.var.tx = c('beta','Tir')
# do three diff categories: all combined, ms/mh, cumc
donotUpdateS1stWave = T # do not allow the filter to update S during first wave
rednUpdateEI = T # do not allow or reduce the level allowed, the filter to update E or I - it takes OEVr all efforts

{
  fn_getWkLowIFR = function(idx.all){
    consecutive =  idx.all[-1] - idx.all[-length(idx.all)]
    i.div = which(consecutive >1)
    if(length(i.div) > 0){
      w.div = idx.all[which(consecutive >1)]
      grps = list()
      if(length(i.div) == 1){
        grps[[1]] = idx.all[1]: w.div[1]
        grps[[2]] =idx.all[i.div[1]+1]: tail(idx.all,1)
      } else {
        for(id in 1: length(i.div)){
          if(id == 1){
            grps[[id]] = idx.all[1]: w.div[id]
          } else if (id == length(i.div)){
            # both before and after
            
            grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]]
            grps[[id+1]] = idx.all[i.div[id]+1]: tail(idx.all,1)
          } else {
            grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]] # (idx.all[i.div[id-1]+1]) : w.div[id]
          }
          
        }
      }
      
      dur = lapply(grps, length) %>% unlist
      igrs = which(dur > 5)
      
      WkLowIFR = NULL
      if(length(igrs) > 0){
        for(igr in igrs){
          WkLowIFR = c(WkLowIFR, grps[[igr]])
        }
      }
    } else {
      WkLowIFR = idx.all
    }
    
    WkLowIFR
  }
}

locs = c('uk', 'sa','br')
for(loc.t in locs){  # 'uk', 'sa','br'
  print(loc.t)
  
  da.t = da[,c('date','year','week',paste0(c('case.','death.','mob.bus.','mob.full.'),loc.t)),with=F] %>% 
    setnames(paste0(c('case.','death.','mob.bus.','mob.full.'),loc.t), c('case','death','mob.bus','mob.full'))
  # only when it is > 2 cases per M
  da.t = da.t[case > 2]
  da.t = da.t[complete.cases(da.t)]
  
  
  if(loc.t == 'uk'){
    IsLargeCountry = F
    excludeLockDown = T
    # 3rd lockdown, starting 1/6/21. some reopening in March
    # https://www.legislation.gov.uk/uksi/2021/8/pdfs/uksi_20210008_en.pdf
    wkLockDown = which(as.Date(da.t$date) %in% seq(as.Date('2021/1/6')+14, as.Date('2021/3/1'), by = 'day'))
    # for vaccination
    vax.start = as.Date('2020/12/08')
    massvax.start = vax.start
    VE1 = .85
    VE2 = .95  # higher b/c we are using mortality data too
    
    redn.priority = runif(num_ens, 1/3, 1/2) # lower risk of infection for priority groups?
    
    seed_max = da.t$case[1] * 100 # 10000
    end1stWave = which(da.t$date == as.Date("2020-08-16")) # 30 - 5; # end of first wave, after which S can be updated by the filter
    main2ndWave = which(da.t$date == as.Date('2020-11-15')) # 38
    
    # identify weeks with lower IFR
    wk.WkLowIFR = which(as.Date(da.t$date) %in% seq(as.Date('2020/6/1'), as.Date('2020/10/1'), by = 'day'))
    # identify weeks with higher IFR
    wk.2strtHighIFR = which(as.Date(da.t$date) %in% seq(as.Date('2020/10/1'), length.out = 7, by = 'day'))
    wk.2strtHigherIFR = which(as.Date(da.t$date) %in% seq(as.Date('2021/1/1'), length.out = 7, by = 'day'))
   
    tm_rednUpdateEI = which(da.t$date == as.Date('2020-12-06')) + 0:3 #  41:45 # main2ndWave + 0: 2 # 
    tm_largerVar = 5 # number of initial weeks to have larger OEV
    pOEV = 1
    
    ifr_bounds = c(.1, 1.5) / 100
    # increasing detection
    alpha_bounds = c(.01, .15);# c(.01, .2); # reporting rate
    alpha_bounds2 = c(.1, .3); # reporting rate, later wave to reflect increase in detection rate
    
    DAalpha_bounds2 = c(.05, .5)  # for s.only, c(.1, .5) is too low
    SRalpha_bounds2 = c(.1, .3) # fall - started to increase
    SRalpha_bounds3 = c(.15, .3)
    wk.summer = which(as.Date(da.t$date) %in% seq(as.Date('2020/06/30'), as.Date('2020/9/1'),by='day'))
    SRalpha_bounds.summer = c(0.05, .2) # summer lower detection
    DAalpha_bounds3 = c(.1, .5)
    
    SRifr_bounds1 = c(.01, 1) / 100 # too low?
    # SRifr_bounds2 = c(.01, 1) / 100  # lower for 2nd wave - may be too high
    SRifr_bounds2 = c(.01, .5) / 100 #  early phase
    SRifr_bounds3 = c(.01, .6) / 100 # later phase
    SRifr_bounds.WkLowIFR = c(.01, .3) / 100 # 
    
    beta_bounds = c(.5, .8) # transmission rate
    p.mob_bounds = c(.5, 1.5); # scaling for mobility
    Td.mean_bounds =  c(5,7) # mean Td: reporting delay
    Td.sd_bounds = c(1,3) # Td, sd: reporting delay sd
    
  } else if(loc.t == 'br'){
    
    IsLargeCountry = T # allow slower changes in S b/c likely traveling wave
    
    # for vaccination
    vax.start = as.Date('2021/1/18')
    massvax.start = as.Date('2021/6/1') # very slow rollout
    VE1 = .45
    VE2 = .55
    
    seed_max = da.t$case[1] * 50 # 10000
    end1stWave = which(da.t$date == as.Date("2020-10-04")) # 30; # end of first wave, after which S can be updated by the filter
    main2ndWave = which(da.t$date == as.Date('2020-11-8')) # 38
    
    c2d = da.t$death[-c(nrow(da.t)-0:1)] / da.t$case[-c(1,2)]
    c2d[is.infinite(c2d)] = max(c2d[!is.infinite(c2d)])
    idx.all = which(c2d[6:main2ndWave] < quantile(c2d[6: (main2ndWave)], probs = .5))
    idx.all = min(idx.all) : max(idx.all)
    # are they consecutive change
    # fn_getWkLowIFR(idx.all)
    wk.WkLowIFR = 5 + fn_getWkLowIFR(idx.all)
    da.t$date[wk.WkLowIFR]
    
    # get the week with increasing IFR during 2nd wave
    # quantile(c2d[main2ndWave: length(c2d)], probs = .5)
    wk.2strtHighIFR = which(as.Date(da.t$date) %in% seq(as.Date('2021/1/1'), length.out = 7, by = 'day'))
    
    # even higher IFR when healthcare systems collapse
    wk.2strtHigherIFR = which(as.Date(da.t$date) %in% seq(as.Date('2021/2/1'), length.out = 7, by = 'day'))
    
    
    tm_rednUpdateEI = main2ndWave + 0: 3 # 35:37
    tm_largerVar = 5 # number of initial weeks to have larger OEV
    pOEV = 2 # LARGER OEV B/C NOISY DATA
    
    ifr_bounds = c(.01, 1) / 100  # over estimate early deaths
    SRifr_bounds1 = c(.01, .7) / 100 # too low?
    SRifr_bounds2 = c(.01, .6) / 100  # early for 2nd wave
    SRifr_bounds3 = c(.01, 1) / 100  # later for 2nd wave
    SRifr_bounds.WkLowIFR = c(.01, .35) / 100
    
    
    # 4/18/21
    alpha_bounds = c(.01, .1); # reporting rate
    alpha_bounds2 = c(.02, .1); # reporting rate, later wave to reflect increase in detection rate
    
    
    DAalpha_bounds2 = c(.02, .3)
    SRalpha_bounds2 = c(.02, .1)
    
    # beta_bounds = c(.5, .8)  - .1 # too low ?
    beta_bounds = c(.5 - .1, .8)  # slightly lower than the UK, performs a bit better
    p.mob_bounds = c(.5, 1.5); # scaling for mobility
    Td.mean_bounds =  c(5,8) # mean Td: reporting delay
    Td.sd_bounds = c(1,3) # Td, sd: reporting delay sd
    
  } else if(loc.t == 'sa'){
    IsLargeCountry = F
    vax.start = as.Date('2021/2/18')  # mid Feb 2021, adding the lag
    massvax.start = as.Date('2021/6/1') # very slow rollout
    VE1 = .6
    VE2 = .7
    
    seed_max = da.t$case[1] * 20 # 10000
    end1stWave = which(da.t$date == as.Date("2020-08-30")) # 25; # end of first wave, after which S can be updated by the filter
    main2ndWave = which(da.t$date == as.Date('2020-11-1')) # 38
    
    c2d = da.t$death[-c(nrow(da.t)-0:1)] / da.t$case[-c(1,2)]
    c2d[is.infinite(c2d)] = max(c2d[!is.infinite(c2d)])
    c2d[is.na(c2d)] = 0
    idx.all = which(c2d[6:main2ndWave] < quantile(c2d[6: (main2ndWave)], probs = .33))
    idx.all = min(idx.all) : max(idx.all)
    # are they consecutive change
    # fn_getWkLowIFR(idx.all)
    wk.WkLowIFR = 5 + fn_getWkLowIFR(idx.all)
    da.t$date[wk.WkLowIFR]
    
    # get the week with increasing IFR during 2nd wave
    wk.2strtHighIFR = which(as.Date(da.t$date) %in% seq(as.Date('2020/12/10'), length.out = 7, by = 'day'))
    
    tm_rednUpdateEI = which(da.t$date == as.Date('2020-11-22')) + 0: 3 # 39:41
    tm_largerVar = 10 # number of initial weeks to have larger OEV - low numbers in the first ~10 wks
    pOEV = 1.5
    
    
    # ifr_bounds = c(.01, 1) / 100  # over estimate early deaths
    ifr_bounds = c(.01, .3) / 100
    SRifr_bounds1 = c(.01, .2) / 100  # lower for 1st wave, later weeks
    SRifr_bounds2 = c(.01, .2) / 100  # lower for 2nd wave
    SRifr_bounds.WkLowIFR = c(.01, .15) / 100
    
    alpha_bounds = c(.01, .06); # reporting rate
    alpha_bounds2 = c(.02, .08); # reporting rate, later wave to reflect increase in detection rate
    
    DAalpha_bounds2 = c(.02, .3)
    SRalpha_bounds2 = c(.02, .08)
    
    beta_bounds = c(.5, .8) - .1  # slightly lower than the UK, performs better
    p.mob_bounds = c(.5, 1.5); # scaling for mobility
    Td.mean_bounds =  c(5,8) # mean Td: reporting delay
    Td.sd_bounds = c(1,3) # Td, sd: reporting delay sd
  }
  

  
  
  # load vx data for each country
  if(loc.t %in% c('uk','sa', 'br')){
    da.vacc = read.csv(paste0(dir_data,'vx.lagged.per1Mpop_',loc.t,'.csv')) %>% data.table()
    da.vacc$date = da.vacc$date %>% as.Date
  }  else {
    da.vacc = data.table(date = as.Date('2021/6/1'), n.v1=0, n.v2=0)
  }
  
  
  
  
  
  obs_i = (da.t$case) %>% as.matrix() 
  obs_vars_i = obs_i
  for(j in 1:num_obs){
    tmp=rep(0,nrow(da.t))
    for (i in 3:nrow(da.t)){
      tmp[i]=mean(obs_i[(i-2):(i-0),j]);
    }
    obs_vars_i[,j]= (c(rep(N/1000,tm_largerVar),rep(N/1e3,nrow(da.t)-tm_largerVar)) + (tmp^2)/50) * pOEV;
    
  }
  # obs_vars_i = ((obs_i + 10) * 4)  %>% as.matrix() 
  
  obs_d = (da.t$death) %>% as.matrix() 
  obs_vars_d = obs_d
  for(j in 1:num_obs){
    tmp=rep(0,nrow(da.t))
    for (i in 3:nrow(da.t)){
      tmp[i]=mean(obs_d[(i-2):(i-0),j]);
    }
    
    if(loc.t == 'uk'){
      
      # lower, b/c more accurate docummentation
      obs_vars_d[,j]= (c(rep(N/1e4,tm_largerVar),rep(N/5e4,nrow(da.t)-tm_largerVar)) + 
                         pmin((tmp^2)/10, tmp*5)
      ) * pOEV;
    } else {
      obs_vars_d[,j]= (c(rep(N/1e4,tm_largerVar),rep(N/1e4,nrow(da.t)-tm_largerVar)) + 
                         pmin((tmp^2)/10, tmp*20)
      ) * pOEV;
    }
  }
  
  # get relative mobility
  # rel.mob = lhs(num_ens, rect = cbind(1 + da.t$mob.bus/100, 1+da.t$mob.full/100))  %>% t # as.matrix() 
  rel.mob = lhs(num_ens, rect = cbind(1 + da.t$mob.bus/100, 1+da.t$mob.bus/100))  %>% t # as.matrix() 
  tmp = cbind(1 + da.t$mob.bus/100, 1+da.t$mob.full/100, (1+pmax(da.t$mob.full * 1.5, da.t$mob.bus)/100))
  
  # get relative R0 for seasonality
  relR0 = fn_getRelR0(loc.t, ref.wk = da.t$week[1], Rwea_parm.bounds=Rwea_parm.bounds)
  
  
  weeks = da.t$week # for seasonality if applicable
  Week.starts = da.t$date
  
  # do 5 runs
  for(ir in tno){
    
    print(paste('run', ir))
    
    
    So=t(lhs(num_ens,rect = rbind(cbind(.99, 1) * N, # S0
                                  cbind(seed_max/20,seed_max/2), # E0
                                  cbind(seed_max/20,seed_max/2), # I0
                                  cbind(0,seed_max/100) # deaths0
    )))
    
    
    S0 = So[1:num_gr,,drop=F]
    E0 = So[1:num_gr+num_gr,,drop=F]
    I0 = So[1:num_gr+num_gr*2,,drop=F]
    D0 = So[1:num_gr+num_gr*3,,drop=F]
    
    newItot = I0; newIobs = I0; 
    rownames(S0)=paste0('S',1:num_gr); 
    rownames(E0)=paste0('E',1:num_gr); 
    rownames(I0)=paste0('I',1:num_gr); 
    rownames(D0)=paste0('death',1:num_gr);
    rownames(newItot)=paste0('newItot',1:num_gr); 
    rownames(newIobs)=paste0('newIobs',1:num_gr); 
    
    
    imm_bounds = c(2, 3) * 365
    
    
    
    
    parm.bounds = rbind(beta_bounds, # beta for all loc's
                        c(2,5), # Tei: time from exposed to infectious: incubation time mean = 4
                        c(2,5), # Tir: time from infectous to not (remOEVd)
                        imm_bounds, # immunity period, Trs
                        Td.mean_bounds, # c(5,7), # mean Td: reporting delay
                        Td.sd_bounds, # c(1,3), # Td, sd: reporting delay sd
                        p.mob_bounds, # scaling for mobility
                        alpha_bounds, # reporting rate
                        ifr_bounds # infection fatality risk
    )
    parm.names = c('beta','Tei','Tir','Trs', 'Td.mean', 'Td.sd', 'p.mob','alpha', 'ifr')
    
    rownames(parm.bounds) = parm.names
    parm.bounds
    
    parm0=t(lhs(num_ens,parm.bounds)); rownames(parm0)=rownames(parm.bounds)
    
    
    STATE0=rbind(S0, E0, I0, D0, newIobs, newItot, parm0)
    state.names=rownames(STATE0)
    idx.obs_i= which(state.names == 'newIobs1')  # the random tests are testing the prevalence of infectious - I
    idx.obs_d= which(state.names == 'death1')  # the random tests are testing the prevalence of infectious - I
    idx.newItot = which(state.names == 'newItot1') 
    idx.e = which(state.names == 'E1') 
    idx.i = which(state.names == 'I1') 
    
    num_state = 4 + 2
    
    DA.bounds = rbind(matrix(c(rep(0,num_state * num_gr), rep(N,num_state)),num_state * num_gr,2),
                      cbind(parm.bounds[,1]*.5, parm.bounds[,2]*1.5)) # cbind(parm.bounds[,1]*.5, parm.bounds[,2]*1.5)
    rownames(DA.bounds)=state.names
    DA.bounds[c('E1','I1'),2] = N * .15 # / 20
    DA.bounds['death1',2] = N / 100 # 200
    DA.bounds['S1',1] = N / 10
    DA.bounds[c('newItot1','newIobs1'),2] = N * .2
    
    DA.bounds['Td.mean',] = parm.bounds['Td.mean',]
    DA.bounds['Trs',1] = 200
    
    
    DA.bounds['beta',2] = parm.bounds['beta',2]*2
    DA.bounds['alpha',1] = parm.bounds['alpha',1]*.75
    DA.bounds['alpha',2] = parm.bounds['alpha',2]*2
    DA.bounds['ifr',] = c(parm.bounds['ifr',1]*.25, parm.bounds['ifr',2]*1.25)
    DA.bounds['ifr',1] = pmin(DA.bounds['ifr',1], 1e-4)
    
    
    SR.bounds = parm.bounds # 
    rownames(SR.bounds)=parm.names
    SR.bounds['alpha',] = parm.bounds['alpha',] 
    SR.bounds['ifr',] = parm.bounds['ifr',] 
    SR.bounds0 = SR.bounds
    
    # FOR TESTING DIFF HYPOTHESES
    SR.bounds.wider = cbind(parm.bounds[,1] * 1.2, parm.bounds[,2]*1.2) # higher than normal, for increase
    # SR.bounds.wider = cbind(parm.bounds[,1] * 1.1, parm.bounds[,2]*1.2) # higher than normal, for increase
    rownames(SR.bounds.wider)=parm.names
    SR.bounds.wider['alpha',] = SR.bounds['alpha',] # * 1.2 #  + .05
    SR.bounds.wider['p.mob',] = SR.bounds['p.mob',] # c(.7, 1.3); # 
    SR.bounds.wider['Tei',] = parm.bounds['Tei',]
    # SR.bounds.wider['Tir',2] = parm.bounds['Tir',2] * 1.2
    SR.bounds.wider['Tir',1] = parm.bounds['Tir',1] #  * 1.2
    SR.bounds.wider['Tir',2] = parm.bounds['Tir',2] * 1.1
    SR.bounds.wider['ifr',] = SR.bounds['ifr',]
    
    SR.bounds.wider2 = cbind(parm.bounds[,1] * 1.3, parm.bounds[,2]*1.4) # higher than normal, for increase
    rownames(SR.bounds.wider2)=parm.names
    SR.bounds.wider2['alpha',] = SR.bounds['alpha',] # * 1.2 # + .05
    SR.bounds.wider2['p.mob',] = SR.bounds['p.mob',] # c(.7, 1.3); # 
    SR.bounds.wider2['Tei',] = parm.bounds['Tei',]
    
    SR.bounds.wider2['Tir',1] = parm.bounds['Tir',1] # * 1.2
    SR.bounds.wider2['Tir',2] = parm.bounds['Tir',2]  * 1.2
    SR.bounds.wider2['ifr',] = SR.bounds['ifr',]
    
    tm.ini=1; tmstep=7; newI.previous = NULL; inflat=1.03; state0=STATE0
    
    severity['death',] = STATE0['ifr',]
    
    res.train = EAKF(epi.model=epi.model, num_ens=num_ens,inflat=1.03, 
                     obs_i=obs_i, obs_vars_i=obs_vars_i, # case
                     obs_d=obs_d, obs_vars_d=obs_vars_d,
                     weeks=weeks,Week.starts=Week.starts,
                     parm.bounds=parm.bounds, DA.bounds=DA.bounds, SR.bounds=SR.bounds, 
                     parm.names = rownames(parm.bounds), rel.mob = rel.mob,
                     state0=STATE0, state.names=rownames(STATE0),
                     severity = severity,
                     tm.ini=1, tmstep=7,
                     newI.previous = NULL,
                     SRparms = SRparms
    )
    
    save(res.train, file = paste0(dir_res, loc.t,'_train_r',ir,'.RData'))
    if(ir == 1)
      save(SRparms, file = paste0(dir_res, 'SRparms.RData'))
  }
  
} # 'uk', 'sa', 'br'

