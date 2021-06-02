# 3/8/21 - to run the synthetic truth


num_runs = 2  # number of runs
tno = 1:num_runs

num_ens = 500 # number of ensemble members
epi.model = 'SEIRSV' # susceptible-exposed-infectious-recovered-susc

seed1stWave = 2
seed2ndWave = 50
alpha.t = .1

alphas = c(.1, .2)

alpha.t = alphas[1] # Detection rate - change to .1 or .2 to get the synthetic data

end1stWave = 20; # end of first wave, after which S can be updated by the filter

stoch = T

# tag.evals = c('eq','rank', 'obs.more', 'obs.most', 'obs.only') # 'eq','rank', 'obs.only' - not good
tag.evals = c('obs.more', 'obs.most','obs.comb')

if(T){
  dir_data = '../../data/'
  dir_code = '../model_inference/'
  
  dir_res = paste0('../../results/syn.test/test_seed1st',seed1stWave,'seed2nd',seed2ndWave,'alpha',alpha.t,'/')
  
  
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

da = read.csv(paste0(dir_data,'da_syn.truth_seed1st',seed1stWave,'_seed2nd', seed2ndWave,'_alpha',alpha.t,'.csv'), stringsAsFactors = F)  %>% data.table()
da[date < as.Date('2020-02-09') & is.na(da)] = 0
tmp = MMWRweek(da$date)
da$year = tmp[,1]
da$week = tmp[,2]
# set population size to 
N = 1e6; # per 1 M
num_gr = num_obs = length(N); # no age structure


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





seed = .05  # default weekly seeding

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
IsLargeCountry = F

loc.t = 'sce5'

for(loc.t in paste0('sce', 1:5)){  # 5 different synthetic test
  print(loc.t)
  
  
  
  da.t = da[,c('date','year','week',paste0(c('case.','death.','mob.','relR0.','wt.suscept.','newitot.','Etot.','Itot.'),loc.t)),with=F] %>% 
    setnames(paste0(c('case.','death.','mob.','relR0.','wt.suscept.','newitot.','Etot.','Itot.'),loc.t), c('case','death','relMob','relR0','S','newitot','E','I'))
  
  tm_largerVar = which(da.t[1:15]$case < 100) %>% length # number of initial weeks to have larger OEV
  
  relR0.t = da.t$relR0[1]
  relR0 = read.csv(paste0(dir_data,'relR0_',relR0.t,'.csv')) %>% data.table()
  relR0 = relR0$rel.R0 %>% matrix(nrow = nrow(relR0), ncol = num_ens)
  
  {
    # for vaccination - not included in the synthtic test here
    vax.start = as.Date('2021/12/08')
    massvax.start = vax.start
    VE1 = .7
    VE2 = .85
    
    seed_max = pmax(mean(da.t$case[1]) / alpha.t, seed1stWave*2) # 10000
    
    # for both alpha - get the weeks seeding initial increases at the begining of 2nd wave
    if(loc.t == 'sce1'){
      tm_rednUpdateEI = 34:37
      main2ndWave = 34
    } else if (loc.t == 'sce2'){
      tm_rednUpdateEI = 32:35
      main2ndWave = 32
    } else if (loc.t == 'sce3'){
      tm_rednUpdateEI = 31:34
      main2ndWave = 31
    } else if (loc.t == 'sce5'){
      tm_rednUpdateEI = 30:33
      main2ndWave = 30
    } else if (loc.t == 'sce4'){
      tm_rednUpdateEI = 34:37
      main2ndWave = 34
    } 
    
    # tm_largerVar = 5 # number of initial weeks to have larger OEV
    pOEV = 1 # LARGER OEV If NOISY DATA
    
    # initial prior range for IFR and other key parameters
    ifr_bounds = c(.3, 1.1) / 100
    # detection could be decreasing during large 2nd wave
    
    
    alpha_bounds = alpha.t + c(-.05, .05) # c(.05, .15); # reporting rate
    alpha_bounds2 = alpha.t + c(-.05, .05) # c(.05, .15); # reporting rate, later wave to reflect increase in detection rate, depending on real situation
    
    DAalpha_bounds2 = alpha.t * c(.5, 2) # c(.05, .2)
    SRalpha_bounds2 = alpha.t + c(-.05, .05) # c(.05, .15);
    
    beta_bounds = c(.5, .8) 
    p.mob_bounds = c(.5, 1.5); # scaling for mobility
    
  } 
  
  
  
  da.vacc = data.table(date = as.Date('2021/6/1'), n.v1=0, n.v2=0)
  
  
  
  
  obs_i = (da.t$case) %>% as.matrix() 
  obs_vars_i = obs_i
  for(j in 1:num_obs){
    tmp=rep(0,nrow(da.t))
    for (i in 3:nrow(da.t)){
      tmp[i]=mean(obs_i[(i-2):(i-0),j]);
    }
    
    obs_vars_i[,j]= (c(rep(N/1000,tm_largerVar),rep(N/1e3,nrow(da.t)-tm_largerVar)) + (tmp^2)/50) * pOEV;
    
  }
  
  obs_d = (da.t$death) %>% as.matrix() 
  obs_vars_d = obs_d
  for(j in 1:num_obs){
    tmp=rep(0,nrow(da.t))
    for (i in 3:nrow(da.t)){
      tmp[i]=mean(obs_d[(i-2):(i-0),j]);
    }
   
    obs_vars_d[,j]= (c(rep(N/1e4,tm_largerVar),rep(N/1e4,nrow(da.t)-tm_largerVar)) + 
                       pmin((tmp^2)/10, tmp*20)
    ) * pOEV;
    
  }
  
  # get relative mobility
  # rel.mob = lhs(num_ens, rect = cbind(1 + da.t$mob.bus/100, 1+da.t$mob.full/100))  %>% t # as.matrix() # add noise if prefer
  rel.mob = matrix(da.t$relMob, nrow(da.t), num_ens) # as.matrix() 
  
  
  weeks = da.t$week # for seasonality if applicable
  Week.starts = da.t$date
  
  
  
  # LOOP
  # do X runs
  for(ir in tno){
    print(paste('run', ir))
    
    
    So=t(lhs(num_ens,rect = rbind(cbind(1, 1) * N, # S0
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
                        c(5,7), # mean Td: reporting delay
                        c(1,3), # Td, sd: reporting delay sd
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
    DA.bounds[c('E1','I1'),2] = N * .15 
    DA.bounds['death1',2] = N / 100 
    DA.bounds['S1',1] = N *.1 # 10  5% too low
    
    DA.bounds['Trs',1] = 200 # allow it to be lower
    
    DA.bounds['beta',2] = parm.bounds['beta',2]*2
    DA.bounds['alpha',1] = pmin(parm.bounds['alpha',1]*.75, .04)
    DA.bounds['ifr',] = c(parm.bounds['ifr',1]*.5, parm.bounds['ifr',2]*1.5)
    
    # SR.bounds = cbind(parm.bounds[,1]*.75, parm.bounds[,2]*1.25)
    SR.bounds = parm.bounds; # cbind(parm.bounds[,1]*.9, parm.bounds[,2]*1.1)
    rownames(SR.bounds)=parm.names
    SR.bounds['alpha',] = parm.bounds['alpha',] 
    SR.bounds['ifr',] = parm.bounds['ifr',] 
    SR.bounds0 = SR.bounds
    
    # FOR TESTING DIFF HYPOTHESES
    SR.bounds.wider = cbind(parm.bounds[,1] * 1.1, parm.bounds[,2]*1.2) # higher than normal, for increase
    rownames(SR.bounds.wider)=parm.names
    SR.bounds.wider['alpha',] = SR.bounds['alpha',] # * 1.2 #  + .05
    SR.bounds.wider['p.mob',] = SR.bounds['p.mob',]
    SR.bounds.wider['ifr',] = SR.bounds['ifr',]
    SR.bounds.wider['Tei',] = parm.bounds['Tei',]
    SR.bounds.wider['Tir',1] = parm.bounds['Tir',1] #  * 1.2
    SR.bounds.wider['Tir',2] = parm.bounds['Tir',2] * 1.1
    
    SR.bounds.wider2 = cbind(parm.bounds[,1] * 1.3, parm.bounds[,2]*1.4) # higher than normal, for increase
    rownames(SR.bounds.wider2)=parm.names
    SR.bounds.wider2['alpha',] = SR.bounds['alpha',] # * 1.2 # + .05
    SR.bounds.wider2['p.mob',] = SR.bounds['p.mob',]
    SR.bounds.wider2['ifr',] = SR.bounds['ifr',]
    SR.bounds.wider2['Tei',] = parm.bounds['Tei',]
    SR.bounds.wider2['Tir',1] = parm.bounds['Tir',1] # * 1.2
    SR.bounds.wider2['Tir',2] = parm.bounds['Tir',2]  * 1.2
    
    
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
  
} # 

