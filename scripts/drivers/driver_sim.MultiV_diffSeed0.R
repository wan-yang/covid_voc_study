# driver to run diff scenarios
# do not include B.1.526 and B.1.429/427
# 4/29/21

loc.t = 'nycish'  # not exactly the same as nyc but with some modification as specified
vec_ini.seed = c('eqB1351P1','moreP1', 'moreB1351'); 
ini.seed.t = 'eqB1351P1'  # CHANGE HERE USING ONE OF THE ABOVE FOR INITIAL SEEDING

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

dir_data = '../../data/'
dir_code = '../simMultiV/'
dir_code1 = '../model_inference/'
dir_code2 =  '../plot_results/' 

dir_res = paste0('../../results/sim_multiV/',ini.seed.t,'/')


source(paste0(dir_code,'multistrainSEIRSV.R'))
source(paste0(dir_code,'get_health.metrics.R'))
source(paste0(dir_code1,'SEIRS.R'))
source(paste0(dir_code2,'getPlot.R'))

if(! file.exists(dir_res))  dir.create(dir_res, recursive = T)

# load initial conditions
load(file =  paste0('../../data/simMultiV_ini.conds.RData'))


senarios_all = c('Current NPI, same VE to new variants',
                 'Current NPI, high VE to new variants',
             'Current NPI, median VE to new variants',
             'Current NPI, low VE to new variants',
             '25% Less NPI, same VE to new variants',
             '25% Less NPI, high VE to new variants',
             '25% Less NPI, median VE to new variants',
             '25% Less NPI, low VE to new variants', # ,
             '50% Less NPI, same VE to new variants',
             '50% Less NPI, high VE to new variants',
             '50% Less NPI, median VE to new variants',
             '50% Less NPI, low VE to new variants',
             'No NPI slow, same VE to new variants',
             'No NPI slow, high VE to new variants',
             'No NPI slow, median VE to new variants',
             'No NPI slow, low VE to new variants',
             'No NPI fast, same VE to new variants',
             'No NPI fast, high VE to new variants',
             'No NPI fast, median VE to new variants',
             'No NPI fast, low VE to new variants'
             )

sce.tags_all = c('curNPI.sameVE',
                 'curNPI.highVE',
             'curNPI.midVE',
             'curNPI.lowVE',
             '25lNPI.sameVE',
             '25lNPI.highVE',
             '25lNPI.midVE',
             '25lNPI.lowVE', #,
             '50lNPI.sameVE',
             '50lNPI.highVE',
             '50lNPI.midVE',
             '50lNPI.lowVE',
             'noNPIslow.sameVE',
             'noNPIslow.highVE',
             'noNPIslow.midVE',
             'noNPIslow.lowVE',
             'noNPIfast.sameVE',
             'noNPIfast.highVE',
             'noNPIfast.midVE',
             'noNPIfast.lowVE'
             )

iit = c(1:3, 5:7, 9:11, 13:15, 17:19)
senarios = senarios_all[iit]
sce.tags = sce.tags_all[iit]


# MODEL SETTING
num_ens = 100  # use small number for test, increase to 1000 for formal analysis
epi.model = 'multistrainSEIRSVage' # susceptible-exposed-infectious-recovered-susc
fn_epi = get(epi.model)

stoch = T # run stochastically
seasonality = T # include seasonality

# Initial condition, seeding, etc - maybe location specific
# based on 4/27/21 nyc data
if(ini.seed.t == 'moreP1'){
  percSeed.lwrs = c(35, .2, 35, .2, 2); 
  percSeed.uprs = c(45, 1, 45, .6, 5); 
  seed.tag = 'moreP1'
} else if (ini.seed.t == 'moreB1351'){
  percSeed.lwrs = c(35, .2, 35, 2, .2); 
  percSeed.uprs = c(45, 1, 45, 5, .6);
  seed.tag = 'moreB1351'
} else if (ini.seed.t == 'eqB1351P1'){
  percSeed.lwrs = c(35, .2, 35, 2, 2); 
  percSeed.uprs = c(45, 1, 45, 5, 5);
  seed.tag = 'eqB1351P1'
} 
names(percSeed.lwrs) = c('b1526','b1427','b117','b1351','p1')
names(percSeed.uprs) = c('b1526','b1427','b117','b1351','p1')
percSeed.mns = colMeans(rbind(percSeed.lwrs, percSeed.uprs))

seeds = c(1/21, 0, 1/21, 1/21, 1/21, 1/21); # keep seeding the same for b1351 and p1 to focus on initial introduction
names(seeds) = c('wt','b1526','b1427','b117','b1351','p1')

#
date.t = '2021/04/25' %>% as.Date  # the prior week, used to set the initial conditions
date.end = '2021/08/30' %>% as.Date  # end of the simulation

parm0 = t(lhs(num_ens,parm.bounds))
rownames(parm0) = rownames(parm.bounds) 
for(ia in 1:8){ # diff by age group
  
  source(paste0(dir_code,'set_tm2event.R'))
  
  source(paste0(dir_code,'set_age.spec.rates.R'))
  
  # update ifr
  severity['death',] = parm0[paste0('ifr.',numbers[ia]),]
  # update edr - ed covid visit rate
  severity['edr',] = parm0[paste0('edr.',numbers[ia]),]
  
  dist_tm.to.detect = NULL
  for(ii in 1:num_ens){
    tmp = generation.time(dist_tm.to.detect.name,c(parm0[paste0('Td.mean.',numbers[ia]),ii],parm0[paste0('Td.sd.',numbers[ia]),ii]),truncate = tm.to.detect.max)
    dist_tm.to.detect=cbind(dist_tm.to.detect,tmp$GT[-1]); 
  }
  
  dist_tm.to.death = NULL  # time to death
  for(ii in 1:num_ens){
    # tmp = generation.time(dist_tm.to.death.name,c(state0['Td.mean',ii]+diff.dd,state0['Td.sd',ii]+diff.sd2),truncate = tm.to.death.max)
    # do not link it to Td
    tmp = generation.time(dist_tm.to.death.name,c(tm.to.outcome.mn['tm.to.death',ii], tm.to.outcome.sd['tm.to.death',ii]),truncate = tm.to.death.max)
    dist_tm.to.death=cbind(dist_tm.to.death,tmp$GT[-1]); 
  }
  if(tm.to.deathFromDiag){ # if the distribution of time to death is from diagnosis, not infectious
    # add time from infectious to diagnosis
    dist_tm.to.death = rbind(matrix(0,round(tm.to.diag,0),num_ens),dist_tm.to.death)
  }
  
  # for other health metrics
  {
    dist_tm.to.hospital = NULL  # time to hospitalization
    for(ii in 1:num_ens){
      # tmp = generation.time(dist_tm.to.hospital.name,c(state0['Td.mean',ii]+diff.hd,state0['Td.sd',ii]),truncate = tm.to.hospital.max)
      # do not link it to Td
      tmp = generation.time(dist_tm.to.hospital.name,c(tm.to.outcome.mn['tm.to.hospital',ii], tm.to.outcome.sd['tm.to.hospital',ii]),truncate = tm.to.hospital.max)
      
      dist_tm.to.hospital=cbind(dist_tm.to.hospital,tmp$GT[-1]); 
    }
    
    dist_tm.to.icu = NULL  # time to icu
    for(ii in 1:num_ens){
      # tmp = generation.time(dist_tm.to.icu.name,c(state0['Td.mean',ii]+diff.id,state0['Td.sd',ii]+diff.sd),truncate = tm.to.icu.max)
      # do not link it to Td
      tmp = generation.time(dist_tm.to.icu.name,c(tm.to.outcome.mn['tm.to.icu',ii], tm.to.outcome.sd['tm.to.icu',ii]),truncate = tm.to.icu.max)
      
      dist_tm.to.icu=cbind(dist_tm.to.icu,tmp$GT[-1]); 
    }
    
    
    # sample retention time for hospitalization
    dist_hospital.stay = rgamma(num_ens,shape = (tm.to.outcome.mn['hospital.stay',]/tm.to.outcome.sd['hospital.stay',])^2,
                                rate=tm.to.outcome.mn['hospital.stay',]/tm.to.outcome.sd['hospital.stay',]^2)
    # sample retention time for icu
    dist_icu.stay = rgamma(num_ens,shape = (tm.to.outcome.mn['icu.stay',]/tm.to.outcome.sd['icu.stay',])^2,
                           rate=tm.to.outcome.mn['icu.stay',]/tm.to.outcome.sd['icu.stay',]^2)
    dist_hospital.stay = round(pmin(dist_hospital.stay,hospital.stay.max),0) # round it to days
    dist_icu.stay = round(pmin(dist_icu.stay,icu.stay.max),0)
    
    # sampe time of ventilator use
    dist_vent.stay = rgamma(num_ens,shape = (tm.to.outcome.mn['vent.stay',]/tm.to.outcome.sd['vent.stay',])^2,
                            rate=tm.to.outcome.mn['vent.stay',]/tm.to.outcome.sd['vent.stay',]^2)
    dist_vent.stay = round(pmin(dist_vent.stay,vent.stay.max),0)
    
    eval(parse(text = paste('dist_tm.to.hospital.',numbers[ia], '=dist_tm.to.hospital', sep='')))
    eval(parse(text = paste('dist_tm.to.icu.',numbers[ia], '=dist_tm.to.icu', sep='')))
    eval(parse(text = paste('dist_hospital.stay.',numbers[ia], '=dist_hospital.stay', sep='')))
    eval(parse(text = paste('dist_icu.stay.',numbers[ia], '=dist_icu.stay', sep='')))
    eval(parse(text = paste('dist_vent.stay.',numbers[ia], '=dist_vent.stay', sep='')))
  }
  
  eval(parse(text = paste('tm.from.inf.to.death.max.',numbers[ia], '=tm.from.inf.to.death.max', sep='')))
  eval(parse(text = paste('dist_tm.to.detect.',numbers[ia], '=dist_tm.to.detect', sep='')))
  eval(parse(text = paste('dist_tm.to.death.',numbers[ia], '=dist_tm.to.death', sep='')))
  eval(parse(text = paste('severity.',numbers[ia], '=severity', sep='')))
}


ivcomb = 1; isce = 1
ivcomb = 2; isce = 10
for(ivcomb in 1: length(variant.combs)){
  
  variants.t = variant.combs[[ivcomb]]
  Ns= length(variants.t) # ncol(obs_i);  # number of strains included
  v.tag = paste(variants.t, collapse = '.')
  
  NaNs = Na * Ns
  
  for(isce in 1:length(senarios)){ # 1: length(senarios)
    
    sce.tag = sce.tags[isce]
    sce.t = senarios[isce]
    
    
    {
      # VE1 = .8; VE2 = .9
      VE1 = .85; VE2 = .95
      if(grepl('same VE to new variants', sce.t)){
        VE1redn = rep(1, 6); 
        VE2redn = rep(1, 6);  
        
      } else if(grepl('high VE to new variants', sce.t)){
        VE1redn = c(1, rep(.95, 3), .8, .85); 
        VE2redn = c(1, rep(1, 3), .9, .95);  
        
      } else if (grepl('median VE to new variants', sce.t)){
        VE1redn = c(1, rep(.95, 3), .7, .8); 
        VE2redn = c(1, rep(.95, 3), .8, .9);  
        
      } else if (grepl('low VE to new variants', sce.t)){
        VE1redn = c(1, rep(.95, 3), .2, .3); 
        VE2redn = c(1, rep(.95, 3), .3, .4);  
      }
      names(VE1redn) = c('wt','b1526','b1427','b117','b1351','p1')
      names(VE2redn) = c('wt','b1526','b1427','b117','b1351','p1')
      
      
      week.starts = seq(as.Date(date.t)+7, date.end, by = 'week')
      num_times = length(week.starts)
      
      for(voc.t in c('b117','b1351','p1')){  # location data were use
        tmp = eval(parse(text = paste('data.table(dRtx.lwr =', 
                                      est.newv[voc == voc.t & state == 'dRtx']$bound.lwr,
                                      ', dRtx.upr =', 
                                      est.newv[voc == voc.t & state == 'dRtx']$bound.upr,
                                      ', dImm.lwr =', 
                                      est.newv[voc == voc.t & state == 'dImm']$bound.lwr,
                                      ', dImm.upr =', 
                                      est.newv[voc == voc.t & state == 'dImm']$bound.upr,
                                      ')',
                                      sep = '')))
        tmp$cimmb_a.lwr = 1 - tmp$dImm.upr
        tmp$cimmb_a.upr = 1 - tmp$dImm.lwr
        
        # infection of new variants, provide protection to wt
        tmp$cimma_b.lwr = 1
        tmp$cimma_b.upr = 1 
        
        
        tmp$percSeed.lwr = percSeed.lwrs[voc.t]
        tmp$percSeed.upr = percSeed.uprs[voc.t]
        eval(parse(text = paste('adjVparm.', voc.t, '= tmp', sep ='')))
      }
      
      for(voc.t in c('b1526','b1427')){
        tmp = eval(parse(text = paste('data.table(dRtx.lwr = (Rt.incr.', 
                                      voc.t, '- 1) * .8',
                                      ', dRtx.upr = (Rt.incr.', voc.t, '- 1) * 1.2',
                                      ', dImm.lwr = 0',
                                      ', dImm.upr = 0.1',
                                      # ', dImm.upr = 0',
                                      ')',
                                      sep = '')))
        tmp$cimmb_a.lwr = 1 - tmp$dImm.upr
        tmp$cimmb_a.upr = 1 - tmp$dImm.lwr
        
        # infection of new variants, provide protection to wt
        tmp$cimma_b.lwr = 1
        tmp$cimma_b.upr = 1 
        
        
        tmp$percSeed.lwr = percSeed.lwrs[voc.t]
        tmp$percSeed.upr = percSeed.uprs[voc.t]
        eval(parse(text = paste('adjVparm.', voc.t, '= tmp', sep ='')))
      }
      
      num_parm=(4)*Ns * Na +(Ns^2-Ns); # number parms: 4 for the SEIRS (beta, Tei, Tir, Trs); (Ns^2-Ns) for the cross immunity matrix
      num_var=(5+4)*Ns * Na +(Ns^2-Ns); # 5 states: S, E, I, (newI, cases);  4 for the SEIRS (beta, Tei, Tir, Trs); 2*n for the cross immunity matrix
      # use the posterior bounds to run simulaitons
      var.state=c(outer(paste0('S.',numbers[1:Na]), letters[1:Ns], FUN = 'paste0') %>% c,
                  outer(paste0('E.',numbers[1:Na]), letters[1:Ns], FUN = 'paste0') %>% c,
                  outer(paste0('I.',numbers[1:Na]), letters[1:Ns], FUN = 'paste0') %>% c,
                  outer(paste0('newItot.',numbers[1:Na]), letters[1:Ns], FUN = 'paste0') %>% c, 
                  outer(paste0('newIobs.',numbers[1:Na]), letters[1:Ns], FUN = 'paste0') %>% c,
                  outer(paste0('death.',numbers[1:Na]), letters[1:Ns], FUN = 'paste0') %>% c
      ) # these are in %
      
      
      # need to assumble the state variables together, accounting for diff variants
      key.sta.names = outer(c('S.','E.','I.'),numbers[1:Na], FUN = 'paste0') %>% c
      key.parm.names = c(outer(c('Tei.','Tir.','Trs.'),numbers[1:Na], FUN = 'paste0') %>% c,
                         beta_age.names)
      key.var.names = c(key.sta.names, key.parm.names)
    }
    
    source(paste0(dir_code, 'Fn_sim_multiVs_age.R'))
  }
  
} # diff variant combo

