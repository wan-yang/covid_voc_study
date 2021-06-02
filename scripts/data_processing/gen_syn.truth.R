# generate synthetic truth for testing
# 3/15/21 - add seasonality



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
dir_code = '../../scripts/model_inference/'
dir_code2 = '../../scripts/plot_results/'
dir_code3 = '../../scripts/simMultiV/'
dir_res = '../../data/'


source(paste0(dir_code,'SEIRS.R'))
source(paste0(dir_code3,'multistrainSEIRSV.R'))
source(paste0(dir_code,'EAKF.R'))
source(paste0(dir_code2,'getPlot.R'))
source(paste0(dir_code,'get_relR0.R'))


if(! file.exists(dir_res))  dir.create(dir_res, recursive = T)

seed1stWave = 2
seed2ndWave = 50
alpha.t = .2
num_ens = 100

seasonality = T

# relR0 = relR0$rel.R0

stoch.t = T 
# no vacc
da.vacc = data.table(date = '2020/03/01', n.v1=0, n.v2=0)

N = 1e6
tm.ini = 1; tmstep = 7
Iexp.t = 1 # non-perfect mixing - maybe not a good idea for the filter

width = c(1,1); # c(.99, 1.01)
# percS0 = 1 - c(.3, .3, .6, .6)
pImm = .5 # prior immunity
percS0 = 1 - rep(pImm, 3)
tx.incr = c(0, .5, .5, .25, .5) # c(0, .5, 0, .5)
imm.loss = c(.8, 0, .8, .4, 0) # c(.8, 0, .8, .5)
sce.lab = c('large prior immunity, no tx incr, 80% imm loss',
            'large prior immunity, 50% tx incr, 0% imm loss', 
            'large prior immunity, 50% tx incr, 80% imm loss',
            'large prior immunity, 25% tx incr, 40% imm loss',
            'minor prior immunity, 50% tx incr, 0% imm loss'
            )
n.sce = length(tx.incr)

percSeed = 10 # 1%
# load estimated tx and imm change
variants = c('wt',paste0('new',1:length(imm.loss)))

seed0 = c(.005, .005, .005, .0003) # diff seeding for the 4 scenarios
seed0 = c(rep(.005,3), .0003) # diff seeding for the 4 scenarios
seed0 = c(rep(.0003,3), .0003) # diff seeding for the 4 scenarios
seeds = rep(.1, length(variants)); names(seeds) = variants
VE1 = .8; VE2 = .95

VE1redn = c(1, imm.loss); names(VE1redn) = variants
VE2redn = c(1, imm.loss); names(VE2redn) = variants

date.t =  '2020/3/1' %>% as.Date
# date.t =  '2020/2/16' %>% as.Date
num_times= 52  # nrow(obs_i);  
week.starts = seq(as.Date(date.t), length.out = num_times, by = 'week')

vax.start = as.Date('2020/12/14')


end1stWave = 25


for(isce in 1:length(imm.loss)){
  tmp = data.table(dRtx.lwr = tx.incr[isce] * width[1], # dRtx
                   dRtx.upr = tx.incr[isce] * width[2], # dRtx
                   dImm.lwr = imm.loss[isce] * width[1],
                   dImm.upr = imm.loss[isce] * width[2],
                   cimm1_2.lwr = 1,
                   cimm1_2.upr = 1,
                   cimm2_1.lwr = 1 - imm.loss[isce] * width[2],
                   cimm2_1.upr = 1 - imm.loss[isce] * width[1],
                   percSeed.lwr = percSeed * width[1],
                   percSeed.upr = percSeed * width[2]
  ) 
  # rownames(adjVparm.new) = c('dRtx','dImm', 'cimm2_1', 'cimm1_2', 'percSeed')
  eval(parse(text = paste('adjVparm.new',isce,'=tmp',sep='')))
}

isce = 1; 
pdf(paste0(dir_res, 'Fig_by.v_seed1st',seed1stWave,'_seed2nd',seed2ndWave, '_alpha',alpha.t,'.pdf'), width = 8, height = 5)
par(mfrow = c(1,1), mar = c(2.5, 2.5, .5, .5), mgp = c(1.5, .5, 0), tck = -.02)
for(isce in 1:5){
  
  if(isce %in% c(1:4)){
    # scenario 1-3:
    relR0.t = 'sa'
    relR0 = read.csv(paste0(dir_data,'relR0_',relR0.t,'.csv')) %>% data.table()
    
    tda.mob = read.csv(paste0(dir_data,'da_case.death.mob_uk.sa.br.csv'), stringsAsFactors = F)  %>% data.table()

    # relMob = 1 + c(relMob$mob.bus.uk[1:17], rowMeans(relMob[18:nrow(relMob),c('mob.full.sa','mob.full.br')])) / 100
    # relMob = 1 + c(relMob$mob.bus.uk[1:30], rowMeans(relMob[31:nrow(relMob),c('mob.full.sa','mob.full.br')])) / 100
    # relMob = 1 + c(rowMeans(relMob[1:nrow(relMob),c('mob.full.uk','mob.full.br')])) / 100
    relMob = 1 + c(tda.mob$mob.bus.sa) / 100
    relMob[is.na(relMob)] = 1
    relMob.full = 1 + c(tda.mob$mob.full.sa) / 100
    # relMob[10: 35] = relMob[10: 35] * .9
    relMob[2: 8] = relMob[2: 8] * .9
    # relMob[11: 18] = relMob[11: 18] * 1.4
    # relMob[11: 18] = relMob.full[11:18]
    p1 = 11:15; p2 = tail(p1,1)+1:4; 
    p3 = tail(p2,1)+1:4; p4 = tail(p3,1)+1:4; p5 = tail(p4,1)+1:4;
    relMob[p1] = (relMob.full[p1])  %>% pmin(.6) # * .9
    relMob[p2] = (relMob[p2] * .7) %>% pmin(.37)
    relMob[p3] = relMob[p3] * .4 # .5
    relMob[p4] = relMob[p4] * .4 # .5
    relMob[p5] = relMob[p5] * .9
    
    if(F){
      relMob[11: 17] = relMob.full[11:17]
      relMob[18: 20] = relMob[18: 20] * .6
      relMob[21: 25] = relMob[21: 25] * .65
      relMob[26: 30] = relMob[26: 30] * .8
      relMob[31: 35] = relMob[31: 35] * .9
    }

    if(isce %in% c(1:4))
      relMob[tail(p5,1): length(relMob)] = (relMob[tail(p5,1): length(relMob)] * 1.1)
    relMob = relMob %>% pmin(1)
    
    relMob %>% plot
  } else if (isce ==5){
    relR0.t = 'uk'
    relR0 = read.csv(paste0(dir_data,'relR0_',relR0.t,'.csv')) %>% data.table()

    # for scenario 4
    relMob = read.csv(paste0(dir_data,'da_case.death.mob_uk.sa.br.csv'), stringsAsFactors = F)  %>% data.table()
    relMob = 1 + c(relMob$mob.bus.uk) / 100
    relMob[is.na(relMob)] = 1
    relMob[5: 25] = relMob[5: 25] * .9
    relMob %>% plot
  }
  
  variants.t = c('wt',paste0('new',isce))
  
  # parms
  parm.bounds = rbind(c(.65, .65),
                      c(3.5, 3.5),
                      c(3.5, 3.5),
                      c(2.5, 2.5) * 365,
                      c(6,6), # mean Td: reporting delay
                      c(2,2), # Td, sd: reporting delay sd
                      alpha.t * width,
                      .7 * width /100
  )
  parm.names = c('beta','Tei','Tir','Trs', 'Td.mean', 'Td.sd','alpha', 'ifr') # , 'p.mob'
  rownames(parm.bounds) = parm.names
  parm.bounds
  
  
  
  source(paste0(dir_code,'set_tm2event.R'))
  
  parm0 = t(lhs(num_ens,parm.bounds))
  rownames(parm0) = parm.names
  dist_tm.to.detect = NULL
  for(ii in 1:num_ens){
    tmp = generation.time(dist_tm.to.detect.name,c(parm0['Td.mean',ii],parm0['Td.sd',ii]),truncate = tm.to.detect.max)
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
  
  
  # wave 1
  {
    # scenarior 1-3
    state0.bounds = rbind(c(1, 1) * N, 
                          # seed0[isce]/100 * width  * N, # E0
                          # seed0[isce]/100 * width  * N
                          seed1stWave, 
                          seed1stWave
                          ) # I0
    rownames(state0.bounds) = c('S','E','I')
    
    num_obs= length(variants.t) # ncol(obs_i);  # number of strains included
    num_parm=(4)*num_obs+(num_obs^2-num_obs); # number parms: 4 for the SEIRS (beta, Tei, Tir, Trs); (num_obs^2-num_obs) for the cross immunity matrix
    num_var=(6+4)*num_obs+(num_obs^2-num_obs); # 5 states: S, E, I, (newI, cases, death);  4 for the SEIRS (beta, Tei, Tir, Trs); 2*n for the cross immunity matrix
    # use the posterior bounds to run simulaitons
    var.state=c(paste0('S',1:num_obs),paste0('E',1:num_obs),paste0('I',1:num_obs),
                paste0('newItot',1:num_obs), paste0('newIobs',1:num_obs),
                paste0('death',1:num_obs)
    ) # these are in %
    
    # 1 = wildtype; 2 = uk; 3 = sa
    
    # need to assumble the state variables together, accounting for diff variants
    key.sta.names = c('S','E','I'); key.parm.names = c('beta', 'Tei','Tir','Trs')
    key.var.names = c(key.sta.names, key.parm.names)
    
    # sample, then adjust
    # adjVparm.t = adjVparm.sa
    VE1redn.t = VE1redn[variants.t] # c(1, .6)
    VE2redn.t = VE2redn[variants.t] # c(1, .8)
    
    base.bounds = rbind(state0.bounds[key.sta.names,], parm.bounds[key.parm.names,])
    So.base=t(lhs(num_ens,base.bounds))
    rownames(So.base)=key.var.names
    
    
    
    So.newi = matrix(0, 3, num_ens)
    rownames(So.newi)= c('newItot', 'newIobs','death')
    So.base = rbind(So.base, So.newi)
    So = NULL
    # add diff variants
    for(iv in 2:num_obs){
      adjVparm.t = get(paste0('adjVparm.', variants.t[iv]))
      So.newv = So.base
      # exclude some gain immunity from vaccination
      cumVx = 0 # sum(da.vacc[date < week.starts[1]]$n.v1) * VE1redn[iv]
      So.newv['S',] = So.base['S',] + (N - So.base['S',] - cumVx) * runif(num_ens, adjVparm.t$dImm.lwr, adjVparm.t$dImm.upr)
      So.newv['E',] = 0 # So.base['E',] / 100 * runif(num_ens, adjVparm.t$percSeed.lwr, adjVparm.t$percSeed.upr)
      So.newv['I',] = 0 # So.base['I',] / 100 * runif(num_ens, adjVparm.t$percSeed.lwr, adjVparm.t$percSeed.upr)
      So.newv['beta',] = So.base['beta',] * (1 + runif(num_ens, adjVparm.t$dRtx.lwr, adjVparm.t$dRtx.upr))
      
      rownames(So.newv) = paste0(rownames(So.newv),iv)
      
      So = rbind(So, So.newv)
    }
    rownames(So.base) = paste0(rownames(So.base),1)
    So = rbind(So.base, So)
    # put together the cross immunity matrix
    icross.lwr = icross.upr = matrix(0, num_obs, num_obs)
    rownames(icross.lwr) = rownames(icross.upr) = variants.t
    colnames(icross.lwr) = colnames(icross.upr) = variants.t
    for(iv in 2:num_obs){
      adjVparm.t = get(paste0('adjVparm.', variants.t[iv]))
      
      icross.lwr['wt', variants.t[iv]] = adjVparm.t$cimm1_2.lwr
      icross.lwr[variants.t[iv], 'wt'] = adjVparm.t$cimm2_1.lwr
      
      icross.upr['wt', variants.t[iv]] = adjVparm.t$cimm1_2.upr
      icross.upr[variants.t[iv], 'wt'] = adjVparm.t$cimm2_1.upr
    }
    
    
    cimm_bounds = matrix(0, num_obs^2 - num_obs, 2) # lwr and upr bounds
    cimm_names = NULL; cnt = 0
    for(i in 1:num_obs){
      for(j in 1:num_obs){
        if(i == j)
          next;
        
        cnt = cnt+1
        cimm.t = paste0('cimm',i,'_',j)
        cimm_names = c(cimm_names, cimm.t)
        cimm_bounds[cnt, ] = c(icross.lwr[variants.t[i], variants.t[j]],
                               icross.upr[variants.t[i], variants.t[j]]
        )
      }
    }
    rownames(cimm_bounds) = cimm_names
    
    So.icross = t(lhs(num_ens,rect = cimm_bounds))
    rownames(So.icross) = cimm_names
    
    So = rbind(So, So.icross)
    
    var.parm=c(paste0('beta',1:num_obs),paste0('Tei',1:num_obs),paste0('Tir',1:num_obs),paste0('Trs',1:num_obs),
               cimm_names
    )
    
    
    x=array(0,c(num_obs*6,num_ens,num_times));  # to store the state variables S, E, I, newItot, newIobs
    dimnames(x)[1]=list(c(paste0('S',1:num_obs),paste0('E',1:num_obs),paste0('I',1:num_obs),
                          paste0('newItot',1:num_obs), paste0('newIobs',1:num_obs), paste0('death',1:num_obs)))
    
    newI.previous = NULL
    
    newI.daily = array(0, c(num_obs, num_ens, num_times * tmstep)) # to save daily newI
    
    for(tt in 1:end1stWave){
      
      seed.t = matrix(c(0, 0), num_obs, num_ens)
      if(tt %in% 1:5){
        seed.t = matrix(c(1/7, 0), num_obs, num_ens)
      } else if (6:10){
        seed.t = matrix(c(1/3, 0), num_obs, num_ens)
      }
        
      
      tcurrent=tm.ini+tmstep*(tt-1);
      
      vdate.t = week.starts[tt]
      # Get the parm for that week use the mean
      parm.tt=So[var.parm,]; row.names(parm.tt)=var.parm
      
      beta.tt=parm.tt[paste0('beta',1:num_obs),] * relMob[tt];
      # adjust for mobility
      cur.wk = MMWRweek(date = week.starts[tt])['MMWRweek'] %>% as.numeric()
      
      if(seasonality){
        beta.tt = beta.tt * unlist(relR0[week == cur.wk,2]) # * rel.mob[tt]
      }
      
      cross=parm.tt[grep('cimm',var.parm),]
      
      if(tt==1){
        S0.tt=So[paste0('S',1:num_obs),];
        E0.tt=So[paste0('E',1:num_obs),];
        I0.tt=So[paste0('I',1:num_obs),];
      } else {
        
        S0.tt=x[paste0('S',1:num_obs),,(tt-1)]
        E0.tt=x[paste0('E',1:num_obs),,(tt-1)]
        I0.tt=x[paste0('I',1:num_obs),,(tt-1)]
        
        # daily
        # S0.tt=x[paste0('S',1:num_obs),,(tt-1)*tmstep]
        # E0.tt=x[paste0('E',1:num_obs),,(tt-1)*tmstep]
        # I0.tt=x[paste0('I',1:num_obs),,(tt-1)*tmstep]
      }
      
      
      percSmax.t = matrix(1, num_obs, num_ens)
      
      tm_strt = tcurrent+1; tm_end=tcurrent+tmstep
      daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
      
      if(nrow(daVacc.t)<1){  # no data yet
        V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
      } else { # yes data
        
        daVacc.t$date = daVacc.t$date %>% as.Date
        
        # make sure it includes a full week
        dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
        daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
        daVacc.t[is.na(daVacc.t)] = 0
        V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
        V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
        
        # print('start vacc!')
        
      }
      
      if(tt ==1){
        newI.previous = NULL
      } else {
        newI.previous = newI.daily[,,1:((tt-1)*tmstep)] # these are at weekly level!
        
      }
      
      
      Sr_tmp=multistrainSEIRSV(tm_strt=tm_strt, tm_end=tm_end, tm_step=1, 
                               S0=S0.tt,E0=E0.tt,I0=I0.tt, 
                               N=N, 
                               Tei = parm.tt[paste0('Tei',1:num_obs),],
                               Tir = parm.tt[paste0('Tir',1:num_obs),], 
                               Trs=parm.tt[paste0('Trs',1:num_obs),],
                               beta=beta.tt, cross=cross,
                               Iexp = Iexp.t,
                               birthrate = 0, 
                               severity = severity,
                               newI.previous = newI.previous,
                               dist_tm.to.detect = dist_tm.to.detect,
                               dist_tm.to.death = dist_tm.to.death,
                               percSmax.t = percSmax.t,
                               V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                               # these are total number of vaccinees with unknown immunity
                               # but pre-ajust for time lag from vaccination to immune protection
                               VE1 = VE1, VE2 = VE2, # Vaccine efficacy, need further adjustment by prior immunity
                               VE1redn = VE1redn.t, VE2redn = VE2redn.t,
                               seed = seed.t,
                               stoch = stoch.t)
      
      
      x[paste0('S',1:num_obs),,tt]=Sr_tmp$S[,,1+tmstep] # tail(Sr_tmp$S,1);
      x[paste0('E',1:num_obs),,tt]=Sr_tmp$E[,,1+tmstep];
      x[paste0('I',1:num_obs),,tt]=Sr_tmp$I[,,1+tmstep];
      x[paste0('newItot',1:num_obs),,tt]=Sr_tmp$cumItot[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      x[paste0('newIobs',1:num_obs),,tt]=Sr_tmp$cumIobs[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      x[paste0('death',1:num_obs),,tt]=Sr_tmp$death[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      
      
      newI.daily[,,1:tmstep+(tt-1)*tmstep]= Sr_tmp$daily.newItot # we want the daily total new cases, without delay, without under-report
      
    } 
    
    
  } # wave 1
  
  
  # wave 2
  {
    # update So
    So['S1',] = x['S1',, end1stWave]
    So['E1',] = x['E1',, end1stWave]
    So['I1',] = x['I1',, end1stWave]
    So['S2',] = So['S1',] + (N - So['S1',] - cumVx) * runif(num_ens, adjVparm.t$dImm.lwr, adjVparm.t$dImm.upr)
    
    if(F){
      So['E2',] = x['E1',, end1stWave] / 100 * runif(num_ens, adjVparm.t$percSeed.lwr, adjVparm.t$percSeed.upr)
      So['I2',] = x['I1',, end1stWave] / 100 * runif(num_ens, adjVparm.t$percSeed.lwr, adjVparm.t$percSeed.upr)
      
      # adjust wt E and I
      So['E1', ] = So['E1', ] * (1 - percSeed[1]/100)
      So['I1', ] = So['I1', ] * (1 - percSeed[1]/100)
    }
    
    So['E2',] = seed2ndWave
    So['I2',] = seed2ndWave

    
    seed.t = matrix(c(0, 0), num_obs, num_ens)
    
    for(tt in (1 + end1stWave):num_times){
      tcurrent=tm.ini+tmstep*(tt-1);
      
      vdate.t = week.starts[tt]
      # Get the parm for that week use the mean
      parm.tt=So[var.parm,]; row.names(parm.tt)=var.parm
      
      beta.tt=parm.tt[paste0('beta',1:num_obs),] * relMob[tt];
      # adjust for mobility
      cur.wk = MMWRweek(date = week.starts[tt])['MMWRweek'] %>% as.numeric()
      
      if(seasonality){
        beta.tt = beta.tt * unlist(relR0[week == cur.wk,2]) # * rel.mob[tt]
      }
      
      cross=parm.tt[grep('cimm',var.parm),]
      
      if(tt==(1 + end1stWave)){
        S0.tt=So[paste0('S',1:num_obs),];
        E0.tt=So[paste0('E',1:num_obs),];
        I0.tt=So[paste0('I',1:num_obs),];
      } else {
        
        S0.tt=x[paste0('S',1:num_obs),,(tt-1)]
        E0.tt=x[paste0('E',1:num_obs),,(tt-1)]
        I0.tt=x[paste0('I',1:num_obs),,(tt-1)]
        
      }
      
      
      percSmax.t = matrix(1, num_obs, num_ens)
      
      tm_strt = tcurrent+1; tm_end=tcurrent+tmstep
      daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
      
      if(nrow(daVacc.t)<1){  # no data yet
        V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, num_ens)
      } else { # yes data
        
        daVacc.t$date = daVacc.t$date %>% as.Date
        
        # make sure it includes a full week
        dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
        daVacc.t = merge(daVacc.t, dates.t, all = T, by = 'date')
        daVacc.t[is.na(daVacc.t)] = 0
        V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
        V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
        
        # print('start vacc!')
        
      }
      
      if(tt ==1){
        newI.previous = NULL
      } else {
        newI.previous = newI.daily[,,1:((tt-1)*tmstep)] # these are at weekly level!
        
      }
      
      
      Sr_tmp=multistrainSEIRSV(tm_strt=tm_strt, tm_end=tm_end, tm_step=1, 
                               S0=S0.tt,E0=E0.tt,I0=I0.tt, 
                               N=N, 
                               Tei = parm.tt[paste0('Tei',1:num_obs),],
                               Tir = parm.tt[paste0('Tir',1:num_obs),], 
                               Trs=parm.tt[paste0('Trs',1:num_obs),],
                               beta=beta.tt, cross=cross,
                               Iexp = Iexp.t,
                               birthrate = 0, 
                               severity = severity,
                               newI.previous = newI.previous,
                               dist_tm.to.detect = dist_tm.to.detect,
                               dist_tm.to.death = dist_tm.to.death,
                               percSmax.t = percSmax.t,
                               V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                               # these are total number of vaccinees with unknown immunity
                               # but pre-ajust for time lag from vaccination to immune protection
                               VE1 = VE1, VE2 = VE2, # Vaccine efficacy, need further adjustment by prior immunity
                               VE1redn = VE1redn.t, VE2redn = VE2redn.t,
                               seed = seed.t,
                               stoch = stoch.t)
      
      
      x[paste0('S',1:num_obs),,tt]=Sr_tmp$S[,,1+tmstep] # tail(Sr_tmp$S,1);
      x[paste0('E',1:num_obs),,tt]=Sr_tmp$E[,,1+tmstep];
      x[paste0('I',1:num_obs),,tt]=Sr_tmp$I[,,1+tmstep];
      x[paste0('newItot',1:num_obs),,tt]=Sr_tmp$cumItot[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      x[paste0('newIobs',1:num_obs),,tt]=Sr_tmp$cumIobs[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      x[paste0('death',1:num_obs),,tt]=Sr_tmp$death[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      
      
      newI.daily[,,1:tmstep+(tt-1)*tmstep]= Sr_tmp$daily.newItot # we want the daily total new cases, without delay, without under-report
      
    }
    
  } # wave 2
  
  
  cols = c('darkgreen','orange','red','brown')
  newItot = x[paste0('newItot',1:num_obs),,] %>% apply(c(1,3), quantile, probs=c(.5, .25, .75)) / N * 100
  newItot[1,,] %>% t %>% matplot(type='l', lty=1, col = cols, xaxt = 'n')
  axis(1, at = 1: ncol(newItot[1,,]), labels = week.starts)
  mtext(isce, line = -1.5, outer = F)
  
  # newIobs = x[paste0('newIobs',1:num_obs),,] %>% apply(c(1,3), quantile, probs=c(.5, .25, .75)) / N * 100
  # newIobs[1,,] %>% t %>% matplot(type='l', lty=1, col = cols, xaxt = 'n')
  # axis(1, at = 1: ncol(newIobs[1,,]), labels = week.starts)
  
  # newD = x[paste0('death',1:num_obs),,] %>% apply(c(1,3), quantile, probs=c(.5, .25, .75)) / N * 100
  # newD[1,,] %>% t %>% matplot(type='l', lty=1, col = cols, xaxt = 'n')
  # axis(1, at = 1: ncol(newD[1,,]), labels = week.starts)
  
  perS = x[paste0('S',1:num_obs),,] %>% apply(c(1,3), quantile, probs=c(.5, .25, .75)) / N * 100
  perS[1,,] %>% t %>% matplot(type='l', lty=1, col = cols, xaxt = 'n')
  axis(1, at = 1: ncol(perS[1,,]), labels = week.starts)
  mtext(isce, line = -1.5, outer = F)
  
  tmp = x[paste0('E',1:num_obs),,] %>% apply(c(1,3), quantile, probs=c(.5))

  # combine both strains and save the time series
  suscept = x[paste0('S',1:num_obs),,] %>% apply(c(1,3),median)  %>% t
  infect = x[paste0('newItot',1:num_obs),,] %>% apply(c(1,3),median) %>% t
  infect = infect / rowSums(infect)
  wt.suscept = (suscept * infect) %>% rowSums() 
  # compute the wt Rtx as well
  Rtx1 = (So['beta1',] * So['Tir1',]) %>% median()
  Rtx2 = (So.newv['beta2',] * So.newv['Tir2',]) %>% median()
  Rtx = c(Rtx1, Rtx2)
  wt.Rtx = (Rtx * t(infect)) %>% colSums()
  tda = data.table(date = week.starts, 
                   case.truth = x[paste0('newIobs',1:num_obs),,] %>% apply(c(2,3),sum) %>% apply(2,median),
                   case = rpois(n = num_times, lambda = x[paste0('newIobs',1:num_obs),,] %>% apply(c(2,3),sum) %>% apply(2,median)),
                   death = x[paste0('death',1:num_obs),,] %>% apply(c(2,3),sum) %>% apply(2,median),
                   mob = relMob[1:num_times],
                   relR0 = relR0.t,
                   wt.suscept = wt.suscept, 
                   wt.Rtx = wt.Rtx,
                   newitot = x[paste0('newItot',1:num_obs),,] %>% apply(c(2,3),sum) %>% apply(2,median),
                   Etot = x[paste0('E',1:num_obs),,] %>% apply(c(2,3),sum) %>% apply(2,median),
                   Itot = x[paste0('I',1:num_obs),,] %>% apply(c(2,3),sum) %>% apply(2,median),
                   x[paste0('S',1:num_obs),,]  %>% apply(c(1,3),median)  %>% t,
                   x[paste0('newItot',1:num_obs),,] %>% apply(c(1,3),median)  %>% t,
                   x[paste0('E',1:num_obs),,] %>% apply(c(1,3),median)  %>% t,
                   x[paste0('I',1:num_obs),,] %>% apply(c(1,3),median)  %>% t
                   )
  
  write.csv(tda, paste0(dir_data, 'syn.truth_sce',isce,'_seed1st',seed1stWave,'_seed2nd',seed2ndWave, '_alpha',alpha.t,'.csv'), row.names = F)
} # end this scenario
dev.off()

variables = c('case.truth', 'case','death','mob','relR0','wt.suscept','wt.Rtx','newitot','Etot','Itot','S1','S2','newItot1', 'newItot2','E1','E2','I1','I2')
da = read.csv(paste0(dir_data, 'syn.truth_sce1_seed1st',seed1stWave,'_seed2nd',seed2ndWave, '_alpha',alpha.t,'.csv')) %>% data.table() %>% 
  setnames(variables, paste0(variables,'.sce1'))
for(isce in 2:length(imm.loss)){
  tda = read.csv(paste0(dir_data, 'syn.truth_sce',isce, '_seed1st',seed1stWave,'_seed2nd',seed2ndWave, '_alpha',alpha.t,'.csv'))%>% data.table() %>% 
    setnames(variables, paste0(variables,'.sce',isce))
  da = merge(da, tda, by = 'date')
}

write.csv(da, paste0(dir_data, 'da_syn.truth_seed1st',seed1stWave,'_seed2nd',seed2ndWave, '_alpha',alpha.t,'.csv'), row.names = F)

