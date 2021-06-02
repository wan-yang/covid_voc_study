# simulate multiple variants

# diff scenarios

if(!exists('withObs.t')) withObs.t = F
if(!exists('obs')) obs = NULL
if(!exists('obs.prev')) obs.prev = NULL
if(!exists('seed.tag')) seed.tag = ''

# sample, then adjust
# adjVparm.t = adjVparm.sa
VE1redn.t = VE1redn[variants.t] # c(1, .6)
VE2redn.t = VE2redn[variants.t] # c(1, .8)

base.bounds = rbind(state0.bounds[key.sta.names,], parm.bounds[key.parm.names,])
So.base=t(lhs(num_ens,base.bounds))
rownames(So.base)=key.var.names
So.newi = matrix(0, 2 * Na, num_ens)
rownames(So.newi)= c(paste0('newItot.',1:Na), paste0('newIobs.',1:Na))
So.base = rbind(So.base, So.newi)
So = NULL
mnImmLoss = rep(0, Na) # for wt
# add diff variants
for(iv in 2:Ns){
  adjVparm.t = get(paste0('adjVparm.', variants.t[iv]))
  # get the cumulative infection rate at the start
  cumIstart.t = get(paste0('cumIstart.', variants.t[iv]))
  adjVparm.t
  sum(cumIstart.t) / N
  
  percSeed.mns.t = percSeed.mns[variants.t[iv]] / 100 # the prevalence for this variant at the start
  
  So.newv = So.base
  
  for(ia in 1: Na){
    # exclude some gain immunity from vaccination
    if(F){
      # if all vacc data were included
      cumVx1 = sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v1) * VE1 
      cumVx2 = sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v2) *  (1 - VE1) * VE2 
      
      cumVx1redn = sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v1) * VE1 * (1 - VE1redn[iv])
      cumVx2redn = sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v2) *  (1 - VE1) * VE2 * (1 - VE2redn[iv])
      cumVx.redn = cumVx1redn + cumVx2redn
      cumVx.redn / Nage[ia] # % of this age group vaccinated but lost immunity
      
    }
    
    # b/c data prior to the simulation were removed, here we add the cumulative prior to the simualtion
    cumVx1 = (sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v1) +
                da.vacc0[age.grp == numbers[ia]]$n.v1) * VE1 
    cumVx2 = (sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v2) +
                da.vacc0[age.grp == numbers[ia]]$n.v2
              ) *  (1 - VE1) * VE2 
    
    
    
    cumVx = cumVx1 + cumVx2
    cumVx / Nage[ia] # % of this age group been fully vaccinated
    
    # those vaccinated lost immunity due to reduced VE
    cumVx1redn = (sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v1) +
                    da.vacc0[age.grp == numbers[ia]]$n.v1) * VE1 * (1 - VE1redn[iv])
    cumVx2redn = (sum(da.vacc[date < week.starts[1] & age.grp == numbers[ia]]$n.v2) +
                    da.vacc0[age.grp == numbers[ia]]$n.v2
    ) *  (1 - VE1) * VE2 * (1 - VE2redn[iv])
    cumVx.redn = cumVx1redn + cumVx2redn
    cumVx.redn / Nage[ia] # % of this age group vaccinated but lost immunity
    
    
    So.newv[paste0('S.',numbers[ia]),] = So.base[paste0('S.',numbers[ia]),] + cumVx.redn + 
       # (Nage[ia] - So.base[paste0('S.',numbers[ia]),] - cumVx) * # this one is close but include some infections due to the newV
        (Nage[ia] - So.base[paste0('S.',numbers[ia]),] - cumVx -
           pmax(0, (cumI0wt[ia] - cumIstart.t[ia])) * percSeed.mns.t /2 * (1 - cumVx/Nage[ia]) # new infections from introduction to now that are due to the newV
           ) *
      runif(num_ens, adjVparm.t$dImm.lwr, adjVparm.t$dImm.upr) # cumulative infection
      
    
    
    # use cumulative infection
    # So.newv[paste0('S.',numbers[ia]),] = So.base[paste0('S.',numbers[ia]),] + cumI0wt[ia] * runif(num_ens, adjVparm.t$dImm.lwr, adjVparm.t$dImm.upr)
    # make sure it's less than N
    So.newv[paste0('S.',numbers[ia]),] = pmin(So.newv[paste0('S.',numbers[ia]),], Nage[ia])
    So.newv[paste0('E.',numbers[ia]),] = So.base[paste0('E.',numbers[ia]),] / 100 * runif(num_ens, adjVparm.t$percSeed.lwr, adjVparm.t$percSeed.upr)
    So.newv[paste0('I.',numbers[ia]),] = So.base[paste0('I.',numbers[ia]),] / 100 * runif(num_ens, adjVparm.t$percSeed.lwr, adjVparm.t$percSeed.upr)
    
    
    mnImmLoss = c(mnImmLoss, mean(adjVparm.t$dImm.lwr, adjVparm.t$dImm.upr))
  } # age group
  
  So.newv[beta_age.names,] = So.base[beta_age.names,] * matrix((1 + runif(num_ens, adjVparm.t$dRtx.lwr, adjVparm.t$dRtx.upr)), length(beta_age.names), num_ens, byrow = T)
  
  
  rownames(So.newv) = paste0(rownames(So.newv),letters[iv])
  
  So = rbind(So, So.newv)
}


rownames(So.base) = paste0(rownames(So.base),letters[1])
# adjust E and I for wt 
So.base[paste0('E.', 1:Na,letters[1]), ] = So.base[paste0('E.', 1:Na,letters[1]), ] * (1 - sum(percSeed.mns)/100)
So.base[paste0('I.', 1:Na,letters[1]), ] = So.base[paste0('I.', 1:Na,letters[1]), ] * (1 - sum(percSeed.mns)/100)
So = rbind(So.base, So)
# put together the cross immunity matrix
# icross.lwr = icross.upr = matrix(0, Ns, Ns)
icross.lwr = icross.upr = matrix(1, Ns, Ns) # make the default 1
diag(icross.lwr) = 0; diag(icross.upr) = 0; 
rownames(icross.lwr) = rownames(icross.upr) = variants.t
colnames(icross.lwr) = colnames(icross.upr) = variants.t
for(iv in 2:Ns){
  adjVparm.t = get(paste0('adjVparm.', variants.t[iv]))
  
  icross.lwr['wt', variants.t[iv]] = adjVparm.t$cimma_b.lwr
  icross.lwr[variants.t[iv], 'wt'] = adjVparm.t$cimmb_a.lwr
  
  icross.upr['wt', variants.t[iv]] = adjVparm.t$cimma_b.upr
  icross.upr[variants.t[iv], 'wt'] = adjVparm.t$cimmb_a.upr
}
# use the same estimate for wt for uk vs. sa and br
if('b1526' %in% variants.t & 'p1' %in% variants.t){
  icross.lwr['b1526','p1'] = icross.lwr['wt','p1']
  icross.lwr['p1','b1526'] = icross.lwr['p1','wt']
  
  icross.upr['b1526','p1'] = icross.upr['wt','p1']
  icross.upr['p1','b1526'] = icross.upr['p1','wt']
}

if('b1526' %in% variants.t & 'b1351' %in% variants.t){
  icross.lwr['b1526','b1351'] = icross.lwr['wt','b1351']
  icross.lwr['b1351','b1526'] = icross.lwr['b1351','wt']
  
  icross.upr['b1526','b1351'] = icross.upr['wt','b1351']
  icross.upr['b1351','b1526'] = icross.upr['b1351','wt']
}
if('b1427' %in% variants.t & 'p1' %in% variants.t){
  icross.lwr['b1427','p1'] = icross.lwr['wt','p1']
  icross.lwr['p1','b1427'] = icross.lwr['p1','wt']
  
  icross.upr['b1427','p1'] = icross.upr['wt','p1']
  icross.upr['p1','b1427'] = icross.upr['p1','wt']
}

if('b1427' %in% variants.t & 'b1351' %in% variants.t){
  icross.lwr['b1427','b1351'] = icross.lwr['wt','b1351']
  icross.lwr['b1351','b1427'] = icross.lwr['b1351','wt']
  
  icross.upr['b1427','b1351'] = icross.upr['wt','b1351']
  icross.upr['b1351','b1427'] = icross.upr['b1351','wt']
}

if('b117' %in% variants.t & 'p1' %in% variants.t){
  icross.lwr['b117','p1'] = icross.lwr['wt','p1']
  icross.lwr['p1','b117'] = icross.lwr['p1','wt']
  
  icross.upr['b117','p1'] = icross.upr['wt','p1']
  icross.upr['p1','b117'] = icross.upr['p1','wt']
}

if('b117' %in% variants.t & 'b1351' %in% variants.t){
  icross.lwr['b117','b1351'] = icross.lwr['wt','b1351']
  icross.lwr['b1351','b117'] = icross.lwr['b1351','wt']
  
  icross.upr['b117','b1351'] = icross.upr['wt','b1351']
  icross.upr['b1351','b117'] = icross.upr['b1351','wt']
}
  
  
if('b1351' %in% variants.t & 'p1' %in% variants.t){
  # study showed very good neutralization of sa serum to br variant; unknown br -> sa
  icross.lwr['p1','b1351'] = .8 # they share quite a lot same mutations
  icross.lwr['b1351','p1'] = .6 # .7 # they share quite a lot same mutations
  
  icross.upr['b1351','p1'] = 1 # they share quite a lot same mutations
  icross.upr['p1','b1351'] = 1 # they share quite a lot same mutations
  
}





cimm_bounds = matrix(0, Ns^2 - Ns, 2) # lwr and upr bounds
cimm_names = NULL; cnt = 0
for(i in 1:Ns){
  for(j in 1:Ns){
    if(i == j)
      next;
    
    cnt = cnt+1
    cimm.t = paste0('cimm',letters[i],'_',letters[j])
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

var.parm=c(outer(key.parm.names, letters[1:Ns], FUN = paste0) %>% c, # paste0('beta',1:Ns),paste0('Tei',1:Ns),paste0('Tir',1:Ns),paste0('Trs',1:Ns),
           cimm_names
)

tm.ini = 1; tmstep = 7
x=array(0,c(Ns * Na *6,num_ens,num_times));  # to store the state variables S, E, I, newItot, newIobs, death
dimnames(x)[1]=list(var.state)

N.nans = NULL
for(iv in 1: Ns) {
  N.nans = c(N.nans, Nage)
}
cumI0 = rep(cumI0wt, Ns) * .75 * (1 - mnImmLoss) # prior infection rate thus far, adjusted for immune loss to new variant; for vx


newI.previous = NULL

# add 2 weeks of prior infection

newI.daily = array(0, c(Ns * Na, num_ens, (num_times + nwk.pre) * tmstep)) # to save daily newI
dimnames(newI.daily)[1]= list(outer(paste0('daily.newI.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
# fill in the first few weeks 
tmp0 = dcast(daily.pre.inf, age.grp ~ date, value.var = 'daily.newI')
tmp0 = tmp0[,-1] %>% as.matrix
tmp = NULL
for(iv in 2:Ns){
  tmp1 = tmp0 * percSeed.lwrs[variants.t[iv]] / 100
  tmp = rbind(tmp, tmp1)
}
tmp0 = tmp0 * (100 - sum(percSeed.lwrs)) /100
tmp = rbind(tmp0, tmp)
for(ii in 1:num_ens){
  newI.daily[,ii,1:(nwk.pre*7)] = tmp
}


seed.t = matrix(outer(Nage / sum(Nage), seeds[variants.t]) %>% c, Na * Ns, num_ens) # age, variant

# relMob = .5

if(grepl('25% Less NPI', sce.t)){
  relMob = seq(1, length.out = num_times, by = .05)
  # make sure it is less than the pre-pandemic level
  relMob = pmin(relMob, pmin(1.25, relMob.noNPI))
} else if(grepl('50% Less NPI', sce.t)){
  relMob = seq(1, length.out = num_times, by = .05)
  # make sure it is less than the pre-pandemic level
  relMob = pmin(relMob, pmin(1.5,relMob.noNPI))
} else if (grepl('Current NPI', sce.t)){
  relMob = rep(1, num_times) 
} else if (grepl('Historic', sce.t)){
  relMob = relMob.hist[date %in% week.starts]$rel.mob
  relMob.ref = relMob.hist[date %in% seq(week.starts[1]-14,week.starts[1],by='week')]$rel.mob %>% mean
  relMob = relMob / relMob.ref  # relative to the first week
  # in case the last week is missing
  if(length(relMob) < num_times){
    relMob = c(relMob, rep(mean(tail(relMob,3)), num_times - length(relMob)))
  }
} else if(grepl('No NPI fast', sce.t)){
  relMob = seq(1, length.out = num_times, by = .1)
  # make sure it is less than the pre-pandemic level
  relMob = pmin(relMob, relMob.noNPI)
} else if(grepl('No NPI slow', sce.t)){
  relMob = seq(1, length.out = num_times, by = .05)
  # make sure it is less than the pre-pandemic level
  relMob = pmin(relMob, relMob.noNPI)
} 



for(tt in 1:num_times){
  tcurrent=tm.ini+tmstep*(tt-1);
  
  vdate.t = week.starts[tt]
  
  # time varying seeding? # for the UK variants etc
  if(exists('tmVarSeed')){
    if(tmVarSeed){
      if(vdate.t >= as.Date('2020/11/1') & vdate.t <= as.Date('2020/11/15')){
        seed.t = matrix(outer(Nage / sum(Nage), seeds.afterDec20[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else if(vdate.t >= as.Date('2020/11/16') & vdate.t <= as.Date('2020/12/31')){
        seed.t = matrix(outer(Nage / sum(Nage), seeds.afterDec20late[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else if(vdate.t >= as.Date('2021/1/1') & vdate.t <= as.Date('2021/1/31')){
        seed.t = matrix(outer(Nage / sum(Nage), seeds.afterJan21[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else if(vdate.t >= as.Date('2021/2/1') & vdate.t <= as.Date('2021/2/28')){
        seed.t = matrix(outer(Nage / sum(Nage), seeds.afterFeb21[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else if(vdate.t >= as.Date('2021/3/1') & vdate.t <= as.Date('2021/3/31')) {
        seed.t = matrix(outer(Nage / sum(Nage), seeds.afterMar21[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else if(vdate.t >= as.Date('2021/4/1') & vdate.t <= as.Date('2021/5/31')){
        seed.t = matrix(outer(Nage / sum(Nage), seeds.afterApr21[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else if(vdate.t >= as.Date('2021/6/1')){
        seed.t = matrix(outer(Nage / sum(Nage), seeds.late[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
      } else {
        seed.t = matrix(outer(Nage / sum(Nage), seeds[variants.t]) %>% c, Na * Ns, num_ens) # age, variant
        
      }
      
    }
  }
  
  
  # Get the parm for that week use the mean
  parm.tt=So[var.parm,]; row.names(parm.tt)=var.parm
  
  # beta.tt=parm.tt[paste0('beta',1:Ns),] * relMob[tt];
  beta.tt=parm.tt[outer(beta_age.names,letters[1:Ns], FUN=paste0) %>% c,] * relMob[tt]
  # adjust for mobility
  cur.wk = MMWRweek(date = week.starts[tt])['MMWRweek'] %>% as.numeric()
  
  if(seasonality){
    beta.tt = beta.tt * unlist(relR0[week == cur.wk,2]) # * rel.mob[tt]
  }
  
  cross=parm.tt[grep('cimm',var.parm),]
  
  if(tt==1){
    S0.tt=So[outer(paste0('S.',1:Na), letters[1:Ns], FUN=paste0) %>% c,];
    E0.tt=So[outer(paste0('E.',1:Na), letters[1:Ns], FUN=paste0) %>% c,];
    I0.tt=So[outer(paste0('I.',1:Na), letters[1:Ns], FUN=paste0) %>% c,];
  } else {
    S0.tt=x[outer(paste0('S.',1:Na), letters[1:Ns], FUN=paste0) %>% c,,tt-1]
    E0.tt=x[outer(paste0('E.',1:Na), letters[1:Ns], FUN=paste0) %>% c,,tt-1]
    I0.tt=x[outer(paste0('I.',1:Na), letters[1:Ns], FUN=paste0) %>% c,,tt-1]
  }
  
  {
    newI.previous = newI.daily[,,1:((tt-1 + nwk.pre)*tmstep)]
    dimnames(newI.previous)[1] = list(outer(paste0('ihist.', 1: Na),letters[1:Ns], FUN = paste0) %>% c)
  }
  
  if(!is.null(newI.previous) & vdate.t >= vax.start){
    dim.t = dim(newI.previous)
    # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
    # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
    # 2/5/21 set t1 to 1 yr given the slow rollout
    t1 = 365
    cumI.t = apply(newI.previous[,,pmax(1,dim.t[3]-t1) : (dim.t[3]),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
    # add cumulative infection before the sumulation
    cumI.t = cumI0 + cumI.t # note cumI0 already adjusted for loss of immunity
    
    
    # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
    tm.t = dim.t[3]
    # higher infection rate for the priority group
    # tm.t = pmax(1, tm.t - t1 + 1) # re-align timing with start of vac
    tm.imm = 365 * 2.5 # assume 3 yr immunity
    p.imm.wane.max = .8; k = .015  # 1/tm.imm  
    p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(tm.t + 60 - tm.imm/2))) # not all infected from day 1
    
    # earlier infections are likely to be in the high priority groups 
    p.imm = 1 *  p.imm.wane # assume 100% prior infection provide immunity, but wane over time
    # and multiple by % excluded if there is prior testing before vax: p.prior.test
    percSmax.t = 1 - cumI.t / N.nans * p.imm    # immune escape? no, these are variant specific
    # no lower than 50%, in case of outliers
    percSmax.t = pmax(percSmax.t, matrix(.5, dim.t[1], dim.t[2]))
    rownames(percSmax.t) = outer(paste0("percSmax.", 1:Na), letters[1:Ns], FUN = paste0) %>% c
    # print(c('cohort %S:',round(summary(apply(percSmax.t, 1,mean)),2)), quote = F)
  } else {
    percSmax.t = matrix(0, Na * Ns, num_ens)
    rownames(percSmax.t) = outer(paste0("percSmax.", 1:Na), letters[1:Ns], FUN = paste0) %>% c
  }
  
  
  
  tm_strt = tcurrent+1; tm_end=tcurrent+tmstep
  daVacc.t = da.vacc[as.Date(date) >= as.Date(vdate.t) & as.Date(date) < as.Date(vdate.t)+tm_end-tm_strt+1] # vaccination data
  
  if(nrow(daVacc.t)<1){  # no data yet
    V1.t = V2.t = matrix(0, tm_end - tm_strt + 1, Na)
  } else { # yes data
    
    daVacc.t$date = daVacc.t$date %>% as.Date
    
    # make sure it includes a full week
    # dates.t = data.table(date = seq(as.Date(vdate.t), length.out = tm_end-tm_strt+1, by='day'))
    # need both days and age groups
    dates.t = data.table(date = rep(seq(as.Date(vdate.t), length.out = (tm_end-tm_strt+1), by='day'), e = Na),
                         age.grp = rep(1:Na, (tm_end-tm_strt+1)))
    daVacc.t = merge(daVacc.t, dates.t, all = T, by = c('date', 'age.grp'))
    daVacc.t[is.na(daVacc.t)] = 0
    # daVacc.t = daVacc.t[complete.cases(daVacc.t)]
    # V1.t = as.matrix(daVacc.t$n.v1, tmstep, num_ens)
    # V2.t = as.matrix(daVacc.t$n.v2, tmstep, num_ens)
    
    # BUG
    V1.t = dcast(daVacc.t, date ~ age.grp, value.var = 'n.v1')
    V2.t = dcast(daVacc.t, date ~ age.grp, value.var = 'n.v2')
    V1.t = V1.t[,-1] %>% as.matrix()
    V2.t = V2.t[,-1] %>% as.matrix()
    # print('start vacc!')
    
  }
  
  
  Sr_tmp=multistrainSEIRSVage(tm_strt=tm_strt, tm_end=tm_end, tm_step=1, 
                S0=S0.tt,E0=E0.tt,I0=I0.tt, 
                Nage = Nage, 
                Tei = parm.tt[outer(paste0('Tei.',1:Na), letters[1:Ns], FUN=paste0) %>% c,],
                Tir = parm.tt[outer(paste0('Tir.',1:Na), letters[1:Ns], FUN=paste0) %>% c,], 
                Trs= parm.tt[outer(paste0('Trs.',1:Na), letters[1:Ns], FUN=paste0) %>% c,],
                beta=beta.tt, cross=cross,
                Iexp = .97,
                birthrate = 0, 
                # severity = severity,
                newI.previous = newI.previous,
                # dist_tm.to.detect = dist_tm.to.detect,
                # dist_tm.to.death = dist_tm.to.death,
                # for vaccination
                percSmax.t = percSmax.t,
                V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                # these are total number of vaccinees with unknown immunity
                # but pre-ajust for time lag from vaccination to immune protection
                VE1 = VE1, VE2 = VE2, # Vaccine efficacy, need further adjustment by prior immunity
                VE1redn = VE1redn.t, VE2redn = VE2redn.t,
                seed = seed.t,
                stoch = T)
  
  x[outer(paste0('S.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,tt]=Sr_tmp$S[,,1+tmstep]  # tail(Sr_tmp$S,1);
  x[outer(paste0('E.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,tt]=Sr_tmp$E[,,1+tmstep];
  x[outer(paste0('I.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,tt]=Sr_tmp$I[,,1+tmstep];
  x[outer(paste0('newItot.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,tt]=Sr_tmp$cumItot[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:Ns),ii])
  x[outer(paste0('newIobs.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,tt]=Sr_tmp$cumIobs[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:Ns),ii])
  x[outer(paste0('death.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,tt]=Sr_tmp$death[,,1+tmstep]; #  %*% diag(So[paste0('sf',1:Ns),ii])
  newI.daily[outer(paste0('daily.newI.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,1:tmstep+(tt-1+nwk.pre)*tmstep]= Sr_tmp$daily.newItot # we want the daily total new cases, without delay, without under-report
  
}

# compute other metrics as well
tmp = fn_get.health.metrics(newI.daily)
stats = tmp$res
perc = tmp$perc
sums.stats = tmp$sums.res
sums.perc = tmp$sums.perc
totals = tmp$sums
cumsums = tmp$cumsums
# save results
save(x, stats, perc, sums.stats, sums.perc, totals, cumsums, file = paste0(dir_res, 'out_',seed.tag, '_',v.tag,'_',sce.tag, '.RData'))

cols = c('darkgreen','orange','red','brown')



# by age group

vname.t = 'newItot'
tda.t = x[outer(paste0(vname.t,'.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,]; 
newInf = fn_get.aggregate.by.v(tda.t, per.capita = T, N = N)$res
p1 = getPlotMultiV(newInf, y.lab = '% Infection')
vname.t = 'S'
tda.t = x[outer(paste0(vname.t,'.',1:Na), letters[1:Ns], FUN = paste0) %>% c,,]; 
percS = fn_get.aggregate.by.v(tda.t, per.capita = T, N = N)$res
p2 = getPlotMultiV(percS, y.lab = '% Susceptible')

newInf = stats[measure == 'New infections' & age.grp !='All']
p3 = getPlotMultiVage.overlay(newInf, title.t = sce.t, y.lab = 'Number Infections')
newcase = stats[measure == 'New cases' & age.grp !='All']
p4 = getPlotMultiVage.overlay(newcase, title.t = sce.t, y.lab = 'Number Cases')
newed = stats[measure == 'New ED visits' & age.grp !='All']
p5 = getPlotMultiVage.overlay(newed, title.t = sce.t, y.lab = 'Number ED Visits')
newhosp = stats[measure == 'New hospitalizations' & age.grp !='All']
p6 = getPlotMultiVage.overlay(newhosp, title.t = sce.t, y.lab = 'Number Hospitalizations')
newdeath = stats[measure == 'New deaths' & age.grp !='All']
p7 = getPlotMultiVage.overlay(newdeath, title.t = sce.t, y.lab = 'Number Deaths')

pdf(paste0(dir_res, 'Fig_byage.sim.',loc.t,'_', seed.tag, '_',v.tag,'_', sce.tag, '.pdf'), width = 8, height = 6)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

tda = stats[measure %in% c('New infections', 'New cases', 'New hospitalizations', 'New deaths')  & age.grp =='All']
tda$measure = factor(tda$measure, levels = c('New infections', 'New cases', 'New hospitalizations', 'New deaths'),
                     labels = c('New infections', 'New cases', 'New hospitalizations', 'New deaths'))
p8 = getPlotMultiVoverlay(tda, title = sce.t, y.lab = 'Estimate', withObs = withObs.t, obs.t = obs)

tda = stats[measure %in% c("Hospital bed needs (mean)", "Non-ICU hospital bed needs (mean)", "ICU bed needs (mean)", "Ventilator needs (mean)")  & age.grp =='All']
tda$measure = factor(tda$measure, levels = c("Hospital bed needs (mean)", "Non-ICU hospital bed needs (mean)", "ICU bed needs (mean)", "Ventilator needs (mean)"),
                     labels = c("Hospital bed needs (mean)", "Non-ICU hospital bed needs (mean)", "ICU bed needs (mean)", "Ventilator needs (mean)"))
p9 = getPlotMultiVoverlay(tda, title = sce.t, y.lab = 'Estimate')

tda = perc[measure %in% c('New infections', 'New cases', 'New hospitalizations', 'New deaths')  & age.grp =='All']
tda$measure = factor(tda$measure, levels = c('New infections', 'New cases', 'New hospitalizations', 'New deaths'),
                     labels = c('New infections', 'New cases', 'New hospitalizations', 'New deaths'))
p10 = getPlotMultiVoverlay(tda, title = sce.t, y.lab = 'Estimate %', withObs = withObs.t, obs.t = obs.prev)

tda = perc[measure %in% c("Hospital bed needs (mean)", "Non-ICU hospital bed needs (mean)", "ICU bed needs (mean)", "Ventilator needs (mean)")  & age.grp =='All']
tda$measure = factor(tda$measure, levels = c("Hospital bed needs (mean)", "Non-ICU hospital bed needs (mean)", "ICU bed needs (mean)", "Ventilator needs (mean)"),
                     labels = c("Hospital bed needs (mean)", "Non-ICU hospital bed needs (mean)", "ICU bed needs (mean)", "Ventilator needs (mean)"))
p11 = getPlotMultiVoverlay(tda, title = sce.t, y.lab = 'Estimate %')


pdf(paste0(dir_res, 'Fig_sim.',loc.t,'_',seed.tag, '_',v.tag,'_', sce.tag, '.pdf'), width = 8, height = 6)
print(p8)
print(p9)
print(p10)
print(p11)
dev.off()
