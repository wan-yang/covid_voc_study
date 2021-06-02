# 3/16/21 - multiple tests using the eakf, 
# each round allow 1 hyp - change beta only, change S only, change both beta and S but more for S, change both but more for beta
# compare result and output best performing est

Fn_checkDA<-function(xnew,bound.low,bound.up){
  b.low=bound.low;
  b.up=bound.up;
  n.var=nrow(xnew); n.ens=ncol(xnew);
  for(vi in 1:n.var){
    #  Corrects if <b.low
    ug=min(xnew[vi,]);
    if (ug<b.low[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]<b.low[vi]){
          # xnew[vi,jj]=b.low[vi];
          xnew[vi,jj]=pmax(b.low[vi],runif(1, min=pmax(b.low[vi],quantile(xnew[vi,],.25)), max=pmax(b.low[vi],quantile(xnew[vi,],.75)))); # biased high
        }
      }
    }
    ug=max(xnew[vi,]);
    if (ug>b.up[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]>b.up[vi]){
          # xnew[vi,jj]=b.up[vi];
          # apply to non-S variables
          if(grepl('S', names(b.up[vi]))){
            xnew[vi,jj]=b.up[vi];
          } else {
            xnew[vi,jj]=runif(1, min=min(b.up[vi]/2,quantile(xnew[vi,],.5)), max=min(quantile(xnew[vi,],.80),b.up[vi]));
          }
          
        }
      }
    }
  }
  xnew;
}

Fn_getR0_SEIR=function(PARMS){
  with(as.list(PARMS), {
    Ro = beta.mean * Tir.mean
    
    Ro
  })
}

Fn_getRt_SEIR=function(PARMS){
  with(as.list(PARMS), {
    Rt = beta.mean * Tir.mean * S/N
    
    Rt
  })
}

fn_getRelMob = function(rel.mob.t, p.mob.t){ # return the scaled moblity for adjusting tx
  (rel.mob.t * p.mob.t) %>% pmin(1) # make sure it's <=1
}

library("truncnorm"); library("tgp"); library('mvtnorm'); # for lhs
library("MASS"); # for multivariate normal distribution

EAKF = function(epi.model=epi.model, num_ens=num_ens,inflat=1.03, 
                obs_i=obs_i, obs_vars_i=obs_vars_i, # case
                obs_d=obs_d, obs_vars_d=obs_vars_d,
                weeks=weeks,Week.starts=Week.starts,
                parm.bounds=parm.bounds, DA.bounds=DA.bounds, SR.bounds=SR.bounds, 
                parm.names = rownames(parm.bounds), rel.mob = rel.mob,
                state0=state0, state.names=rownames(state0),
                severity = severity,
                tm.ini=1, tmstep=7,
                newI.previous = NULL,
                SRparms = NULL
){
  # parm.bounds: prior bounds for the parms
  # x.prior: priors for all state variables, and model parameters
  # obs_i: observations
  # idx.obs: idx for obseverations to filter through
  
  if(!exists('wk.summer')) wk.summer = NULL
  if(!exists('IsLargeCountry')) IsLargeCountry = F
  if(!exists('wk.WkLowIFR')) wk.WkLowIFR = NULL
  if(!exists('excludeLockDown')) excludeLockDown = F
  
  # save inital condition passed in
  state00 = state0
  DA.bounds00=DA.bounds;
  SR.bounds00=SR.bounds;
  newI.previous00 = newI.previous;
  doSR00 = doSR
  
  p.major.cut = .25
  p.median.cut = .1
  
  if(!exists('redn.priority')) redn.priority = 1
  
  if(is.null(dim(obs_i))){ # only 1 observation (no subgroup)
    num_times_i=length(obs_i);  
  } else {
    num_times_i = nrow(obs_i)
  }
  if(is.null(dim(obs_d))){ # only 1 observation (no subgroup)
    num_times_d=length(obs_d);  
  } else {
    num_times_d = nrow(obs_d)
  }
  num_var = nrow(state0)
  
  cumlike=NULL; # cummulative likelihood
  

  # diff hypotheses
  
  # 'noSR', - no good
  hyps = c('beta.only.minor','beta.only.major',
           's.only.minor','s.only.major',
           'both.minor','both.major',
           'both.minorb.majors','both.majorb.minors',
           'both.flex',
           # for slow changes in s
           's.only.slow','both.minorb.slows', 'both.majorb.slows',
           'both.minorb.slows.minor', 'both.majorb.slows.minor'  # for BR
           ) 
  
  if(loc.t !='br'){
    hyps = hyps[!grepl('slows.minor', hyps)]  # to save time
  }
  
  rrmse = NULL
  
  for(ihyp in 1:length(hyps)){
    
    hyp.t = hyps[ihyp]
    
    print('', quote = F)
    print(hyp.t)
    
    # number parms allow the probe
    # n.parm.sr = ifelse(grepl('only',hyp.t), 1, 2)
    
    # slightly low for beta as it's more challenging and more restricted
    if (grepl('beta.only',hyp.t)){
      n.parm.sr = 1
      
    } else if(grepl('s.only',hyp.t)){
      n.parm.sr = 1 # 1
      
    } else if (hyp.t %in% c('both.major','both.flex')){
      n.parm.sr = 1.2 # 1.5 not so good # 1.3 # 4
    } else if (hyp.t %in% c('both.minor')){
      n.parm.sr = 1.2 # 1.5 # 2
    } else if (hyp.t %in% c('both.minorb.majors','both.majorb.minors')){
      n.parm.sr = 1.2 # 1.5 # 1.3 # 3
    } else if (hyp.t %in% c('both.majorb.slows', 'both.minorb.slows', 
                            'both.minorb.slows.minor', 'both.majorb.slows.minor')){
      # put all those with slows together, as the changes on s is minor - for large places like BR
      
      if(IsLargeCountry){
        # if it's a large country, slow transmission and traveling wave is more possible
        # otherwise, not much of good idea
        n.parm.sr = 1.1 # 1.25 # 1.2
      } else {
        n.parm.sr = 2 # 2
      }
    }
    
    wave2.1st.major = end1stWave # place holder for the first time during 2nd wave needing major SR
    S2ndwave0 = NULL; # record the susceptibility at the beginging of 2nd wave
    
    Smajor = 0; Smajor.cnt = 0
    # Smajor.cut = SRparms$Smajor.cut # 1.5 
    # Is it a good idea to change immediately? 

    # count number of major SR in the 2nd wave
    cntSRwave2 = 0 
    
    doSR = doSR00
    if(grepl('noSR', hyp.t))
      doSR = F # only allow minor adjustment to S
    
    
    donotUpdateSallWks = F; wks.srS4beta.only = 1:num_times_i;
    if(grepl('beta.only', hyp.t)){
      donotUpdateSallWks = T # only allow minor adjustment to S
      # allow it to update S after wave 1 and before wave 2
      wks.srS4beta.only = (end1stWave+0:7) %>% pmin(main2ndWave)
    }
      
    
    state0 = state00 # get the inital copy
    DA.bounds=DA.bounds00;
    SR.bounds=SR.bounds00;
    newI.previous = newI.previous00;
    
    # get the bounds for this hypothesis
    # for beta
    if(hyp.t %in% c('beta.only.minor','both.minor','both.minorb.majors','both.minorb.slows','both.minorb.slows.minor')){ # minor change for beta
      SR.bounds.wider.t = SR.bounds.wider
      SR.bounds.wider2.t = SR.bounds.wider  # only 1 lower level for beta
    } else if (hyp.t %in% c('beta.only.major','both.major')){ # minor change for beta
      SR.bounds.wider.t = SR.bounds.wider2 # both very high from the start
      SR.bounds.wider2.t = SR.bounds.wider2  # only 1 lower level for beta
    } else if (hyp.t %in% c('beta.only.grad','both.flex','both.majorb.slows','both.majorb.minors','both.majorb.slows.minor')){
      SR.bounds.wider.t = SR.bounds.wider
      SR.bounds.wider2.t = SR.bounds.wider2 
    } else {
      SR.bounds.wider.t = SR.bounds
      SR.bounds.wider2.t = SR.bounds
    }
    # for S - adjust the number and level of adjustment
    if(hyp.t %in% c('s.only.minor','both.minor','both.majorb.minors')){ # minor change for s
      s.upr.t = SRparms$s.upr.t.minor # upper bound for level of change in %S
      s.lwr.t = SRparms$s.lwr.t.minor # make the distribution wide
      p.cum.dS.large = SRparms$p.cum.dS.large.minor; 
      p.cum.dS.median = SRparms$p.cum.dS.median.minor; 
      p.cum.dS.small = SRparms$p.cum.dS.small.minor; 
      
      Smajor.cut.lwr = SRparms$Smajor.cut.lwr.minor  # only start at a later point?
      Smajor.cut.upr = SRparms$Smajor.cut.upr.minor
      Smajor.cnt.cut = SRparms$Smajor.cnt.cut.minor
      cntSadj_large.tot.t = SRparms$cntSadj_large.tot.minor # 2 #  # 0
      cntSadj_median.tot.t = SRparms$cntSadj_median.tot.minor # 3
      cntSadj_small.tot.t = SRparms$cntSadj_small.tot.minor # 4
      cntSadjtot.cut = SRparms$cntSadjtot.cut.minor
      cum.dS.cut = SRparms$cum.dS.cut.minor
      
      Spb_large.t = SRparms$Spb_large.minor # .25
      Spb_median.t = SRparms$Spb_median.minor  # .1
      Spb_small.t = SRparms$Spb_small.minor # .05
    } else if (hyp.t %in% c('both.minorb.majors','s.only.major')){ # 'both.major', minor change for beta
      s.upr.t = SRparms$s.upr.t.major # upper bound for level of change in %S
      s.lwr.t = SRparms$s.lwr.t.major # make the distribution wide
      p.cum.dS.large = SRparms$p.cum.dS.large.major; 
      p.cum.dS.median = SRparms$p.cum.dS.median.major; 
      p.cum.dS.small = SRparms$p.cum.dS.small.major; 
      
      Smajor.cut.lwr = SRparms$Smajor.cut.lwr.major  # only start at a later point?
      Smajor.cut.upr = SRparms$Smajor.cut.upr.major
      Smajor.cnt.cut = SRparms$Smajor.cnt.cut.major
      cntSadj_large.tot.t = SRparms$cntSadj_large.tot.major # 8 # 4 + 2 # 2
      cntSadj_median.tot.t = SRparms$cntSadj_median.tot.major # 8 # 6 + 2 # 3
      cntSadj_small.tot.t = SRparms$cntSadj_small.tot.major # 8 + 2 # 4
      cntSadjtot.cut = SRparms$cntSadjtot.cut.major
      Spb_large.t = SRparms$Spb_large.major # .25
      Spb_median.t = SRparms$Spb_median.major  # .1
      Spb_small.t = SRparms$Spb_small.major # .05
      cum.dS.cut = SRparms$cum.dS.cut.major
    } else if (hyp.t %in% c('both.major')){ # 'both.major', minor change for beta
      s.upr.t = SRparms$s.upr.t.max # upper bound for level of change in %S
      s.lwr.t = SRparms$s.lwr.t.max # make the distribution wide
      p.cum.dS.large = SRparms$p.cum.dS.large.max; 
      p.cum.dS.median = SRparms$p.cum.dS.median.max; 
      p.cum.dS.small = SRparms$p.cum.dS.small.max; 
      
      Smajor.cut.lwr = SRparms$Smax.cut.lwr.max  # only start at a later point?
      Smajor.cut.upr = SRparms$Smax.cut.upr.max
      Smajor.cnt.cut = SRparms$Smax.cnt.cut.max
      cntSadj_large.tot.t = SRparms$cntSadj_large.tot.max # 8 # 4 + 2 # 2
      cntSadj_median.tot.t = SRparms$cntSadj_median.tot.max # 8 # 6 + 2 # 3
      cntSadj_small.tot.t = SRparms$cntSadj_small.tot.max # 8 + 2 # 4
      cntSadjtot.cut = SRparms$cntSadjtot.cut.max
      Spb_large.t = SRparms$Spb_large.max # .25
      Spb_median.t = SRparms$Spb_median.max  # .1
      Spb_small.t = SRparms$Spb_small.max # .05
      cum.dS.cut = SRparms$cum.dS.cut.max
    } else if (hyp.t %in% c('both.flex')){ # minor change for beta
      s.upr.t = SRparms$s.upr.t.flex # upper bound for level of change in %S
      s.lwr.t = SRparms$s.lwr.t.flex # make the distribution wide
      p.cum.dS.large = SRparms$p.cum.dS.large.flex; 
      p.cum.dS.median = SRparms$p.cum.dS.median.flex; 
      p.cum.dS.small = SRparms$p.cum.dS.small.flex; 
      
      Smajor.cut.lwr = SRparms$Smajor.cut.lwr.flex  # only start at a later point?
      Smajor.cut.upr = SRparms$Smajor.cut.upr.flex
      Smajor.cnt.cut = SRparms$Smajor.cnt.cut.flex
      cntSadj_large.tot.t = SRparms$cntSadj_large.tot.flex # 8 # 4 + 2 # 2
      cntSadj_median.tot.t = SRparms$cntSadj_median.tot.flex # 8 # 6 + 2 # 3
      cntSadj_small.tot.t = SRparms$cntSadj_small.tot.flex # 8 + 2 # 4
      cntSadjtot.cut = SRparms$cntSadjtot.cut.flex
      Spb_large.t = SRparms$Spb_large.flex # .25
      Spb_median.t = SRparms$Spb_median.flex  # .1
      Spb_small.t = SRparms$Spb_small.flex # .05
      cum.dS.cut = SRparms$cum.dS.cut.flex
    } else if (hyp.t %in% c('s.only.slow','both.minorb.slows', 'both.majorb.slows')){ # minor change for beta
      s.upr.t = SRparms$s.upr.t.slow # upper bound for level of change in %S
      s.lwr.t = SRparms$s.lwr.t.slow # make the distribution wide
      p.cum.dS.large = SRparms$p.cum.dS.large.slow; 
      p.cum.dS.median = SRparms$p.cum.dS.median.slow; 
      p.cum.dS.small = SRparms$p.cum.dS.small.slow; 
      
      Smajor.cut.lwr = SRparms$Smajor.cut.lwr.slow  # only start at a later point?
      Smajor.cut.upr = SRparms$Smajor.cut.upr.slow
      Smajor.cnt.cut = SRparms$Smajor.cnt.cut.slow
      cntSadj_large.tot.t = SRparms$cntSadj_large.tot.slow # 8 # 4 + 2 # 2
      cntSadj_median.tot.t = SRparms$cntSadj_median.tot.slow # 8 # 6 + 2 # 3
      cntSadj_small.tot.t = SRparms$cntSadj_small.tot.slow # 8 + 2 # 4
      cntSadjtot.cut = SRparms$cntSadjtot.cut.slow
      Spb_large.t = SRparms$Spb_large.slow # .25
      Spb_median.t = SRparms$Spb_median.slow  # .1
      Spb_small.t = SRparms$Spb_small.slow # .05
      cum.dS.cut = SRparms$cum.dS.cut.slow
    } else if (hyp.t %in% c('both.minorb.slows.minor','both.majorb.slows.minor')){ # minor change for beta
      s.upr.t = SRparms$s.upr.t.slow.minor # upper bound for level of change in %S
      s.lwr.t = SRparms$s.lwr.t.slow.minor # make the distribution wide
      p.cum.dS.large = SRparms$p.cum.dS.large.slow.minor; 
      p.cum.dS.median = SRparms$p.cum.dS.median.slow.minor; 
      p.cum.dS.small = SRparms$p.cum.dS.small.slow.minor; 
      
      Smajor.cut.lwr = SRparms$Smajor.cut.lwr.slow.minor  # only start at a later point?
      Smajor.cut.upr = SRparms$Smajor.cut.upr.slow.minor
      Smajor.cnt.cut = SRparms$Smajor.cnt.cut.slow.minor
      cntSadj_large.tot.t = SRparms$cntSadj_large.tot.slow.minor # 8 # 4 + 2 # 2
      cntSadj_median.tot.t = SRparms$cntSadj_median.tot.slow.minor # 8 # 6 + 2 # 3
      cntSadj_small.tot.t = SRparms$cntSadj_small.tot.slow.minor # 8 + 2 # 4
      cntSadjtot.cut = SRparms$cntSadjtot.cut.slow.minor
      Spb_large.t = SRparms$Spb_large.slow.minor # .25
      Spb_median.t = SRparms$Spb_median.slow.minor  # .1
      Spb_small.t = SRparms$Spb_small.slow.minor # .05
      cum.dS.cut = SRparms$cum.dS.cut.slow.minor
    } else {
      s.upr.t = s.lwr.t = Smajor.cut.lwr = Smajor.cut.upr = Smajor.cnt.cut = 
        cntSadj_large.tot.t = cntSadj_median.tot.t = cntSadj_small.tot.t = cntSadjtot.cut = 
        cum.dS.cut = Spb_large.t = Spb_median.t = Spb_small.t = 0
      p.cum.dS.large = p.cum.dS.median = p.cum.dS.small = 0
    }
    
    # Spb.adj_upr = 1.5; Spb.adj_lwr = .5
    # Spb.adj_upr = 1; Spb.adj_lwr = 1
    # Spb.adj_upr = 1.25; Spb.adj_lwr = .75 # likely help
    Spb.adj_upr = SRparms$Spb.adj_upr; Spb.adj_lwr = SRparms$Spb.adj_lwr

    # 4/18/21 need to monitor the cumulative change in S
    cum.dS.perc = 0; # percentage, relative to total from first wave
    cum.dS.po = 0; # posterior, only include weeks with probing on s
    stat.sr.S = NULL # record the details, wk, level of sr, etc
    wk.sr.S = rep(0, num_times_i)
    over.adjS = F; # record wether it over adjust S when doing SR
    wk.srS.last = num_times_i + 1
    stop.srS = F
    res.imm = NULL
    # p.S2beta = 3; # test 12 the penalty level be higher for S than beta
    p.S2beta = 1; # test 13 - yes, better!
    # initialize system
    xprior=array(0,c(num_var,num_ens,num_times_i+1));
    xpost=array(0,c(num_var,num_ens,num_times_i));
    dimnames(xpost)[1] = list(state.names)
    
    # to estimate daily cases
    xprior.daily = matrix(0,(num_times_i+1)*tmstep, num_ens);
    xpost.daily = matrix(0,(num_times_i)*tmstep, num_ens);
    Eprior = Iprior = numeric(num_times_i+1)
    Epost = Ipost = Spost = Itotpost = numeric(num_times_i)
    
    # to record SR for each step
    SRflag = data.table(time = 1: nrow(obs_i), SR = 'minor')
    # cntSadj = 0 # count number times adjust S
    cntSadj_large = cntSadj_small = cntSadj_median = 0; cntSadjtot = 0
    cntUpdateSRbound = 0
    update2widerSR1 = update2widerSR2 = F
    wk.update2widerSR2 = num_times_i
    
    
     
    
    {
      dist_tm.to.detect = NULL
      for(ii in 1:num_ens){
        tmp = generation.time(dist_tm.to.detect.name,c(state0['Td.mean',ii],state0['Td.sd',ii]),truncate = tm.to.detect.max)
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
    }
    
    # integrate forward 1 step to get the prior
    tm_strt = tm.ini+1; tm_end = tm.ini + tmstep;
    cur.wk = weeks[1];
    seed = seed
    vdate.t = Week.starts[1]
    
    severity.t = severity
    
    # get prior for week 1
    {
      if(!is.null(newI.previous) & vdate.t >= vax.start){
        tm.t = nrow(newI.previous)
        # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
        # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
        # 2/5/21 set t1 to 1 yr given the slow rollout
        t1 = 365
        cumI.t = colSums(newI.previous[pmax(1,tm.t-t1) : (tm.t),]) #  %>% apply(1, median) # excl last two weeks and get the median
        # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
        # higher infection rate for the priority group
        # tm.t = pmax(1, tm.t - t1 + 1) # re-aline timing with start of vac
        tm.imm = 365*2.5 # assume 3 yr immunity
        p.imm.wane.max = .8; k = .015  # 1/tm.imm  
        p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(tm.t + 60 - tm.imm/2))) # not all infected from day 1
        # earlier infections are likely to be in the high priority groups 
        p.imm = 1 *  p.imm.wane * redn.priority # assume 100% prior infection provide immunity, but wane over time
        # and multiple by % excluded if there is prior testing before vax: p.prior.test
        percSmax.t = 1 - cumI.t / N * p.imm
        # no lower than 50%, in case of outliers
        percSmax.t = pmax(percSmax.t, .5)
        # print(c('cohort %S:',round(summary(mean(percSmax.t)),2)), quote = F)
      } else {
        percSmax.t = 0
        # print('no vax yet')
      }
      
      beta_tt = state0['beta',] * fn_getRelMob(rel.mob[1,], state0['p.mob',])
      if(seasonality) {
        # beta_tt = state0['beta',] * seasonal.cycle[seasonal.cycle$week==cur.wk,]$relR0 
        beta_tt = beta_tt * relR0[cur.wk,] 
      } 
      
      if(epi.model == 'SEIRS'){
        
        simEpi=SEIRS(tm_strt, tm_end, tm_step=1, # 1 day time-step
                     tmstep = tmstep,
                     state0 = state0,
                     S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                     I0=state0[paste0('I',1:num_gr),], 
                     beta=beta_tt, 
                     Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                     seed=seed, stoch=stoch, 
                     severity = severity.t,
                     newI.previous = newI.previous,
                     dist_tm.to.detect = dist_tm.to.detect,
                     dist_tm.to.death = dist_tm.to.death) # include the delay reporting etc.
        
      } else if(epi.model == 'SEIRSV') {
        
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
        simEpi=SEIRSV(tm_strt, tm_end, tm_step=1, # 1 day time-step
                      tmstep = tmstep,
                      state0 = state0,
                      S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                      I0=state0[paste0('I',1:num_gr),], 
                      beta=beta_tt, 
                      Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                      seed=seed, stoch=stoch, 
                      severity = severity.t,
                      newI.previous = newI.previous,
                      dist_tm.to.detect = dist_tm.to.detect,
                      dist_tm.to.death = dist_tm.to.death,
                      percSmax.t = percSmax.t,
                      V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                      # these are total number of vaccinees with unknown immunity
                      # but pre-ajust for time lag from vaccination to immune protection
                      VE1 = VE1, VE2=VE2 # Vaccine efficacy, need further adjustment by prior immunity 
        ) # include the delay reporting etc.
      }
      
      # re-assemble to the same order as the prior: state0
      n.end = tm_end - tm_strt + 2
      state.new = NULL
      for(i in 1:(length(simEpi)-1)){
        tmp = simEpi[[i]][n.end,,drop=F]; 
        rownames(tmp)=gsub('cumI','newI',paste0(names(simEpi)[i],1:num_gr))
        state.new = rbind(state.new,tmp)
      }
      
      state.new = rbind(state.new, state0[parm.names,])
      state.new = state.new[rownames(state0),] # make sure the order is the same
      
      
      xprior[,,1]= state.new
      xprior.daily[1:tmstep,] = simEpi$daily.newItot
      # save prior for E and I
      Eprior[1] = state.new['E1',] %>% mean
      Iprior[1] = state.new['I1',] %>% mean
    }
    
    dx.t.newitot = dx.t.newiobs = dx.t.e = dx.t.i = numeric(num_times_i)
    #### Begin looping through observations
    for (tt in 1:num_times_i){ # num_times_i
      
      # update detection rate
      if(tt == 5 & loc.t == 'uk') {
        DA.bounds['alpha',1] = .03
        
        SR.bounds['alpha',1] = .1
        SR.bounds['alpha',2] = SRalpha_bounds2[2] * .8
      } 
      
      
      
      if((tt == 5 & loc.t == 'br') | (tt == 5 & loc.t == 'sa')) {
        DA.bounds['alpha',1] = .02  # reduce uncertainty
        DA.bounds['alpha',2] = pmax(DA.bounds['alpha',2], .2)
        SR.bounds['alpha',2] = pmax(SRalpha_bounds2[2] * .8, alpha_bounds2[2])
        
        SR.bounds['ifr',] = SRifr_bounds1
      }
      
      if(tt == end1stWave  & loc.t %in% c('sa', 'br')){  # only for the uk
        DA.bounds['alpha',] = DAalpha_bounds2
        SR.bounds['alpha',] = SRalpha_bounds2
        
        SR.bounds.wider.t['alpha',] = SRalpha_bounds2
        SR.bounds.wider2.t['alpha',] = SRalpha_bounds2
        
      }
      
      if(loc.t == 'uk'){
        if(tt %in% wk.summer){  # only for the uk
          SR.bounds['alpha',] = SRalpha_bounds.summer
        } else if(tt == max(wk.summer)+1){  # only for the uk
          DA.bounds['alpha',] = DAalpha_bounds2
          SR.bounds['alpha',] = SRalpha_bounds2
          
          SR.bounds.wider.t['alpha',] = SRalpha_bounds2
          SR.bounds.wider2.t['alpha',] = SRalpha_bounds2
          
        }
      }
      
      
      
      
      if(tt == end1stWave + 10 & loc.t == 'uk') {
        # DA.bounds['alpha',1] = .15
        SR.bounds['alpha',] = SRalpha_bounds3
        DA.bounds['alpha',] = DAalpha_bounds3
        SR.bounds.wider.t['alpha',] = SRalpha_bounds3
        SR.bounds.wider2.t['alpha',] = SRalpha_bounds3
        
      }
      
      # update IFR bounds
      if(loc.t %in% c('uk','sa', 'br')){
        if(tt %in% wk.WkLowIFR){  # WkLowIFR low ifr
          SR.bounds['ifr',] = SRifr_bounds.WkLowIFR
          SR.bounds.wider.t['ifr',] = SRifr_bounds.WkLowIFR
          SR.bounds.wider2.t['ifr',] = SRifr_bounds.WkLowIFR
        } else if(tt == wk.2strtHighIFR){ # make it later as the IFR increases usually lag
          # else if(tt == max(wk.WkLowIFR)+1){  #  & loc.t == 'uk'
          SR.bounds['ifr',] = SRifr_bounds2
          SR.bounds.wider.t['ifr',] = SRifr_bounds2
          SR.bounds.wider2.t['ifr',] = SRifr_bounds2
        }
        if(loc.t %in% c('uk', 'br')){
          if(tt == wk.2strtHigherIFR){ # even higher IFR which the system collapse
            SR.bounds['ifr',] = SRifr_bounds3
            SR.bounds.wider.t['ifr',] = SRifr_bounds3
            SR.bounds.wider2.t['ifr',] = SRifr_bounds3
          }
        }
      }

      
      # if(tt == 20 & loc.t == 'uk'){ # update ifr bounds for the WkLowIFR - loss of young people get infected & low ifr
      #   SR.bounds['ifr',] = SRifr_bounds2
      # }
      
      # prevent it shifting off
      if(F){
        if(tt == 10 & loc.t %in% paste0('sce',1:4)) {
          SR.bounds['alpha',1] = SR.bounds['alpha',1] * 1.2  # reduce uncertainty
          SR.bounds['alpha',2] = SR.bounds['alpha',2] * 1.2
        }
      }
      
      
      
      # print(tt)
      # only inflate the observed
      # inflat=diag(x=c(1,1,1,1,1,1,lambda),7,7);
      # inflat all states and parameters
      xmn=rowMeans(xprior[,,tt]); 
      xnew=inflat*(xprior[,,tt]-xmn%*%matrix(1,1,num_ens))+xmn%*%matrix(1,1,num_ens);
      
      for(H in idx.obs_i){  # LOOP THROUGH OBSERVATIONS
        ####  Get the variance of the ensemble
        ####  Define the mapping operator H
        ####  HH is the number of variables and parameters in state space
        ####  H is only the number of columns being observed
        io = H - min(idx.obs_i) + 1
        obs_var = obs_vars_i[tt,io];
        prior_var = pmax(var(xnew[H,]),1);
        # post_var = 1/(1/prior_var + 1/obs_var);
        post_var = prior_var*obs_var/(prior_var+obs_var);
        
        prior_mean = mean(xnew[H,]);
        post_mean = post_var*(prior_mean/prior_var + obs_i[tt,io]/obs_var);
        
        #### Compute alpha and adjust distribution to conform to posterior moments
        alp = sqrt(obs_var/(obs_var+prior_var));
        
        dy = post_mean + alp*(xnew[H,]-prior_mean)-xnew[H,];
        ###  Getting the covariance of the prior state space and
        ###  observations  (which could be part of state space, e.g. infections)
        rr=NULL;
        for (j in 1:num_var){
          C=cov(xnew[j,],xnew[H,])/prior_var;  # covariance/varance of x.obs
          if((donotUpdateS1stWave & tt < end1stWave) |
             (donotUpdateSallWks & ! tt %in% wks.srS4beta.only)  | # for 'beta.only'
             # cum.dS.perc > cum.dS.cut # after reaching the threshold for probing S
             tt > wk.srS.last # yes, this seems to help
             ){
            if(grepl('S',state.names[j])) C = C / 10
          }
          if(rednUpdateEI  & tt %in% tm_rednUpdateEI){
            if((grepl('E',state.names[j]) | grepl('I',state.names[j])) &
               !grepl('new',state.names[j])){
              C = C / 10 # reduce the level of update by a factor of 10
            }
          }
          # restrict probing if S is the only change
          if(grepl('s.only',hyp.t) & tt> end1stWave + 5){ # & tt> end1stWave
            if(state.names[j] %in% SR.var.tx){
              C = C / 10 # reduce the level of update by a factor of 10
            }
          }
          rr=append(rr,C);
        }
        dx=rr%*%t(dy);
        
        # record the level of adjustment
        dx.t = rowMeans(dx); names(dx.t) = state.names
        dx.t.newitot[tt] = dx.t['newItot1']
        dx.t.newiobs[tt] = dx.t['newIobs1']
        dx.t.e[tt] = dx.t['E1']
        dx.t.i[tt] = dx.t['I1']
        
        ###  Get the new ensemble and save prior and posterior
        xnew = xnew + dx;  # ADJUSTED NUM_OBS TIMES, ONCE BY EACH OBSERVATION
      }
      
      # filter using the mortability data
      if(tt %in% 1:num_times_d){
        for(H in idx.obs_d){  # LOOP THROUGH OBSERVATIONS
          ####  Get the variance of the ensemble
          ####  Define the mapping operator H
          ####  HH is the number of variables and parameters in state space
          ####  H is only the number of columns being observed
          io = H - min(idx.obs_d) + 1
          obs_var = obs_vars_d[tt,io];
          prior_var = pmax(var(xnew[H,]),1);
          # post_var = 1/(1/prior_var + 1/obs_var);
          post_var = prior_var*obs_var/(prior_var+obs_var);
          
          prior_mean = mean(xnew[H,]);
          post_mean = post_var*(prior_mean/prior_var + obs_d[tt,io]/obs_var);
          
          #### Compute alpha and adjust distribution to conform to posterior moments
          alp = sqrt(obs_var/(obs_var+prior_var));
          
          dy = post_mean + alp*(xnew[H,]-prior_mean)-xnew[H,];
          ###  Getting the covariance of the prior state space and
          ###  observations  (which could be part of state space, e.g. infections)
          rr=NULL;
          for (j in 1:num_var){
            C=cov(xnew[j,],xnew[H,])/prior_var;  # covariance/varance of x.obs
            if((donotUpdateS1stWave & tt < end1stWave) |
               (donotUpdateSallWks & ! tt %in% wks.srS4beta.only)  # | # for 'beta.only'
               # cum.dS.perc > cum.dS.cut # after reaching the threshold for probing S
               # tt > wk.srS.last
               ){
              if(grepl('S',state.names[j])) C = C / 10
            }
            if(rednUpdateEI & tt %in% tm_rednUpdateEI){ # & tt> end1stWave
              if((grepl('E',state.names[j]) | grepl('I',state.names[j])) &
                 !grepl('new',state.names[j])){
                C = C / 10 # reduce the level of update by a factor of 10
              }
            }
            
            # restrict probing if S is the only change
            if(grepl('s.only',hyp.t) & tt> end1stWave + 5){ # & tt> end1stWave
              if(state.names[j] %in% SR.var.tx){
                C = C / 10 # reduce the level of update by a factor of 10
              }
            }
            
            rr=append(rr,C);
          }
          dx=rr%*%t(dy);
          
          ###  Get the new ensemble and save prior and posterior
          xnew = xnew + dx;  # ADJUSTED NUM_OBS TIMES, ONCE BY EACH OBSERVATION
        }
      }
      
      #  Corrections to DA produced aphysicalities
      xnew=Fn_checkDA(xnew, bound.low = DA.bounds[,1], bound.up =  DA.bounds[,2]);
      row.names(xnew)=state.names
      xpost[,,tt]=xnew;
      state0 = xnew
      
      # save post for E and I
      Epost[tt] = state0['E1',] %>% mean
      Ipost[tt] = state0['I1',] %>% mean
      Spost[tt] = state0['S1',] %>% mean
      Itotpost[tt] = state0['newItot1',] %>% mean
      
      # update daily estimates as well
      {
        
        # avoid dividing by 0 here!
        if(tmstep>1){
          f.adj = state0[idx.newItot,]/colSums(xprior.daily[1:tmstep+(tt-1)*tmstep,])
        } else {
          f.adj = state0[idx.newItot,]/(xprior.daily[1:tmstep+(tt-1)*tmstep,])
        }
        
        
        {
          i0=which(is.na(f.adj) | is.infinite(f.adj)); 
          xx=1:num_ens; 
          i.non0=xx[!(xx %in% i0)]
          f.adj[i0] = sample(f.adj[i.non0],size=length(i0),replace = T)
        }
        f.adj_all = matrix(f.adj, tmstep, num_ens, byrow = T)
        # if obs=0 for the whole week, set all days to 0 - that is taken care of by default
        # so all the problem comes from obs!=0, but the prior says all days have 0 cases
        # in that case, distribute the cases evenly
        xpost.daily[1:tmstep+(tt-1)*tmstep,]=xprior.daily[(1:tmstep)+(tt-1)*tmstep,]*f.adj_all
        
      } # update daily estimates
      #  Integrate forward one time step
      tcurrent = tm.ini+tmstep*tt;
      tm_strt = tcurrent+1; tm_end = tcurrent + tmstep
      cur.wk = weeks[tt+1];
      if(is.na(cur.wk))
        cur.wk = weeks[tt];
      # print(paste('Week:',cur.wk),quote=F)
      vdate.t = Week.starts[pmin(tt+1,length(Week.starts))]
      
      newI.previous = xpost.daily[1:(tt*tmstep),] # these are at weekly level!
      
      
      # check the level of adjust for newItot to gauge how the filter has been working
      mn.pr.p = xprior[idx.newItot,,tt-1] %>% mean
      mn.po.p = xpost[idx.newItot,,tt-1] %>% mean
      mn.pr = xprior[idx.newItot,,tt] %>% mean
      mn.po = xpost[idx.newItot,,tt] %>% mean
      e.po = xpost[idx.e,,tt] %>% mean
      i.po = xpost[idx.i,,tt] %>% mean
      e.pr = xprior[idx.e,,tt] %>% mean
      i.pr = xprior[idx.i,,tt] %>% mean
      
      # SR - depending on hypthesis
      if(doSR){
        p.sr.beta = 1
        SR.var.local.t = SR.var.local
        extraSRbeta = F # whether extra SR on beta is needed
        
        # for updating beta and SR bounds for probing
        npre.tot1 = 2; npre.tot2 = 4
        n.cut1 = 2; n.cut2 = ifelse(hyp.t == "beta.only.major", 4, 4)
        if(grepl('slow', hyp.t)){ # if it's a traveling wave, make update of beta slower as well
          if(grepl('minorb',hyp.t)){
            n.cut1 = 10; n.cut2 = 16
            npre.tot1 = 16; npre.tot2 = 20
          } else if (grepl('majorb',hyp.t)) {
            n.cut1 = 10; n.cut2 = 16
            npre.tot1 = 16; npre.tot2 = 20
          }
          
        }
        
        if(tt > 2){
          # didnotovershot = ((Eprior[tt-1] - Epost[tt-1])/Epost[tt-1] < .25 & mn.po > N * .2/100) &
          #   ((Eprior[tt] - Epost[tt])/Epost[tt] < .25 & mn.po > N * .2/100)
          
          overshot = ((Eprior[tt-1] - Epost[tt-1])/Epost[tt-1] > .25 & mn.po > N * .2/100) &
            ((Eprior[tt] - Epost[tt])/Epost[tt] > .25 & mn.po > N * .2/100) & 
            (mn.pr.p - mn.po.p)/mn.po.p > .1  & 
            (mn.pr - mn.po)/mn.po > .1  
          
          didnotovershot = !overshot
          
          ascend = (mn.po.p - mn.pr.p)/mn.po > .1 & 
            (mn.po - mn.pr)/mn.po > .1 &  
            (obs_i[tt-1] - obs_i[tt-2])/obs_i[tt-2] > .1 & 
            (obs_i[tt] - obs_i[tt-1])/obs_i[tt-1] > .1  # filter is trying hard to get the numbers up
          descend = (mn.pr.p - mn.po.p)/mn.po.p > .1  
            (mn.pr - mn.po)/mn.po > .1 & 
            (obs_i[tt-2] - obs_i[tt-1])/obs_i[tt-1] > .1 & 
            (obs_i[tt-1] - obs_i[tt])/obs_i[tt] > .1 # filter is trying hard to get the numbers down
          
        }
        
        # same for wave 1 for all
        if(tt <= end1stWave){
          if(((mn.po - mn.pr) / mn.pr > p.major.cut  | (mn.po - mn.pr) > N * .2/100
          ) & mn.po > N * .2/100){ # only when epi is increasing
            
            
            SR.perc.full.t = pmax(SR.perc.full * 2, percSRmajor) # do not increase % to too large
            SR.perc.local.t = pmax(SR.perc.local * 2, .1)
            
            SR.var.full.t = c(SR.var.full,'alpha', 'ifr', 'p.mob') %>% unique()
            SRbounds.full.t = SR.bounds[SR.var.full.t,]
            
            SRflag[time == tt]$SR = 'major'
            
            print(paste(tt, 'need major SR'), quote = F)
            
            
          } else if (abs(mn.pr - mn.po) / mn.pr > p.median.cut & mn.po > N * .2/100){
            
            # 
            # if(tt < end1stWave | (mn.po - mn.pr) / mn.pr > .3)
            #  SR.perc.full.t = pmax(SR.perc.full * 2, .05) # not good for 2nd wave, only increase it for wave 1 or when it's increasing?
            
            SR.perc.full.t = pmax(SR.perc.full * 2, .05) # not good for 2nd wave
            SR.perc.local.t = pmax(SR.perc.local * 2, .1)
            
            
            SR.var.full.t = c(SR.var.full,'p.mob', 'ifr')  %>% unique()  # ,'alpha', 'ifr'
            SRbounds.full.t = SR.bounds[SR.var.full.t,]
            # SRbounds.full.t['alpha', 1] = alpha_bounds[1]
            
            SRflag[time == tt]$SR = 'median'
            print(paste(tt, 'need additional SR'), quote = F)
            
          } else {
            
            SR.perc.full.t = SR.perc.full
            SR.perc.local.t = SR.perc.local
            
            SR.var.full.t = SR.var.full
            SRbounds.full.t = SR.bounds[SR.var.full.t,]
          }
          
          
        } # END ALL FIRST WAVE
        
        if(tt == end1stWave){
          state0.2ndWave.strt = state0
          cumIwave1 = N - mean(state0.2ndWave.strt['S1',])
          if(cumIwave1 < N * .16){
            # if the first wave is small, restrict the time major SR on S can take as too many times increase uncertainty
            cntSadjtot.cut = pmin(5, cntSadjtot.cut)
          }
          
          # if the prior immunity is low, chances of increase due to immune evasion is low
          if (cumIwave1 < N * .16 & grepl('beta.only',hyp.t)){
            n.parm.sr = n.parm.sr * .5
            # cntSadjtot.cut = pmin(5, cntSadjtot.cut)
          } 
          
          if(is.null(S2ndwave0)) 
            S2ndwave0 = Spost[end1stWave]
          
        } 
          
        
        # FOR 2ND WAVE
        if(tt > end1stWave){
          
          # compute the posterior of increase in S - immune evasion
          if(wk.sr.S[tt-1] # |
            # (stop.srS & tt <= wk.srS.last)
             ){ # last week it probe s
            
            dS.t = xpost['S1',,tt] - xpost['S1',, tt-1] + xpost['newItot1',,tt]  # not good b/c the eakf uses the mean for indiv variables
            # dS.t = mean(xpost['S1',,tt]) - mean(xpost['S1',, tt-1]) + mean(xpost['newItot1',,tt])
            
            xx1 = xpost['S1',,tt]; xx1mn = xx1 %>% mean; 
            xx2 = xpost['S1',, tt-1]; xx2mn = xx2 %>% mean; 
            xx3=xpost['newItot1',,tt]; xx3mn = xx3 %>% mean; 
            # var1 = var(xx1); var2 = var(xx2);  var3=var(xx3)
            m.cov = cov(x=cbind(xx1, xx2, xx3))
            
            dS.sd = deltamethod(~ x1 - x2 + x3, mean = c(xx1mn, xx2mn, xx3mn), 
                                  cov = m.cov
            ) %>% sqrt
            
            dS.mn = dS.t  %>% mean
            
            stat.sr.S = rbind(stat.sr.S, 
                              data.table(week = tt-1, mean = dS.mn,
                                         sd = dS.sd,
                                         dS.t %>% t
                                         ))
            
            if(wk.sr.S[tt-1]){ # only do this if last step pr on S
              cum.dS.po = cum.dS.po + (dS.mn %>% pmax(0))
              res.imm = (cumIwave1 - sum(stat.sr.S$mean)) %>% pmax(0)  # resisual immunity
              
              # 4/25/21 further consider waning immunity - those have return to S already if a long time goes by
              limm = mean(state0['Trs',])
              res.imm = res.imm * exp(-(tt - end1stWave)*7 /limm)
            }
            
          }
          
          
          # the tigher SR bound for beta, etc
          {
            
            sr.local.mean = rowMeans(state0[parm.names,]) # local SR bounds
            sr.local.lwr = apply(state0[parm.names,], 1, quantile, .1) # .2
            sr.local.upr = apply(state0[parm.names,], 1, quantile, .9) # .7
            sr.local.lwr = pmin(sr.local.mean * .7, sr.local.lwr) # .75
            sr.local.upr = pmax(sr.local.mean * 1.3, sr.local.upr) # 1.25
            sr.local.lwr = pmax(sr.local.lwr, DA.bounds[parm.names,1]) # .75
            sr.local.upr = pmin(sr.local.upr, DA.bounds[parm.names,2]) # 1.25
            
            
            SR.bounds.tight =  cbind(sr.local.lwr, sr.local.upr)
            # this allow beta to drift, so not working well
            # replace them with the original parm bound
            SR.bounds.tight[SR.var.full,] = parm.bounds[SR.var.full,]
            
            
            # update ifr
            if(loc.t %in% c('uk','br','sa')){
              if(tt %in% (wk.WkLowIFR)){
                SR.bounds.tight['ifr',1] = pmin(SR.bounds.tight['ifr',1], SR.bounds['ifr',1])
                SR.bounds.tight['ifr',2] = pmin(SR.bounds.tight['ifr',2], SR.bounds['ifr',2])
              } else {
                SR.bounds.tight['ifr',] = SR.bounds['ifr',]
              }
            }
            
              
            # lower the upper bound for beta and Tir?
            # SR.bounds.tight[c('beta'),2] = parm.bounds[c('beta'),2] * .8
            # SR.bounds.tight[c('Tir'),2] = parm.bounds[c('Tir'),2] * .8
            # SR.bounds.tight[c('Tei'),1] = 3
          }
          
          
          # for others parameters
          if(((mn.po - mn.pr) / mn.pr > p.major.cut | (mn.po - mn.pr) > N * .2/100 |
              dx.t.newitot[tt] > N * .1/100
          ) & mn.po > N * .2/100){ # only when epi is increasing
            
            
            p.sr.beta = ifelse(tt< main2ndWave, 1, ((mn.po - mn.pr) / mn.pr / .1) %>% pmin(4))
            if(! grepl('s.only', hyp.t) &
               ! grepl('slow', hyp.t)  # don't do it if it's changing slowly
               ){
              perc.extra.sr.beta = (SR.perc.full * 2 * p.sr.beta - percSRmajor) %>% pmin(.1)
              if(perc.extra.sr.beta > .05){
                extraSRbeta = T
                
              }
              
              
            }
            
            
            
            
            # (xpost_mean$newItot1 - xprior_mean$newItot1)/mn.pr
            # (xpost_mean$newItot1 - xprior_mean$newItot1)
            
            # if the level of change is extremely large, adjust SRflag to make sure it gets updated earlier
            if((((mn.po - mn.pr) / mn.pr >= p.major.cut * 1.5 & (mn.po - mn.pr) / mn.pr < p.major.cut * 2) | 
                ((mn.po - mn.pr) >= N * .3/100 & (mn.po - mn.pr) < N * .7/100) |
                dx.t.newitot[tt] > N * .15/100
            ) & mn.po > N * .4/100 & 
            tt >= main2ndWave  # only do this after the major 2nd wave
            ){
              
              if(update2widerSR1){  # already updated beta bound once
                # SRflag[time %in% c(tt-1, tt-2, tt-3)]$SR = 'major'
                # may be to aggressive
                SRflag[time %in% c(tt-1, tt-2)]$SR = 'major'
              } else {
                SRflag[time %in% (tt-1)]$SR = 'major'
              }
              
              
            } else if (((mn.po - mn.pr) / mn.pr >= p.major.cut * 2 | (mn.po - mn.pr) >= N * .7/100 |
                        # (dx.t.newitot[tt] > N * .15/100 & dx.t.newitot[tt] > N * .3/100) ?
                        (dx.t.newitot[tt-1] > N * .15/100 & dx.t.newitot[tt] > N * .3/100)
            ) & mn.po > N * .4/100 & tt >= main2ndWave # only do this after the major 2nd wave
            ){
              # may be to aggressive
              # SRflag[time %in% c(tt-1, tt-2, tt-3)]$SR = 'major'
              SRflag[time %in% c(tt-1, tt-2)]$SR = 'major'
            } 

            
            
            if(cntSRwave2 == 0)
              wave2.1st.major = tt
            
            
              
            
            # SR.perc.full.t = pmax(SR.perc.full * 2 * p.sr.beta, percSRmajor) # only do it for beta 
            SR.perc.full.t = pmax(SR.perc.full * 2, percSRmajor)
            SR.perc.local.t = pmax(SR.perc.local * 2, .1)
            
            SR.var.full.t = c(SR.var.full,'alpha', 'ifr', 'p.mob') %>% unique()
            
            
            SRflag[time == tt]$SR = 'major'
            
            if(grepl('s.only', hyp.t)){
              cntSRwave2 = cntSRwave2 + .5 * length(SR.var.tx) / length(SR.var.full.t) # SR.perc.full.t / .1 * .2 # .2
              print(paste('cntSRwave2 = ', round(cntSRwave2, 2)))
            } else {
              cntSRwave2 = cntSRwave2 + 1 # SR.perc.full.t / .1
              print(paste('cntSRwave2 = ', round(cntSRwave2, 2)))
            }
            
            print(paste(tt, 'need major SR'), quote = F)
            
            # check if last time step also need major SR
            # first time need major update, no update on beta yet
            if(! update2widerSR1){
              SRbounds.full.t = cbind(SR.bounds[SR.var.full.t,1], SR.bounds.wider.t[SR.var.full.t,2])
            } else if (! update2widerSR2){ # 1st update done, 2nd not yet
              SRbounds.full.t = cbind(SR.bounds.wider.t[SR.var.full.t,1], SR.bounds.wider2.t[SR.var.full.t,1])
            } else {
              SRbounds.full.t = SR.bounds
            }
            
            
            if(p.sr.beta > 1.5 |
               dx.t.newitot[tt] > N * .1/100
               ){
              
              SRbounds.full.t['beta',1] = (SRbounds.full.t['beta',1] * (p.sr.beta %>% pmin(1.6))) %>% pmin(SRbounds.full.t['beta',2]*.6)
              SRbounds.full.t['beta',2] = SR.bounds.wider2.t['beta',2]  # in case it's after the first update but before the 2nd
              
            }
            
            
          } else if (abs(mn.pr - mn.po) / mn.pr > p.median.cut & mn.po > N * .1/100){
            
            if(grepl('s.only', hyp.t)){
              cntSRwave2 = cntSRwave2 + .5/2 * length(SR.var.tx) / length(SR.var.full.t) # .1
              print(paste('s.only, median SR, cntSRwave2 = ', round(cntSRwave2, 2)))
            } else {
              cntSRwave2 = cntSRwave2 + .5
              print(paste('beta, median SR, cntSRwave2 = ', round(cntSRwave2, 2)))
            }
            
            SR.perc.full.t = pmax(SR.perc.full * 2, .05) # not good for 2nd wave
            SR.perc.local.t = pmax(SR.perc.local * 2, .1)  # .1
            
            
            SR.var.full.t = c(SR.var.full,'p.mob','alpha')  %>% unique()  # ,'alpha', 'ifr'
            
            SRbounds.full.t = SR.bounds
            
            SRflag[time == tt]$SR = 'median'
            print(paste(tt, 'need additional SR'), quote = F)
            
          } else {
            
            SR.perc.full.t = SR.perc.full
            SR.perc.local.t = SR.perc.local
            
            SR.var.full.t = SR.var.full
            SRbounds.full.t = SR.bounds
            
          }
          
          
          # for S - determimine if probe on S is needed
          Sflag = F
          if(grepl('beta.only', hyp.t)){
            Sflag = F
          } else {
            
            # if (tt > pmax(end1stWave, main2ndWave - 6)) # don't allow major update of S until later?
            # not good 

            Sflag = F
            
            # cum.dS.perc = (Spost[tt] - Spost[end1stWave] + sum(Itotpost[end1stWave:tt])) / sum(Itotpost[1:end1stWave])
           
            
            # cum.dS.perc = cum.dS.po / (N - S2ndwave0)
            
            # for the cum.dS.po need to futher account for those lost to time
            limm = mean(state0['Trs',])
            cumIwave1loss = (cumIwave1 - cum.dS.po) * (1-exp(-(tt - end1stWave)*7 /limm))
            cum.dS.perc = (cum.dS.po + cumIwave1loss) / cumIwave1
            
            # only do it once
            if((cum.dS.perc > cum.dS.cut) & !stop.srS){  # tendency to over adjust s?
              # overshot - no good
              # cntSRwave2 = cntSRwave2 + abs(dS.mn) / cumIwave1 * 10
              stop.srS = T
              wk.srS.last = tail(stat.sr.S$week,1)
              
              if(F){
                if(length(which(stat.sr.S$mean < 0))>=2){
                  over.adjS = T
                  wk.srS.last = tail(stat.sr.S$week,1) + 2 # did not help, some further incr
                } else {
                  # over.adjS = T
                  wk.srS.last = tail(stat.sr.S$week,1)
                }
              }
              
              
              # extend it by 2 weeks if so
            }
            
            if(cum.dS.perc >= cum.dS.cut & ! grepl('beta.only', hyp.t)) {
              print(paste(tt, 'cum.dS >=cum.dS.cut'), quote = F)
            } else {
              # after 4/1/21
              p.small = 1.01; p.small.cb = 1.02
              p.median = 1.02; p.median.cb = 1.04; 
              p.large = 1.03; p.large.cb = 1.06;
              
              if(grepl('s.only', hyp.t) | 
                 hyp.t %in% c('both.major','both.minorb.majors') |
                 # update2widerSR2 |  # already updated to the highest beta, so it'd been on S more
                 update2widerSR1 |  # already updated to the highest beta, so it'd been on S more
                 grepl('slows', hyp.t)
                 ){
                Sdrift_small = (mean(xpost['S1',,tt-1]) > (mean(xpost['S1',,tt-2]) - mean(xpost['newItot1',,tt-1])) * p.small & 
                                  mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.small) |
                  (mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.small.cb) |
                  SRflag[time == tt]$SR == 'major'  # make it more senstive to the need to change 
                
                
                Sdrift_median = (mean(xpost['S1',,tt-1]) > (mean(xpost['S1',,tt-2]) - mean(xpost['newItot1',,tt-1])) * p.small & 
                                   mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.median) |
                  (mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.median.cb) |
                  (SRflag[time == tt-1]$SR %in% c('median', 'major')) & SRflag[time == tt]$SR == 'major'  # make it more senstive to the need to change 
                
                Sdrift_large = (mean(xpost['S1',,tt-1]) > (mean(xpost['S1',,tt-2]) - mean(xpost['newItot1',,tt-1])) * p.median & 
                                  mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.large) |
                  mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.large.cb |
                  (SRflag[time == tt-1]$SR %in% c('major')) & SRflag[time == tt]$SR == 'major'  # make it more senstive to the need to change
                
              } else {
                
                Sdrift_small = (mean(xpost['S1',,tt-1]) > (mean(xpost['S1',,tt-2]) - mean(xpost['newItot1',,tt-1])) * p.small & 
                                  mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.small) |
                  (mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.small.cb) |
                  (SRflag[time == tt-1]$SR %in% c('median', 'major')) & SRflag[time == tt]$SR == 'major'   # make it more senstive to the need to change 
                
                
                Sdrift_median = (mean(xpost['S1',,tt-1]) > (mean(xpost['S1',,tt-2]) - mean(xpost['newItot1',,tt-1])) * p.small & 
                                   mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.median) |
                  (mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.median.cb) |
                  (SRflag[time == tt-1]$SR %in% c('major')) & SRflag[time == tt]$SR == 'major'  # make it more senstive to the need to change 
                
                Sdrift_large = (mean(xpost['S1',,tt-1]) > (mean(xpost['S1',,tt-2]) - mean(xpost['newItot1',,tt-1])) * p.median & 
                                  mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.large) |
                  mean(xpost['S1',,tt]) > (mean(xpost['S1',,tt-1]) - mean(xpost['newItot1',,tt])) * p.large.cb |
                  (SRflag[time == tt-2]$SR %in% c('major','median')) & SRflag[time == tt]$SR == 'major' & (SRflag[time == tt-1]$SR %in% c('major')) & SRflag[time == tt]$SR == 'major'  # make it more senstive to the need to change
                
              }
              
              
              # can this be feed back to the level of changes?
              Sdrift.t = (mean(xpost['S1',,tt]) - mean(xpost['S1',,tt-1]) + mean(xpost['newItot1',,tt]))/mean(xpost['S1',,tt-1])
              pSdrift.large.t = Sdrift.t / .04 # compare to a 10% increase as baseline
              pSdrift.median.t = Sdrift.t / .02 # compare to a 5% increase as baseline
              
              # make sure it did not over shot
              Sdrift_small = Sdrift_small # & didnotovershot
              Sdrift_median = Sdrift_median # & didnotovershot
              Sdrift_large = Sdrift_large # & didnotovershot
              
              if(Sdrift_small |
                 Sdrift_median |
                 Sdrift_large # S drifting up
              ){
                if(! grepl('beta.only', hyp.t))
                  print('increase S helps')
                Sflag = T
              }
            }
            
          } #  determimine if probe on S is needed
          
          # determine bounds based on hyp
          if(grepl('s.only',hyp.t)){ #   & ((tt > end1stWave + 8 | tt >= wave2.1st.major) & tt < pmax(wave2.1st.major+15, num_times_i - 5))
            SR.perc.full.t = SR.perc.full
            SR.perc.local.t = SR.perc.local  # reduce level of local probing
            SR.var.full.t = SR.var.full.t[!(SR.var.full.t %in% SR.var.tx)] # exclude beta and Tir
            SR.var.local.t = SR.var.local.t [!(SR.var.local.t %in% SR.var.tx)] 
            SRbounds.full.t = SR.bounds.tight # restrict probing on beta/Tir, etc
            SRbounds.full.t['alpha',] = SR.bounds['alpha',]
          } else {
            SRbounds.full.t = SRbounds.full.t
          }
          
         
          
          # SR on S if needed
          if(Sflag & !grepl('beta.only',hyp.t) & ! descend &
             (SRflag[time == tt]$SR == 'major' | SRflag[time == tt]$SR == 'median') &
             cum.dS.perc < cum.dS.cut & cntSadjtot < cntSadjtot.cut & 
             (as.Date(vdate.t) < massvax.start + 7) & 
             (cntSadj_large < cntSadj_large.tot.t) &  
             cntSadj_median <cntSadj_median.tot.t & 
             cntSadj_small < cntSadj_small.tot.t){
            
            # record fist week needing SR on S?
            if(cntSadjtot==0) # first time
              S2ndwave0 = Spost[tt - 3] # back to 3wk ago before the likely surge
            # record the weeks needing SR on S
            wk.sr.S[tt] = 1
            
            
            # smp specifically those in the lowest half
            if(Sdrift_large){
              
              # Spb = (Spb_large.t - .01 * pmax(0, (tt - end1stWave - 15))) %>% pmax(Spb_large.t /3)
              # further adjust by the level needed to change for this case
              Spb = Spb_large.t  # restricted below with the remaining %S
              Spb.adj.t = pmin(Spb.adj_upr, pmax(Spb.adj_lwr, pSdrift.large.t)) # make the adjustment b/w .5 - 1.5
              Spb = Spb * Spb.adj.t
              
              print(paste0('prob S, large; Spb.adj.t=',Spb.adj.t %>% round(2)), quote = F)
              
              cntSadj_large = cntSadj_large + 1
              cntSadjtot = cntSadjtot + 1
              # cntSRwave2 = cntSRwave2 + Spb * 10 * p.S2beta
              # print(paste('SR on S, cntSRwave2 = ', round(cntSRwave2, 2)))
              
              if(cum.dS.perc < cum.dS.cut * p.cum.dS.large & Smajor >= Smajor.cut.lwr
                # cntSadj_large<=2 & # Smajor < Smajor.cut
                #  Smajor >= Smajor.cut.lwr & Smajor < Smajor.cut.upr & 
                #  Smajor.cnt < Smajor.cnt.cut
              ){
                Slow = which(state0['S1',] < quantile(state0['S1',], probs = .33 * Spb.adj.t))
                print(paste('Smajor=',Smajor, '; lowest tercile'), quote = F)
                Smajor = Smajor + Spb_large.t * 2
                Smajor.cnt = Smajor.cnt + 1
              } else {
                Slow = 1:num_ens
                print(paste('Smajor=',Smajor, '; random'), quote = F)
                Smajor = Smajor + Spb_large.t * 2
                
              }
              
              
            } else if (Sdrift_median){
              
              
              
              # Spb = (Spb_median.t - .005 * pmax(0, (tt - end1stWave - 15))) %>% pmax(Spb_median.t /3)
              # .25
              Spb = Spb_median.t  # restricted below with the remaining %S
              Spb.adj.t = pmin(Spb.adj_upr, pmax(Spb.adj_lwr, pSdrift.median.t)) # make the adjustment b/w .5 - 1.5
              Spb = Spb * Spb.adj.t
              print(paste0('prob S, median; Spb.adj.t=',Spb.adj.t %>% round(2)), quote = F)
              
              
              cntSadj_median = cntSadj_median + 1
              cntSadjtot = cntSadjtot + 1
              # cntSRwave2 = cntSRwave2 + Spb * 10 * p.S2beta
              # print(paste('SR on S, cntSRwave2 = ', round(cntSRwave2, 2)))
              
              if(cum.dS.perc < cum.dS.cut * p.cum.dS.median & Smajor >= Smajor.cut.lwr
                # cntSadj_median <=2  & # Smajor < Smajor.cut
                #  Smajor >= Smajor.cut.lwr & Smajor < Smajor.cut.upr & 
                #  Smajor.cnt < Smajor.cnt.cut
              ){
                Slow = which(state0['S1',] < quantile(state0['S1',], probs = .33))
                print(paste('Smajor=',Smajor, '; lowest tercile'), quote = F)
                Smajor = Smajor + Spb_median.t * 2
                Smajor.cnt = Smajor.cnt + .5
              } else {
                Slow = 1:num_ens
                print(paste('Smajor=',Smajor, '; random'), quote = F)
                Smajor = Smajor + Spb_median.t * 2
                
              }
              
            }  else if (Sdrift_small){
              
              # Spb = (Spb_small.t - .003 * pmax(0, (tt - end1stWave - 15))) %>% pmax(Spb_small.t /3)
              # .15
              Spb = Spb_small.t  # restricted below with the remaining %S
              cntSadj_small = cntSadj_small + 1
              cntSadjtot = cntSadjtot + .5
              
              # cntSRwave2 = cntSRwave2 + Spb * 10 * p.S2beta
              # print(paste('SR on S, cntSRwave2 = ', round(cntSRwave2, 2)))
              
              print('prob S, small')
              
              if(cum.dS.perc < cum.dS.cut * p.cum.dS.small & Smajor >= Smajor.cut.lwr
                # cntSadj_small <=2 & # Smajor < Smajor.cut
                #  Smajor >= Smajor.cut.lwr & Smajor < Smajor.cut.upr & 
                #  Smajor.cnt < Smajor.cnt.cut
              ){
                Slow = which(state0['S1',] < quantile(state0['S1',], probs = .5))
                print(paste('Smajor=',Smajor, '; lower half'), quote = F)
                Smajor = Smajor + Spb_small.t * 2
                Smajor.cnt = Smajor.cnt + .25
              } else {
                Slow = 1:num_ens
                print(paste('Smajor=',Smajor, '; random'), quote = F)
                Smajor = Smajor + Spb_small.t * 2
                
              }
            } 
            
            # Slow = which(state0['S1',] < quantile(state0['S1',], probs = Spb))
            # Sidx = Slow 
            # Slow = which(state0['S1',] < quantile(state0['S1',], probs = .5))
            # Slow = 1:num_ens
            Sidx =  sample(Slow, round(Spb * num_ens, 0), replace = F)
            # state0['S1', Sidx] = state0['S1', Sidx] + (N - state0['S1', Sidx]) * runif(length(Sidx), .75, 1)
            # s.lwr.t = (s.upr.t * .8 - .03 * pmax(0, (tt - end1stWave - 10))) %>% pmax(s.upr.t*.5)
            
            # 4/18/21 - later on, S may be too low and overshot!
            # state0['S1', Sidx] = state0['S1', Sidx] + (N - state0['S1', Sidx]) * runif(length(Sidx), s.lwr.t, s.upr.t)
            
            # 4/18/21
            # state0['S1', Sidx] = state0['S1', Sidx] + (N - state0.2ndWave.strt['S1', Sidx]) * runif(length(Sidx), s.lwr.t, s.upr.t)
            
            # 4/19/21 - adjust to the resisual immunity
            if(!is.null(res.imm)){
              tmp.mn = res.imm / ((N - state0['S1', Sidx]) %>% mean)  # res.imm
            } else {
              tmp.mn = 1
            }
            
            state0['S1', Sidx] = state0['S1', Sidx] + (N - state0['S1', Sidx]) * tmp.mn * runif(length(Sidx), s.lwr.t, s.upr.t)
            
            cntSRwave2 = cntSRwave2 + Spb * 10 * p.S2beta * tmp.mn
            print(paste('SR on S, cntSRwave2 = ', round(cntSRwave2, 2)))
            
          } # END PROBE S ONLY
          
          
          
          # UPDATE BETA BOUNDS?
          if(tt >= main2ndWave # only do this after the major 2nd wave
             ){
            if(cntUpdateSRbound < 2 & tt > (end1stWave+2)){  # add the 2 week lag
              
              
              # go back 2 weeks at most
              cntSR3wk = as.data.table(table(SRflag[time %in% pmax(main2ndWave-1, tt-c(0:npre.tot1))]$SR))
              cntSR5wk = as.data.table(table(SRflag[time %in% pmax(main2ndWave-2, tt-c(0:npre.tot2))]$SR))
              
              cntSRfull = data.table(V1 = c('minor','median','major'))
              cntSR3wk = merge(cntSRfull, cntSR3wk, all = T, by = 'V1')
              cntSR3wk[is.na(cntSR3wk)] = 0
              cntSR5wk = merge(cntSRfull, cntSR5wk, all = T, by = 'V1')
              cntSR5wk[is.na(cntSR5wk)] = 0
              
              
              # make sure it did not over shot
              
              if(!update2widerSR1 & # haven't updated beta bound yet
                 ((cntSR3wk[V1 == 'major']$N >= n.cut1 & cntSR5wk[V1 == 'major']$N < n.cut2) | # & didnotovershot
                 (! grepl('slow', hyp.t) & (dx.t.newitot[tt-1] > N * .1/100 & dx.t.newitot[tt] > N * .1/100))
                 ) # or when the filter needing to up-adjust newItot a lot
              ){  # 2 out of 3 weeks needing major SR
                # only do it once
                if(SR.bounds['beta', 1] < SR.bounds.wider.t['beta',1] | SR.bounds['beta', 2] < SR.bounds.wider.t['beta',2]){
                  SR.bounds[SR.var.tx,] = SR.bounds.wider.t[SR.var.tx,]
                  cntUpdateSRbound = cntUpdateSRbound + 1
                  update2widerSR1 = T
                  
                  print(paste0(tt, ':', n.cut1, 'times need major changes, update tx SR bounds 1'), quote = F)
                  
                  if(grepl('minorb', hyp.t) |
                     SR.bounds['beta', 2] == SR.bounds.wider2.t['beta',2]
                     ){
                    wk.update2widerSR2 = tt
                    update2widerSR2 = T
                    print('final update on beta')
                  } # if it's only 1 level, stop increase SR earlier
                    
                  
                  
                }
                
              } 
              
              if ((update2widerSR1 & !update2widerSR2) & 
                  (cntSR5wk[V1 == 'major']$N >= n.cut2 | #  & didnotovershot
                   (! grepl('slow', hyp.t) &  (dx.t.newitot[tt-1] > N * .2/100 & dx.t.newitot[tt] > N * .2/100))
                   )
              ){ # 4 out of 5 weeks needing major SR
                if(SR.bounds['beta', 1] < SR.bounds.wider2.t['beta',1] | SR.bounds['beta', 2] < SR.bounds.wider2.t['beta',2]){
                  SR.bounds[SR.var.tx,] = SR.bounds.wider2.t[SR.var.tx,]
                  cntUpdateSRbound = cntUpdateSRbound + 1
                  update2widerSR2 = T
                  wk.update2widerSR2 = tt
                  print(paste0(tt,':',n.cut2, 'times need major changes, update tx SR bounds 2'), quote = F)
                }
                
              }
              
            }
          }
            
        } # end 2nd wave
        
        num_SR.t = round(num_ens*(SR.perc.local.t + SR.perc.full.t), 0)
        num_SR.local.t = round(num_ens*SR.perc.local.t, 0)
        num_SR.full.t = round(num_ens*SR.perc.full.t, 0)
        SR.idx.t = sample(1:num_ens, num_SR.t, replace = F)
        SR.idx.local.t = head(SR.idx.t, num_SR.local.t)
        SR.idx.full.t = tail(SR.idx.t, num_SR.full.t)
        
        sr.local.mean = rowMeans(state0[SR.var.local,]) # local SR bounds
        sr.local.lwr = apply(state0[SR.var.local,], 1, quantile, .1) # .2
        sr.local.upr = apply(state0[SR.var.local,], 1, quantile, .9) # .7
        sr.local.lwr = pmin(sr.local.mean * .7, sr.local.lwr) # .75
        sr.local.upr = pmax(sr.local.mean * 1.3, sr.local.upr) # 1.25
        sr.local.lwr = pmax(sr.local.lwr, DA.bounds[SR.var.local,1]) # .75
        sr.local.upr = pmin(sr.local.upr, DA.bounds[SR.var.local,2]) # 1.25
        
        # if it's not summer in the UK
        
        if(tt %in% wk.summer & loc.t == 'uk'){ # low detection
          sr.local.lwr['alpha'] = pmin(sr.local.lwr['alpha'], SR.bounds['alpha',1])
          sr.local.upr['alpha'] = pmin(sr.local.upr['alpha'], SR.bounds['alpha',2])
          print('uk, summer!')
        } else {
          sr.local.lwr['alpha'] = pmax(sr.local.lwr['alpha'], SR.bounds['alpha',1])
          sr.local.upr['alpha'] = pmin(sr.local.upr['alpha'], SR.bounds['alpha',2])
        }
        
        
        # update ifr
        if(loc.t %in% c('uk','br','sa')){
          if(tt %in% (wk.WkLowIFR)){
          sr.local.lwr['ifr'] = pmin(sr.local.lwr['ifr'], SR.bounds['ifr',1])
          sr.local.upr['ifr'] = pmin(sr.local.upr['ifr'], SR.bounds['ifr',2])
          }
          
        }
        
        
        
        SR.bounds.local.t = cbind(sr.local.lwr, sr.local.upr)
        state0[SR.var.local.t,SR.idx.local.t] = t(lhs(num_SR.local.t, rect =  SR.bounds.local.t[SR.var.local.t,]))
        state0[SR.var.full.t,SR.idx.full.t] = t(lhs(num_SR.full.t, rect = SRbounds.full.t[SR.var.full.t,]))
        
        # addition sr on ifr for the uk enterting the summer
        if(loc.t == 'uk' & tt %in% wk.WkLowIFR[c(1,3)]){  # do it twice?
          idx = sample(1:num_ens, round(num_ens * .1, 0), replace = F)
          state0['ifr', idx] = runif(length(idx), min = SRifr_bounds.WkLowIFR[1], max = SRifr_bounds.WkLowIFR[2])
        }
        
        # additional SR on beta if needed
        if(extraSRbeta & (! grepl('s.only', hyp.t)) &
           tt>= main2ndWave & 
           # update2widerSR1 & # only do it after updated on beta
           tt <= wk.update2widerSR2  # ifelse(grepl('minor', hyp.t),2, 0)
           ){
          
          # samps = 1:num_ens; samps = samps[!(samps %in% SR.idx.full.t)]
          
          # make it more effecient, sample the lower half
          samps = which(state0['beta',] < quantile(state0['beta',], .33))
          # samps = which(state0['beta',] < quantile(state0['beta',], .5))
          samps = samps[!(samps %in% SR.idx.full.t)]
          
          SR.idx.extrabeta.t = sample(samps, round(num_ens * perc.extra.sr.beta, 0), replace = F)
          SR.idx.full.t = c(SR.idx.full.t, SR.idx.extrabeta.t)
          state0['beta',SR.idx.extrabeta.t] = runif(length(SR.idx.extrabeta.t), 
                                                    pmax(SR.bounds.wider.t['beta',1]*ifelse(ascend, 1.2, 1.1), SRbounds.full.t['beta',1]), 
                                                    pmax(SR.bounds.wider.t['beta',2]*ifelse(ascend, 1.1, 1.05), SRbounds.full.t['beta',2]))
          cntSRwave2 = cntSRwave2 + perc.extra.sr.beta / .1
          print(paste('extraSRbeta, cntSRwave2 = ', round(cntSRwave2, 2)))
          print(paste('extraSRbeta', perc.extra.sr.beta %>% round(2)))
        }
        
        # extra SR on beta
        # is the filter massively adjusting newItot up?
        if(dx.t.newitot[tt] > N * .15/100 & tt>= main2ndWave & (! grepl('s.only', hyp.t)) &
           tt <= pmin(wk.update2widerSR2, main2ndWave+8)  # ifelse(grepl('minor', hyp.t),2, 0)
           ){
          
          p.tmp.raw = (dx.t.newitot[tt]/N / (.15/100))
          p.tmp = p.tmp.raw %>% pmin(2)
          # make sure it's not the same ones probed
          # samps = 1:num_ens; samps = samps[!(samps %in% SR.idx.full.t)]
          # make it more effecient, sample the lower half
          samps = which(state0['beta',] < quantile(state0['beta',], .5))
          samps = samps[!(samps %in% SR.idx.full.t)]
          
          # SR.idx.t = sample(samps, round(num_ens * .1 * p.tmp, 0), replace = F)
          SR.idx.t = sample(samps, pmin(round(num_ens * .1 * p.tmp, 0), length(samps)), replace = F)
          
          
          if((! grepl('s.only', hyp.t))){
            state0['beta',SR.idx.t] = runif(round(num_ens * .1  * p.tmp, 0), 
                                            pmax(SR.bounds.wider.t['beta',1]*ifelse(ascend, 1.2, 1.1), SRbounds.full.t['beta',1]), 
                                            pmax(SR.bounds.wider.t['beta',2]*ifelse(ascend, 1.1, 1.05), SRbounds.full.t['beta',2])
                                            # pmax(SR.bounds.wider.t['beta',1]*1.2, SRbounds.full.t['beta',1]), 
                                            # pmax(SR.bounds.wider.t['beta',2]*1.05, SRbounds.full.t['beta',2])
                                            # SRbounds.full.t['beta',1], 
                                            # pmax(SR.bounds.wider2.t['beta',2]*.9, SRbounds.full.t['beta',2])
                                            )
            
          
            # also penalize fitting based on up-adjust newItot
            cntSRwave2 = cntSRwave2 + p.tmp # * p.tmp2 # may not be a good idea
            print(paste('extra SR on beta, cntSRwave2 = ', round(cntSRwave2, 2)))
            
          } 
          
          
          print(paste(tt, 'filter upadjusts newItot:', round(p.tmp, 2)))
        }
        
        # look at level of adjustment for both ups and downs, for all
        if(tt>= main2ndWave & mn.po > N * .1/100){
          
         
          # more penelty for up-adjustment
          p.up2down = 1.5; 
          pp = 2 # 5 too large
          cntSRwave2 = cntSRwave2 + (abs(dx.t.newitot[tt]) / pmin(mn.po, mn.pr, N * .3/100)  * ifelse(sign(dx.t.newitot[tt])>0, p.up2down, 1) + 
            abs(dx.t.newiobs[tt]) / pmin(obs_i[tt,1], mean(obs_i)*.25)  * ifelse(sign(dx.t.newiobs[tt])>0, p.up2down, 1) + 
            abs(dx.t.e[tt]) / pmin(e.po, e.pr, N * .15/100) * ifelse(sign(dx.t.e[tt])>0, p.up2down, 1) + 
            abs(dx.t.i[tt]) / pmin(i.po, i.pr, N * .25/100) * ifelse(sign(dx.t.i[tt])>0, p.up2down, 1)) * pp
          print(paste('dx, cntSRwave2 = ', round(cntSRwave2, 2)))
        }
        
        # extra SR on IFR
        # if(loc.t == 'uk' & tt %in% c(12:35)){
        #   SR.idx.local.t = sample(1:num_ens, round(num_ens * .15, 0), replace = F)
        #   state0['ifr',SR.idx.local.t] = runif(round(num_ens * .15, 0), sr.local.lwr['ifr'], sr.local.upr['ifr'])
        #   print(paste('SR on IFR', (c(sr.local.lwr['ifr'], sr.local.upr['ifr']) * 100) %>% round(2)))
        #   
        # }
        
      } # end doSR
      
      
      # get prior
      {
        if(!is.null(newI.previous) & vdate.t >= vax.start){
          tm.t = nrow(newI.previous)
          # cumI.t = apply(newI.previous00[,,1:(dim.t[3]-14),drop=F],c(1,2),sum) #  %>% apply(1, median) # excl last two weeks and get the median
          # t1 = (as.Date('2020/12/14') - as.Date('2020/3/1')) %>% as.numeric() # 1st day of vaccination
          # 2/5/21 set t1 to 1 yr given the slow rollout
          t1 = 365
          cumI.t = colSums(newI.previous[pmax(1,tm.t-t1) : (tm.t),]) #  %>% apply(1, median) # excl last two weeks and get the median
          # only count the last 12 months? so as the epidemic unfold, you don't over count cum infect?
          # higher infection rate for the priority group
          # tm.t = pmax(1, tm.t - t1 + 1) # re-aline timing with start of vac
          tm.imm = 365* 2.5 # assume 3 yr immunity
          p.imm.wane.max = .8; k = .015  # 1/tm.imm  
          p.imm.wane = 1 - p.imm.wane.max / (1+exp(-k*(tm.t + 60 - tm.imm/2))) # not all infected from day 1
          # earlier infections are likely to be in the high priority groups 
          p.imm = 1 *  p.imm.wane * redn.priority # assume 100% prior infection provide immunity, but wane over time
          # and multiple by % excluded if there is prior testing before vax: p.prior.test
          percSmax.t = 1 - cumI.t / N * p.imm
          # no lower than 50%, in case of outliers
          percSmax.t = pmax(percSmax.t, .5)
          # print(c('cohort %S:',round(summary(mean(percSmax.t)),2)), quote = F)
        } else {
          percSmax.t = 0
          # print('no vax yet')
        }
        
        beta_tt = state0['beta',] * fn_getRelMob(rel.mob[pmin(tt+1, nrow(rel.mob)),],state0['p.mob',]) 
        if(seasonality) {
          # beta_tt = state0['beta',] * seasonal.cycle[seasonal.cycle$week==cur.wk,]$relR0 # * rel.mob[tt]
          beta_tt = beta_tt * relR0[cur.wk,] # * rel.mob[tt]
        } 
        
        severity.t = severity
        severity.t['death',] = state0['ifr',]
        
        # [2/23/21] check this if problem!
        dist_tm.to.detect = NULL
        for(ii in 1:num_ens){
          tmp = generation.time(dist_tm.to.detect.name,c(state0['Td.mean',ii],state0['Td.sd',ii]),truncate = tm.to.detect.max)
          dist_tm.to.detect=cbind(dist_tm.to.detect,tmp$GT[-1]); 
        }
        
        if(epi.model == 'SEIRS'){
          
          simEpi=SEIRS(tm_strt, tm_end, tm_step=1, # 1 day time-step
                       tmstep = tmstep,
                       state0 = state0,
                       S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                       I0=state0[paste0('I',1:num_gr),], 
                       beta=beta_tt, 
                       Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                       seed=seed, stoch=stoch,
                       severity = severity.t,
                       newI.previous = newI.previous,
                       dist_tm.to.detect = dist_tm.to.detect,
                       dist_tm.to.death = dist_tm.to.death)
          
        } else if(epi.model == 'SEIRSV'){
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
          simEpi=SEIRSV(tm_strt, tm_end, tm_step=1, # 1 day time-step
                        tmstep = tmstep,
                        state0 = state0,
                        S0=state0[paste0('S',1:num_gr),], E0=state0[paste0('E',1:num_gr),], 
                        I0=state0[paste0('I',1:num_gr),], 
                        beta=beta_tt, 
                        Tei=state0['Tei',], Tir=state0['Tir',], Trs = state0['Trs',],
                        seed=seed, stoch=stoch, 
                        severity = severity.t,
                        newI.previous = newI.previous,
                        dist_tm.to.detect = dist_tm.to.detect,
                        dist_tm.to.death = dist_tm.to.death,
                        percSmax.t = percSmax.t,
                        V1 = V1.t, V2 = V2.t, # add vaccination for dose 1 and dose 2 -
                        # these are total number of vaccinees with unknown immunity
                        # but pre-ajust for time lag from vaccination to immune protection
                        VE1 = VE1, VE2=VE2 # Vaccine efficacy, need further adjustment by prior immunity 
          ) # include the delay reporting etc.
        }
        
        # re-assemble to the same order as the prior: state0
        n.end = tm_end - tm_strt + 2
        state.new = NULL
        for(i in 1:(length(simEpi)-1)){
          tmp = simEpi[[i]][n.end,,drop=F]; 
          rownames(tmp)=gsub('cumI','newI',paste0(names(simEpi)[i],1:num_gr))
          state.new = rbind(state.new,tmp)
        }
        
        state.new = rbind(state.new, state0[parm.names,])
        state.new = state.new[rownames(state0),] # make sure the order is the same
        
        
        xprior[,,tt+1]= state.new
        xprior.daily[1:tmstep+tt*tmstep, ]=simEpi$daily.newItot # we want the daily total new cases, without delay, without under-report
        
        # save prior for E and I
        Eprior[tt+1] = state.new['E1',] %>% mean
        Iprior[tt+1] = state.new['I1',] %>% mean
      } # end get prior
      
    } # end for-loop 
    
    # calculate the mean of ensemble
    xprior_mean=xpost_mean=matrix(0,num_times_i,num_var)
    for (tt in 1:num_times_i){
      xprior_mean[tt,]=apply(xprior[,,tt],1,mean, na.rm=T)
      xpost_mean[tt,]=apply(xpost[,,tt],1,mean, na.rm=T)
      
    }
    colnames(xprior_mean)=colnames(xpost_mean)=state.names
    xpost_mean = data.table(Week.start = Week.starts, xpost_mean)
    xprior_mean = data.table(Week.start = Week.starts, xprior_mean)
    rrmse = rbind(rrmse,data.table(hyp = hyp.t, n.parm.sr = n.parm.sr, 
                                   n.parm.sr.tot = n.parm.sr * length(which(SRflag[time >= main2ndWave]$SR == 'major')),
                                   cntSRmajor = length(which(SRflag[time >= main2ndWave]$SR == 'major')),
                                   cntSRwave2 = cntSRwave2, 
      E = sqrt(mean((xprior_mean$E1 - xpost_mean$E1)[end1stWave:num_times_i]^2)) / mean(xpost_mean$E1[end1stWave:num_times_i]),
      I = sqrt(mean((xprior_mean$I1 - xpost_mean$I1)[end1stWave:num_times_i]^2)) / mean(xpost_mean$I1[end1stWave:num_times_i]),
      newItot = sqrt(mean((xprior_mean$newItot1 - xpost_mean$newItot1)[end1stWave:num_times_i]^2)) / mean(xpost_mean$newItot1[end1stWave:num_times_i]),
      newIobs = sqrt(mean((xprior_mean$newIobs1 - xpost_mean$newIobs1)[end1stWave:num_times_i]^2)) / mean(xpost_mean$newIobs1[end1stWave:num_times_i]),
      newD = sqrt(mean((xprior_mean$death1 - xpost_mean$death1)[end1stWave:num_times_i]^2)) / mean(xpost_mean$death1[end1stWave:num_times_i]),
      obs.pr = sqrt(mean((xprior_mean$newIobs1 - obs_i)[end1stWave:num_times_i]^2)) / mean(obs_i[end1stWave:num_times_i]),
      death.pr = sqrt(mean((xprior_mean$death1 - obs_d)[end1stWave:num_times_i]^2)) / mean(obs_d[end1stWave:num_times_i]),
      obs.po = sqrt(mean((xpost_mean$newIobs1 - obs_i)[end1stWave:num_times_i]^2)) / mean(obs_i[end1stWave:num_times_i]),
      death.po = sqrt(mean((xpost_mean$death1 - obs_d)[end1stWave:num_times_i]^2)) / mean(obs_d[end1stWave:num_times_i])
    ))
    
    # save results
    {
      eval(parse(text = paste0('xpost',ihyp,'=xpost', sep='')))
      eval(parse(text = paste0('xprior',ihyp,'=xprior', sep='')))
      eval(parse(text = paste0('mpost',ihyp,'=xpost_mean', sep='')))
      eval(parse(text = paste0('mprior',ihyp,'=xprior_mean', sep='')))
      eval(parse(text = paste0('stat.sr.S',ihyp,'=stat.sr.S', sep='')))
      # eval(parse(text = paste0('cum.dS.po',ihyp,'=cum.dS.po', sep='')))
      # eval(parse(text = paste0('cum.dS.perc',ihyp,'=cum.dS.perc', sep='')))
      
      if(!is.null(stat.sr.S)){
        dS.tot.mn = sum(stat.sr.S$mean)
        # did it wait a long time to do sr?
        # go back 3 weeks ago get a baseline before the surge due to the newV
        wk.t = pmax(stat.sr.S$week[1] - 3, end1stWave+3) # %>% round(0)
        
        S0.2nd.strt = xpost['S1',,wk.t]
        S0.2nd.mn = S0.2nd.strt %>% mean; 
        S0.2nd.var = S0.2nd.strt %>% var
        # dS.ens = colSums(stat.sr.S[,3+1:num_ens,with=F])  # probably not a good way
        Imm.ens = N - S0.2nd.strt; Imm.mn = Imm.ens %>% mean; 
        Imm.var = Imm.ens %>% var; Imm.sd = Imm.ens %>% sd
        
        # compute the combined SD for dS.tot.mn
        if(nrow(stat.sr.S) > 1){
          xx.mean = numeric(nrow(stat.sr.S))
          formula.t = 'x1'
          for(i in 1:nrow(stat.sr.S)){
            xx.t = stat.sr.S[i,3+1:num_ens,with=F] %>% unlist
            xx.t.mn =xx.t %>% mean; 
            # xx.t.var = var(xx.t); 
            xx.mean[i] = xx.t.mn
            # eval(parse(text=paste('xx',i,'=xx.t.mn',sep='')))
            if(i > 1)
              formula.t = paste0(formula.t,'+x',i)
          }
          
          m.cov = cov(x=stat.sr.S[,3+1:num_ens,with=F]%>% t)
          eval(parse(text = paste('tmp = deltamethod(~', formula.t, ', mean = xx.mean, cov = m.cov)'))) 
          dS.tot.var = tmp; dS.tot.sd = dS.tot.var %>% sqrt
        } else {
          dS.tot.sd = stat.sr.S$sd; dS.tot.var = dS.tot.sd^2
        }
        cov.hat = dS.tot.sd * Imm.sd
        m.cov = matrix(c(dS.tot.var, cov.hat, cov.hat, Imm.var),2,2)
        
        # m.cov = cov(x=cbind(dS.ens, Imm.ens)) # probably not a good way
        perc.dImm.mn = dS.tot.mn / Imm.mn * 100
        perc.dImm.sd = (deltamethod(~x1/x2, mean = c(dS.tot.mn, Imm.mn), cov = m.cov) %>% sqrt) * 100
        
        # for changes in tx
        
        
        newVstat = data.table(cum.dS = cum.dS.po, cum.Imm = N - S2ndwave0, 
                              perc.dImm.mn = perc.dImm.mn, perc.dImm.sd = perc.dImm.sd)
      } else {
        newVstat = data.table(cum.dS = cum.dS.po, cum.Imm = N - S2ndwave0, 
                              perc.dImm.mn = 0, perc.dImm.sd = 0)
      }
      
      
      
      eval(parse(text = paste0('newVstat',ihyp,'=newVstat', sep='')))
      
      # eval(parse(text = paste0('cntSR',ihyp,'=cntSRwave2', sep='')))
      
      
    }
    
  } # diff hypothesis
  
  # for SR, relative to the max
  rrmse.raw = rrmse
  
  rrmse$cntSRwave2scaled = rrmse$cntSRwave2 / max(rrmse$cntSRwave2)
  rrmse$n.parm.sr.scaled = rrmse$n.parm.sr / max(rrmse$n.parm.sr)
  rrmse$n.parm.sr.tot.scaled = rrmse$n.parm.sr.tot / max(rrmse$n.parm.sr.tot)
  
  vars.tt = c('n.parm.sr.tot.scaled', 'n.parm.sr.scaled', 'cntSRwave2scaled','E','I','newItot','newIobs','obs.pr','death.pr','obs.po','death.po')
  
  # not good
  # rrmse = rrmse[,lapply(.SD, FUN = function(x){x/max(x)}), .SDcols = vars.tt]
  # rrmse$hyp = rrmse.raw$hyp
  wt.nparm = 1
  wt.cntSRwave2 = 1 # may be too strong
  # wt.cntSRwave2 = .5
  # rrmse$cb = rowMeans(rrmse[,c('E','I','newItot','newIobs'),with=F])
  # with probing on beta only, the changes in E and I tend to be milder and ending up on the top
  # so perhaps a weighted measure would be better
  # rrmse$cb = rowMeans(rrmse[,c('newItot','newIobs','obs'),with=F])
  # rrmse0 = rrmse
  fn_wt.rrmse = function(hyp.t, vars.t, wts.t, rrmse.t){
    sum(rrmse.t[hyp == hyp.t,vars.t,with=F] * wts.t)
  }
  
  # vars.tt = c('E','I','newItot','newIobs','newD','obs','death') # this double count the observations
  # wts.tt = c(rep(1,2),rep(2,5)); wts.tt = wts.tt/sum(wts.tt)
  # vars.tt = c('cntSRwave2','E','I','newItot','newIobs','obs.pr','death.pr','obs.po','death.po')
  wts.tt.eq = c(wt.nparm, wt.nparm, wt.cntSRwave2, rep(1,3),rep(1,5)); wts.tt.eq = wts.tt.eq/sum(wts.tt.eq)
  
  # wts.tt.obs.more = c(wt.cntSRwave2, rep(2,2), rep(2,2),rep(3,2), rep(4,2)); wts.tt.obs.more = wts.tt.obs.more/sum(wts.tt.obs.more)
  # wts.tt.obs.more = c(wt.cntSRwave2, rep(.5,2), rep(1,2),rep(2,2), rep(4,2)); wts.tt.obs.more = wts.tt.obs.more/sum(wts.tt.obs.more)
  pmax.others = mean(rrmse[,c('newItot','newIobs','obs.pr','death.pr','obs.po','death.po'),with=F] %>% unlist)
  wts.tt.obs.more = c(0, wt.nparm, wt.cntSRwave2, rep(0,2), rep(.5,2),rep(2,2), rep(4,2));  wts.tt.obs.more = wts.tt.obs.more/sum(wts.tt.obs.more)
  names(wts.tt.obs.more) = vars.tt
  # less for death?
  # wts.tt.obs.more = c(rep(1,2), rep(2,2),c(3,1.5), c(4,2)); wts.tt.obs.more = wts.tt.obs.more/sum(wts.tt.obs.more)
  
  # wts.tt.obs.most = c(rep(.5,2), rep(2,2),rep(3,2), rep(4,2)); wts.tt.obs.most = wts.tt.obs.most/sum(wts.tt.obs.most) # not good
  wts.tt.obs.most = c(wt.nparm, 0, wt.cntSRwave2, rep(0,2), rep(.5,2),rep(2,2), rep(4,2));  wts.tt.obs.most = wts.tt.obs.most/sum(wts.tt.obs.most) # tend to be better?
  wts.tt.obs.most2 = c(wt.nparm, 0, wt.cntSRwave2 * .5, rep(0,2), rep(.5,2),rep(2,2), rep(4,2));  wts.tt.obs.most2 = wts.tt.obs.most2/sum(wts.tt.obs.most2) # tend to be better?
  
  # somewhere in bw?
  # wts.tt.obs.comb = c(wt.nparm * .75, wt.nparm * .75, wt.cntSRwave2, rep(0,2), rep(.5,2),rep(2,2), rep(4,2));  wts.tt.obs.comb = wts.tt.obs.comb/sum(wts.tt.obs.comb) # tend to be better?
  wts.tt.obs.comb = c(wt.nparm * 2/3, wt.nparm * 2/3, wt.cntSRwave2 * 2/3, rep(0,2), rep(.5,2),rep(2,2), rep(4,2));  wts.tt.obs.comb = wts.tt.obs.comb/sum(wts.tt.obs.comb) # tend to be better?
  
  wts.tt.obs.only = c(wt.nparm, wt.cntSRwave2, rep(0,2), rep(0,2),rep(0,2), rep(1,2)); wts.tt.obs.only = wts.tt.obs.only/sum(wts.tt.obs.only)
  tmp.eq = rrmse[,list(cb.eq=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.eq, rrmse.t=rrmse)),
                by = 'hyp']
  tmp.obs.more = rrmse[,list(cb.obs.more=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.obs.more, rrmse.t=rrmse)),
                 by = 'hyp']
  tmp.obs.most = rrmse[,list(cb.obs.most=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.obs.most, rrmse.t=rrmse)),
                       by = 'hyp']
  tmp.obs.most2 = rrmse[,list(cb.obs.most2=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.obs.most2, rrmse.t=rrmse)),
                       by = 'hyp']
  
  tmp.obs.comb = rrmse[,list(cb.obs.comb=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.obs.comb, rrmse.t=rrmse)),
                       by = 'hyp']
  tmp.obs.only = rrmse[,list(cb.obs.only=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.obs.only, rrmse.t=rrmse)),
                       by = 'hyp']
  # rank the rrmse
  ranks = rrmse[,lapply(.SD, rank), .SDcols = vars.tt]
  ranks$hyp = rrmse$hyp
  wts.tt.rank = c(wt.nparm,wt.nparm, wt.cntSRwave2, rep(1,2), rep(1,2),rep(1,2), rep(1,2)); wts.tt.rank = wts.tt.rank/sum(wts.tt.rank)
  tmp.rank = ranks[,list(cb.rank=fn_wt.rrmse(hyp.t = hyp, vars.t = vars.tt, wts.t = wts.tt.rank, rrmse.t=ranks)),
                       by = 'hyp']
  
  rrmse = merge(rrmse, tmp.eq, by = 'hyp')
  rrmse = merge(rrmse, tmp.obs.more, by = 'hyp')
  rrmse = merge(rrmse, tmp.obs.most, by = 'hyp')
  rrmse = merge(rrmse, tmp.obs.most2, by = 'hyp')
  rrmse = merge(rrmse, tmp.obs.comb, by = 'hyp')
  rrmse = merge(rrmse, tmp.obs.only, by = 'hyp')
  rrmse = merge(rrmse, tmp.rank, by = 'hyp')
  rrmse$hyp = factor(rrmse$hyp, levels = hyps, labels = hyps)
  rrmse = rrmse[order(hyp)]
  
  
  # i.best = which.min(rrmse$cb.obs.only)
  
  # output for all eval method
  Rt_stats_all = R0_stats_all = Rtx_stats_all = Rtx_ens_all = states_stats_all = xpost_mean_all = 
    xpost_sd_all = xprior_mean_all = xprior_sd_all = cumIperc_stats_all = cumIperc_ens_all = 
    Susceptibility_stats_all = Susceptibility_ens_all = xpost.last_all = rrmse_all = NULL
  hyp.best_all = NULL; newVstat_all = NULL
  
  for(tag.eval in tag.evals){
    
    i.best = eval(parse(text = paste('which.min(rrmse$cb.', tag.eval,')', sep='')))
    
    hyp.best = rrmse$hyp[i.best]
    print(paste0('best for ',loc.t,': ',hyp.best), quote = F)
    eval(parse(text = paste0('xpost=xpost',i.best, sep='')))
    eval(parse(text = paste0('xprior=xprior',i.best, sep='')))
    eval(parse(text = paste0('newVstat=newVstat',i.best, sep='')))
    # newVstat
    
    hyp.best_all = rbind(hyp.best_all, data.table(eval = tag.eval, hyp.test = hyp.best))
    
    # calculate the mean of ensemble
    xprior_mean=xpost_mean=xprior_sd=xpost_sd=matrix(0,num_times_i,num_var)
    for (tt in 1:num_times_i){
      xprior_mean[tt,]=apply(xprior[,,tt],1,mean, na.rm=T)
      xprior_sd[tt,]=apply(xprior[,,tt],1,sd, na.rm=T)
      xpost_mean[tt,]=apply(xpost[,,tt],1,mean, na.rm=T)
      xpost_sd[tt,]=apply(xpost[,,tt],1,sd, na.rm=T)
      
    }
    colnames(xprior_mean)=colnames(xpost_mean)=colnames(xpost_sd)=colnames(xprior_sd)=state.names
    xpost_mean = data.table(Week.start = Week.starts, xpost_mean)
    xprior_mean = data.table(Week.start = Week.starts, xprior_mean)
    
    # get stats instead
    dimnames(xpost)[1] = list(state.names)
    states_stats = NULL
    for(var in state.names){
      tmp = xpost[var,,]
      tmp = tmp %>% apply(2, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
      colnames(tmp) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
      
      tmp.mean = xpost[var,,] %>% apply(2, mean)
      states_stats = rbind(states_stats, 
                           data.table(Week.start = Week.starts, state = var, mean = tmp.mean, tmp))
    }
    
    
    # compute Rt
    Rt_ens = matrix(0, num_times_i, num_ens)
    R0_ens = matrix(0, num_times_i, num_ens)
    Rtx_ens = matrix(0, num_times_i, num_ens)
    for(tt in 1:num_times_i){
      
      cur.wk = weeks[tt]
      
      for(ii in 1:num_ens){
        
        beta.mean = xpost['beta',ii,tt] * fn_getRelMob(rel.mob[tt,ii],xpost['p.mob',ii,tt]) 
        # add seasonality? 
        if(seasonality) {
          # beta.mean =beta.mean * seasonal.cycle[seasonal.cycle$week==weeks[tt],]$relR0
          beta.mean = beta.mean * relR0[cur.wk,ii] # * rel.mob[tt]
        }
        PARMS = list(beta.mean = beta.mean,  Tir.mean= xpost['Tir',ii,tt], S = xpost['S1',ii,tt], N = N)
        Rt_ens[tt, ii] = Fn_getRt_SEIR(PARMS)
        R0_ens[tt, ii] = Fn_getR0_SEIR(PARMS)
        Rtx_ens[tt, ii] = xpost['beta',ii,tt] * xpost['Tir',ii,tt]
      } # ens
    }
    
    Rt_stats = Rt_ens %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
    colnames(Rt_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
    
    Rt_stats = data.table(Week.start = Week.starts, 
                          mean = Rt_ens %>% apply(1, mean),
                          sd = Rt_ens %>% apply(1, sd),
                          Rt_stats)
    
    R0_stats = R0_ens %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
    colnames(R0_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
    
    R0_stats = data.table(Week.start = Week.starts, 
                          mean = R0_ens %>% apply(1, mean),
                          sd = R0_ens %>% apply(1, sd),
                          R0_stats)
    
    Rtx_stats = Rtx_ens %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
    colnames(Rtx_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
    
    Rtx_stats = data.table(Week.start = Week.starts, 
                           mean = Rtx_ens %>% apply(1, mean),
                           sd = Rtx_ens %>% apply(1, sd),
                           Rtx_stats)
    
    xpost.last=xpost[,,tt] # save all the ens members
    row.names(xpost.last)=state.names
    
    # get the change in tx
    idxW1main = 5: (end1stWave)  # exclude 1st and last 3 weeks
    # identify the week it get to a high level during 2nd wave
    if (hyp.t == 'both.minorb.slows'){
      i2 = which.max(Rtx_stats[main2ndWave+0:pmin(npre.tot1, num_times_i - main2ndWave-8)]$mean) + main2ndWave -1  # give it 2 months to see
    } else if (hyp.t == 'both.majorb.slows') {
      i2 = which.max(Rtx_stats[main2ndWave+0:pmin(npre.tot2-3, num_times_i - main2ndWave-8)]$mean) + main2ndWave -1  # give it 2 months to see
    } else {
      i2 = which.max(Rtx_stats[main2ndWave+0:9]$mean) + main2ndWave -1  # give it 2 months to see
    }
    
    
    idxW2main = i2 : num_times_i # (main2ndWave + 5) : num_times_i
    
    if(loc.t == 'uk' & excludeLockDown){
      # exclude the week with lockdown
      idxW2main = idxW2main[!(idxW2main %in% wkLockDown)]
    }
    
    Rtx1 = Rtx_stats[idxW1main, c("mean", 'sd', "median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr"), with=F] %>% colMeans()
    Rtx2 = Rtx_stats[idxW2main, c("mean", 'sd',"median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr"), with=F] %>% colMeans()
    
    tmp1 =  Rtx_ens[idxW1main,];
    tmp2 = Rtx_ens[idxW2main,];
    formula.t = paste0('x',1:nrow(tmp1), collapse = '+')
    mean.t = rowMeans(tmp1); cov.t = cov(t(tmp1))
    tmp1sd = eval(parse(text=paste('deltamethod(~', formula.t, ', mean= mean.t, cov=cov.t)', sep=''))) %>% sqrt
    formula.t = paste0('x',1:nrow(tmp2), collapse = '+')
    mean.t = rowMeans(tmp2); cov.t = cov(t(tmp2))
    tmp2sd = eval(parse(text=paste('deltamethod(~', formula.t, ', mean= mean.t, cov=cov.t)', sep=''))) %>% sqrt
    cov.t = matrix(c(tmp1sd^2, tmp1sd*tmp2sd,tmp1sd*tmp2sd, tmp2sd^2),2,2)
    dRtx.sd = (deltamethod(~(x2 - x1)/x1, mean = c(Rtx1['mean'], Rtx2['mean']) %>% unname, cov = cov.t) %>% sqrt) * 100
    
    dRtx = data.table(perc.dRtx.mean = (Rtx2['mean'] - Rtx1['mean'])/Rtx1['mean'] * 100, 
                      perc.dRtx.median = (Rtx2['median'] - Rtx1['median'])/Rtx1['median'] * 100,
                      perc.dRtx.sd = dRtx.sd)
    newVstat = cbind(newVstat, dRtx)
    # eval(parse(text = paste0('newVstat',i.best,'=newVstat', sep='')))
    
    if(T){ # don't need this
      # cumulative infection
      cumIperc = xpost['newItot1',,] %>% t %>% apply(2, cumsum)
      cumIperc = cumIperc / N * 100 # %
      cumIperc_stats = cumIperc %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
      colnames(cumIperc_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
      cumIperc_stats = data.table(Week.start = Week.starts, 
                                  mean = cumIperc %>% apply(1, mean),
                                  sd = cumIperc %>% apply(1, sd),
                                  cumIperc_stats)
      # for susceptible
      Susceptibility = xpost['S1',,] %>% t 
      Susceptibility = Susceptibility / N * 100 # %
      Susceptibility_stats = Susceptibility %>% apply(1, quantile, prob = c(.5, .25, .75, .025, .975)) %>% t
      colnames(Susceptibility_stats) = c('median', 'iqr.lwr','iqr.upr','ci95.lwr','ci95.upr')
      Susceptibility_stats = data.table(Week.start = Week.starts, 
                                        mean = Susceptibility %>% apply(1, mean),
                                        sd = Susceptibility %>% apply(1, sd),
                                        Susceptibility_stats)
    }
    
    newVstat_all = rbind(newVstat_all,data.table(eval = tag.eval, newVstat))
    Rt_stats_all = rbind(Rt_stats_all,data.table(eval = tag.eval, Rt_stats))
    R0_stats_all = rbind(R0_stats_all, data.table(eval = tag.eval, R0_stats))
    Rtx_stats_all = rbind(Rtx_stats_all, data.table(eval = tag.eval, Rtx_stats))
    Rtx_ens_all = rbind(Rtx_ens_all, data.table(eval = tag.eval, Rtx_ens))
    states_stats_all = rbind(states_stats_all, data.table(eval = tag.eval, states_stats))
    xpost_mean_all = rbind(xpost_mean_all, data.table(eval = tag.eval, xpost_mean))
    xpost_sd_all = rbind(xpost_sd_all, data.table(eval = tag.eval, xpost_sd))
    xprior_mean_all = rbind(xprior_mean_all, data.table(eval = tag.eval, xprior_mean))
    xprior_sd_all = rbind(xprior_sd_all, data.table(eval = tag.eval, xprior_sd))
    cumIperc_stats_all = rbind(cumIperc_stats_all, data.table(eval = tag.eval, cumIperc_stats))
    cumIperc_ens_all = rbind(cumIperc_ens_all, data.table(eval = tag.eval, cumIperc))
    Susceptibility_stats_all = rbind(Susceptibility_stats_all, data.table(eval = tag.eval, Susceptibility_stats))
    Susceptibility_ens_all = rbind(Susceptibility_ens_all, data.table(eval = tag.eval, Susceptibility))
    xpost.last_all = rbind(xpost.last_all, data.table(eval = tag.eval, xpost.last))
    # rrmse_all = rbind(rrmse_all, data.table(eval = tag.eval, rrmse))
    
    
  } # end eval
  

  
  return(list(Rt_stats = Rt_stats_all, R0_stats = R0_stats_all, 
              Rtx_stats = Rtx_stats_all, 
              Rtx_ens = Rtx_ens_all,
              states_stats = states_stats_all, 
              xpost_mean = xpost_mean_all, xpost_sd = xpost_sd_all, 
              xprior_mean = xprior_mean_all, xprior_sd = xprior_sd_all,
              cumIperc_stats = cumIperc_stats_all,
              cumIperc_ens = cumIperc_ens_all,
              Susceptibility_stats = Susceptibility_stats_all,
              Susceptibility_ens = Susceptibility_ens_all,
              xpost.last=xpost.last_all,
              newVstat = newVstat_all,
              rrmse = rrmse,
              hyp.best_all = hyp.best_all))
  
}
