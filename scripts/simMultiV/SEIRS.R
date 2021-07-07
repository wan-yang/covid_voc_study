obs.delayedSEIR = function(newI_ts, dist_tm.to.detect){
  # newI_ts: time series of new cases
  # distribution of time to detection
  n.days = nrow(newI_ts); n.ens = ncol(newI_ts)
  tm.max = tm.to.detect.max # longest delay
  est.obs = matrix(0,n.days+tm.max, n.ens);
  for(i in 1:n.days){
    tmp = dist_tm.to.detect %*% diag(newI_ts[i,]) # distribute the number to the possible onset dates
    est.obs[i-1+1:tm.max,] = est.obs[i-1+1:tm.max,] + tmp; # cummulative
  }
  est.obs
}

# for any delay in general
outcome.delayedSEIR = function(newI_ts, dist_tm.to.outcome, tm.to.outcome.max){
  # newI_ts: time series of new cases
  # distribution of time to detection
  n.days = nrow(newI_ts); n.ens = ncol(newI_ts)
  tm.max = tm.to.outcome.max # longest delay
  est.outcome = matrix(0,n.days+tm.max, n.ens);
  for(i in 1:n.days){
    tmp = dist_tm.to.outcome %*% diag(newI_ts[i,]) # distribute the number to the possible onset dates
    est.outcome[i-1+1:tm.max,] = est.outcome[i-1+1:tm.max,] + tmp; # cummulative
  }
  est.outcome
}

outcomeInsOuts.delayedSEIR = function(newI_ts, dist_tm.to.outcomeIns, 
                                      dist_tm.to.outcomeOuts,
                                      tm.to.outcomeIns.max, 
                                      tm.to.outcomeOuts.max){
  # newI_ts: time series of new cases
  # distribution of time to detection
  # also include when those will leave based on retention time
  n.days = nrow(newI_ts); n.ens = ncol(newI_ts)
  tmIns.max = tm.to.outcomeIns.max # longest delay
  tmOuts.max = tm.to.outcomeOuts.max
  est.outcomeIns = matrix(0,n.days+tmIns.max, n.ens);
  est.outcomeOuts = matrix(0,n.days+tmOuts.max, n.ens);
  est.outcomePrev = matrix(0,n.days+tmOuts.max, n.ens);
  for(i in 1:n.days){
    tmpIns = dist_tm.to.outcomeIns %*% diag(newI_ts[i,]) # distribute the number to the possible onset dates
    # shift the timing for outs, but keep the numbers
    est.outcomeIns[i-1+1:tmIns.max,] = est.outcomeIns[i-1+1:tmIns.max,] + tmpIns; # cummulative
    for(ii in 1:n.ens){
      est.outcomeOuts[i-1+1:tmIns.max+dist_tm.to.outcomeOuts[ii],ii]= est.outcomeOuts[i-1+1:tmIns.max+dist_tm.to.outcomeOuts[ii],ii] + tmpIns[,ii]
    }
  }
  # compute net - prevalance
  tmpIns = rbind(est.outcomeIns,matrix(0,tmOuts.max-tmIns.max,n.ens))
  cumIns = apply(tmpIns, 2, cumsum)
  cumOuts = apply(est.outcomeOuts, 2, cumsum)
  est.outcomePrev = cumIns - cumOuts
  
  # plot(est.outcomePrev[,2]); lines(est.outcomeIns[,2],col='red')
  
  return(list(est.outcomeIns=est.outcomeIns, est.outcomeOuts=est.outcomeOuts, est.outcomePrev=est.outcomePrev))
}



SEIRS = function(tm_strt, tm_end, tm_step=1, # 1 day time-step
                 tmstep = tmstep,
                 state0 = state0,
                 S0, E0, I0, # H0, # N0,
                 beta, Tei, Tir, Trs,
                 seed=0, stoch=T,
                 severity = severity,
                 newI.previous = newI.previous,
                 dist_tm.to.detect = dist_tm.to.detect,
                 dist_tm.to.death = dist_tm.to.death,
                 birth.rate = 0 # set to 0 here as
                 # in the context of the campus, birth.rate is like recruitment of new students/staff
){
  
  newI.previous00 = newI.previous;
  
  # RK 4 step algorithm, with random grab from Poisson dist
  cnt=1;
  # beta stores only data during the time used for the truth
  
  tm_vec=seq(tm_strt, tm_end, by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  np=length(S0); # 
  # if integrating parallelly, S0, and I0 should be passed in as a vector (Np particles)
  # but if so, can't do multiple age groups
  S=E=I=cumItot=matrix(0,tm_sz,np)
  S[1,]=S0; E[1,]=E0; I[1,]=I0; 
  cumIobs=matrix(0,tm_sz,np); cumIobs[1,]=0; 
  cumItot[1,]=0; # N[1,]=N0; 
  
  for(tt in tm_vec){
    cnt=cnt+1
    
    # step 1
    mu.infect = tm_step * (beta * I[cnt-1,]) * S[cnt-1,] / N  # new infection
    mu.ei = tm_step * E[cnt-1,] / Tei # E->I
    mu.ir = tm_step * I[cnt-1,] / Tir # I -> R community transmission
    mu.rs = tm_step * (N - S[cnt-1,] - E[cnt-1,] - I[cnt-1,]) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, S[cnt-1,]) # new cases < S
    mu.ei = pmin(mu.ei, E[cnt-1,])
    mu.ir = pmin(mu.ir, I[cnt-1,])
    mu.rs = pmin(mu.rs, (N - S[cnt-1,] - E[cnt-1,] - I[cnt-1,]))
    
    sk1 = - mu.infect + mu.rs # - seed 
    ek1 = mu.infect - mu.ei # + seed
    ik1 = mu.ei - mu.ir 
    ik1t = mu.ei
    ik1o = mu.ir
    
    Ts1 = S[cnt-1,] + round(sk1/2)
    Te1 = E[cnt-1,] + round(ek1/2)
    Ti1 = I[cnt-1,] + round(ik1/2)
    
    # step 2
    mu.infect = tm_step * (beta * Ti1) * Ts1 / N  # new infection
    mu.ei = tm_step * Te1 / Tei # E->I
    mu.ir = tm_step * Ti1 / Tir # I -> R community transmission
    mu.rs = tm_step * (N - Ts1 - Te1 - Ti1) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, Ts1) # new cases < S
    mu.ei = pmin(mu.ei, Te1)
    mu.ir = pmin(mu.ir, Ti1)
    mu.rs = pmin(mu.rs, (N - Ts1 - Te1 - Ti1))
    
    sk2 = - mu.infect + mu.rs # - seed 
    ek2 = mu.infect - mu.ei # + seed
    ik2 = mu.ei - mu.ir 
    ik2t = mu.ei
    ik2o = mu.ir
    
    Ts2 = S[cnt-1,] + round(sk2/2)
    Te2 = E[cnt-1,] + round(ek2/2)
    Ti2 = I[cnt-1,] + round(ik2/2)
    
    # step 3
    mu.infect = tm_step * (beta * Ti2) * Ts2 / N  # new infection
    mu.ei = tm_step * Te2 / Tei # E->I
    mu.ir = tm_step * Ti2 / Tir # I -> R community transmission
    mu.rs = tm_step * (N - Ts2 - Te2 - Ti2) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, Ts2) # new cases < S
    mu.ei = pmin(mu.ei, Te2)
    mu.ir = pmin(mu.ir, Ti2)
    mu.rs = pmin(mu.rs, (N - Ts2 - Te2 - Ti2))
    
    sk3 = - mu.infect + mu.rs # - seed 
    ek3 = mu.infect - mu.ei # + seed
    ik3 = mu.ei - mu.ir 
    ik3t = mu.ei
    ik3o = mu.ir
    
    Ts3 = S[cnt-1,] + round(sk3)
    Te3 = E[cnt-1,] + round(ek3)
    Ti3 = I[cnt-1,] + round(ik3)
    
    # step 4
    mu.infect = tm_step * (beta * Ti3) * Ts3 / N  # new infection
    mu.ei = tm_step * Te3 / Tei # E->I
    mu.ir = tm_step * Ti3 / Tir # I -> R community transmission
    mu.rs = tm_step * (N - Ts3 - Te3 - Ti3) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, Ts3) # new cases < S
    mu.ei = pmin(mu.ei, Te3)
    mu.ir = pmin(mu.ir, Ti3)
    mu.rs = pmin(mu.rs, (N - Ts3 - Te3 - Ti3))
    
    sk4 = - mu.infect + mu.rs # - seed 
    ek4 = mu.infect - mu.ei # + seed
    ik4 = mu.ei - mu.ir 
    ik4t = mu.ei
    ik4o = mu.ir
    
    # seed = seed # seeding
    # seed = runif(np, seed*p.home.lwr, seed*p.home.upr) # city baseline from outside, p.home.lwr - p.home.upr from home 
    seed = seed # * p.home 
    
    S[cnt,]=S[cnt-1,]+round(sk1/6+sk2/3+sk3/3+sk4/6,0)-seed + birth.rate * (N - S[cnt-1,])
    E[cnt,]=E[cnt-1,]+round(ek1/6+ek2/3+ek3/3+ek4/6,0)+seed - birth.rate * E[cnt-1,]
    I[cnt,]=I[cnt-1,]+round(ik1/6+ik2/3+ik3/3+ik4/6,0) - birth.rate * I[cnt-1,]
    cumItot[cnt,]=cumItot[cnt-1,]+round(ik1t/6+ik2t/3+ik3t/3+ik4t/6,0)
    cumIobs[cnt,]=cumIobs[cnt-1,]+round(ik1o/6+ik2o/3+ik3o/3+ik4o/6,0)
  } # end time step
  
  # include the delayed reporting etc.
  # account for delay in case reporting and death/reporting
  death=cumIobs # place holder
  {
    if(! is.null(newI.previous00)){
      newI.previous00.t = newI.previous00
    } else {
      newI.previous00.t = NULL
    }
    
    
    newI_ts = cumItot # these are cummulative cases
    newI_ts = newI_ts[-1,,drop=F] - newI_ts[-nrow(newI_ts),,drop=F] # total cases, without delay or under-reporting
    
    # need to include previous cases (not yet detected as well)
    newI.previous = tail(newI.previous00.t,tm.to.detect.max) 
    newI.combined = rbind(newI.previous, newI_ts)
    
    
    # est.obs = obs.delayedSEIR(newI_ts, dist_tm.to.detect)
    est.daily.tot = obs.delayedSEIR(newI.combined, dist_tm.to.detect)
    idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
    est.daily.tot = est.daily.tot[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
    
    # now account for report rate
    est.obs.daily = est.daily.tot * matrix(state0[grep('alpha',rownames(state0)),],tmstep,num_ens,byrow = T)
    
    # death:
    newI.previous = tail(newI.previous00.t,tm.from.inf.to.death.max)
    newI.combined = rbind(newI.previous, newI_ts)
    est.daily.death = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.death, tm.to.outcome.max=tm.from.inf.to.death.max)
    idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
    est.daily.death =est.daily.death[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
    # now account for severity
    est.daily.death = est.daily.death * matrix(severity[grep('death',rownames(severity)),],tmstep,num_ens,byrow = T)
    
    # back to cumIobs:
    est.obs.this=rbind(0, apply(est.obs.daily,2,cumsum))
    
    # cummulative deaths
    est.death.this = rbind(0, apply(est.daily.death,2,cumsum))
    
    cumIobs = est.obs.this # re-assign observed
    
    death = est.death.this # re-assign observed
    
    
  }
  
  # return
  return(list(S=S,E=E,I=I,cumIobs=cumIobs,cumItot=cumItot, death = death, daily.newItot = newI_ts))
}

# include vaccination
SEIRSV = function(tm_strt, tm_end, tm_step=1, # 1 day time-step
                 tmstep = tmstep,
                 state0 = state0,
                 S0, E0, I0, # H0, # N0,
                 beta, Tei, Tir, Trs,
                 seed=0, stoch=T,
                 severity = severity,
                 newI.previous = newI.previous,
                 dist_tm.to.detect = dist_tm.to.detect,
                 dist_tm.to.death = dist_tm.to.death,
                 birth.rate = 0, # set to 0 here as
                 # in the context of the campus, birth.rate is like recruitment of new students/staff
                 # for vaccination
                 percSmax.t = percSmax.t,
                 V1, V2, # add vaccination for dose 1 and dose 2 -
                 # these are total number of vaccinees with unknown immunity
                 # but pre-ajust for time lag from vaccination to immune protection
                 VE1, VE2 # Vaccine efficacy, need further adjustment by prior immunity 
){
  
  newI.previous00 = newI.previous;
  
  # RK 4 step algorithm, with random grab from Poisson dist
  cnt=1;
  # beta stores only data during the time used for the truth
  
  tm_vec=seq(tm_strt, tm_end, by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  np=length(S0); # 
  # if integrating parallelly, S0, and I0 should be passed in as a vector (Np particles)
  # but if so, can't do multiple age groups
  S=E=I=cumItot=matrix(0,tm_sz,np)
  S[1,]=S0; E[1,]=E0; I[1,]=I0; 
  cumIobs=matrix(0,tm_sz,np); cumIobs[1,]=0; 
  cumItot[1,]=0; # N[1,]=N0; 
  
  for(tt in tm_vec){
    cnt=cnt+1
    
    # step 1
    mu.infect = tm_step * (beta * I[cnt-1,]) * S[cnt-1,] / N  # new infection
    mu.ei = tm_step * E[cnt-1,] / Tei # E->I
    mu.ir = tm_step * I[cnt-1,] / Tir # I -> R community transmission
    mu.rs = tm_step * (N - S[cnt-1,] - E[cnt-1,] - I[cnt-1,]) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, S[cnt-1,]) # new cases < S
    mu.ei = pmin(mu.ei, E[cnt-1,])
    mu.ir = pmin(mu.ir, I[cnt-1,])
    mu.rs = pmin(mu.rs, (N - S[cnt-1,] - E[cnt-1,] - I[cnt-1,]))
    
    sk1 = - mu.infect + mu.rs # - seed 
    ek1 = mu.infect - mu.ei # + seed
    ik1 = mu.ei - mu.ir 
    ik1t = mu.ei
    ik1o = mu.ir
    
    Ts1 = S[cnt-1,] + round(sk1/2)
    Te1 = E[cnt-1,] + round(ek1/2)
    Ti1 = I[cnt-1,] + round(ik1/2)
    
    # step 2
    mu.infect = tm_step * (beta * Ti1) * Ts1 / N  # new infection
    mu.ei = tm_step * Te1 / Tei # E->I
    mu.ir = tm_step * Ti1 / Tir # I -> R community transmission
    mu.rs = tm_step * (N - Ts1 - Te1 - Ti1) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, Ts1) # new cases < S
    mu.ei = pmin(mu.ei, Te1)
    mu.ir = pmin(mu.ir, Ti1)
    mu.rs = pmin(mu.rs, (N - Ts1 - Te1 - Ti1))
    
    sk2 = - mu.infect + mu.rs # - seed 
    ek2 = mu.infect - mu.ei # + seed
    ik2 = mu.ei - mu.ir 
    ik2t = mu.ei
    ik2o = mu.ir
    
    Ts2 = S[cnt-1,] + round(sk2/2)
    Te2 = E[cnt-1,] + round(ek2/2)
    Ti2 = I[cnt-1,] + round(ik2/2)
    
    # step 3
    mu.infect = tm_step * (beta * Ti2) * Ts2 / N  # new infection
    mu.ei = tm_step * Te2 / Tei # E->I
    mu.ir = tm_step * Ti2 / Tir # I -> R community transmission
    mu.rs = tm_step * (N - Ts2 - Te2 - Ti2) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, Ts2) # new cases < S
    mu.ei = pmin(mu.ei, Te2)
    mu.ir = pmin(mu.ir, Ti2)
    mu.rs = pmin(mu.rs, (N - Ts2 - Te2 - Ti2))
    
    sk3 = - mu.infect + mu.rs # - seed 
    ek3 = mu.infect - mu.ei # + seed
    ik3 = mu.ei - mu.ir 
    ik3t = mu.ei
    ik3o = mu.ir
    
    Ts3 = S[cnt-1,] + round(sk3)
    Te3 = E[cnt-1,] + round(ek3)
    Ti3 = I[cnt-1,] + round(ik3)
    
    # step 4
    mu.infect = tm_step * (beta * Ti3) * Ts3 / N  # new infection
    mu.ei = tm_step * Te3 / Tei # E->I
    mu.ir = tm_step * Ti3 / Tir # I -> R community transmission
    mu.rs = tm_step * (N - Ts3 - Te3 - Ti3) / Trs
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = rpois(np, mu.infect)
      mu.ei = rpois(np,mu.ei)
      mu.ir = rpois(np,mu.ir)
      mu.rs = rpois(np,mu.rs)
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, Ts3) # new cases < S
    mu.ei = pmin(mu.ei, Te3)
    mu.ir = pmin(mu.ir, Ti3)
    mu.rs = pmin(mu.rs, (N - Ts3 - Te3 - Ti3))
    
    sk4 = - mu.infect + mu.rs # - seed 
    ek4 = mu.infect - mu.ei # + seed
    ik4 = mu.ei - mu.ir 
    ik4t = mu.ei
    ik4o = mu.ir
    
    # add vaccination: adjust for ve and prior infection 
    # but later on, with vaccinees making up for the majority of immune, 
    # this would no longer reflective of the ture prob of vaccinee's prior immunity
    # should put a lower bound for suscept/upper bound for prior immunity
    # percS = pmax(S[,,cnt-1] / Npops, matrix(percSmax, Ngr, Np))
    # use the cumulative infection rate instead, 
    percS = percSmax.t # matrix: num_gr x num_ens
    vacc.d1 = percS * V1[cnt-1,] *  VE1  # number ppl vaccinated after first dose of vaccine
    # dose 2: should we account for prior immunity as well? if so it should be ~1 month ago
    # for simplicity, use the same percS
    # adjust for additional VE for the 2nd dose
    VE2acc = 1-(1-VE2)/(1-VE1)
    vacc.d2 = percS * V2[cnt-1,] * (1 - VE1) * VE2acc 
    
    # seed = seed # seeding
    # seed = runif(np, seed*p.home.lwr, seed*p.home.upr) # city baseline from outside, p.home.lwr - p.home.upr from home 
    seed = seed # * p.home 
    
    S[cnt,]=S[cnt-1,]+round(sk1/6+sk2/3+sk3/3+sk4/6,0)-seed + birth.rate * (N - S[cnt-1,]) - vacc.d1 - vacc.d2
    E[cnt,]=E[cnt-1,]+round(ek1/6+ek2/3+ek3/3+ek4/6,0)+seed - birth.rate * E[cnt-1,]
    I[cnt,]=I[cnt-1,]+round(ik1/6+ik2/3+ik3/3+ik4/6,0) - birth.rate * I[cnt-1,]
    cumItot[cnt,]=cumItot[cnt-1,]+round(ik1t/6+ik2t/3+ik3t/3+ik4t/6,0)
    cumIobs[cnt,]=cumIobs[cnt-1,]+round(ik1o/6+ik2o/3+ik3o/3+ik4o/6,0)
  } # end time step
  
  # include the delayed reporting etc.
  # account for delay in case reporting and death/reporting
  death=cumIobs # place holder
  {
    if(! is.null(newI.previous00)){
      newI.previous00.t = newI.previous00
    } else {
      newI.previous00.t = NULL
    }
    
    
    newI_ts = cumItot # these are cummulative cases
    newI_ts = newI_ts[-1,,drop=F] - newI_ts[-nrow(newI_ts),,drop=F] # total cases, without delay or under-reporting
    
    # need to include previous cases (not yet detected as well)
    newI.previous = tail(newI.previous00.t,tm.to.detect.max) 
    newI.combined = rbind(newI.previous, newI_ts)
    
    
    # est.obs = obs.delayedSEIR(newI_ts, dist_tm.to.detect)
    est.daily.tot = obs.delayedSEIR(newI.combined, dist_tm.to.detect)
    idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
    est.daily.tot = est.daily.tot[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
    
    # now account for report rate
    est.obs.daily = est.daily.tot * matrix(state0[grep('alpha',rownames(state0)),],tmstep,num_ens,byrow = T)
    
    # death:
    newI.previous = tail(newI.previous00.t,tm.from.inf.to.death.max)
    newI.combined = rbind(newI.previous, newI_ts)
    est.daily.death = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.death, tm.to.outcome.max=tm.from.inf.to.death.max)
    idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
    est.daily.death =est.daily.death[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
    # now account for severity
    est.daily.death = est.daily.death * matrix(severity[grep('death',rownames(severity)),],tmstep,num_ens,byrow = T)
    
    # back to cumIobs:
    est.obs.this=rbind(0, apply(est.obs.daily,2,cumsum))
    
    # cummulative deaths
    est.death.this = rbind(0, apply(est.daily.death,2,cumsum))
    
    cumIobs = est.obs.this # re-assign observed
    
    death = est.death.this # re-assign observed
    
    
  }
  
  # return
  return(list(S=S,E=E,I=I,cumIobs=cumIobs,cumItot=cumItot, death = death, daily.newItot = newI_ts))
}
