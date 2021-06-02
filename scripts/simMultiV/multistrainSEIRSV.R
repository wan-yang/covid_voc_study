# Epidemiology models
# model: multi-strain SIRS

## cross immunity: slighly diff than the Gog model, sirs model fast version
multistrainSEIRSV <-function(tm_strt, tm_end, tm_step, S0, E0, I0, N, 
                             Tei, Tir, Trs, beta, cross, 
                             Iexp = 1, # exponent for I to reduce the exponential growth rate
                            birthrate = 0, # set to 0 here as
                            # in the context of the campus, birth.rate is like recruitment of new students/staff
                            severity = severity,
                            newI.previous = newI.previous,
                            dist_tm.to.detect = dist_tm.to.detect,
                            dist_tm.to.death = dist_tm.to.death,
                            # for vaccination
                            percSmax.t = percSmax.t,
                            V1, V2, # add vaccination for dose 1 and dose 2 -
                            # these are total number of vaccinees with unknown immunity
                            # but pre-ajust for time lag from vaccination to immune protection
                            VE1, VE2, # Vaccine efficacy, need further adjustment by prior immunity 
                            VE1redn, VE2redn,
                            seed,
                            stoch = T
                            ){
  # function to integrate to the next time step
  # use multistrain SIRS model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         Tir: infection period, day; matrix(num_obs,Np)
  #         L: immune period, day; matrix(num_obs,Np)
  #         alpha: rate from exposed to infectious; 
  #         beta: transmission matrix at time t; matrix(num_obs,Np)
  #         cross: cross immunity matrix, c11=c22=c33=0; 
  #         c12:infected by strain 2, provide cross-protection to strain 1
  #         NOTE: C12 may be equal to C21, 0<cij<1
  #         CROSS=cross: matrix[diff cij, Np]
  #         birthrate: per day per population
  #
  # output: S, I for all time steps
  
  seed0 = seed;
  newI.previous00 = newI.previous;
  
  cnt=1;

  tm_vec=seq(tm_strt, tm_end, by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  
  Ns=dim(I0)[1]; # number of strains
  Np=dim(I0)[2]; # number of ensemble members
  # if integrating parallelly, S0, and I0 should be passed in as a vector (Np particles)
  
  # S=E=I=newI=array(0,c(Ns,Np,tm_sz))
  # S[,,1]=S0; E[,,1]=E0;  I[,,1]=I0; # R[,1]=N-S0-I0;
  # newI[,,1]=0;
  
  S=E=I=cumIobs=cumItot=array(0,c(Ns,Np,tm_sz)); # matrix(0,tm_sz,np)
  S[,,1]=S0; E[,,1]=E0; I[,,1]=I0; 
  cumIobs[,,1]=0; cumItot[,,1]=0; # N[,,1]=N0; 
  # also track cummulative cross-imm 
  
  cross.newI=array(0,c(Ns,Np,tm_sz))
  
  # cross immunity: slighly diff than the Gog model
  # Si= -beta_i*Si*Ii/N - sum(beta_i * cimm_i_j * Sj * Ij/N) 
  # note the infection term by other strain here is Sj*Ij/N (cp Si*Ij/N in Gog & Grenfell)
  
  cross.sums=rep('0',Ns);
  for(si in 1:Ns){  # compute cross imm for each strain
    for(sj in 1:Ns){
      if(si==sj) next;
      # cross.sums[si]=paste0(cross.sums[si],"+cross['cimm",si,'_',sj,"',]*Einf[",sj,',]')
      cross.sums[si]=paste0(cross.sums[si],"+cross['cimm",si,'_',sj,"',]*mu.infect[",sj,',]')
    }
  }
  
  for (t in 1:length(tm_vec)){
    cnt=cnt+1;
    
    # step 1
    N.cur=N; 
    
    # step 1
    S.cur = S[,,cnt-1]; E.cur = E[,,cnt-1]; I.cur = I[,,cnt-1]; 
    R.cur = matrix(N.cur, Ns, Np) - S.cur - E.cur - I.cur
    # I.cur = (I.cur^Iexp) %>% apply(2, round, 0) # to account for imperfect mixing
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    
    mu.infect=tm_step*(beta*I.cur*S.cur/N)
    # cross immunity
    mu.icross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      mu.icross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / matrix(Tei,Ns,Np,byrow = T) # E->I
    mu.ir = tm_step * I.cur / matrix(Tir,Ns,Np,byrow = T) # I -> R community transmission
    mu.rs = tm_step * R.cur / matrix(Trs,Ns,Np,byrow = T) # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*Ns,mu.infect),Ns,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*Ns,mu.ei),Ns,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*Ns,mu.ir),Ns,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*Ns,mu.rs),Ns,Np);
      mu.icross=matrix(rpois(Np*Ns,mu.icross),Ns,Np);
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)

    
    # update the transmission part
    
    sk1 = - mu.infect + mu.rs - mu.icross # - seed 
    ek1 = mu.infect - mu.ei # + seed
    ik1 = mu.ei - mu.ir 
    ik1t = mu.ei
    ik1o = mu.ir
    
    Ts1 = S[,,cnt-1] + apply(sk1/2,2,round,0)
    Te1 = E[,,cnt-1] + apply(ek1/2,2,round,0)
    Ti1 = I[,,cnt-1] + apply(ik1/2,2,round,0)
    Ts1[Ts1<0]=0; Te1[Te1<0]=0; Ti1[Ti1<0]=0; 
    
    # step 2
    S.cur = Ts1; E.cur = Te1; I.cur = Ti1; 
    R.cur = matrix(N.cur, Ns, Np) - S.cur - E.cur - I.cur
    # I.cur = (I.cur^Iexp) %>% apply(2, round, 0) # to account for imperfect mixing
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    
    mu.infect=tm_step*(beta*I.cur*S.cur/N)
    # cross immunity
    mu.icross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      mu.icross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / matrix(Tei,Ns,Np,byrow = T) # E->I
    mu.ir = tm_step * I.cur / matrix(Tir,Ns,Np,byrow = T) # I -> R community transmission
    mu.rs = tm_step * R.cur / matrix(Trs,Ns,Np,byrow = T) # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*Ns,mu.infect),Ns,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*Ns,mu.ei),Ns,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*Ns,mu.ir),Ns,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*Ns,mu.rs),Ns,Np);
      mu.icross=matrix(rpois(Np*Ns,mu.icross),Ns,Np);
    }
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    # update the transmission part
    
    sk2 = - mu.infect + mu.rs - mu.icross # - seed 
    ek2 = mu.infect - mu.ei # + seed
    ik2 = mu.ei - mu.ir 
    ik2t = mu.ei
    ik2o = mu.ir
    
    Ts2 = S[,,cnt-1] + apply(sk2/2,2,round,0)
    Te2 = E[,,cnt-1] + apply(ek2/2,2,round,0)
    Ti2 = I[,,cnt-1] + apply(ik2/2,2,round,0)
    Ts2[Ts2<0]=0; Te2[Te2<0]=0; Ti2[Ti2<0]=0; 
    
    # step 3
    S.cur = Ts2; E.cur = Te2; I.cur = Ti2; 
    R.cur = matrix(N.cur, Ns, Np) - S.cur - E.cur - I.cur
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    
    mu.infect=tm_step*(beta*I.cur*S.cur/N)
    # cross immunity
    mu.icross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      mu.icross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / matrix(Tei,Ns,Np,byrow = T) # E->I
    mu.ir = tm_step * I.cur / matrix(Tir,Ns,Np,byrow = T) # I -> R community transmission
    mu.rs = tm_step * R.cur / matrix(Trs,Ns,Np,byrow = T) # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*Ns,mu.infect),Ns,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*Ns,mu.ei),Ns,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*Ns,mu.ir),Ns,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*Ns,mu.rs),Ns,Np);
      mu.icross=matrix(rpois(Np*Ns,mu.icross),Ns,Np);
    }
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    # update the transmission part
    
    sk3 = - mu.infect + mu.rs - mu.icross # - seed 
    ek3 = mu.infect - mu.ei # + seed
    ik3 = mu.ei - mu.ir 
    ik3t = mu.ei
    ik3o = mu.ir
    
    Ts3 = S[,,cnt-1] + apply(sk3,2,round,0)
    Te3 = E[,,cnt-1] + apply(ek3,2,round,0)
    Ti3 = I[,,cnt-1] + apply(ik3,2,round,0)
    Ts3[Ts3<0]=0; Te3[Te3<0]=0; Ti3[Ti3<0]=0; 
    
    # step 4
    S.cur = Ts3; E.cur = Te3; I.cur = Ti3; 
    R.cur = matrix(N.cur, Ns, Np) - S.cur - E.cur - I.cur
    # I.cur = I.cur^Iexp %>% apply(2, round, 0) # to account for imperfect mixing
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    
    mu.infect=tm_step*(beta*I.cur*S.cur/N)
    # cross immunity
    mu.icross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      mu.icross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / matrix(Tei,Ns,Np,byrow = T) # E->I
    mu.ir = tm_step * I.cur / matrix(Tir,Ns,Np,byrow = T) # I -> R community transmission
    mu.rs = tm_step * R.cur / matrix(Trs,Ns,Np,byrow = T) # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*Ns,mu.infect),Ns,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*Ns,mu.ei),Ns,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*Ns,mu.ir),Ns,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*Ns,mu.rs),Ns,Np);
      mu.icross=matrix(rpois(Np*Ns,mu.icross),Ns,Np);
    }
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    # update the transmission part
    
    sk4 = - mu.infect + mu.rs - mu.icross # - seed 
    ek4 = mu.infect - mu.ei # + seed
    ik4 = mu.ei - mu.ir 
    ik4t = mu.ei
    ik4o = mu.ir
    
    
    # add vaccination
    # add vaccination: adjust for ve and prior infection 
    # but later on, with vaccinees making up for the majority of immune, 
    # this would no longer reflective of the ture prob of vaccinee's prior immunity
    # should put a lower bound for suscept/upper bound for prior immunity
    # percS = pmax(S[,,cnt-1] / Npops, matrix(percSmax, Ns, Np))
    # use the cumulative infection rate instead, 
    vacc.d1 = vacc.d2 = matrix(0, Ns, Np)
    
    for(iv in 1:Ns){  # different effectiveness for diff variants
      percS = percSmax.t[iv,] # matrix: Ns x Np
      vacc.d1[iv,] = percS * V1[cnt-1,] *  VE1 * VE1redn[iv] # number ppl vaccinated after first dose of vaccine
      # dose 2: should we account for prior immunity as well? if so it should be ~1 month ago
      # for simplicity, use the same percS
      vacc.d2[iv,] = percS * V2[cnt-1,] * (1 - VE1) * VE2  * VE2redn[iv]
    }
    
    
    # seed = seed # 
    
    seed = matrix(rpois(Np*Ns, seed0),Ns,Np);
    
    S[,,cnt]=S[,,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + tm_step*birthrate*(N-S[,,cnt-1]) - vacc.d1 - vacc.d2 # add births/deaths
    E[,,cnt] = E[,,cnt-1] + apply(ek1/6+ek2/3+ek3/3+ek4/6,2,round,0) + seed
    I[,,cnt] = I[,,cnt-1] + apply(ik1/6+ik2/3+ik3/3+ik4/6,2,round,0)
    cumItot[,,cnt] = cumItot[,,cnt-1] + apply(ik1t/6+ik2t/3+ik3t/3+ik4t/6,2,round,0)
    cumIobs[,,cnt] = cumIobs[,,cnt-1] + apply(ik1o/6+ik2o/3+ik3o/3+ik4o/6,2,round,0)
    
    S[S<0]=0; E[E<0]=0; I[I<0]=0; cumItot[cumItot<0]=0; cumIobs[cumIobs<0]=0; 

  }
  
  # include the delayed reporting etc.
  # account for delay in case reporting and death/reporting
  death=cumIobs # place holder
  {
    
    # for all groups
    est.daily.tot_all= array(0, c(Ns, Np, tmstep )) # this account for delay in reporting, but not under-reporing
    # not the newIobs in state.new account for both delay in reporting and under-reporing
    est.daily.death_all = array(0, c(Ns, Np, tmstep ))
    newI_ts_all = array(0, c(Ns, Np, tmstep ))
    
    # NEED TO DO IT FOR EACH GROUP
    for(ig in 1:Ns){
      
      if(! is.null(newI.previous00)){
        newI.previous00.t = t(newI.previous00[ig,,])
      } else {
        newI.previous00.t = NULL
      }
      
      
      newI_ts = t(cumItot[ig,,]) # these are cummulative cases
      newI_ts = newI_ts[-1,,drop=F] - newI_ts[-nrow(newI_ts),,drop=F] # total cases, without delay or under-reporting
      
      # need to include previous cases (not yet detected as well)
      newI.previous = tail(newI.previous00.t,tm.to.detect.max) 
      newI.combined = rbind(newI.previous, newI_ts)
      # est.obs = obs.delayedSEIR(newI_ts, dist_tm.to.detect)
      est.daily.tot = obs.delayedSEIR(newI.combined, dist_tm.to.detect)
      idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
      est.daily.tot = est.daily.tot[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
      
      # now account for report rate
      est.obs.daily = est.daily.tot * matrix(parm0[grep('alpha',rownames(parm0)),],tmstep,Np,byrow = T)
      
      # death:
      newI.previous = tail(newI.previous00.t,tm.from.inf.to.death.max)
      newI.combined = rbind(newI.previous, newI_ts)
      est.daily.death = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.death, tm.to.outcome.max=tm.from.inf.to.death.max)
      idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
      est.daily.death =est.daily.death[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
      # now account for severity
      est.daily.death = est.daily.death * matrix(severity[grep('death',rownames(severity)),],tmstep,Np,byrow = T)
      
      # back to cumIobs:
      est.obs.this=rbind(0, apply(est.obs.daily,2,cumsum))
      
      # cummulative deaths
      est.death.this = rbind(0, apply(est.daily.death,2,cumsum))
      
      cumIobs[ig,,] = t(est.obs.this) # re-assign observed
      
      death[ig,,] = t(est.death.this) # re-assign observed
      
      # save for this group
      est.daily.tot_all[ig,,] = t(est.daily.tot)  # this account for delay in reporting, but not under-reporing
      newI_ts_all[ig,,] = t(newI_ts)
      
      
    } # group
    
  } # get case and death, with delay and under reporting
  
  # return
  return(list(S=S,E=E,I=I,cumIobs=cumIobs,cumItot=cumItot, death = death, daily.newItot = newI_ts_all))
}

# add age structure
multistrainSEIRSVage <-function(tm_strt, tm_end, tm_step, S0, E0, I0, Nage, 
                             Tei, Tir, Trs, beta, cross, 
                             Iexp = 1, # exponent for I to reduce the exponential growth rate
                             birthrate = 0, # set to 0 here as
                             # in the context of the campus, birth.rate is like recruitment of new students/staff
                             # severity = severity,
                             newI.previous = newI.previous,
                             # dist_tm.to.detect = dist_tm.to.detect,
                             # dist_tm.to.death = dist_tm.to.death,
                             # for vaccination
                             percSmax.t = percSmax.t,
                             V1, V2, # add vaccination for dose 1 and dose 2 -
                             # these are total number of vaccinees with unknown immunity
                             # but pre-ajust for time lag from vaccination to immune protection
                             VE1, VE2, # Vaccine efficacy, need further adjustment by prior immunity 
                             VE1redn, VE2redn,
                             seed,
                             stoch = T
){
  # function to integrate to the next time step
  # use multistrain SIRS model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         Tir: infection period, day; matrix(num_obs,Np)
  #         L: immune period, day; matrix(num_obs,Np)
  #         alpha: rate from exposed to infectious; 
  #         beta: transmission matrix at time t; matrix(num_obs,Np)
  #         cross: cross immunity matrix, c11=c22=c33=0; 
  #         c12:infected by strain 2, provide cross-protection to strain 1
  #         NOTE: C12 may be equal to C21, 0<cij<1
  #         CROSS=cross: matrix[diff cij, Np]
  #         birthrate: per day per population
  #
  # output: S, I for all time steps
  
  seed0 = seed
  newI.previous00 = newI.previous;
  
  cnt=1;
  
  tm_vec=seq(tm_strt, tm_end, by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  
  Np=dim(I0)[2]; # number of ensemble members
  NaNs=dim(I0)[1]; # number of strains x number of age group
  
  
  # if integrating parallelly, S0, and I0 should be passed in as a vector (Np particles)
  
  # S=E=I=newI=array(0,c(Ns,Np,tm_sz))
  # S[,,1]=S0; E[,,1]=E0;  I[,,1]=I0; # R[,1]=N-S0-I0;
  # newI[,,1]=0;
  N.nans = NULL
  for(iv in 1: Ns) {
    N.nans = c(N.nans, Nage)
  }
  
  S=E=I=cumIobs=cumItot=array(0,c(NaNs,Np,tm_sz)); # matrix(0,tm_sz,np)
  S[,,1]=S0; E[,,1]=E0; I[,,1]=I0; 
  cumIobs[,,1]=0; cumItot[,,1]=0; # N[,,1]=N0; 
  # also track cummulative cross-imm 
  
  
  # cross immunity: slighly diff than the Gog model
  # Si= -beta_i*Si*Ii/N - sum(beta_i * cimm_i_j * Sj * Ij/N) 
  # note the infection term by other strain here is Sj*Ij/N (cp Si*Ij/N in Gog & Grenfell)
  
  cross.sums.t=rep('0',Ns);
  for(si in 1:Ns){  # compute cross imm for each strain
    for(sj in 1:Ns){
      if(si==sj) next;
      # cross.sums[si]=paste0(cross.sums[si],"+cross['cimm",si,'_',sj,"',]*Einf[",sj,',]')
      cross.sums.t[si]=paste0(cross.sums.t[si],"+cross['cimm",letters[si],'_',letters[sj],"',]*mu.infect.t[",sj,',]')
    }
  }
  
  
  # get the command assembling the bits for FOI: I_i/N_i
  fn_get.str.mu.infect = function(ig){
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur 
    str.mu.infect = paste0('mu.infect.t[',ig,',]=tm_step * (')
    for(igr in 1:(Na-1)){
      str.mu.infect=paste0(str.mu.infect,"rtBETA.t[,",ig,",",igr,"]*I.cur.t[",igr,",]/N.cur[",igr,"]+")
    }
    str.mu.infect=paste0(str.mu.infect,"rtBETA.t[,",ig,",",Na,"]*I.cur.t[",Na,",]/N.cur[",Na,"])*S.cur.t[",ig,",]")
    str.mu.infect
  }
  
  for(iv in 1:Ns){
    rtBETA.t = beta[paste0(beta_age.names,letters[iv]),]
    rownames(rtBETA.t) = beta_age.names
    # re-arrange
    rtBETA.arr = array(0, c(Np, Na, Na))
    for(ai in 1:Na){
      for(aj in 1:Na){
        rtBETA.arr[,ai,aj] = rtBETA.t[paste0('beta.',ai, aj), ]
      }
    }
    eval(parse(text = paste('rtBETA.',letters[iv],'=rtBETA.arr', sep='')))
  }
  
  for (t in 1:length(tm_vec)){
    cnt=cnt+1;
    
    # step 1
    N.cur=N.nans; 
    
    # step 1
    S.cur = S[,,cnt-1]; E.cur = E[,,cnt-1]; I.cur = I[,,cnt-1]; 
    R.cur = matrix(N.cur, NaNs, Np) - S.cur - E.cur - I.cur
    # I.cur = (I.cur^Iexp) %>% apply(2, round, 0) # to account for imperfect mixing
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    rownames(S.cur) = rownames(S0); rownames(E.cur) = rownames(E0); rownames(I.cur) = rownames(I0)
    
    # mu.infect=tm_step*(beta*I.cur*S.cur/N)
    
    # compute for each strain/variant
    mu.infect = NULL # matrix(0,NaNs,Np);
    for(iv in 1:Ns){
      rtBETA.t = get(paste0('rtBETA.', letters[iv]))
      S.cur.t = S.cur[paste0('S.', 1:Na, letters[iv]),]
      E.cur.t = E.cur[paste0('E.', 1:Na, letters[iv]),]
      I.cur.t = I.cur[paste0('I.', 1:Na, letters[iv]),]
      
      mu.infect.t =matrix(0,Na,Np);
      for (ig in 1:Na){
        # force of infection for each age group
        eval(parse(text=fn_get.str.mu.infect(ig)))
      }
      rownames(mu.infect.t) = paste0('newi.', 1:Na, letters[iv])
      mu.infect = rbind(mu.infect, mu.infect.t)
    }
    mu.infect = mu.infect[gsub('S','newi', rownames(S.cur)),]
    
    # cross immunity
    mu.icross = NULL
    for(ia in 1:Na){
      mu.icross.t=matrix(0,Ns,Np)
      mu.infect.t = mu.infect[paste0('newi.',ia, letters[1:Ns]),]
      for(si in 1:Ns){
        mu.icross.t[si,]=eval(parse(text=eval(parse(text=paste('cross.sums.t[',si,']',sep='')))))
      }
      rownames(mu.icross.t) = paste0('icross.',ia,letters[1:Ns])
      mu.icross = rbind(mu.icross, mu.icross.t)
    }
    # re-arrage to match order
    mu.icross = mu.icross[gsub('S','icross', rownames(S.cur)),]
    
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / Tei # E->I
    mu.ir = tm_step * I.cur / Tir # I -> R community transmission
    mu.rs = tm_step * R.cur / Trs # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*NaNs,mu.infect),NaNs,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*NaNs,mu.ei),NaNs,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*NaNs,mu.ir),NaNs,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*NaNs,mu.rs),NaNs,Np);
      mu.icross=matrix(rpois(Np*NaNs,mu.icross),NaNs,Np);
    }
    
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    
    # update the transmission part
    
    sk1 = - mu.infect + mu.rs - mu.icross # - seed 
    ek1 = mu.infect - mu.ei # + seed
    ik1 = mu.ei - mu.ir 
    ik1t = mu.ei
    ik1o = mu.ir
    
    Ts1 = S[,,cnt-1] + apply(sk1/2,2,round,0)
    Te1 = E[,,cnt-1] + apply(ek1/2,2,round,0)
    Ti1 = I[,,cnt-1] + apply(ik1/2,2,round,0)
    Ts1[Ts1<0]=0; Te1[Te1<0]=0; Ti1[Ti1<0]=0; 
    
    # step 2
    S.cur = Ts1; E.cur = Te1; I.cur = Ti1; 
    R.cur = matrix(N.cur, NaNs, Np) - S.cur - E.cur - I.cur
    # I.cur = (I.cur^Iexp) %>% apply(2, round, 0) # to account for imperfect mixing
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    rownames(S.cur) = rownames(S0); rownames(E.cur) = rownames(E0); rownames(I.cur) = rownames(I0)
    
    # mu.infect=tm_step*(beta*I.cur*S.cur/N)
    
    # compute for each strain/variant
    mu.infect = NULL # matrix(0,NaNs,Np);
    for(iv in 1:Ns){
      rtBETA.t = get(paste0('rtBETA.', letters[iv]))
      S.cur.t = S.cur[paste0('S.', 1:Na, letters[iv]),]
      E.cur.t = E.cur[paste0('E.', 1:Na, letters[iv]),]
      I.cur.t = I.cur[paste0('I.', 1:Na, letters[iv]),]
      
      mu.infect.t =matrix(0,Na,Np);
      for (ig in 1:Na){
        # force of infection for each age group
        eval(parse(text=fn_get.str.mu.infect(ig)))
      }
      rownames(mu.infect.t) = paste0('newi.', 1:Na, letters[iv])
      mu.infect = rbind(mu.infect, mu.infect.t)
    }
    mu.infect = mu.infect[gsub('S','newi', rownames(S.cur)),]
    
    # cross immunity
    mu.icross = NULL
    for(ia in 1:Na){
      mu.icross.t=matrix(0,Ns,Np)
      mu.infect.t = mu.infect[paste0('newi.',ia, letters[1:Ns]),]
      for(si in 1:Ns){
        mu.icross.t[si,]=eval(parse(text=eval(parse(text=paste('cross.sums.t[',si,']',sep='')))))
      }
      rownames(mu.icross.t) = paste0('icross.',ia,letters[1:Ns])
      mu.icross = rbind(mu.icross, mu.icross.t)
    }
    # re-arrage to match order
    mu.icross = mu.icross[gsub('S','icross', rownames(S.cur)),]
    
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / Tei # E->I
    mu.ir = tm_step * I.cur / Tir # I -> R community transmission
    mu.rs = tm_step * R.cur / Trs # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*NaNs,mu.infect),NaNs,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*NaNs,mu.ei),NaNs,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*NaNs,mu.ir),NaNs,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*NaNs,mu.rs),NaNs,Np);
      mu.icross=matrix(rpois(Np*NaNs,mu.icross),NaNs,Np);
    }
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    # update the transmission part
    
    sk2 = - mu.infect + mu.rs - mu.icross # - seed 
    ek2 = mu.infect - mu.ei # + seed
    ik2 = mu.ei - mu.ir 
    ik2t = mu.ei
    ik2o = mu.ir
    
    Ts2 = S[,,cnt-1] + apply(sk2/2,2,round,0)
    Te2 = E[,,cnt-1] + apply(ek2/2,2,round,0)
    Ti2 = I[,,cnt-1] + apply(ik2/2,2,round,0)
    Ts2[Ts2<0]=0; Te2[Te2<0]=0; Ti2[Ti2<0]=0; 
    
    # step 3
    S.cur = Ts2; E.cur = Te2; I.cur = Ti2; 
    R.cur = matrix(N.cur, NaNs, Np) - S.cur - E.cur - I.cur
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    rownames(S.cur) = rownames(S0); rownames(E.cur) = rownames(E0); rownames(I.cur) = rownames(I0)
    
    # mu.infect=tm_step*(beta*I.cur*S.cur/N)
    
    # compute for each strain/variant
    mu.infect = NULL # matrix(0,NaNs,Np);
    for(iv in 1:Ns){
      rtBETA.t = get(paste0('rtBETA.', letters[iv]))
      S.cur.t = S.cur[paste0('S.', 1:Na, letters[iv]),]
      E.cur.t = E.cur[paste0('E.', 1:Na, letters[iv]),]
      I.cur.t = I.cur[paste0('I.', 1:Na, letters[iv]),]
      
      mu.infect.t =matrix(0,Na,Np);
      for (ig in 1:Na){
        # force of infection for each age group
        eval(parse(text=fn_get.str.mu.infect(ig)))
      }
      rownames(mu.infect.t) = paste0('newi.', 1:Na, letters[iv])
      mu.infect = rbind(mu.infect, mu.infect.t)
    }
    mu.infect = mu.infect[gsub('S','newi', rownames(S.cur)),]
    
    # cross immunity
    mu.icross = NULL
    for(ia in 1:Na){
      mu.icross.t=matrix(0,Ns,Np)
      mu.infect.t = mu.infect[paste0('newi.',ia, letters[1:Ns]),]
      for(si in 1:Ns){
        mu.icross.t[si,]=eval(parse(text=eval(parse(text=paste('cross.sums.t[',si,']',sep='')))))
      }
      rownames(mu.icross.t) = paste0('icross.',ia,letters[1:Ns])
      mu.icross = rbind(mu.icross, mu.icross.t)
    }
    # re-arrage to match order
    mu.icross = mu.icross[gsub('S','icross', rownames(S.cur)),]
    
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / Tei # E->I
    mu.ir = tm_step * I.cur / Tir # I -> R community transmission
    mu.rs = tm_step * R.cur / Trs # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*NaNs,mu.infect),NaNs,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*NaNs,mu.ei),NaNs,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*NaNs,mu.ir),NaNs,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*NaNs,mu.rs),NaNs,Np);
      mu.icross=matrix(rpois(Np*NaNs,mu.icross),NaNs,Np);
    }
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    # update the transmission part
    
    sk3 = - mu.infect + mu.rs - mu.icross # - seed 
    ek3 = mu.infect - mu.ei # + seed
    ik3 = mu.ei - mu.ir 
    ik3t = mu.ei
    ik3o = mu.ir
    
    Ts3 = S[,,cnt-1] + apply(sk3,2,round,0)
    Te3 = E[,,cnt-1] + apply(ek3,2,round,0)
    Ti3 = I[,,cnt-1] + apply(ik3,2,round,0)
    Ts3[Ts3<0]=0; Te3[Te3<0]=0; Ti3[Ti3<0]=0; 
    
    # step 4
    S.cur = Ts3; E.cur = Te3; I.cur = Ti3; 
    R.cur = matrix(N.cur, NaNs, Np) - S.cur - E.cur - I.cur
    # I.cur = I.cur^Iexp %>% apply(2, round, 0) # to account for imperfect mixing
    I.cur = (pmin(I.cur^Iexp, I.cur)) %>% apply(2, round, 0) # to account for imperfect mixing
    rownames(S.cur) = rownames(S0); rownames(E.cur) = rownames(E0); rownames(I.cur) = rownames(I0)
    
    # mu.infect=tm_step*(beta*I.cur*S.cur/N)
    
    # compute for each strain/variant
    mu.infect = NULL # matrix(0,NaNs,Np);
    for(iv in 1:Ns){
      rtBETA.t = get(paste0('rtBETA.', letters[iv]))
      S.cur.t = S.cur[paste0('S.', 1:Na, letters[iv]),]
      E.cur.t = E.cur[paste0('E.', 1:Na, letters[iv]),]
      I.cur.t = I.cur[paste0('I.', 1:Na, letters[iv]),]
      
      mu.infect.t =matrix(0,Na,Np);
      for (ig in 1:Na){
        # force of infection for each age group
        eval(parse(text=fn_get.str.mu.infect(ig)))
      }
      rownames(mu.infect.t) = paste0('newi.', 1:Na, letters[iv])
      mu.infect = rbind(mu.infect, mu.infect.t)
    }
    mu.infect = mu.infect[gsub('S','newi', rownames(S.cur)),]
    
    # cross immunity
    mu.icross = NULL
    for(ia in 1:Na){
      mu.icross.t=matrix(0,Ns,Np)
      mu.infect.t = mu.infect[paste0('newi.',ia, letters[1:Ns]),]
      for(si in 1:Ns){
        mu.icross.t[si,]=eval(parse(text=eval(parse(text=paste('cross.sums.t[',si,']',sep='')))))
      }
      rownames(mu.icross.t) = paste0('icross.',ia,letters[1:Ns])
      mu.icross = rbind(mu.icross, mu.icross.t)
    }
    # re-arrage to match order
    mu.icross = mu.icross[gsub('S','icross', rownames(S.cur)),]
    
    
    # mu.infect = tm_step * (beta.i * I.cur + beta.h * H.cur) * S.cur / N.cur  # new infection
    mu.ei = tm_step * E.cur / Tei # E->I
    mu.ir = tm_step * I.cur / Tir # I -> R community transmission
    mu.rs = tm_step * R.cur / Trs # R -> S loss immunity
    
    mu.infect[mu.infect<0]=0; mu.ei[mu.ei<0]=0; mu.ir[mu.ir<0]=0; mu.rs[mu.rs<0]=0
    mu.icross[mu.icross<0]=0;
    
    
    # stochastisity - random grape from a Poisson distribution
    if(stoch){
      mu.infect = matrix(rpois(Np*NaNs,mu.infect),NaNs,Np); # rpois(Np, mu.infect)
      mu.ei = matrix(rpois(Np*NaNs,mu.ei),NaNs,Np); # rpois(Np,mu.ei)
      mu.ir = matrix(rpois(Np*NaNs,mu.ir),NaNs,Np); # rpois(Np,mu.ir)
      mu.rs = matrix(rpois(Np*NaNs,mu.rs),NaNs,Np);
      mu.icross=matrix(rpois(Np*NaNs,mu.icross),NaNs,Np);
    }
    # check DA physicality
    mu.infect = pmin(mu.infect, S.cur) # new cases < S
    mu.ei = pmin(mu.ei, E.cur)
    mu.ir = pmin(mu.ir, I.cur)
    mu.rs = pmin(mu.rs, R.cur)
    mu.icross = pmin(mu.icross, S.cur)
    
    # update the transmission part
    
    sk4 = - mu.infect + mu.rs - mu.icross # - seed 
    ek4 = mu.infect - mu.ei # + seed
    ik4 = mu.ei - mu.ir 
    ik4t = mu.ei
    ik4o = mu.ir
    
    
    # add vaccination
    # add vaccination: adjust for ve and prior infection 
    # but later on, with vaccinees making up for the majority of immune, 
    # this would no longer reflective of the ture prob of vaccinee's prior immunity
    # should put a lower bound for suscept/upper bound for prior immunity
    # percS = pmax(S[,,cnt-1] / Npops, matrix(percSmax, Ns, Np))
    # use the cumulative infection rate instead, 
    vacc.d1 = vacc.d2 = NULL # matrix(0, NaNs, Np)
    
    # V1: vacc by date (row) for each age group (col)
    # assemble for different variants, account for reduction
    for(iv in 1:Ns){  # different effectiveness for diff variants
      percS = percSmax.t[paste0('percSmax.',1:Na,letters[iv]),] # matrix: Ns x Np
      vacc.d1.t = percS * (V1[cnt-1,] *  VE1 * VE1redn[iv]) # for diff age groups # number ppl vaccinated after first dose of vaccine
      # dose 2: should we account for prior immunity as well? if so it should be ~1 month ago
      # for simplicity, use the same percS
      vacc.d2.t = percS * (V2[cnt-1,] * (1 - VE1) * VE2  * VE2redn[iv])
      
      rownames(vacc.d1.t) = paste0('vx.',1:Na,letters[iv])
      rownames(vacc.d2.t) = paste0('vx.',1:Na,letters[iv])
      
      vacc.d1 = rbind(vacc.d1, vacc.d1.t)
      vacc.d2 = rbind(vacc.d2, vacc.d2.t)
    }
    vacc.d1 = vacc.d1[gsub('S','vx', rownames(S.cur)),]
    vacc.d2 = vacc.d2[gsub('S','vx', rownames(S.cur)),]
    # seed = seed # 
    
    seed.t = matrix(rpois(Np*NaNs, seed0),NaNs,Np);
    
    S[,,cnt]=S[,,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + tm_step*birthrate*(N.nans-S[,,cnt-1]) - vacc.d1 - vacc.d2 # add births/deaths
    E[,,cnt] = E[,,cnt-1] + apply(ek1/6+ek2/3+ek3/3+ek4/6,2,round,0) + seed.t
    I[,,cnt] = I[,,cnt-1] + apply(ik1/6+ik2/3+ik3/3+ik4/6,2,round,0)
    cumItot[,,cnt] = cumItot[,,cnt-1] + apply(ik1t/6+ik2t/3+ik3t/3+ik4t/6,2,round,0)
    cumIobs[,,cnt] = cumIobs[,,cnt-1] + apply(ik1o/6+ik2o/3+ik3o/3+ik4o/6,2,round,0)
    
    S[S<0]=0; E[E<0]=0; I[I<0]=0; cumItot[cumItot<0]=0; cumIobs[cumIobs<0]=0; 
    
    dimnames(cumItot)[1] = list(gsub('S','cumItot',rownames(S.cur)))
    dimnames(cumIobs)[1] = list(gsub('S','cumIobs',rownames(S.cur)))
  }
  
  # include the delayed reporting etc.
  # account for delay in case reporting and death/reporting
  death=cumIobs # place holder
  dimnames(death)[1] = list(gsub('S','death',rownames(S.cur)))
  
  {
    
    # for diff age groups
    newI_ts_all = array(0, c(NaNs, Np, tmstep ))
    dimnames(newI_ts_all)[1] = list(gsub('S','daily.newItot',rownames(S.cur)))
      
    
    for(ia in 1:Na){
      severity.t = get(paste0('severity.',ia))
      alpha.t = parm0[paste0('alpha.',ia),]
      dist_tm.to.detect.t = get(paste0('dist_tm.to.detect.',ia))
      dist_tm.to.death.t = get(paste0('dist_tm.to.death.',ia))
      tm.from.inf.to.death.max.t = get(paste0('tm.from.inf.to.death.max.',ia))
      
      # est.daily.tot_all= array(0, c(Ns, Np, tmstep )) # this account for delay in reporting, but not under-reporing
      # not the newIobs in state.new account for both delay in reporting and under-reporing
      # est.daily.death_all = array(0, c(Ns, Np, tmstep ))
      
      
      # NEED TO DO IT FOR EACH GROUP
      for(ig in 1:Ns){
        
        if(! is.null(newI.previous00)){
          newI.previous00.t = t(newI.previous00[paste0('ihist.',ia,letters[ig]),,])
        } else {
          newI.previous00.t = NULL
        }
        
        
        newI_ts = t(cumItot[paste0('cumItot.',ia,letters[ig]),,]) # these are cummulative cases
        newI_ts = newI_ts[-1,,drop=F] - newI_ts[-nrow(newI_ts),,drop=F] # total cases, without delay or under-reporting
        
        # need to include previous cases (not yet detected as well)
        newI.previous = tail(newI.previous00.t,tm.to.detect.max) 
        newI.combined = rbind(newI.previous, newI_ts)
        # est.obs = obs.delayedSEIR(newI_ts, dist_tm.to.detect)
        est.daily.tot = obs.delayedSEIR(newI.combined, dist_tm.to.detect.t)
        idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
        est.daily.tot = est.daily.tot[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
        
        # now account for report rate
        est.obs.daily = est.daily.tot * matrix(alpha.t,tmstep,Np,byrow = T)
        
        # death:
        newI.previous = tail(newI.previous00.t,tm.from.inf.to.death.max)
        newI.combined = rbind(newI.previous, newI_ts)
        est.daily.death = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.death.t, tm.to.outcome.max=tm.from.inf.to.death.max.t)
        idx.start = nrow(newI.previous); if(is.null(idx.start)) idx.start=0
        est.daily.death =est.daily.death[idx.start+1:tmstep,] # exclude previous week and after, only include days for this time step
        # now account for severity
        est.daily.death = est.daily.death * matrix(severity.t[grep('death',rownames(severity.t)),],tmstep,Np,byrow = T)
        
        # back to cumIobs:
        est.obs.this=rbind(0, apply(est.obs.daily,2,cumsum))
        
        # cummulative deaths
        est.death.this = rbind(0, apply(est.daily.death,2,cumsum))
        
        cumIobs[paste0('cumIobs.',ia,letters[ig]),,] = t(est.obs.this) # re-assign observed
        
        death[paste0('death.',ia,letters[ig]),,] = t(est.death.this) # re-assign observed
        
        # save for this group
        # est.daily.tot_all[paste0('dailyItot.',ia,letters[ig]),,] = t(est.daily.tot)  # this account for delay in reporting, but not under-reporing
        newI_ts_all[paste0('daily.newItot.',ia,letters[ig]),,] = t(newI_ts)
        
        
      } # group
    } # age
    
    
    
  } # get case and death, with delay and under reporting
  
  # return
  return(list(S=S,E=E,I=I,cumIobs=cumIobs,cumItot=cumItot, death = death, daily.newItot = newI_ts_all))
}





