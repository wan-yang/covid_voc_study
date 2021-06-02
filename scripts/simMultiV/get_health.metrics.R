# to compute hospitalization, ed, healthcare demand etc.
# get weekly summary

# aggregate across all ages for each variant
fn_get.aggregate.by.v = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  res = NULL; tmp.ens = array(0, c(Ns, num_ens, num_times))
  for(iv in 1:Ns){
    
    tmp0 = tda[paste0('hm.', 1:Na, letters[iv]),,] %>% apply(c(2,3), sum, na.rm=T) # sum over age groups
    # save individual ens
    tmp.ens[iv,,] = tmp0
    tmp = cbind(tmp0 %>% apply(2, mean, na.rm=T),
                tmp0 %>% apply(2, sd, na.rm=T),
                tmp0 %>% apply(2, quantile, probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
    
    if(per.capita)
      tmp = tmp / N * 100
    colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
    res = rbind(res, data.table(variant = variants.t[iv], age.grp = 'All', Week.start = week.starts, tmp))
  }
  # res$variant = factor(res$variant, levels = c('All','wt', 'b117', 'b1351', 'p1'), 
  #                     labels = c('All', 'Wildtype','UK variant', 'SA variant', 'BR variant'))
  res$variant = factor(res$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                       labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  
  # also compute the percentage by each variant
  tmp.perc = array(0, c(num_ens, num_times, Ns))
  for(ii in 1:num_ens){
    tmp = tmp.ens[,ii,] %>% t
    tmp = tmp / rowSums(tmp) * 100
    tmp.perc[ii,,] = tmp
  }
  perc = NULL
  for(iv in 1:Ns){
    
    tmp0 = tmp.perc[,,iv]
    
    tmp = cbind(tmp0 %>% apply(2, mean, na.rm=T),
                tmp0 %>% apply(2, sd, na.rm=T),
                tmp0 %>% apply(2, quantile, probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
    
    colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
    perc = rbind(perc, data.table(variant = variants.t[iv], age.grp = 'All', Week.start = week.starts, tmp))
  }
  perc$variant = factor(perc$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                       labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  perc$age.grp = factor(perc$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  
  
  return(list(res = res, perc = perc))

}

fn_get.prec.v = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  res = NULL
  for(iv in 1:Ns){
    
    tmp0 = tda[paste0('hm.', 1:Na, letters[iv]),,] %>% apply(c(2,3), sum, na.rm=T) # sum over age groups
    tmp = cbind(tmp0 %>% apply(2, mean, na.rm=T),
                tmp0 %>% apply(2, sd, na.rm=T),
                tmp0 %>% apply(2, quantile, probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
    
    if(per.capita)
      tmp = tmp / N * 100
    colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
    res = rbind(res, data.table(variant = variants.t[iv], age.grp = 'All', Week.start = week.starts, tmp))
  }
  # res$variant = factor(res$variant, levels = c('All','wt', 'b117', 'b1351', 'p1'), 
  #                     labels = c('All', 'Wildtype','UK variant', 'SA variant', 'BR variant'))
  res$variant = factor(res$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                       labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  
  # also compute the percentage by each variant
  tmp = dcast(res, Week.start ~ variant, value.var = 'v.mean')
  perc = tmp[,-1] / rowSums(tmp[,-1]) * 100
  
  perc
}

# get aggregate, over all variants, ages
fn_get.aggregate = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  tmp0 = tda %>% apply(c(2,3), sum, na.rm=T) # sum over age groups
  tmp = cbind(tmp0 %>% apply(2, mean, na.rm=T),
              tmp0 %>% apply(2, sd, na.rm=T),
              tmp0 %>% apply(2, quantile, probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
  
  if(per.capita)
    tmp = tmp / N * 100
  colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
  
  res = data.table(variant = 'All', age.grp = 'All', Week.start = week.starts, tmp)
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  res
}

# get aggregate, over all variants, ages, weeks
fn_get.sum = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  tmp0 = tda %>% apply(c(2), sum, na.rm=T) # sum over age groups c(2,3)
  tmp = cbind(tmp0 %>% mean(na.rm=T),
              tmp0 %>% sd(na.rm=T),
              tmp0 %>% quantile(probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
  
  if(per.capita)
    tmp = tmp / N * 100
  colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
  
  res = data.table(variant = 'All', age.grp = 'All', tmp)
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  res
}

fn_get.cumsum = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  # first sum over age groups and variant, then cumsum
  tmp00 = tda %>% apply(c(2, 3), sum, na.rm=T)
  tmp0 = matrix(0, num_ens, ncol(tmp00))
  for(ii in 1: num_ens){
    tmp0[ii,] = tmp00[ii,] %>% cumsum
  }
  
  tmp = cbind(tmp0 %>% apply(2, mean, na.rm=T),
              tmp0 %>% apply(2, sd, na.rm=T),
              tmp0 %>% apply(2, quantile, probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
  
  
  if(per.capita)
    tmp = tmp / N * 100
  colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
  
  res = data.table(variant = 'All', age.grp = 'All', Week.start = week.starts, tmp)
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  res
}

fn_get.sum.by.v = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  res = NULL; tmp.ens = matrix(0, Ns, num_ens)
  for(iv in 1:Ns){
    
    tmp0 = tda[paste0('hm.', 1:Na, letters[iv]),,] %>% apply(c(2), sum, na.rm=T) # sum over age groups
    # save individual ens
    tmp.ens[iv,] = tmp0
    tmp = cbind(tmp0 %>% mean(na.rm=T),
                tmp0 %>% sd(na.rm=T),
                tmp0 %>% quantile(probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
    
    if(per.capita)
      tmp = tmp / N * 100
    colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
    res = rbind(res, data.table(variant = variants.t[iv], age.grp = 'All',  tmp))
  }
  # res$variant = factor(res$variant, levels = c('All','wt', 'b117', 'b1351', 'p1'), 
  #                     labels = c('All', 'Wildtype','UK variant', 'SA variant', 'BR variant'))
  res$variant = factor(res$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                       labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  
  # also compute the percentage by each variant
  tmp.perc =  t(tmp.ens) / colSums(tmp.ens) * 100
  
  perc = NULL
  for(iv in 1:Ns){
    
    tmp0 = tmp.perc[,iv]
    
    tmp = cbind(tmp0 %>% mean(na.rm=T),
                tmp0 %>% sd(na.rm=T),
                tmp0 %>% quantile(probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
    
    colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
    perc = rbind(perc, data.table(variant = variants.t[iv], age.grp = 'All',  tmp))
  }
  perc$variant = factor(perc$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                        labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  perc$age.grp = factor(perc$age.grp, levels = c('All', numbers[1:Na]), 
                        labels = c('All', age.grp.names))
  
  
  return(list(res = res, perc = perc))
  
}

fn_get.cumsum.by.v = function(tda, per.capita = F, N = 1e6){
  
  dimnames(tda)[1] = list(outer(paste0('hm.', numbers[1:Na]), letters[1:Ns], FUN = paste0) %>% c)
  
  res = NULL; tmp.ens = array(0, c(Ns, num_ens, num_times))
  for(iv in 1:Ns){
    
    tmp00 = tda[paste0('hm.', 1:Na, letters[iv]),,] %>% 
              apply(c(2, 3), sum, na.rm=T)
      
    tmp0 = matrix(0, num_ens, ncol(tmp00))
    for(ii in 1: num_ens){
      tmp0[ii,] = tmp00[ii,] %>% cumsum
    }
    
    # save individual ens
    tmp.ens[iv,,] = tmp0
    tmp = cbind(tmp0 %>% apply(2, mean, na.rm=T),
                tmp0 %>% apply(2, sd, na.rm=T),
                tmp0 %>% apply(2, quantile, probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
    
    
    if(per.capita)
      tmp = tmp / N * 100
    colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
    
    res = rbind(res, data.table(variant = variants.t[iv], age.grp = 'All', Week.start = week.starts, tmp))
    
  }
  # res$variant = factor(res$variant, levels = c('All','wt', 'b117', 'b1351', 'p1'), 
  #                     labels = c('All', 'Wildtype','UK variant', 'SA variant', 'BR variant'))
  res$variant = factor(res$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                       labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  
  # also compute the percentage by each variant
  if(F){
    perc = NULL
    for(tt in 1:num_times){
      tmp.perc =  t(tmp.ens[,,tt]) / colSums(tmp.ens[,,tt]) * 100
      for(iv in 1:Ns){
        tmp0 = tmp.perc[,iv]
        tmp = cbind(tmp0 %>% mean(na.rm=T),
                    tmp0 %>% sd(na.rm=T),
                    tmp0 %>% quantile(probs=c(.5, .25, .75, .025, .975), na.rm=T) %>% t)
        
        colnames(tmp) = c('v.mean', 'v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')
        perc = rbind(perc, data.table(variant = variants.t[iv], age.grp = 'All', Week.start = week.starts[tt], tmp))
      }
    }
    # not useful
  }
  
  res
  
}

fn_get.stats.v.age = function(x, NaNs){  # by variant, by age,
  res = data.frame(variant=rep(rep(variants.t,e=Na), length(week.starts)), 
             age.grp = rep(numbers[1:Na], length(variants.t) * length(week.starts)),
             Week.start=rep(week.starts, e = NaNs),
             v.mean=c(apply(x,c(1,3),mean, na.rm=T)),
             v.sd=c(apply(x,c(1,3),sd, na.rm=T)), 
             v.median=c(apply(x,c(1,3),quantile,.5, na.rm=T)),
             iqr.lwr=c(apply(x,c(1,3),quantile,.25, na.rm=T)),
             iqr.upr=c(apply(x,c(1,3),quantile,.75, na.rm=T)),
             # v.90lwr=c(apply(x,c(1,3),quantile,.05)),
             # v.90upr=c(apply(x,c(1,3),quantile,.95)),
             ci95.lwr=c(apply(x,c(1,3),quantile,.025, na.rm=T)),
             ci95.upr=c(apply(x,c(1,3),quantile,.975, na.rm=T)))
  
  # res$variant = factor(res$variant, levels = c('All','wt', 'b117', 'b1351', 'p1'), 
  #                     labels = c('All', 'Wildtype','UK variant', 'SA variant', 'BR variant'))
  res$variant = factor(res$variant, levels = c('All','wt', 'b1526', 'b1427', 'b117', 'b1351', 'p1'), 
                       labels = c('All', 'Wildtype','B.1.526','B.1.427/9', 'B.1.1.7', 'B.1.351','P.1'))
  res$age.grp = factor(res$age.grp, levels = c('All', numbers[1:Na]), 
                       labels = c('All', age.grp.names))
  res
}


fn_get.health.metrics = function(newI.daily # daily number of infections
                                 ){
  
  newI.previous00 = newI.daily
  dim.t = dim(newI.previous00)
  NaNs = dim.t[1]; Np = dim.t[2]
  num_days = dim.t[3] - nwk.pre * 7 
  
  xpost.daily.inf= array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.case= array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.death= array(0, c(NaNs, num_ens, tmstep * num_times))
  
  xpost.daily.ed= array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.hospital= array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.icu = array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.nonicuhosp = array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.vent = array(0, c(NaNs, num_ens, tmstep * num_times))
  
  xpost.daily.hospitalPrev=array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.icuPrev = array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.ventPrev = array(0, c(NaNs, num_ens, tmstep * num_times))
  xpost.daily.nonicuhospPrev = array(0, c(NaNs, num_ens, tmstep * num_times))
  
  dimnames(xpost.daily.inf)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.case)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.death)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  
  
  dimnames(xpost.daily.ed)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.hospital)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.icu)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.vent)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.nonicuhosp)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.hospitalPrev)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.icuPrev)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.ventPrev)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  dimnames(xpost.daily.nonicuhospPrev)[1] = list(outer(paste0('hm.',1:Na), letters[1:Ns], FUN = paste0) %>% c)
  
  {
    
    # for diff age groups
    
    for(ia in 1:Na){
      severity.t = get(paste0('severity.',ia))
      alpha.t = parm0[paste0('alpha.',ia),]
      dist_tm.to.detect.t = get(paste0('dist_tm.to.detect.',ia))
      dist_tm.to.death.t = get(paste0('dist_tm.to.death.',ia))
      tm.from.inf.to.death.max.t = get(paste0('tm.from.inf.to.death.max.',ia))
      
      
      dist_tm.to.hospital.t = get(paste0('dist_tm.to.hospital.',ia))
      dist_tm.to.icu.t = get(paste0('dist_tm.to.icu.',ia))
      dist_hospital.stay.t = get(paste0('dist_hospital.stay.',ia))
      dist_icu.stay.t = get(paste0('dist_icu.stay.',ia))
      dist_vent.stay.t = get(paste0('dist_vent.stay.',ia))
      
      
      # NEED TO DO IT FOR EACH GROUP
      for(ig in 1:Ns){
        
        if(! is.null(newI.previous00)){
          newI.previous00.t = t(newI.previous00[paste0('daily.newI.',ia,letters[ig]),,])
        } else {
          newI.previous00.t = NULL
        }
        
        # need to include previous cases (not yet detected as well)
        newI.combined = newI.previous00.t # both training and wkly
        n.tot = nrow(newI.combined)
        idx.start = nwk.pre * 7 # from previousl week
        
        
        tmp = outcomeInsOuts.delayedSEIR(newI_ts=newI.combined, 
                                         dist_tm.to.outcomeIns=dist_tm.to.hospital.t, 
                                         dist_tm.to.outcomeOuts=dist_hospital.stay.t,
                                         tm.to.outcomeIns.max=tm.to.hospital.max, 
                                         tm.to.outcomeOuts.max=tm.to.hospital.max + hospital.stay.max)
        tmpIns = tmp$est.outcomeIns
        tmpPrev = tmp$est.outcomePrev
        tmpOuts = tmp$est.outcomeOuts
        
        this.daily.hospital =tmpIns[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        this.daily.hospitalPrev =tmpPrev[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        # this.daily.hosp.discharge =tmpOuts[(idx.start+1):n.tot,] 
        # now account for severity
        this.daily.hospital = this.daily.hospital * matrix(severity.t[grep('hospital',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T)
        this.daily.hospitalPrev = this.daily.hospitalPrev * matrix(severity.t[grep('hospital',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T)
        
        # icu:
        # newI.previous = tail(newI.previous00.t,tm.to.icu.max) # do not truncate, need all for discharge
        # this.daily.icu = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.icu, tm.to.outcome.max=tm.to.icu.max)
        tmp = outcomeInsOuts.delayedSEIR(newI_ts=newI.combined, 
                                         dist_tm.to.outcomeIns=dist_tm.to.icu.t, 
                                         dist_tm.to.outcomeOuts=dist_icu.stay.t,
                                         tm.to.outcomeIns.max=tm.to.icu.max, 
                                         tm.to.outcomeOuts.max=tm.to.icu.max + icu.stay.max)
        tmpIns = tmp$est.outcomeIns
        tmpPrev = tmp$est.outcomePrev
        tmpOuts = tmp$est.outcomeOuts
        this.daily.icu =tmpIns[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        this.daily.icuPrev =tmpPrev[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        
        # now account for severity
        this.daily.icu = this.daily.icu * matrix(severity.t[grep('icu',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T)
        this.daily.icuPrev = this.daily.icuPrev * matrix(severity.t[grep('icu',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T)
        
        # ADD EST FOR VENTILATOR
        # same as icu
        tmp = outcomeInsOuts.delayedSEIR(newI_ts=newI.combined, 
                                         dist_tm.to.outcomeIns=dist_tm.to.icu.t, # new icu rate, need further scale but same rate
                                         dist_tm.to.outcomeOuts=dist_vent.stay.t, # going out, based on length of vent use
                                         tm.to.outcomeIns.max=tm.to.icu.max, 
                                         tm.to.outcomeOuts.max=tm.to.icu.max + vent.stay.max)
        tmpIns = tmp$est.outcomeIns
        tmpPrev = tmp$est.outcomePrev
        tmpOuts = tmp$est.outcomeOuts
        this.daily.vent =tmpIns[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        this.daily.ventPrev =tmpPrev[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        # now account for severity
        this.daily.vent = this.daily.vent * matrix(severity.t[grep('icu',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T) * 
          matrix(severity.t[grep('vent',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T)
        this.daily.ventPrev = this.daily.ventPrev * matrix(severity.t[grep('icu',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T) *
          matrix(severity.t[grep('vent',rownames(severity.t)),],tmstep * num_times,num_ens,byrow = T)
        
        # non icu hospitalization
        this.daily.nonicuhosp = this.daily.hospital - this.daily.icu
        this.daily.nonicuhospPrev = this.daily.hospitalPrev - this.daily.icuPrev
        this.daily.nonicuhosp[this.daily.nonicuhosp<0] = 0
        this.daily.nonicuhospPrev[this.daily.nonicuhospPrev<0] = 0
        
        # case
        this.daily.case = obs.delayedSEIR(newI.combined, dist_tm.to.detect.t)
        this.daily.case = this.daily.case[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        # now account for report rate
        this.daily.case = this.daily.case * matrix(alpha.t,tmstep * num_times, num_ens, byrow = T)
        
        # death:
        this.daily.death = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.death.t, tm.to.outcome.max=tm.from.inf.to.death.max.t)
        this.daily.death =this.daily.death[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        # now account for severity
        this.daily.death = this.daily.death * matrix(severity.t[grep('death',rownames(severity.t)),],tmstep * num_times, num_ens,byrow = T)
        
        # ed.visit:
        # est.daily.hospital = outcome.delayedSEIR(newI_ts=newI.combined, dist_tm.to.outcome=dist_tm.to.hospital, tm.to.outcome.max=tm.to.hospital.max)
        this.daily.ed = outcome.delayedSEIR(newI_ts=newI.combined, 
                                           dist_tm.to.outcome=dist_tm.to.hospital.t, 
                                           tm.to.outcome.max=tm.to.hospital.max)
        
        this.daily.ed = this.daily.ed[(idx.start+1):n.tot,] # exclude previous week and after, only include days for this time step
        # now account for severity
        this.daily.ed = this.daily.ed * matrix(severity.t[grep('edr',rownames(severity.t)),],tmstep * num_times, num_ens,byrow = T)
        
        # infection
        this.daily.inf = newI.combined[(idx.start+1):n.tot,]
        
        # save for this group
        xpost.daily.inf[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.inf)
        xpost.daily.case[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.case)
        xpost.daily.death[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.death)
        xpost.daily.ed[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.ed)
        xpost.daily.hospital[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.hospital)
        xpost.daily.icu[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.icu)
        xpost.daily.nonicuhosp[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.nonicuhosp)
        xpost.daily.vent[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.vent)
        xpost.daily.hospitalPrev[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.hospitalPrev)
        xpost.daily.icuPrev[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.icuPrev)
        xpost.daily.ventPrev[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.ventPrev)
        xpost.daily.nonicuhospPrev[paste0('hm.', numbers[ia],letters[ig]),,] = t(this.daily.nonicuhospPrev)
        
        
      } # variant
    } # age
    
    # get the weekly est
    wkly.inf = wkly.case = wkly.death =
      wkly.ed = wkly.hospital = wkly.icu = wkly.nonicuhosp = wkly.vent =
      wkly.hospitalPrev.mean = wkly.icuPrev.mean = wkly.nonicuhospPrev.mean = wkly.ventPrev.mean =
      wkly.hospitalPrev.max = wkly.icuPrev.max = wkly.nonicuhospPrev.max = wkly.ventPrev.max =
      array(0, c(NaNs, num_ens, num_times))
    
    for(iwk in 1:num_times){
      idx = (iwk - 1) * tmstep + 1:tmstep
      
      wkly.inf[,,iwk]=apply(xpost.daily.inf[,,idx],c(1,2), sum, na.rm=T)
      wkly.case[,,iwk]=apply(xpost.daily.case[,,idx],c(1,2), sum, na.rm=T)
      wkly.death[,,iwk]=apply(xpost.daily.death[,,idx],c(1,2), sum, na.rm=T)
      
      wkly.ed[,,iwk]=apply(xpost.daily.ed[,,idx],c(1,2), sum, na.rm=T)
      wkly.hospital[,,iwk]=apply(xpost.daily.hospital[,,idx],c(1,2), sum, na.rm=T)
      wkly.icu[,,iwk]=apply(xpost.daily.icu[,,idx],c(1,2), sum, na.rm=T)
      wkly.nonicuhosp[,,iwk]=apply(xpost.daily.nonicuhosp[,,idx],c(1,2), sum, na.rm=T)
      wkly.vent[,,iwk]=apply(xpost.daily.vent[,,idx],c(1,2), sum, na.rm=T)
      
      wkly.hospitalPrev.mean[,,iwk]=apply(xpost.daily.hospitalPrev[,,idx],c(1,2), mean, na.rm=T)
      wkly.icuPrev.mean[,,iwk]=apply(xpost.daily.icuPrev[,,idx],c(1,2),mean, na.rm=T)
      wkly.nonicuhospPrev.mean[,,iwk]=apply(xpost.daily.nonicuhospPrev[,,idx],c(1,2),mean, na.rm=T)
      wkly.ventPrev.mean[,,iwk]=apply(xpost.daily.ventPrev[,,idx],c(1,2),mean, na.rm=T)
      
      wkly.hospitalPrev.max[,,iwk]=apply(xpost.daily.hospitalPrev[,,idx],c(1,2), max, na.rm=T)
      wkly.icuPrev.max[,,iwk]=apply(xpost.daily.icuPrev[,,idx],c(1,2),max, na.rm=T)
      wkly.nonicuhospPrev.max[,,iwk]=apply(xpost.daily.nonicuhospPrev[,,idx],c(1,2),max, na.rm=T)
      wkly.ventPrev.max[,,iwk]=apply(xpost.daily.ventPrev[,,idx],c(1,2),max, na.rm=T)
    }
  } # get health metrics, with delay and under reporting
  
  metrics = paste0('wkly.',c('inf', 'case', 'death','ed','hospital','icu','nonicuhosp','vent',
              'hospitalPrev.mean','icuPrev.mean','nonicuhospPrev.mean','ventPrev.mean',
              'hospitalPrev.max','icuPrev.max','nonicuhospPrev.max','ventPrev.max')
  )
  measures = c('New infections', 'New cases', 'New deaths',
               'New ED visits', 'New hospitalizations','New ICU admissions', 'New non-ICU hospitalizations', 'New intubations',
               'Hospital bed needs (mean)','ICU bed needs (mean)','Non-ICU hospital bed needs (mean)','Ventilator needs (mean)',
               'Hospital bed needs (max)','ICU bed needs (max)','Non-ICU hospital bed needs (max)','Ventilator needs (max)')
  
  res = NULL # by.age.by.variants
  
  for(im in 1: length(metrics)){
    mm = metrics[im]
    mm = get(mm)
    eval(parse(text = paste('tmp = fn_get.stats.v.age(mm, NaNs = NaNs)',sep=''))) 
    
    res = rbind(res, data.table(measure = measures[im], tmp))
  }
  
  # combine all age groups, by variant
  perc = NULL
  for(im in 1: length(metrics)){
    mm = metrics[im]
    mm = get(mm)
    eval(parse(text = paste('tmp = fn_get.aggregate.by.v(tda = mm)',sep=''))) 
    
    res = rbind(res, data.table(measure = measures[im], tmp$res))
    perc = rbind(perc, data.table(measure = measures[im], tmp$perc))
  }
  
  # combine all age groups, and all variants
  for(im in 1: length(metrics)){
    mm = metrics[im]
    mm = get(mm)
    eval(parse(text = paste('tmp = fn_get.aggregate(tda = mm)',sep=''))) 
    
    res = rbind(res, data.table(measure = measures[im], tmp))
  }
  
  # combine all age groups, and all variants, all weeks
  sums = NULL
  for(im in 1: length(metrics)){
    mm = metrics[im]
    mm = get(mm)
    eval(parse(text = paste('tmp = fn_get.sum(tda = mm)',sep=''))) 
    
    sums = rbind(sums, data.table(measure = measures[im], tmp))
  }
  
  sums.perc = NULL; sums.res = NULL
  for(im in 1: length(metrics)){
    mm = metrics[im]
    mm = get(mm)
    eval(parse(text = paste('tmp = fn_get.sum.by.v(tda = mm)',sep=''))) 
    
    sums.res = rbind(sums.res, data.table(measure = measures[im], tmp$res))
    sums.perc = rbind(sums.perc, data.table(measure = measures[im], tmp$perc))
  }
  
  # get cumulative sums
  cumsums = NULL
  for(im in 1: length(metrics)){
    mm = metrics[im]
    mm = get(mm)
    eval(parse(text = paste('tmp = fn_get.cumsum(tda = mm)',sep=''))) 
    cumsums = rbind(cumsums, data.table(measure = measures[im], tmp))
    
    eval(parse(text = paste('tmp = fn_get.cumsum.by.v(tda = mm)',sep=''))) 
    cumsums = rbind(cumsums, data.table(measure = measures[im], tmp))
  }
  
  return(list(res=res, perc = perc, sums.res =  sums.res, sums.perc = sums.perc, sums = sums, cumsums = cumsums))
}
