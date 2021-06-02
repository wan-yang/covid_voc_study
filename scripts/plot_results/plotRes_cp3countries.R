# to plot and compare the overall trends in the three countries
# for the MS

library(data.table)
library(magrittr)
library(ggplot2)
library(xlsx)
library(MMWRweek)
library(lemon)
library('readr')
library('readxl')
library('writexl')
library('stringi')
# library(tidyverse)
library(tgp)

N = 1e6; # population size is set to 1M

eval.t = 'obs.more'
dir_data = '../../data/'
dir_code = '../model_inference/'
dir_code2 =  './' 

dir_res = paste0('../../results/voc/')
dir_plot = paste0('../../results/plots/')

if(!file.exists(dir_plot)) dir.create(dir_plot)

source(paste0(dir_code2,'fn_util.R'))
source(paste0(dir_code2,'getPlot.R'))
source(paste0(dir_code,'get_relR0.R'))


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

num_ens = 2
da = read.csv(paste0(dir_data,'da_case.death.mob_uk.sa.br.csv'), stringsAsFactors = F)  %>% data.table()
da[date < as.Date('2020-02-09') & is.na(da)] = 0
da = da[date >= as.Date('2020-03-01')]
da = da[complete.cases(da)]

cnty = c('(A) United Kingdom', "(B) South Africa", '(C) Brazil')
cnt = 0
pdf(paste0(dir_plot, 'FigS2_cp.sn.npi.trends.pdf'), width = 6, height = 6)
par(mfrow = c(3, 1), mar=c(0, 1.5, 0, 1.5), cex = .8, cex.lab=.95, cex.axis = .95, oma = c(2, 1.5, .5, 1.5), mgp = c(1.3, .3, 0), tck = -.02)
for(loc.t in c('uk','sa', 'br')){
  cnt = cnt + 1
  da.t = da[,c('date','year','week',paste0(c('case.','death.','mob.bus.','mob.full.'),loc.t)),with=F] %>% 
    setnames(paste0(c('case.','death.','mob.bus.','mob.full.'),loc.t), c('case','death','mob.bus','mob.full'))
  # only when it is > 10 cases
  # da.t = da.t[case > 1]
  # da.t = da.t[complete.cases(da.t)]
  da.t$date = da.t$date %>% as.Date
  
  rel.mob = 1 + da.t$mob.bus/100# as.matrix() 
  
  # get relative R0 for seasonality
  relR0 = data.table(week = 1:53, rel.R0 = fn_getRelR0(loc.t, ref.wk = da.t$week[1], Rwea_parm.bounds=Rwea_parm.bounds) %>% rowMeans ) 
  
  da.t$rel.mob = rel.mob
  da.t = merge(da.t, relR0, by = 'week')
  da.t = da.t[order(date)]
  da2.t = da.t; da2.t$death = da2.t$death * 10
  da2.t = melt(da2.t[,c('date','case','death')], id.vars = 'date')
    
  ymin = min(da.t$rel.mob, da.t$rel.R0) * .9; ymax = max(da.t$rel.mob, da.t$rel.R0) * 1.1
  xx = barplot(value ~ variable + date, ylim = c(0, max(da.t$case) * 1.1), data = da2.t, xaxt = 'n', ylab = '', beside = T, border = 'transparent', col = c('grey50', 'red'))
  # xx = barplot(da.t$case, ylim = c(0, max(da.t$case) * 1.1), ylab = '', border = 'grey50', col = 'transparent')
  # barplot(da.t$death, ylim = c(0, max(da.t$case) * 1.1), add =T, beside = T, col = 'red', border = NA)
  mtext(cnty[cnt], side = 3, line = -1.2, cex = .8, adj = 0.01 , outer = F)
  
  par(new = T)
  xx = colMeans(xx)
  plot(xx, da.t$rel.mob, ylim = c(ymin, ymax), type = 'l', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = 'blue', lwd = 1.5)
  lines(xx, da.t$rel.R0, ylim = c(ymin, ymax), col = 'orange', lwd = 1.5)
  abline(h = 1, col = 'grey50', lty=2)
  axis(4)
  
  
  if(cnt == 1)
    legend('top', legend = c('case', 'death', 'mobility', 'seasonality'), pch = c(15, 15, NA, NA), seg.len = .8, ncol=2,
           lty = c(NA, NA, 1,1), col = c('grey50', 'red', 'blue', 'orange'), cex = 1, lwd = 2, bty = 'n')
  
  if(cnt == 3){
    axis(1, at = xx, labels = da.t$date %>% format('%m/%d/%y'))
    mtext('Case per million / Death per 100,000', side = 2, cex = .9, outer = T, line = 0)
    mtext('Relative mobility / Estimated Seasonal trend', side = 4, cex = .9, outer = T, line = 0)
  }
  
}
dev.off()

# plot modeling results for the three countries
# load results
load(paste0(dir_res, "res.summary.RData"))
res.train = res.train[eval == eval.t]
tda = res.train[loc == 'uk']
tda1 = tda[state == 'Infectious']
locs = c('uk', 'sa', 'br')
# plot model fit and overlay with estimated cumulative infection rate?



p.titles = c('(A) UK: Estimated transmissibility', '(B) UK: Estimated susceptibility',
             '(C) South Africa: Estimated transmissibility', '(D) South Africa: Estimated susceptibility',
             '(E) Brazil: Estimated transmissibility', '(F) Brazil: Estimated susceptibility'
)


p.titles2 = c('(A) UK: Rt v. infection rate', '(B) UK: Transmissibility', '(C) UK: Susceptibility', # v. infection rate
              '(D) South Africa: Rt v. infection rate', '(E) South Africa: Transmissibility', '(F) South Africa: Susceptibility',
              '(G) Brazil: Rt v. infection rate', '(H) Brazil: Transmissibility', '(I) Brazil: Susceptibility'
)

# plot model fit vs Rt/Infection rate
p.titles3 = c('(A) Model-fit: UK', '(B) Model-fit: South Africa','(C) Model-fit: Brazil', 
              '(D) Estimates: UK',
              '(E) Estimates: South Africa',
              '(F) Estimates: Brazil'
)

events = read.csv(paste0(dir_data,'events.csv'), stringsAsFactors = F) %>% data.table()
events$start = events$start %>% as.Date(format = '%m/%d/%y')
events$end = events$end %>% as.Date(format = '%m/%d/%y')

p.titles4 = c('(A) UK: Model-fit', 
              '(C) South Africa: Model-fit', 
              '(E) Brazil: Model-fit', 
              '(B) UK: Estimated Rt and infection rate',
              '(D) South Africa: Estimated Rt and infection rate',
              '(F) Brazil: Estimated Rt and infection rate'
)
# validation using independent data
# UK: REACT study, 10 rounds of testing
vda.uk = read.csv(paste0(dir_data, 'uk_react.study.csv'), stringsAsFactors = F, header = F) %>% data.table()
vda.uk = vda.uk[c(1:3,11),] %>% t  %>% data.table()
colnames(vda.uk)= vda.uk[1] %>% unlist
vda.uk = vda.uk[-1,]
setnames(vda.uk, c("First sample", "Last sample", "Weighted prevalence (95% CI)"), c('start','end','wt.prev'))
vda.uk$start = vda.uk$start %>% as.Date(format = '%m/%d/%y')
vda.uk$end = vda.uk$end %>% as.Date(format = '%m/%d/%y')
tmp = vda.uk$wt.prev %>% as.character() %>% strsplit('%') %>% unlist %>% matrix(nrow=nrow(vda.uk), ncol=4, byrow = T)
tmp[,1] = tmp[,1] %>% trimws()
tmp[,2] = gsub('\\(','',tmp[,2]) %>% trimws()
tmp[,3] = gsub(',','',tmp[,3]) %>% trimws()
tmp = tmp[,1:3]; mode(tmp) = 'numeric'
vda.uk$mean = tmp[,1]; vda.uk$ci95.lwr = tmp[,2]; vda.uk$ci95.upr = tmp[,3];

est.uk = res.train[loc == 'uk' & state == 'Infectious'] %>% setnames(c('variable','Week.start'),c('stat','date'))
est.uk$value = est.uk$value / N * 100 # %
est.uk = dcast(est.uk, date ~ stat, value.var = 'value')
dates.uk.t = unique(est.uk$date) %>% as.Date %>% sort
# asign the data to the closest week as the model
obs.uk = data.table(date = dates.uk.t, mean = -1, ci95.lwr = -1, ci95.upr = -1)
for(i in 1:nrow(obs.uk)){
  # if it's with the sampling period, assign that obs
  dt = obs.uk[i]$date; # week start
  for(j in 1:nrow(vda.uk)){
    if(dt %in% seq(as.Date(vda.uk[j]$start), as.Date(vda.uk[j]$end - 6), by = 'day')){
      obs.uk[i, c('mean','ci95.lwr','ci95.upr')] = vda.uk[j, c('mean','ci95.lwr','ci95.upr'),with=F]
    }
  }
}
obs.uk[obs.uk<0] = NA

# for brazil
vda.br = read.csv(paste0(dir_data, 'br_serosurveys.csv'), stringsAsFactors = F) %>% data.table()
vda.br$start = vda.br$start %>% as.Date(format = '%m/%d/%y')
vda.br$end = vda.br$end %>% as.Date(format = '%m/%d/%y')
tmp = vda.br$wt.prev %>% as.character() %>% strsplit('%') %>% unlist %>% matrix(nrow=nrow(vda.br), ncol=4, byrow = T)
tmp[,1] = tmp[,1] %>% trimws()
tmp[,2] = gsub('\\(','',tmp[,2]) %>% trimws()
tmp[,3] = gsub(',','',tmp[,3]) %>% trimws()
tmp = tmp[,1:3]; mode(tmp) = 'numeric'
vda.br$mean = tmp[,1]; vda.br$ci95.lwr = tmp[,2]; vda.br$ci95.upr = tmp[,3];

# first wave only
est.br = res.train[loc == 'br' & Week.start <= as.Date('2020/11/1') & state == 'Cumulative infection rate'] %>% setnames(c('variable','Week.start'),c('stat','date'))
est.br = dcast(est.br, date ~ stat, value.var = 'value')
dates.br.t = unique(est.br$date) %>% as.Date %>% sort
# asign the data to the closest week as the model
obs.br = data.table(date = dates.br.t, mean = -1, ci95.lwr = -1, ci95.upr = -1)
for(i in 1:nrow(obs.br)){
  # if it's with the sampling period, assign that obs
  dt = obs.br[i]$date; # week start
  for(j in 1:nrow(vda.br)){
    # account for delay in antibody generation, and the model date is week start
    if(dt %in% seq(as.Date(vda.br[j]$start-14 -6), as.Date(vda.br[j]$end -14 - 6), by = 'day')){
      obs.br[i, c('mean','ci95.lwr','ci95.upr')] = vda.br[j, c('mean','ci95.lwr','ci95.upr'),with=F]
    }
  }
}
obs.br[obs.br<0] = NA

# for SA
# vda.sa = read.csv(paste0(dir_plot, 'sa_serosurveys.csv'), stringsAsFactors = F) %>% data.table()
vda.sa = read_xlsx(paste0(dir_data, 'sa_serosurveys.xlsx'), sheet = 1) %>% data.table()
vda.sa$start = vda.sa$start %>% as.Date(format = '%m/%d/%y')
vda.sa$end = vda.sa$end %>% as.Date(format = '%m/%d/%y')
tmp = vda.sa$wt.prev %>% as.character() %>% strsplit('%') %>% unlist %>% matrix(nrow=nrow(vda.sa), ncol=4, byrow = T)
tmp[,1] = tmp[,1] %>% trimws()
tmp[,2] = gsub('\\(','',tmp[,2]) %>% trimws()
tmp[,3] = gsub(',','',tmp[,3]) %>% trimws()
tmp = tmp[,1:3,drop=F]; mode(tmp) = 'numeric'
vda.sa$mean = tmp[,1]; vda.sa$ci95.lwr = tmp[,2]; vda.sa$ci95.upr = tmp[,3];

# first wave only
est.sa = res.train[loc == 'sa' & Week.start <= as.Date('2020/9/20') & state == 'Cumulative infection rate'] %>% setnames(c('variable','Week.start'),c('stat','date'))
est.sa = dcast(est.sa, date ~ stat, value.var = 'value')
dates.sa.t = unique(est.sa$date) %>% as.Date %>% sort
# asign the data to the closest week as the model
obs.sa = data.table(date = dates.sa.t, mean = -1, ci95.lwr = -1, ci95.upr = -1)
for(i in 1:nrow(obs.sa)){
  # if it's with the sampling period, assign that obs
  dt = obs.sa[i]$date; # week start
  for(j in 1:nrow(vda.sa)){
    # account for delay in antibody generation, and the model date is week start
    if(dt %in% seq(as.Date(vda.sa[j]$start-14 -6), as.Date(vda.sa[j]$end -14 - 6), by = 'day')){
      # obs.sa[i, c('mean','ci95.lwr','ci95.upr')] = vda.sa[j, c('mean','ci95.lwr','ci95.upr'),with=F]
      obs.sa[i]$mean = vda.sa[j]$mean
    }
  }
}
obs.sa[obs.sa<0] = NA
# seperate the two datasets
obs.sa1 = obs.sa[date <= as.Date(vda.sa$end[1])-14]
obs.sa2 = obs.sa[date > as.Date(vda.sa$end[1])-14]


# put model fit and validation together
pdf(paste0(dir_plot, 'Fig2_modelfit_validation.pdf'), width = 8, height = 7)
par(mfcol=c(3,2),mar=c(1.7,2.5,1.5,2),oma=c(0.1,.1, 0.1,0),mgp=c(1.0,.1,0),cex=.8,cex.axis=.85,cex.lab=.9,tck=-.015)
cnt=0
for(loc.t in locs){
  cnt=cnt+1
  
  obs = da[,c('date','year','week',paste0(c('case.','death.','mob.bus.','mob.full.'),loc.t)),with=F] %>% 
    setnames(paste0(c('case.','death.','mob.bus.','mob.full.'),loc.t), c('case','death','mob.bus','mob.full'))
  
  obs$date = obs$date %>% as.Date
  # setnames(obs,'case', 'value')
  obs = obs[case > 2]
  
  # ymax.case = max(da[,paste0('case.',locs),with=F],na.rm = T) * 1.2
  # ymax.death = max(da[,paste0('death.',locs),with=F],na.rm = T) * 1.2
  
  dates.t = unique(obs$date)  %>% as.Date %>% sort
  x=seq(1, length.out = length(dates.t)*2, by = 1)
  
  mm = res.train[loc == loc.t & state == 'case'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  stats = matrix(0, 5, length(dates.t))
  for(id in 1:length(dates.t)){
    tmp = mm[date==dates.t[id]]
    stats[,id] = c(tmp[stat=='ci95.lwr']$value, tmp[stat=='iqr.lwr']$value, tmp[stat=='mean']$value, tmp[stat=='iqr.upr']$value, tmp[stat=='ci95.upr']$value)
  }
  colnames(stats) = dates.t
  
  # include both case and death
  mm.case = res.train[loc == loc.t & state == 'case'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  mm.death = res.train[loc == loc.t & state == 'death'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  mm.death$value = mm.death$value * 10
  ymax.case = max(obs$case, mm.case$value, mm.death$value, na.rm = T)*1.05
  x=seq(1, length.out = length(dates.t)*2, by = 1)
  stats = matrix(0, 5, length(dates.t)*2)
  id2 = 0
  for(id in seq(1, length.out = length(dates.t), by = 2)){
    id2 = id2 + 1
    tmp1 = mm.case[date==dates.t[id2]]
    tmp2 = mm.death[date==dates.t[id2]]
    stats[,id] = c(tmp1[stat=='ci95.lwr']$value, tmp1[stat=='iqr.lwr']$value, tmp1[stat=='mean']$value, tmp1[stat=='iqr.upr']$value, tmp1[stat=='ci95.upr']$value)
    stats[,id+1] = c(tmp2[stat=='ci95.lwr']$value, tmp2[stat=='iqr.lwr']$value, tmp2[stat=='mean']$value, tmp2[stat=='iqr.upr']$value, tmp2[stat=='ci95.upr']$value)
  }
  # colnames(stats) = dates.t
  x=seq(1, length.out = length(dates.t)*2, by = 1)
  summarydata=list(stats=stats,n=rep(dates.t, e=2),names=rep('',length(dates.t)*2))
  bxp(summarydata, box.width = .1, lwd=.5, ylab='', xaxt='n', yaxt = 'n', at=x, ylim = c(0, ymax.case), 
      xlim=c(0.5,length(dates.t)*2+.5), border = c('blue','red'), fill='transparent') # 
  points(x[seq(1,length.out = length(dates.t),by=2)], obs$case, pch = 8, cex = .4, col='blue')
  points(x[seq(2,length.out = length(dates.t),by=2)], obs$death*10, pch = 8, cex = .4, col='red')
  axis(1,at=x[seq(1.5,length.out = length(dates.t),by=2)],labels = format(dates.t,'%m/%d/%y'),mgp=c(1.0,.1,0),cex.axis=.85)
  axis(2,mgp=c(1.0,.1,0),cex.axis=.85, col.ticks = 'blue', col.lab = 'blue', col.axis='blue')
  mtext('Cases per million', side=2, outer = F, line = .9, cex=.75, col = 'blue')
  axis(4,mgp=c(1.0,.1,0),cex.axis=.85, col.ticks = 'red', col.lab = 'red', col.axis='red')
  mtext('Deaths per 100,000', side=4, outer = F, line = .9, cex=.75, col = 'red')
  
  
  # add time lines
  ymax = ymax.case
  events.t = events[loc == loc.t]
  events.t[is.na(end)]$end = max(dates.t)
  offset= -.8*2; y.offset = ymax*.006
  z=1; z.width = .02; 
  e.t = events.t[type == 'npi']
  if(nrow(e.t)>0){
    for(i in 1:nrow(e.t)){
      d.t1 = x[which.min(abs(as.Date(dates.t) - e.t[i]$start))*2]; 
      d.t2 = x[which.min(abs(as.Date(dates.t) - e.t[i]$end))*2];
      d.t3 = mean(c(d.t1, d.t2))
      tx.t = e.t[i]$event
      rect(xleft=d.t1, ybottom=-1000, xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey',.2), border = 'transparent') # ymax*(z-z.width*2)
      # arrows(x0=d.t, y0=ymax*(z-z.width*(i*2+1)), x1=max(x), length=.05)
      text(d.t3,ymax-y.offset, tx.t, adj = .5, cex=.75,srt=0, font = 1); # pos=4,
    }
  }
  e.t = events.t[type == 'virus']
  if(nrow(e.t)>0){
    for(i in 1:nrow(e.t)){
      d.t1 = x[which.min(abs(as.Date(dates.t) - e.t[i]$start))*2]; 
      d.t2 = x[which.min(abs(as.Date(dates.t) - e.t[i]$end))*2];
      d.t3 = mean(c(d.t1, d.t2))
      tx.t = e.t[i]$event
      # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
      if(loc.t == 'uk'){
        arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
        text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.75,srt=0, font = 1); # pos=4,
      } else {
        arrows(x0=d.t1, y0=ymax*(z-z.width*2), x1=d.t2, length=.05)
        text(d.t1+offset*2,ymax*(z-z.width*.5), tx.t,  pos=4, cex=.75,srt=0, font = 1); # pos=4,
      }
      
    }
  }
  e.t = events.t[type == 'vx']
  if(nrow(e.t)>0){
    for(i in 1:nrow(e.t)){
      d.t1 = x[which.min(abs(as.Date(dates.t) - e.t[i]$start))*2]; 
      d.t2 = x[which.min(abs(as.Date(dates.t) - e.t[i]$end))*2];
      d.t3 = mean(c(d.t1, d.t2))
      tx.t = e.t[i]$event
      # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
      if(loc.t == 'uk'){
        arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
        text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.75,srt=0, font = 1); # pos=4,
      } else {
        arrows(x0=d.t1, y0=ymax*(z-z.width*(4+2)), x1=d.t2, length=.05)
        text(d.t1+offset*2,ymax*(z-z.width*4.5), tx.t,  pos=4, cex=.75,srt=0, font = 1); # pos=4,
      }
      
    }
  }
  
  mtext(p.titles4[cnt],cex=.85,side=3,outer = F,line=.1,adj=0)
}
# validation
{
  # for the uk
  x = 1: length(dates.uk.t); ymax = max(est.uk$ci95.upr, obs.uk$ci95.upr, na.rm = T)*1.05
  plot(x,est.uk$mean,ylab='Weekly prevalence of infection (%)', ylim=c(0, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  lines(x,est.uk$mean, col='blue',lwd=2)
  polygon(c(x,rev(x)),c(est.uk$ci95.lwr,rev(est.uk$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(est.uk$iqr.lwr,rev(est.uk$iqr.upr)),col=alpha('blue',.3),border='transparent')
  points(x, obs.uk$mean, pch = "*", lwd = 2, cex = 1, col='red')
  arrows(x, obs.uk$ci95.lwr, x, obs.uk$ci95.upr, length=0.03, angle=90, code=3, col = 'red', lwd = 1)
  axis(1,at=x,labels = format(dates.uk.t,'%m/%d/%y'),mgp=c(1.0,.1,0),cex.axis=.85)
  legend('topleft', c('Model estimates','Data from the REACT study\n(Riley et al. 2021)'), seg.len = .8, cex = .8, lty = 1, lwd = 2, col = c('blue','red'), bty = 'n')
  mtext('(B) UK: Model validation', line = .1, adj = 0, cex = .85)
  # for SA
  x1 = pmatch(obs.sa1$date, dates.sa.t)
  x2 = pmatch(obs.sa2$date, dates.sa.t)
  x = 1: length(dates.sa.t); ymax = max(est.sa$ci95.upr, obs.sa$ci95.upr, na.rm = T)*1.05
  plot(x,est.sa$mean,ylab='Cumulative infection rate (%)', ylim=c(0, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  lines(x,est.sa$mean, col='blue',lwd=2)
  polygon(c(x,rev(x)),c(est.sa$ci95.lwr,rev(est.sa$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(est.sa$iqr.lwr,rev(est.sa$iqr.upr)),col=alpha('blue',.3),border='transparent')
  points(x1, obs.sa1$mean, pch = '*', lwd = 2, cex = 2, col='red')
  # for the trial
  points(x2, obs.sa2$mean, pch = '*', lwd = 2, cex = 2, col='brown')
  # arrows(x, obs.sa$ci95.lwr, x, obs.sa$ci95.upr, length=0.03, angle=90, code=3, col = 'red', lwd = 1)
  axis(1,at=x,labels = format(dates.sa.t,'%m/%d/%y'),mgp=c(1.0,.1,0),cex.axis=.85)
  mtext('(D) South Africa: Model validation', line = .1, adj = 0, cex = .85)
  legend('topleft', c('Model estimates','Data from a serosurvey in Cape Town (Shaw et al. 2021)',
                      'Data from a vaccine trial (Shinde et al. 2021)'), 
         seg.len = .8, cex = .8, lty = c(1, NA, NA), pch = c(NA, '*', '*'), pt.cex = 2, lwd = 2, col = c('blue','red','brown'), bty = 'n')
  # for br
  x = 1: length(dates.br.t); ymax = max(est.br$ci95.upr, obs.br$ci95.upr, na.rm = T)*1.05
  plot(x,est.br$mean,ylab='Cumulative infection rate (%)', ylim=c(0, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  lines(x,est.br$mean, col='blue',lwd=2)
  polygon(c(x,rev(x)),c(est.br$ci95.lwr,rev(est.br$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(est.br$iqr.lwr,rev(est.br$iqr.upr)),col=alpha('blue',.3),border='transparent')
  points(x, obs.br$mean, pch = '*', lwd = 2, cex = 1, col='red')
  arrows(x, obs.br$ci95.lwr, x, obs.br$ci95.upr, length=0.03, angle=90, code=3, col = 'red', lwd = 1)
  axis(1,at=x,labels = format(dates.br.t,'%m/%d/%y'),mgp=c(1.0,.1,0),cex.axis=.85)
  mtext('(F) Brazil: Model validation', line = .1, adj = 0, cex = .85)
  legend('topleft', c('Model estimates','Data from nationwide serosurveys\n(Hallal et al. 2020)'), seg.len = .8, cex = .8, lty = 1, lwd = 2, col = c('blue','red'), bty = 'n')
  
}
dev.off()

# put all key estimates together

pdf(paste0(dir_plot, 'Fig3_key.est.pdf'), width = 10.5, height = 7)
par(mfrow=c(3,3),mar=c(1.7,2.1,1.3,2.1),oma=c(0,.1, 0,.1),mgp=c(.9,.1,0),cex=.8,cex.axis=.85,cex.lab=.9,tck=-.01)
cnt=0
for(loc.t in locs){
  
  # do estimated infections and Rt
  cnt=cnt+1
  
  mm = res.train[loc == loc.t & state == 'infection'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  dates.t = unique(as.Date(mm$date)) %>% sort
  x=1:length(dates.t)
  stats = matrix(0, 5, length(dates.t))
  for(id in 1:length(dates.t)){
    tmp = mm[date==dates.t[id]]
    stats[,id] = c(tmp[stat=='ci95.lwr']$value, tmp[stat=='iqr.lwr']$value, tmp[stat=='mean']$value, tmp[stat=='iqr.upr']$value, tmp[stat=='ci95.upr']$value)
  }
  colnames(stats) = dates.t
  ymax=max(mm$value, na.rm = T)*1.1;
  summarydata=list(stats=stats,n=dates.t,names=rep('',length(dates.t)))
  bxp(summarydata, ylab='', yaxt='n', xaxt='n', at=x, ylim = c(0, ymax), border = 'grey50', fill='transparent') # xlim=c(0.5,length(dates.t)+.5),
  axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  axis(4,mgp=c(.9,.1,0),cex.axis=.85, col.ticks = 'grey30', col.lab = 'grey30', col.axis='grey30')
  mtext('Estimated infections per million', col = 'grey30',side=4,line = .9, outer = F,cex=.75)
  
  if(T){
    # add time lines
    events.t = events[loc == loc.t]
    events.t[is.na(end)]$end = max(dates.t)
    offset= -2; y.offset = ymax*.006
    z=1; z.width = .02; 
    e.t = events.t[type == 'npi']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        rect(xleft=d.t1, ybottom=-1000, xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey',.2), border = 'transparent') # ymax*(z-z.width*2)
        # arrows(x0=d.t, y0=ymax*(z-z.width*(i*2+1)), x1=max(x), length=.05)
        text(d.t3,ymax-y.offset, tx.t, adj = .5, cex=.7,srt=0, font = 1); # pos=4,
      }
    }
    e.t = events.t[type == 'virus']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*2), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
    e.t = events.t[type == 'vx']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+2)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
  }
  
  # overlay with cumulative infection rate
  par(new=T)
  tda=res.train[loc == loc.t & state == 'Rt'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  tda=dcast(tda, date ~ stat, value.var = 'value')
  ymax=max(tda$ci95.upr)*1.1; 
  
  plot(x,tda$mean,ylab='', yaxt='n',ylim=c(0, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  abline(h = 1, lty = 1, col = 'grey')
  lines(x,tda$mean, col='blue',lwd=2)
  polygon(c(x,rev(x)),c(tda$ci95.lwr,rev(tda$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(tda$iqr.lwr,rev(tda$iqr.upr)),col=alpha('blue',.3),border='transparent')
  axis(2,mgp=c(1.0,.1,0),cex.axis=.85, col.ticks = 'blue', col.lab = 'blue', col.axis='blue')
  mtext('Rt', side=2, outer = F, line = .9, cex=.75, col = 'blue')
  mtext(p.titles2[cnt],cex=.85,side=3,outer = F,line=.1,adj=0)
  
  
  # midle panel
  # do estimated infections and susceptibility
  cnt=cnt+1
  x=1:length(dates.t)
  mm = res.train[loc == loc.t & state == 'infection'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  stats = matrix(0, 5, length(dates.t))
  for(id in 1:length(dates.t)){
    tmp = mm[date==dates.t[id]]
    stats[,id] = c(tmp[stat=='ci95.lwr']$value, tmp[stat=='iqr.lwr']$value, tmp[stat=='mean']$value, tmp[stat=='iqr.upr']$value, tmp[stat=='ci95.upr']$value)
  }
  colnames(stats) = dates.t
  ymax=max(mm$value, na.rm = T)*1.05;
  summarydata=list(stats=stats,n=dates.t,names=rep('',length(dates.t)))
  bxp(summarydata, ylab='', yaxt='n', xaxt='n', at=x, ylim = c(0, ymax), border = 'grey50', fill='transparent') # xlim=c(0.5,length(dates.t)+.5),
  axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  axis(4,mgp=c(.9,.1,0),cex.axis=.85, col.ticks = 'grey30', col.lab = 'grey30', col.axis='grey30')
  mtext('Estimated infections per million', col = 'grey30',side=4,line = .9, outer = F,cex=.75)
  if(T){
    # add time lines
    events.t = events[loc == loc.t]
    events.t[is.na(end)]$end = max(dates.t)
    offset= -2; y.offset = ymax*.006
    z=1; z.width = .02; 
    e.t = events.t[type == 'npi']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        rect(xleft=d.t1, ybottom=-1000, xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey',.2), border = 'transparent') # ymax*(z-z.width*2)
        # arrows(x0=d.t, y0=ymax*(z-z.width*(i*2+1)), x1=max(x), length=.05)
        text(d.t3,ymax-y.offset, tx.t, adj = .5, cex=.7,srt=0, font = 1); # pos=4,
      }
    }
    e.t = events.t[type == 'virus']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*2), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
    e.t = events.t[type == 'vx']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+2)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
  }
  
  # overlay with cumulative infection rate
  par(new=T)
  tda=res.train[loc == loc.t & state == 'Rtx'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  tda=dcast(tda, date ~ stat, value.var = 'value')
  ymax=max(tda$ci95.upr)*1.05; 
  
  plot(x,tda$mean,ylab='Estimated transmissibility',ylim=c(-.5, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  polygon(c(x,rev(x)),c(tda$ci95.lwr,rev(tda$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(tda$iqr.lwr,rev(tda$iqr.upr)),col=alpha('blue',.3),border='transparent')
  # axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  
  mtext(p.titles2[cnt],cex=.85,side=3,outer = F,line=.1,adj=0)
  
  # right panel
  # do estimated infections and susceptibility
  cnt=cnt+1
  x=1:length(dates.t)
  mm = res.train[loc == loc.t & state == 'infection'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  stats = matrix(0, 5, length(dates.t))
  for(id in 1:length(dates.t)){
    tmp = mm[date==dates.t[id]]
    stats[,id] = c(tmp[stat=='ci95.lwr']$value, tmp[stat=='iqr.lwr']$value, tmp[stat=='mean']$value, tmp[stat=='iqr.upr']$value, tmp[stat=='ci95.upr']$value)
  }
  colnames(stats) = dates.t
  ymax=max(mm$value, na.rm = T)*1.05;
  summarydata=list(stats=stats,n=dates.t,names=rep('',length(dates.t)))
  bxp(summarydata, ylab='', yaxt='n', xaxt='n', at=x, ylim = c(0, ymax), border = 'grey50', fill='transparent') # xlim=c(0.5,length(dates.t)+.5),
  axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  axis(4,mgp=c(.9,.1,0),cex.axis=.85, col.ticks = 'grey30', col.lab = 'grey30', col.axis='grey30')
  mtext('Estimated infections per million', col = 'grey30',side=4,line = .9, outer = F,cex=.75)
  if(T){
    # add time lines
    events.t = events[loc == loc.t]
    events.t[is.na(end)]$end = max(dates.t)
    offset= -2; y.offset = ymax*.006
    z=1; z.width = .02; 
    e.t = events.t[type == 'npi']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        rect(xleft=d.t1, ybottom=-1000, xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey',.2), border = 'transparent') # ymax*(z-z.width*2)
        # arrows(x0=d.t, y0=ymax*(z-z.width*(i*2+1)), x1=max(x), length=.05)
        text(d.t3,ymax-y.offset, tx.t, adj = .5, cex=.7,srt=0, font = 1); # pos=4,
      }
    }
    e.t = events.t[type == 'virus']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*2), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
    e.t = events.t[type == 'vx']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+2)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
  }
  
  # overlay with cumulative infection rate
  par(new=T)
  tda=res.train[loc == loc.t & state == 'Susceptibility'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  tda=dcast(tda, date ~ stat, value.var = 'value')
  ymax=max(tda$ci95.upr)*1.05; 
  
  plot(x,tda$mean,ylab='Estimated susceptibility (%)',ylim=c(0, 100),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  polygon(c(x,rev(x)),c(tda$ci95.lwr,rev(tda$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(tda$iqr.lwr,rev(tda$iqr.upr)),col=alpha('blue',.3),border='transparent')
  # axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  
  mtext(p.titles2[cnt],cex=.85,side=3,outer = F,line=.1,adj=0)
}
dev.off()


# OTHER KEY PARAM
p.titles5 = c('(A) UK: Infection-detection rate', '(B) UK: Infection-fatality risk', # v. infection rate
              '(C) South Africa: Infection-detection rate', '(D) South Africa: Infection-fatality risk',
              '(E) Brazil: Infection-detection rate', '(F) Brazil: Infection-fatality risk'
)
pdf(paste0(dir_plot, 'FigS3_other.parm.est.pdf'), width = 7, height = 7)
par(mfrow=c(3,2),mar=c(1.7,2.1,1.3,2.1),oma=c(0,.1, 0,.1),mgp=c(.9,.1,0),cex=.8,cex.axis=.85,cex.lab=.9,tck=-.01)
cnt=0
for(loc.t in locs){
  
  # do estimated infections and Rt
  cnt=cnt+1
  
  mm = res.train[loc == loc.t & state == 'infection'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  dates.t = unique(as.Date(mm$date)) %>% sort
  x=1:length(dates.t)
  stats = matrix(0, 5, length(dates.t))
  for(id in 1:length(dates.t)){
    tmp = mm[date==dates.t[id]]
    stats[,id] = c(tmp[stat=='ci95.lwr']$value, tmp[stat=='iqr.lwr']$value, tmp[stat=='mean']$value, tmp[stat=='iqr.upr']$value, tmp[stat=='ci95.upr']$value)
  }
  colnames(stats) = dates.t
  ymax=max(mm$value, na.rm = T)*1.1;
  summarydata=list(stats=stats,n=dates.t,names=rep('',length(dates.t)))
  bxp(summarydata, ylab='', yaxt='n', xaxt='n', at=x, ylim = c(0, ymax), border = 'grey50', fill='transparent') # xlim=c(0.5,length(dates.t)+.5),
  axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  axis(4,mgp=c(.9,.1,0),cex.axis=.85, col.ticks = 'grey30', col.lab = 'grey30', col.axis='grey30')
  mtext('Estimated infections per million', col = 'grey30',side=4,line = .9, outer = F,cex=.75)
  
  if(T){
    # add time lines
    events.t = events[loc == loc.t]
    events.t[is.na(end)]$end = max(dates.t)
    offset= -2; y.offset = ymax*.006
    z=1; z.width = .02; 
    e.t = events.t[type == 'npi']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        rect(xleft=d.t1, ybottom=-1000, xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey',.2), border = 'transparent') # ymax*(z-z.width*2)
        # arrows(x0=d.t, y0=ymax*(z-z.width*(i*2+1)), x1=max(x), length=.05)
        text(d.t3,ymax-y.offset, tx.t, adj = .5, cex=.7,srt=0, font = 1); # pos=4,
      }
    }
    e.t = events.t[type == 'virus']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*2), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
    e.t = events.t[type == 'vx']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+2)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
  }
  
  # overlay with cumulative infection rate
  par(new=T)
  tda=res.train[loc == loc.t & state == 'infection detection rate'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  tda$value = tda$value * 100
  tda=dcast(tda, date ~ stat, value.var = 'value')
  ymax=max(tda$ci95.upr)*1.1; 
  
  plot(x,tda$mean,ylab='', yaxt='n',ylim=c(0, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  abline(h = 1, lty = 1, col = 'grey')
  lines(x,tda$mean, col='blue',lwd=2)
  polygon(c(x,rev(x)),c(tda$ci95.lwr,rev(tda$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(tda$iqr.lwr,rev(tda$iqr.upr)),col=alpha('blue',.3),border='transparent')
  axis(2,mgp=c(1.0,.1,0),cex.axis=.85, col.ticks = 'blue', col.lab = 'blue', col.axis='blue')
  mtext('Infection-detection rate (%)', side=2, outer = F, line = .9, cex=.75, col = 'blue')
  mtext(p.titles5[cnt],cex=.85,side=3,outer = F,line=.1,adj=0)
  
  
  # midle panel
  # do estimated infections and susceptibility
  cnt=cnt+1
  x=1:length(dates.t)
  mm = res.train[loc == loc.t & state == 'infection'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  stats = matrix(0, 5, length(dates.t))
  for(id in 1:length(dates.t)){
    tmp = mm[date==dates.t[id]]
    stats[,id] = c(tmp[stat=='ci95.lwr']$value, tmp[stat=='iqr.lwr']$value, tmp[stat=='mean']$value, tmp[stat=='iqr.upr']$value, tmp[stat=='ci95.upr']$value)
  }
  colnames(stats) = dates.t
  ymax=max(mm$value, na.rm = T)*1.05;
  summarydata=list(stats=stats,n=dates.t,names=rep('',length(dates.t)))
  bxp(summarydata, ylab='', yaxt='n', xaxt='n', at=x, ylim = c(0, ymax), border = 'grey50', fill='transparent') # xlim=c(0.5,length(dates.t)+.5),
  axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  axis(4,mgp=c(.9,.1,0),cex.axis=.85, col.ticks = 'grey30', col.lab = 'grey30', col.axis='grey30')
  mtext('Estimated infections per million', col = 'grey30',side=4,line = .9, outer = F,cex=.75)
  if(T){
    # add time lines
    events.t = events[loc == loc.t]
    events.t[is.na(end)]$end = max(dates.t)
    offset= -2; y.offset = ymax*.006
    z=1; z.width = .02; 
    e.t = events.t[type == 'npi']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        rect(xleft=d.t1, ybottom=-1000, xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey',.2), border = 'transparent') # ymax*(z-z.width*2)
        # arrows(x0=d.t, y0=ymax*(z-z.width*(i*2+1)), x1=max(x), length=.05)
        text(d.t3,ymax-y.offset, tx.t, adj = .5, cex=.7,srt=0, font = 1); # pos=4,
      }
    }
    e.t = events.t[type == 'virus']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*2), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
    e.t = events.t[type == 'vx']
    if(nrow(e.t)>0){
      for(i in 1:nrow(e.t)){
        d.t1 = which.min(abs(as.Date(dates.t) - e.t[i]$start)); 
        d.t2 = which.min(abs(as.Date(dates.t) - e.t[i]$end));
        d.t3 = mean(c(d.t1, d.t2))
        tx.t = e.t[i]$event
        # rect(xleft=d.t1, ybottom=ymax*(z-z.width*2), xright=d.t2, ytop=ymax*(z+z.width*2), angle = 45,col = alpha('grey50',.3), border = 'transparent')
        # arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
        # text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        if(loc.t == 'uk'){
          arrows(x0=d.t1, y0=ymax*(z-z.width*(8+1.5)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*8), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        } else {
          arrows(x0=d.t1, y0=ymax*(z-z.width*(4+2)), x1=d.t2, length=.05)
          text(d.t1+offset*2,ymax*(z-z.width*4.5), tx.t,  pos=4, cex=.7,srt=0, font = 1); # pos=4,
        }
      }
    }
  }
  
  # overlay with cumulative infection rate
  par(new=T)
  tda=res.train[loc == loc.t & state == 'IFR'] %>% setnames(c('variable','Week.start'),c('stat','date'))
  tda$value = tda$value * 100
  tda=dcast(tda, date ~ stat, value.var = 'value')
  ymax=max(tda$ci95.upr)*1.05; 
  
  plot(x,tda$mean,ylab='Infection-fatality risk (%)',ylim=c(0, ymax),type='l',col='blue',lwd=2,xlab='',xaxt='n')
  polygon(c(x,rev(x)),c(tda$ci95.lwr,rev(tda$ci95.upr)),col=alpha('blue',.15),border='transparent')
  polygon(c(x,rev(x)),c(tda$iqr.lwr,rev(tda$iqr.upr)),col=alpha('blue',.3),border='transparent')
  # axis(1,at=x,labels = format(dates.t,'%m/%d/%y'),mgp=c(.9,.1,0),cex.axis=.85)
  
  mtext(p.titles5[cnt],cex=.85,side=3,outer = F,line=.1,adj=0)
  
}
dev.off()

# generate table
tab = read.csv(paste0(dir_res, "tab_voc.summary.est.csv")) %>% data.table()
tab = tab[eval == eval.t]
tab$cb = apply(tab[,c('mean', 'ci95.lwr', 'ci95.upr')], 1, FUN = fn_format, roundigt=1)
tab2 = dcast(tab, loc ~ state, value.var = 'cb')
setnames(tab2, c('loc', 'dRtx', 'dImm'), c('Location','Changes in Transmissibility (%)','Immune evasion (%)'))
tab2$Variant = factor(tab2$Location, levels = c('uk','sa','br'), labels = c('B.1.1.7', "B.1.351", 'P.1'))
tab2$Location = factor(tab2$Location, levels = c('uk','sa','br'), labels = c('United Kingdom', "South Africa", 'Brazil'))

tab2 = tab2[order(Location),  c('Location', 'Variant','Changes in Transmissibility (%)','Immune evasion (%)')]
write.csv(tab2, file = paste0(dir_plot, 'Tab_est.dRtx.dImm.csv'), row.names = F)


# plot weather variables in the three countries
wea1 = read.csv(paste0(dir_data,'wea.by.week_UK.csv')) %>% data.table()
wea2 = read.csv(paste0(dir_data,'wea.by.week_SA.csv')) %>% data.table()
wea3 = read.csv(paste0(dir_data,'wea.by.week_BR.csv')) %>% data.table()
wea1$Country = 'UK'
wea2$Country = 'South Africa'
wea3$Country = 'Brazil'
wea = rbind(wea1, wea2, wea3)
wea$Country = factor(wea$Country, levels = c('UK', 'South Africa', 'Brazil'))
wea = melt(wea[, c('Country','week','temp', 'spec.hum')], id.vars = c('Country','week'))
wea$variable = factor(wea$variable, levels = c('temp', 'spec.hum'), labels = c('Temperature (C)', 'Specific humidity (kg/kg)'))
wea = wea[week != 53]

theme.t = theme(plot.title = element_text(v=0, size = 11, margin=margin(0,0,2,0)), 
                strip.placement = "outside", strip.text = element_text(size = 11, margin=margin(2,0,2,0)),
                axis.title = element_text(size =10, margin=margin(0,0.2,0,0)), 
                axis.text.y = element_text(size=9, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=9,angle = 0),
                plot.margin=unit(c(c(.3, 1, 1, .5)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=9), legend.text=element_text(size=9),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                legend.key.size = unit(.3, 'cm'), #change legend key size
                legend.key.height = unit(.3, 'cm'), #change legend key height
                legend.key.width = unit(.3, 'cm')) #change legend key width)
p = ggplot(wea) + 
  geom_line(aes(x = week, y = value, color = Country), size = 1) +  # no ctrl
  facet_rep_wrap(~ variable, scales = 'free_y', 
                 repeat.tick.labels = T, ncol = 2) + 
  labs(x = 'Week of the Year', y = 'Weekly average') +
  theme_minimal() + theme.t

pdf(paste0(dir_plot, 'FigS_weather3countires.pdf'), width = 7, height = 3.5)
print(p)
dev.off()



