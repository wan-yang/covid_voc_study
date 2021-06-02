# script to process runs
# NOTE: WE DO NOT INCLUDE THE ORGINAL MODLE OUTPUTS HERE, B/C FILES ARE TOO LARGE
# PLEASE RUN THE DIRVER CODE ON YOUR LOCAL MACHINE AND THEN USE THIS FOR COMPILATION
# OTHERWISE, IT WON'T RUN
# HOWEVER, THE SUMMARY RESULTS ARE INCLUDED IN THE RESUTLS FOLDER

# end1stWave = 30; # end of first wave, after which S can be updated by the filter
num_runs = 100
stoch = T

tag.model = ''
# tag.evals = c('eq','rank', 'obs.more', 'obs.most', 'obs.only') # 'eq','rank', 'obs.only' - not good

# tag.eval = 'rank'



# for(tag.eval in tag.evals){
#  print(paste('start', tag.eval, Sys.time()))

if(T){
  dir_data = '../../data/'
  dir_code = './'
  dir_code1 = '../model_inference/'
  dir_code2 =  '../plot_results/' 
  
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
  
  # write.table(1:50, paste0(dir_code,'RUN_NO'), col.names = F, row.names = F, sep='/n')
}


# scripts

source(paste0(dir_code1,'SEIRS.R'))
source(paste0(dir_code1,'EAKF.R'))
source(paste0(dir_code2,'getPlot.R'))
source(paste0(dir_code1,'get_relR0.R'))


if(! file.exists(dir_res))  dir.create(dir_res,recursive = T)

da = read.csv(paste0(dir_data,'da_case.death.mob_uk.sa.br.csv'), stringsAsFactors = F)  %>% data.table()
da[date < as.Date('2020-02-09') & is.na(da)] = 0
# set population size to 20k? make it smaller for stochasticity
N = 1e6; # per 1 M
loc.t = 'uk'
num_gr = num_obs = length(N); # no age structure
num_ens = 500
epi.model = 'SEIRSV' # susceptible-exposed-infectious-recovered-susc
# for seeding: p.home.lwr - p.home.upr from home


num_runs = 100
stoch = T


seasonality = T; 

# parms from flu hk paper, AH/T model
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


source(paste0(dir_code1,'set_tm2event.R'))


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

locs = c('uk', 'sa','br')
loc.names = c('United Kingdom', 'South Africa', 'Brazil')

# get the combined results from multiple runs
source(paste0(dir_code,'getCombinedRes_diff.eval.R'))
res.train00 = res.train
save(newVstat, res.train, file = paste0(dir_res, 'res.summary.RData'))

for(tag.eval in tag.evals){
  
  res.train = res.train00[eval==tag.eval] %>% data.table()
  
  res.train = dcast(res.train, loc + state + Week.start ~ variable, value.var = 'value')
  res.train$Week.start = res.train$Week.start %>% as.Date
  
  Rt = res.train[state %in% c('Rt','Rtx')]
  Rt$obs = NA; Rt$threshold = 1;
  
  fits = res.train[state %in% c('case','death') ]
  fits$loc = factor(fits$loc, levels = locs, labels = loc.names)
  
  tda1 = melt(da[,c('date',paste0('case.',locs)),with=F], id.vars = 'date')
  tda1$variable = factor(tda1$variable, levels = paste0('case.',locs), labels = loc.names)
  setnames(tda1, c('date','variable','value'), c("Week.start",'loc','obs'))
  tda1$state = 'case'
  tda2 = melt(da[,c('date',paste0('death.',locs)),with=F], id.vars = 'date')
  tda2$variable = factor(tda2$variable, levels = paste0('death.',locs), labels = loc.names)
  setnames(tda2, c('date','variable','value'), c("Week.start",'loc','obs'))
  tda2$state = 'death'
  tda = rbind(tda1, tda2)
  tda$Week.start = tda$Week.start %>% as.Date
  fits = merge(fits, tda, by = c('loc', 'state', "Week.start"))
  # fits$threshold = NA;
  
  tda = fits
  tda$state = factor(tda$state, levels = c('case', 'death'), labels = c('case', 'death'))
  # tda$loc = factor(tda$loc, levels = paste0('sce',1:n.sce), labels = c('truth 1','truth 2','truth 3', 'truth n.sce'))
  pdf(paste0(dir_res, 'Fig_voc_model_fits', tail(da$Week.start,1),'_',tag.eval,'.pdf'), width = 8, height =length(locs) * 2)
  p = getPlot(tda)
  print(p)
  dev.off()
  
  
  res.train$loc = factor(res.train$loc, levels = locs, labels = loc.names)
  
  pdf(paste0(dir_res, 'Fig_voc_parms_', tail(da$Week.start,1),'_',tag.eval,'.pdf'), width = 12, height = 8)
  p = getPlotStates(res.train, var = c('Rt', 'Rtx','infection detection rate','IFR','Susceptibility'))
  print(p)
  dev.off()
  
  pdf(paste0(dir_res, 'Fig_voc_tx.suscept_', tail(da$Week.start,1),'_',tag.eval,'.pdf'), width = 12, height = length(locs) * 2)
  p = getPlotStates(res.train, var = c('Susceptibility','Rtx'))
  print(p)
  dev.off()
  
  
  # compute dRtx and dImm for individual runs and then combine
  
  res = newVstat[eval == tag.eval] %>% data.table()
  res = melt(res, id.vars = c('loc','run','hyp.test','eval')) 
  res = res[variable %in% c('perc.dImm.mn','perc.dRtx.mean')]
  res$state = factor(res$variable, levels = c('perc.dRtx.mean','perc.dImm.mn'), labels =  c('dRtx', 'dImm'))
  res$loc = factor(res$loc, levels = locs)
  res = res[order(loc, run, state)]

  res = res %>% data.table()

 # get the mean and plot it
  means = res[,list(value = mean(value)), by = c('loc','state')]
  
  n <- prod(dim(with(res, table(state, loc))))
  VEC <- seq(1, n/2, length.out=n)*2 - c(0, .2)
  pdf(paste0(dir_res, 'Fig_voc_est_tx_imm_change_',tag.eval,'.pdf'),width = 6, height = 3)
  par(mfrow=c(1,1), mar=c(2, 2.5, 1.5, .5), mgp = c(1.3, .4, 0), tck = -.02)
  boxplot(value ~ loc:state, data =res, ylab = 'Relative change (%)',
          boxwex=0.5, col=alpha(c("orange", "yellow"),0.2),
          xlab="", xaxt = 'n', outcex = 0.2,
          sep=":", lex.order=TRUE, ylim=c(-10, 110), yaxs="i",
          cex.axis=.8, at=VEC)
  x0s <- VEC - 0.2; x1s <- VEC + 0.2
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- means$value
  segments(x0 = x0s, x1 = x1s, y0 = y0s, lwd=1.5, col = "red")
  # points(VEC, means$value, pch = '----', cex=2, col = 'red')
  axis(1, at = VEC, labels = outer(c('tx','S'), locs, FUN = paste, sep=': ') %>% c)
  dev.off()
} # tag.eval

# put all together
pdf(paste0(dir_res, 'Fig_voc_est_tx_imm_change.pdf'),width = 6, height = length(tag.evals) * 2)
par(mfrow=c(length(tag.evals),1), mar=c(2, 2.5, 1.5, .5), mgp = c(1.2, .3, 0), cex=1, tck = -.03)
for(tag.eval in tag.evals){
  res = newVstat[eval == tag.eval] %>% data.table()
  res = melt(res, id.vars = c('loc','run','hyp.test','eval')) 
  res = res[variable %in% c('perc.dImm.mn','perc.dRtx.mean')]
  res$state = factor(res$variable, levels = c('perc.dRtx.mean','perc.dImm.mn'), labels =  c('dRtx', 'dImm'))
  res$loc = factor(res$loc, levels = locs)
  res = res[order(loc, run, state)]
  
  
  n <- prod(dim(with(res, table(state, loc))))
  VEC <- seq(1, n/2, length.out=n)*2 - c(0, .2)
  
  boxplot(value ~ loc:state, data =res, ylab = 'Relative change (%)',
          boxwex=0.5, col=alpha(c("orange", "yellow"),0.2),
          xlab="", xaxt = 'n', outcex = 0.2,
          sep=":", lex.order=TRUE, ylim=c(-10, 110), yaxs="i",
          cex.axis=.95, cex.lab = 1, at=VEC)
  x0s <- VEC - 0.2; x1s <- VEC + 0.2
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- means$value
  segments(x0 = x0s, x1 = x1s, y0 = y0s, lwd=1.5, col = "red")
  axis(1, at = VEC, cex=.95, labels = outer(c('tx','S'), locs, FUN = paste, sep=': ') %>% c)
  mtext(tag.eval,line = .1, adj = .05)
}
dev.off()

# output in a table
tab = NULL
for(tag.eval in tag.evals){
  res = newVstat[eval == tag.eval] %>% data.table()
  res = melt(res, id.vars = c('loc','run','hyp.test','eval')) 
  res = res[variable %in% c('perc.dImm.mn','perc.dRtx.mean')]
  res$state = factor(res$variable, levels = c('perc.dRtx.mean','perc.dImm.mn'), labels =  c('dRtx', 'dImm'))
  res$loc = factor(res$loc, levels = locs)
  res = res[order(loc, run, state)]
  
  for(loc.t in locs){
    tmp = res[loc == loc.t]
    tab = rbind(tab,
                data.table(eval = tag.eval, loc = loc.t, state = 'dRtx', 
                           mean = mean(tmp[state=='dRtx']$value) %>% round(2), sd = sd(tmp[state=='dRtx']$value)%>% round(2),
                           quantile(tmp[state=='dRtx']$value, probs = c(.5, .25, .75, .025, .975, .1, .9, .05, .95)) %>% round(2) %>% t),
                data.table(eval = tag.eval, loc = loc.t, state = 'dImm', 
                           mean = mean(tmp[state=='dImm']$value) %>% round(2), sd = sd(tmp[state=='dImm']$value) %>% round(2),
                           quantile(tmp[state=='dImm']$value, probs = c(.5, .25, .75, .025, .975, .1, .9, .05, .95)) %>% round(2) %>% t)
                )
  }
}
setnames(tab, c("50%", "25%","75%","2.5%","97.5%","10%", '90%',"5%", '95%'), c('median','iqr.lwr','iqr.upr','ci95.lwr','ci95.upr','ci80.lwr','ci80.upr','ci90.lwr','ci90.upr'))
write.csv(tab, paste0(dir_res, 'tab_voc.summary.est.csv'), row.names = F)

loc.t = 'uk'
tmp = newVstat[eval == 'obs.most' & loc == loc.t]
table(tmp$hyp.test)
tmp2 = newVstat[eval == 'obs.more' & loc == loc.t]
table(tmp2$hyp.test)

loc.t = 'sa'
tmp = newVstat[eval == 'obs.most' & loc == loc.t]
table(tmp$hyp.test)
tmp2 = newVstat[eval == 'obs.more' & loc == loc.t]
table(tmp2$hyp.test)

loc.t = 'br'
tmp = newVstat[eval == 'obs.most' & loc == loc.t]
table(tmp$hyp.test)
tmp = newVstat[eval == 'obs.more' & loc == loc.t]
table(tmp$hyp.test)
tmp = newVstat[eval == 'obs.comb' & loc == loc.t]
table(tmp$hyp.test)

hist(newVstat[eval == 'obs.more' & loc == loc.t]$perc.dImm.mn)
tmp1 = newVstat[eval == 'obs.more' & loc == loc.t]
tmp2 = newVstat[eval == 'obs.more' & loc == loc.t & perc.dImm.mn != 0]
summary(tmp2$perc.dImm.mn)
summary(tmp2$perc.dRtx.mean)
summary(tmp1$perc.dRtx.mean)

# get cumI
source(paste0(dir_code,'getCombinedRes_diff.eval_cumI.R'))
save(cumIperc_ens, file = paste0(dir_res, 'res.cumIperc_ens.RData'))


