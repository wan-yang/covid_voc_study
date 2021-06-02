# script to process runs
# NOTE: WE DO NOT INCLUDE THE ORGINAL MODLE OUTPUTS HERE, B/C FILES ARE TOO LARGE
# PLEASE RUN THE DIRVER CODE ON YOUR LOCAL MACHINE AND THEN USE THIS FOR COMPILATION
# OTHERWISE, IT WON'T RUN
# HOWEVER, THE SUMMARY RESULTS ARE INCLUDED IN THE RESUTLS FOLDER



num_runs = 100
n.sce = 5
seed1stWave = 2
seed2ndWave = 50

alpha.t = .2

# alphas = c(.1, .2)

# alpha.t = alphas[dummy_ialpha] # low detection rate

end1stWave = 20; # end of first wave, after which S can be updated by the filter

stoch = T

tag.model = ''
# tag.evals = c('eq','rank', 'obs.more', 'obs.most', 'obs.only') # 'eq','rank', 'obs.only' - not good
tag.evals = c('obs.more','obs.most','rank')
# tag.eval = 'rank'



# for(tag.eval in tag.evals){
#  print(paste('start', tag.eval, Sys.time()))

if(T){
  dir_data = '../../data/'
  dir_code = './'
  dir_code1 = '../model_inference/'
  dir_code2 =  '../plot_results/'  
  dir_res = paste0('../../results/syn.test/test_',seed1stWave,'seed2nd',seed2ndWave,'alpha',alpha.t,'/')
  
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



# scripts

source(paste0(dir_code1,'SEIRS.R'))
source(paste0(dir_code1,'EAKF.R'))
source(paste0(dir_code2,'getPlot.R'))
source(paste0(dir_code1,'get_relR0.R'))


if(! file.exists(dir_res))  dir.create(dir_res,recursive = T)


da = read.csv(paste0(dir_data,'da_syn.truth_seed1st',seed1stWave,'_seed2nd', seed2ndWave,'_alpha',alpha.t,'.csv'), stringsAsFactors = F)  %>% data.table()
da[date < as.Date('2020-02-09') & is.na(da)] = 0
tmp = MMWRweek(da$date)
da$year = tmp[,1]
da$week = tmp[,2]
# set population size to 20k? make it smaller for stochasticity
N = 1e6; # per 1 M
loc.t = 'uk'
num_gr = num_obs = length(N); # no age structure
num_ens = 500
epi.model = 'SEIRSV' # susceptible-exposed-infectious-recovered-susc
# for seeding: p.home.lwr - p.home.upr from home


seasonality = T; # 
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


seed = .05

# SR
doSR = T
percSRmajor = .08
SR.perc = .05
SR.perc.local = .03;
SR.perc.full = 0.03
SR.var.local= c('beta','Tei','Tir','Trs','Td.mean','Td.sd','p.mob','alpha','ifr')
SR.var.full= c('beta','Tir','Trs','ifr') # ,'alpha','ifr'
SR.var.tx = c('beta','Tir')
# do three diff categories: all combined, ms/mh, cumc
donotUpdateS1stWave = T # do not allow the filter to update S during first wave
rednUpdateEI = T # do not allow or reduce the level allowed, the filter to update E or I - it takes OEVr all efforts

locs = paste0('sce', 1:n.sce)
loc.names = paste('truth', 1:n.sce)

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
  fits$loc = factor(fits$loc, levels = paste0('sce',1:n.sce), labels = paste('truth',1:n.sce))
  
  tda1 = melt(da[,c('date',paste0('case.sce',1:n.sce)),with=F], id.vars = 'date')
  tda1$variable = factor(tda1$variable, levels = paste0('case.sce',1:n.sce), labels = paste('truth',1:n.sce))
  setnames(tda1, c('date','variable','value'), c("Week.start",'loc','obs'))
  tda1$state = 'case'
  tda2 = melt(da[,c('date',paste0('death.sce',1:n.sce)),with=F], id.vars = 'date')
  tda2$variable = factor(tda2$variable, levels = paste0('death.sce',1:n.sce), labels = paste('truth',1:n.sce))
  setnames(tda2, c('date','variable','value'), c("Week.start",'loc','obs'))
  tda2$state = 'death'
  tda = rbind(tda1, tda2)
  tda$Week.start = tda$Week.start %>% as.Date
  fits = merge(fits, tda, by = c('loc', 'state', "Week.start"))
  # fits$threshold = NA;
  
  tda = fits
  tda$state = factor(tda$state, levels = c('case', 'death'), labels = c('case', 'death'))
  # tda$loc = factor(tda$loc, levels = paste0('sce',1:n.sce), labels = c('truth 1','truth 2','truth 3', 'truth n.sce'))
  pdf(paste0(dir_res, 'a',itest, 'Fig_syn.test_model_fits_alpha',alpha.t,'_',tag.eval,'.pdf'), width = 8, height = n.sce * 2)
  p = getPlot(tda)
  print(p)
  dev.off()
  
  truth1 = melt(da[,c('date', paste0('wt.suscept.sce', 1:n.sce)),with=F], id.vars = 'date') %>% setnames('variable','loc')
  truth1$state = 'Susceptibility'; truth1$value = truth1$value / N * 100
  truth1$loc = factor(truth1$loc, levels = paste0('wt.suscept.sce',1:n.sce), labels = paste('truth',1:n.sce))
  truth2 = melt(da[,c('date', paste0('wt.Rtx.sce', 1:n.sce)),with=F], id.vars = 'date') %>% setnames('variable','loc')
  truth2$state = 'Rtx'; 
  truth2$loc = factor(truth2$loc, levels = paste0('wt.Rtx.sce',1:n.sce), labels = paste('truth',1:n.sce))
  truth = rbind(truth1, truth2) %>% setnames(c('date','value'),c('Week.start', 'obs'))
  truth$Week.start = truth$Week.start %>% as.Date
  # res.train$Week.start = res.train$Week.start %>% as.Date
  res.train$loc = factor(res.train$loc, levels = paste0('sce',1:n.sce), labels = paste('truth',1:n.sce))
  res.train = merge(res.train, truth, all = T, by = c('loc','state','Week.start'))
  
  
  pdf(paste0(dir_res, 'a',itest, 'Fig_syn.test_parms_alpha',alpha.t,'_',tag.eval,'.pdf'), width = 12, height = n.sce * 3)
  p = getPlotStates(res.train, var = c('Rt', 'Rtx','infection detection rate','IFR','Susceptibility'))
  print(p)
  dev.off()
  
  pdf(paste0(dir_res, 'a',itest, 'Fig_syn.test_tx.suscept_alpha',alpha.t,'_',tag.eval,'.pdf'), width = 12, height = n.sce * 3)
  p = getPlot(res.train[state %in% c('Susceptibility','Rtx')])
  print(p)
  dev.off()
  
  
  # compute dRtx and dImm for individual runs and then combine

  res = newVstat[eval == tag.eval] %>% data.table()
  res = melt(res, id.vars = c('loc','run','hyp.test','eval')) 
  res = res[variable %in% c('perc.dImm.mn','perc.dRtx.mean')]
  res$state = factor(res$variable, levels = c('perc.dRtx.mean','perc.dImm.mn'), labels =  c('dRtx', 'dImm'))
  res = res[order(loc, run, state)]
  truths = c(0, 80, 50, 0, 50, 80, 25, 40, 50, 0)
  
  n <- prod(dim(with(res, table(state, loc))))
  VEC <- seq(1, n/2, length.out=n)*2 - c(0, .2)
  pdf(paste0(dir_res, 'a',itest, 'Fig_syn.test_est_tx_imm_change_alpha',alpha.t,tag.eval,'.pdf'),width = 6, height = 3)
  par(mfrow=c(1,1), mar=c(2, 2.5, 1.5, .5), mgp = c(1.3, .4, 0), tck = -.02)
  boxplot(value ~ loc:state, data =res, ylab = 'Relative change (%)',
          boxwex=0.5, col=alpha(c("orange", "yellow"),0.2),
          xlab="", xaxt = 'n', outcex = 0.2,
          sep=":", lex.order=TRUE, ylim=c(-10, 110), yaxs="i",
          cex.axis=.8, at=VEC)
  points(VEC, truths, pch = '*', col='red', cex=2)
  axis(1, at = VEC, labels = outer(c('tx','S'), 1:n.sce, FUN = paste0) %>% c)
  dev.off()
} # tag.eval

# put all together
pdf(paste0(dir_res, 'a',itest, 'Fig_syn.test_est_tx_imm_change_alpha',alpha.t,'.pdf'),width = 6, height = length(tag.evals) * 2)
par(mfrow=c(length(tag.evals),1), mar=c(2, 2.5, 1.5, .5), mgp = c(1.2, .3, 0), cex=1, tck = -.03)
for(tag.eval in tag.evals){
  res = newVstat[eval == tag.eval] %>% data.table()
  res = melt(res, id.vars = c('loc','run','hyp.test','eval')) 
  res = res[variable %in% c('perc.dImm.mn','perc.dRtx.mean')]
  res$state = factor(res$variable, levels = c('perc.dRtx.mean','perc.dImm.mn'), labels =  c('dRtx', 'dImm'))
  res = res[order(loc, run, state)]
  truths = c(0, 80, 50, 0, 50, 80, 25, 40, 50, 0)
  
  n <- prod(dim(with(res, table(state, loc))))
  VEC <- seq(1, n/2, length.out=n)*2 - c(0, .2)
  
  boxplot(value ~ loc:state, data =res, ylab = 'Relative change (%)',
          boxwex=0.5, col=alpha(c("orange", "yellow"),0.2),
          xlab="", xaxt = 'n', outcex = 0.2,
          sep=":", lex.order=TRUE, ylim=c(-10, 110), yaxs="i",
          cex.axis=.95, at=VEC)
  points(VEC, truths, pch = '*', col='red', cex=2)
  axis(1, at = VEC, cex=.95, labels = outer(c('tx','S'), 1:n.sce, FUN = paste0) %>% c)
  mtext(tag.eval,cex=.95, line = .1, adj = .05)
}
dev.off()


loc.t = 'sce5'
tmp = newVstat[eval == 'obs.most' & loc == loc.t]
table(tmp$hyp.test)
tmp2 = newVstat[eval == 'obs.more' & loc == loc.t]
table(tmp2$hyp.test)
tmp2 = newVstat[eval == 'obs.comb' & loc == loc.t]
table(tmp2$hyp.test)
