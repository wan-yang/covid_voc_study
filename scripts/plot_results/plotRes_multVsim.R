# plot results from the multi-variant simulation runs 

dir_data = '../../data/'
dir_code = '../model_inference/'
dir_code2 =  './' 

dir_res.up = paste0('../../results/sim_multiV/')
dir_res1 = paste0(dir_res.up,'eqB1351P1/')
dir_res2 = paste0(dir_res.up,'moreP1/')
dir_res3 = paste0(dir_res.up,'moreB1351/')
dir_plot = paste0('../../results/plots/')

seed.tags = c('eqB1351P1', 'moreP1', 'moreB1351')
dir.tags = c('Equal seeding','More P.1','More B.1.351')
  
if(!file.exists(dir_plot)) dir.create(dir_plot)

library(data.table); library(magrittr); library(ggplot2); library(MMWRweek); library(lemon)
library(gridExtra)

# functions
source(paste0(dir_code2,'Fn_util.R'))
{
  theme.t = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,3,0)), 
                  strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                  axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                  axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=8,angle = 0),
                  plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                  legend.title = element_text(size=8), legend.text=element_text(size=8),
                  legend.margin=margin(0,0,0,0),
                  legend.box.margin=margin(-10,-10,-10,-10),
                  legend.key.size = unit(.2, 'cm'), #change legend key size
                  legend.key.height = unit(.2, 'cm'), #change legend key height
                  legend.key.width = unit(.2, 'cm')) #change legend key width)
  
  theme.t2 = theme(plot.title = element_text(v=0, size = 10, margin=margin(2,0,3,0)), 
                   strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                   axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                   axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                   axis.text.x = element_text(size=8,angle = 45, hjust = .9),
                   plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                   legend.title = element_text(size=8), legend.text=element_text(size=8),
                   legend.margin=margin(0,0,0,0),
                   legend.box.margin=margin(-10,-10,-10,-10),
                   legend.key.size = unit(.2, 'cm'), #change legend key size
                   legend.key.height = unit(.2, 'cm'), #change legend key height
                   legend.key.width = unit(.2, 'cm')) #change legend key width)
  theme.t3 = theme(plot.title = element_text(v=0, size = 10, margin=margin(2,0,3,0)), 
                   strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1.5,0)),
                   axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                   axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                   axis.text.x = element_text(size=8,angle = 20, hjust = .9),
                   plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                   legend.title = element_text(size=8), legend.text=element_text(size=8),
                   legend.margin=margin(0,0,0,0),
                   legend.box.margin=margin(-10,-10,-10,-10),
                   legend.key.size = unit(.2, 'cm'), #change legend key size
                   legend.key.height = unit(.2, 'cm'), #change legend key height
                   legend.key.width = unit(.2, 'cm')) #change legend key width)
  getPlotMultiVoverlay = function(tda, title.t, y.lab = 'Estimate', ncol.t = 3, withObs = F, col.set){
    
    dates.t = unique(tda$Week.start) %>% as.Date
    p = ggplot(tda) + ggtitle(title.t) +
      geom_line(aes(x = Week.start, y = v.median, color = variant), size = .7) +  # no ctrl
      # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
      geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = variant), alpha = .3) +
      facet_rep_wrap(~ seed.ve, scales = 'free_y',  # sce.cb + seeding
                     repeat.tick.labels = T, ncol = ncol.t) +  
      # scale_fill_brewer(palette = col.set) +
      # scale_color_brewer(palette = col.set) +
      labs(x = '', y = paste(y.lab, '(median, IQR)')) +
      scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                   labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
      theme_minimal() + theme.t
    
    if(withObs){
      p = p + geom_point(data = obs[date %in% dates.t], aes(x = date, y = value))
    }
    p
  }
  
  getPlotTotals = function(tda, title.t, y.lab = 'Projections, per 1 M popuplation (median)', ncol.t = 9, col.set, theme.tt = theme.t2){
    p = ggplot(tda, aes(fill=variant, y=v.median, x=sce.npi)) + 
      geom_bar(position="stack", stat="identity", alpha = .6) +
      ggtitle(title.t) + labs(x = '', y = y.lab) +
      facet_wrap(~seeding + sce.ve, ncol = ncol.t) + 
      # facet_wrap(~seed.ve, ncol = ncol.t) + 
      # scale_fill_brewer(palette = col.set) +
      theme_minimal() + theme.tt
    
    p
  }
  
  
  
}

Nnyc.tot = 8.4e6
N = 1e6
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

variant.combs = list(c('wt','b117','b1351','p1'))

##########################################################################
# NOTE: DUE TO THE LARGE FILE SIZE, WE DO NOT INCLUDE THE SIMULATION RESULTS
# THUS, THIS SECTION OF CODE WOULD NOT RUN 
# PLEASE RUN THE SIMULATIONS ON YOU LOCAL MACHINE AND THEN COMPILE
PERC = STATS = SUMS.PERC = SUMS.STATS = TOTALS = NULL
for(idir in 1:3){
  dir.t = get(paste0('dir_res',idir))
  dir.tag.t = dir.tags[idir]
  seed.tag = seed.tags[idir]
  
  for(ivcomb in 1: length(variant.combs)){
    
    variants.t = variant.combs[[ivcomb]]
    Ns= length(variants.t) # ncol(obs_i);  # number of strains included
    v.tag = paste(variants.t, collapse = '.')
    
    
    for(isce in 1:length(senarios)){ # 1: length(senarios)
      
      sce.tag = sce.tags[isce]
      sce.t = senarios[isce]
      
      load(paste0(dir.t, 'out_', seed.tag, '_', v.tag,'_',sce.tag, '.RData'))
      perc$seeding = dir.tag.t
      stats$seeding = dir.tag.t
      sums.perc$seeding = dir.tag.t
      sums.stats$seeding = dir.tag.t
      totals$seeding = dir.tag.t
      
      perc$sce.cb = sce.tag
      stats$sce.cb = sce.tag
      sums.perc$sce.cb = sce.tag
      sums.stats$sce.cb = sce.tag
      totals$sce.cb = sce.tag
      
      perc$vcomb = v.tag
      stats$vcomb = v.tag
      sums.perc$vcomb = v.tag
      sums.stats$vcomb = v.tag
      totals$vcomb = v.tag

      
      PERC = rbind(PERC, perc)
      STATS = rbind(STATS, stats)
      SUMS.PERC = rbind(SUMS.PERC, sums.perc)
      SUMS.STATS = rbind(SUMS.STATS, sums.stats)
      TOTALS = rbind(TOTALS, totals)
      
      rm(x)
    }
  }
} # diff seeding
##########################################################################

# PLOT SIMULATED WEEKLY INFECTION RATE
mea.t = 'New infections'
ex.sce = c('noNPIslow.sameVE', 'noNPIslow.highVE')
ex.lab = c('Same VE','High VE')
tda = STATS[sce.cb %in% ex.sce & measure == mea.t & age.grp =='All' & variant != 'All']
tda$sce.cb = factor(tda$sce.cb, levels = ex.sce, labels = ex.lab)  
tda$seeding = factor(tda$seeding, levels = dir.tags)  
tda$variant = factor(tda$variant, levels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'), labels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'))
tda$seed.ve = apply(tda[,c('sce.cb','seeding'),with =F], 1, paste, collapse = ', ')
tda$seed.ve = factor(tda$seed.ve, levels = c("Same VE, Equal seeding", "Same VE, More P.1", "Same VE, More B.1.351",
                                             "High VE, Equal seeding","High VE, More P.1", "High VE, More B.1.351"), 
                      labels = c("Same VE, equal seeding", "Same VE, more P.1", "Same VE, more B.1.351",
                                 "High VE, equal seeding","High VE, more P.1", "High VE, more B.1.351"))

# scale to 1 M pop
tda[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')] = 
  tda[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr'), with =F] / Nnyc.tot * N


# compare fast v slow reopening
ex.sce = c('noNPIslow.sameVE', 'noNPIfast.sameVE')
ex.lab = c('Slow full reopen','Fast full reopen')
tda1 = STATS[sce.cb %in% ex.sce & measure == mea.t & age.grp =='All' & variant != 'All']
tda1$sce.cb = factor(tda1$sce.cb, levels = ex.sce, labels = ex.lab)  
tda1$seeding = factor(tda1$seeding, levels = dir.tags)  
tda1$variant = factor(tda1$variant, levels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'), labels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'))
tda1$seed.ve = apply(tda1[,c('sce.cb','seeding'),with =F], 1, paste, collapse = ', ')
tda1$seed.ve = factor(tda1$seed.ve, levels = c("Slow full reopen, Equal seeding", "Slow full reopen, More P.1", "Slow full reopen, More B.1.351",
                                             "Fast full reopen, Equal seeding","Fast full reopen, More P.1", "Fast full reopen, More B.1.351"), 
                     labels = c("Slow full reopen, equal seeding", "Slow full reopen, more P.1", "Slow full reopen, more B.1.351",
                                "Fast full reopen, equal seeding","Fast full reopen, more P.1", "Fast full reopen, more B.1.351"))

# scale to 1 M pop
tda1[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')] = 
  tda1[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr'), with =F] / Nnyc.tot * N

  
tda2 = SUMS.STATS[measure == mea.t]  #  & seeding == 'Equal seeding'
tda2$sce.npi = factor(tda2$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                              '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                              '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                              'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                              'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                      labels = rep(c('Current NPI','25% less','50% less','Fully open, slow','Fully open, fast'), e = 3))
tda2$sce.ve = factor(tda2$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                              '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                              '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                              'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                              'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                      labels = rep(c('Same VE','High VE','Median VE'), 5))
tda2$seeding = factor(tda2$seeding, levels = dir.tags)  
tda2$seed.ve = apply(tda2[,c('seeding','sce.ve'),with =F], 1, paste, collapse = ', ')
tda2$seed.ve = factor(tda2$seed.ve, levels = c("Equal seeding, Same VE", "Equal seeding, High VE","Equal seeding, Median VE", 
                                               "More P.1, Same VE","More P.1, High VE","More P.1, Median VE",
                                               "More B.1.351, Same VE", "More B.1.351, High VE", "More B.1.351, Median VE"), 
                      labels = c("Equal seeding, same VE", "Equal seeding, high VE","Equal seeding, median VE", 
                                 "More P.1, same VE","More P.1, high VE","More P.1, median VE",
                                 "More B.1.351, same VE", "More B.1.351, high VE", "More B.1.351, median VE"))
# # scale to 1 M pop
tda2[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')] = 
  tda2[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr'), with =F] / Nnyc.tot * N

col.set.t = 'Set3'
p1 = getPlotMultiVoverlay(tda, title.t = '(A) Example projections: Slow full reopen', 
                          y.lab = paste0('Projected weekly number of ',gsub('New ','', mea.t), ' per million'), ncol.t = 3, withObs = F, col.set = col.set.t)
# print(p1)
# plot summary with multiple comparison
p2 = getPlotTotals(tda2, title.t = paste0('(B) Projected cumulative ',gsub('New ','', mea.t,': Assuming mRNA vaccines used')), 
                   y.lab = 'Projected totals per million (median)', ncol.t = 9, col.set = col.set.t)
# print(p2)
p3 = getPlotMultiVoverlay(tda1, title.t = '(A) Example projections: Assuming mRNA vaccines used and same VE as for the wildtype virus', 
                          y.lab = paste0('Projected weekly number of ',gsub('New ','', mea.t), ' per million'), ncol.t = 3, withObs = F, col.set = col.set.t)

pdf(paste0(dir_plot,'Fig4v0_proj_', mea.t, date.tag,'.pdf'),width = 9, height = 7)
grid.arrange(
  grobs = list(p1, p2),
  layout_matrix = rbind(c(1, 1),
                        c(1, 1),
                        c(1, 1),
                        c(2, 2),
                        c(2, 2))
)
dev.off()

pdf(paste0(dir_plot,'Fig4_proj_', mea.t, date.tag,'.pdf'),width = 9, height = 7)
grid.arrange(
  grobs = list(p3, p2),
  layout_matrix = rbind(c(1, 1),
                        c(1, 1),
                        c(1, 1),
                        c(2, 2),
                        c(2, 2))
)
dev.off()

tda3 = SUMS.PERC[measure == mea.t]  #  & seeding == 'Equal seeding'
tda3$sce.npi = factor(tda3$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                              '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                              '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                              'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                              'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                      labels = rep(c('Current NPI','25% less','50% less','Fully open, slow','Fully open, fast'), e = 3))
tda3$sce.ve = factor(tda3$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                             '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                             '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                             'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                             'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                     labels = rep(c('Same VE','High VE','Median VE'), 5))

tda3$cb = apply(tda3[,c('v.median','iqr.lwr','iqr.upr'),with=F], 1, FUN = fn_format, roundigt = 1)
tda4 = dcast(tda3, seeding + sce.npi + sce.ve ~ variant, value.var = 'cb')
tda4$seeding = factor(tda4$seeding, levels = dir.tags)
tda4 = tda4[order(seeding, sce.npi, sce.ve)]
write.csv(tda4, paste0(dir_plot, 'TableS_sim.totals_', mea.t, '_perc.prevalence.csv'), row.names = F)



# for mortality, b/c uncertainty in IFR among breakthrough infections of vaccinated, only show the 'same VE' scenario
mea.t = 'New deaths'
ex.sce = c('noNPIslow.sameVE', 'noNPIslow.highVE')
ex.lab = c('Same VE','High VE')
tda = STATS[sce.cb %in% ex.sce & measure == mea.t & age.grp =='All' & variant != 'All']
tda$sce.cb = factor(tda$sce.cb, levels = ex.sce, labels = ex.lab)  
tda$seeding = factor(tda$seeding, levels = dir.tags)  
tda$variant = factor(tda$variant, levels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'), labels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'))
tda$seed.ve = apply(tda[,c('sce.cb','seeding'),with =F], 1, paste, collapse = ', ')
tda$seed.ve = factor(tda$seed.ve, levels = c("Same VE, Equal seeding", "Same VE, More P.1", "Same VE, More B.1.351",
                                             "High VE, Equal seeding","High VE, More P.1", "High VE, More B.1.351"), 
                     labels = c("Same VE, equal seeding", "Same VE, more P.1", "Same VE, more B.1.351",
                                "High VE, equal seeding","High VE, more P.1", "High VE, more B.1.351"))

# scale to 1 M pop
tda[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')] = 
  tda[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr'), with =F] / Nnyc.tot * N


# compare fast v slow reopening
ex.sce = c('noNPIslow.sameVE', 'noNPIfast.sameVE')
ex.lab = c('Slow full reopen','Fast full reopen')
tda1 = STATS[sce.cb %in% ex.sce & measure == mea.t & age.grp =='All' & variant != 'All']
tda1$sce.cb = factor(tda1$sce.cb, levels = ex.sce, labels = ex.lab)  
tda1$seeding = factor(tda1$seeding, levels = dir.tags)  
tda1$variant = factor(tda1$variant, levels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'), labels = c('Wildtype','B.1.1.7', 'B.1.351', 'P.1','All'))
tda1$seed.ve = apply(tda1[,c('sce.cb','seeding'),with =F], 1, paste, collapse = ', ')
tda1$seed.ve = factor(tda1$seed.ve, levels = c("Slow full reopen, Equal seeding", "Slow full reopen, More P.1", "Slow full reopen, More B.1.351",
                                               "Fast full reopen, Equal seeding","Fast full reopen, More P.1", "Fast full reopen, More B.1.351"), 
                      labels = c("Slow full reopen, equal seeding", "Slow full reopen, more P.1", "Slow full reopen, more B.1.351",
                                 "Fast full reopen, equal seeding","Fast full reopen, more P.1", "Fast full reopen, more B.1.351"))

# scale to 1 M pop
tda1[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')] = 
  tda1[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr'), with =F] / Nnyc.tot * N


tda2 = SUMS.STATS[measure == mea.t & sce.cb %in% c('curNPI.sameVE','25lNPI.sameVE','50lNPI.sameVE','noNPIslow.sameVE','noNPIfast.sameVE')]  #  & seeding == 'Equal seeding'
tda2$sce.npi = factor(tda2$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                              '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                              '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                              'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                              'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                      labels = rep(c('Current NPI','25% less','50% less','Fully open, slow','Fully open, fast'), e = 3))
tda2$sce.ve = factor(tda2$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                             '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                             '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                             'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                             'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                     labels = rep(c('Same VE','High VE','Median VE'), 5))
tda2$seeding = factor(tda2$seeding, levels = dir.tags)  
tda2$seed.ve = apply(tda2[,c('seeding','sce.ve'),with =F], 1, paste, collapse = ', ')
tda2$seed.ve = factor(tda2$seed.ve, levels = c("Equal seeding, Same VE", "Equal seeding, High VE","Equal seeding, Median VE", 
                                               "More P.1, Same VE","More P.1, High VE","More P.1, Median VE",
                                               "More B.1.351, Same VE", "More B.1.351, High VE", "More B.1.351, Median VE"), 
                      labels = c("Equal seeding, same VE", "Equal seeding, high VE","Equal seeding, median VE", 
                                 "More P.1, same VE","More P.1, high VE","More P.1, median VE",
                                 "More B.1.351, same VE", "More B.1.351, high VE", "More B.1.351, median VE"))
# # scale to 1 M pop
tda2[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr')] = 
  tda2[,c('v.mean','v.sd', 'v.median', 'iqr.lwr', 'iqr.upr', 'ci95.lwr','ci95.upr'), with =F] / Nnyc.tot * N

col.set.t = 'Set3'
p1 = getPlotMultiVoverlay(tda, title.t = '(A) Example projections: Slow full reopen', 
                          y.lab = paste0('Projected weekly number of ',gsub('New ','', mea.t), ' per million'), ncol.t = 3, withObs = F, col.set = col.set.t)
# print(p1)
# plot summary with multiple comparison
p2 = getPlotTotals(tda2, title.t = paste0('(B) Projected cumulative ',gsub('New ','', mea.t,': Assuming mRNA vaccines used')), 
                   y.lab = 'Projected totals per million (median)', ncol.t = 3, col.set = col.set.t, theme.tt = theme.t3)
# print(p2)
p3 = getPlotMultiVoverlay(tda1, title.t = '(A) Example projections: Assuming mRNA vaccines used and same VE as for the wildtype virus', 
                          y.lab = paste0('Projected weekly number of ',gsub('New ','', mea.t), ' per million'), ncol.t = 3, withObs = F, col.set = col.set.t)



pdf(paste0(dir_plot,'FigS4_proj_sameVE_', mea.t, date.tag,'.pdf'),width = 9, height = 7)
grid.arrange(
  grobs = list(p3, p2),
  layout_matrix = rbind(c(1, 1),
                        c(1, 1),
                        c(1, 1),
                        c(2, 2),
                        c(2, 2))
)
dev.off()


tda3 = SUMS.PERC[measure == mea.t]  #  & seeding == 'Equal seeding'
tda3$sce.npi = factor(tda3$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                              '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                              '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                              'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                              'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                      labels = rep(c('Current NPI','25% less','50% less','Fully open, slow','Fully open, fast'), e = 3))
tda3$sce.ve = factor(tda3$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                             '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                             '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                             'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                             'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                     labels = rep(c('Same VE','High VE','Median VE'), 5))

tda3$cb = apply(tda3[,c('v.median','iqr.lwr','iqr.upr'),with=F], 1, FUN = fn_format, roundigt = 1)
tda4 = dcast(tda3, seeding + sce.npi + sce.ve ~ variant, value.var = 'cb')
tda4$seeding = factor(tda4$seeding, levels = dir.tags)
tda4 = tda4[order(seeding, sce.npi, sce.ve)]
write.csv(tda4, paste0(dir_plot, 'TableS_sim.totals_', mea.t, '_perc.prevalence.csv'), row.names = F)

tda5 = TOTALS[measure == mea.t]
tda5$sce.npi = factor(tda5$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                              '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                              '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                              'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                              'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                      labels = rep(c('Current NPI','25% less','50% less','Fully open, slow','Fully open, fast'), e = 3))
tda5$sce.ve = factor(tda5$sce.cb, levels = c('curNPI.sameVE','curNPI.highVE','curNPI.midVE',
                                             '25lNPI.sameVE','25lNPI.highVE','25lNPI.midVE',
                                             '50lNPI.sameVE','50lNPI.highVE','50lNPI.midVE',
                                             'noNPIslow.sameVE','noNPIslow.highVE','noNPIslow.midVE',
                                             'noNPIfast.sameVE','noNPIfast.highVE','noNPIfast.midVE'),
                     labels = rep(c('Same VE','High VE','Median VE'), 5))

tda5$cb = apply(tda5[,c('v.median','iqr.lwr','iqr.upr'),with=F], 1, FUN = fn_format, roundigt = 1)
tda6 = dcast(tda5, seeding + sce.ve + sce.npi ~ variant, value.var = 'cb')
