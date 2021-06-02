# plot results from the synthetic test
# 5/4/21


alpha.t = .2
eval.t = 'obs.more'
dir_data = '../../data/'
dir_code = '../model_inference/'
dir_code2 =  './' 

dir_res.up = paste0('../../results/syn.test/')
dir_plot = paste0('../../results/plots/')

dir_res.sub = paste0(dir_res.up, 'test_seed1st2seed2nd50alpha', alpha.t,'/')


if(!file.exists(dir_plot)) dir.create(dir_plot)

library(data.table); library(magrittr); library(ggplot2); library(MMWRweek); library(lemon)
library(gridExtra)
# read in data and model estimates
seed1stWave = 2
seed2ndWave = 50
n.sce = 5
da = read.csv(paste0(dir_data,'da_syn.truth_seed1st',seed1stWave,'_seed2nd', seed2ndWave,'_alpha',alpha.t,'.csv'), stringsAsFactors = F)  %>% data.table()
da[date < as.Date('2020-02-09') & is.na(da)] = 0
tmp = MMWRweek(da$date)
da$year = tmp[,1]
da$week = tmp[,2]
# set population size to 20k? make it smaller for stochasticity
N = 1e6; # per 1 M

# load results
load(paste0(dir_res.sub, "res.summary.RData"))
res.train = res.train[eval == eval.t]

theme.t = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,2,0)), 
                strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1,0)),
                axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=8,angle = 30),
                plot.margin=unit(c(c(.3, 1, -1, .5)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=8), legend.text=element_text(size=8),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                legend.key.size = unit(.2, 'cm'), #change legend key size
                legend.key.height = unit(.2, 'cm'), #change legend key height
                legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t2 = theme(plot.title = element_text(v=0, size = 10, margin=margin(0,0,2,0)), 
                strip.placement = "outside", strip.text = element_text(size = 9, margin=margin(1.5,0,1,0)),
                axis.title = element_text(size =9, margin=margin(0,0.2,0,0)), 
                axis.text.y = element_text(size=8, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=8,angle = 30),
                plot.margin=unit(c(c(.3, 1, -1, -.25)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=8), legend.text=element_text(size=8),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,-10,-10,-10),
                legend.key.size = unit(.2, 'cm'), #change legend key size
                legend.key.height = unit(.2, 'cm'), #change legend key height
                legend.key.width = unit(.2, 'cm')) #change legend key width)

getPlotMultiLoc = function(tda, ytitle = 'Estimate (median, IQR, 95% CI)', ptitle = '', ncol.t = 2){
  p = ggplot(tda) +
    geom_line(aes(x = date, y = median,color = loc)) +  # no ctrl
    geom_ribbon(aes(x = date, ymin = ci95.lwr, ymax = ci95.upr, fill = loc),  alpha = .15) +
    geom_ribbon(aes(x = date, ymin = iqr.lwr, ymax = iqr.upr, fill = loc), alpha = .35) +
    geom_point(mapping = aes(x = date, y=obs, color = loc), size = .6) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_rep_wrap(~state, scales = 'free_y', repeat.tick.labels = T, ncol = ncol.t) + 
    labs(x = '', y = ytitle) + ggtitle(ptitle) +
    guides(color = guide_legend(title = '', override.aes = list(size = 0.2)), fill = F) + 
    scale_x_date(breaks = seq(min(tda$date), max(tda$date), by = '2 months'),
                 labels = format(seq(min(tda$date), max(tda$date), by = '2 months'),'%Y%b')) +
    theme_minimal() + theme.t 
  # truths.t = tda
  p = ggplot(tda) +
    geom_line(aes(x = date, y = median,color = loc)) +  # no ctrl
    geom_ribbon(aes(x = date, ymin = ci95.lwr, ymax = ci95.upr, fill = loc),  alpha = .15) +
    geom_ribbon(aes(x = date, ymin = iqr.lwr, ymax = iqr.upr, fill = loc), alpha = .35) +
    geom_point(mapping = aes(x = date, y=obs, color = loc, shape = loc),  size = .6) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    # scale_color_brewer(palette="Paired") + 
    facet_rep_wrap(~state, scales = 'free_y', repeat.tick.labels = T, ncol = ncol.t) + 
    labs(x = '', y = ytitle) + ggtitle(ptitle) +
    guides(color = guide_legend(title = 'Fitted', override.aes = list(size = 1, shape = NA)), fill = F,
           shape = guide_legend(title = 'Truth', override.aes = list(size = 0.5))) + 
    scale_x_date(breaks = seq(min(tda$date), max(tda$date), by = '2 months'),
                 labels = format(seq(min(tda$date), max(tda$date), by = '2 months'),'%Y%b')) +
    theme_minimal() + theme.t 
  p
}

getPlotIndivLoc = function(tda, ytitle = 'Estimate (median, IQR, 95% CI)', ptitle = '',ncol.t = 5){
  # p = ggplot(tda) +
  #   geom_line(aes(x = date, y = median), color = 'blue') +  # no ctrl
  #   geom_ribbon(aes(x = date, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue',  alpha = .15) +
  #   geom_ribbon(aes(x = date, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue',  alpha = .35) +
  #   geom_point(mapping = aes(x = date, y=obs),  size = .6) + 
  #   # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
  #   facet_grid(state ~ loc, scales = 'free_y', switch = 'y') + 
  #   labs(x= '', y = '') + ggtitle(ptitle) + 
  #   guides(guide_legend(override.aes = list(size = 0.2))) + 
  #   scale_x_date(breaks = seq(min(tda$date), max(tda$date), by = '3 months'),
  #                labels = format(seq(min(tda$date), max(tda$date), by = '3 months'),'%Y %b')) +
  #   theme_minimal() + theme(plot.title = element_text(v=0, size = 10), strip.placement = "outside", strip.text = element_text(size = 9), axis.title = element_text(size = 9), axis.text.y = element_text(size=8), axis.text.x = element_text(size=8,angle = 20))
  # 
  p = ggplot(tda) +
    geom_line(aes(x = date, y = median, color = 'blue')) +  # no ctrl
    geom_ribbon(aes(x = date, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue',  alpha = .15) +
    geom_ribbon(aes(x = date, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue',  alpha = .35) +
    geom_point(mapping = aes(x = date, y=obs, color = 'black'),  size = .6) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_grid(state ~ loc, scales = 'free_y', switch = 'y') + 
    labs(x= '', y = '') + ggtitle(ptitle) + 
    scale_color_identity(name = "",
                         breaks = c( "black","blue"),
                         labels = c("Truth","Fitted"),
                         guide = guide_legend(override.aes = list(linetype = c(NA,1), 
                                                                  shape = c(20,NA))))+
    guides(guide_legend(override.aes = list(size = 0.2))) + 
    scale_x_date(breaks = seq(min(tda$date), max(tda$date), by = '2 months'),
                 labels = format(seq(min(tda$date), max(tda$date), by = '2 months'),'%Y%b')) +
    theme_minimal() + theme.t2
  
  p
}




getPlotOverall = function(newVstat, eval.t, ptitle, truths = truths.t){
  res = newVstat[eval == eval.t] %>% data.table()
  res = melt(res, id.vars = c('loc','run','hyp.test','eval')) 
  res = res[variable %in% c('perc.dImm.mn','perc.dRtx.mean')]
  res$state = factor(res$variable, levels = c('perc.dRtx.mean','perc.dImm.mn'), labels =  c('Tx', 'S'))
  res$loc = factor(res$loc, levels = paste0('sce', 1:n.sce), labels = paste('Truth', 1:n.sce))
  res = res[order(loc, run, state)]
  res = merge(res, truths, by = c('loc', 'state'))
  res$state = factor(res$state, levels = c('Tx','S'))
  truths$state = factor(truths$state, levels = c('Tx','S'))
  
  # p<-ggplot(res, aes(x=loc, y=value, fill = as.factor(state))) +
  #   geom_boxplot(outlier.size= .5, alpha = .2, lwd = .3, color = 'grey30') + 
  #   labs(fill = '', x='', y='Relative change (%)') + ggtitle(ptitle) + 
  #   geom_point(aes(x = loc, y = obs, color = as.factor(state)), size = .3, alpha =1, position = position_jitterdodge()) +
  #   guides(color = FALSE, fill = guide_legend(override.aes = list(size = 0.2, alpha =1))) + 
  #   theme_minimal() + theme.t
  
  p<-ggplot(res, aes(x=loc, y=value, fill = as.factor(state))) +
    geom_boxplot(outlier.size= .5, alpha = .3, lwd = .3, color = 'grey30') + 
    labs(fill = '', x='', y='Relative change (%)') + ggtitle(ptitle) + 
    geom_point(data=truths, aes(x = loc, y = obs, shape = as.factor(state)), color = 'brown', size = 1.2, alpha =1, position = position_dodge(width=.75)) +
    scale_shape_manual(values=c(4,8))+
    scale_fill_brewer(palette="Paired") + 
    # scale_fill_manual(values=c('navy','red'))+
    guides(shape = guide_legend(title = 'Truth', override.aes = list(size = 0.8, alpha =1)), fill = guide_legend(title = 'Fitted',override.aes = list(size = 0.2, color = 'black', alpha =.5))) + 
    theme_minimal() + theme.t
  p
}
# states.t = c('case','death','wt.suscept','wt.Rtx')
# states_names.t = c('case','death','Susceptibility','Rtx')

states.t = c('case','death')
states_names.t = c('case','death')
state_labs = c('Incidence','Mortality')
loc_labs = paste('Truth', 1:n.sce)
states_all = outer(states.t, paste0('.sce',1:n.sce), FUN = paste0) %>% c
da1 = melt(da[,c('date', states_all), with=F], id.vars = 'date')
da1$loc = factor(da1$variable, levels = states_all, labels = rep(paste0('sce',1:n.sce),e=length(states.t)))
da1$state = factor(da1$variable, levels = states_all, labels = rep(states_names.t, n.sce))
# da1[state == 'Susceptibility']$value = da1[state == 'Susceptibility']$value / N * 100
da1$date = da1$date %>% as.Date
setnames(da1,'value','obs')
est1 = dcast(res.train, loc + state + Week.start ~ variable, value.var = 'value')
est1$date = est1$Week.start %>% as.Date
tda1 = merge(da1, est1[state %in% states_names.t], by = c('loc','state','date'))
tda1$loc = factor(tda1$loc, levels = paste0('sce',1:n.sce), labels = loc_labs)
tda1$state  = factor(tda1$state, levels = states_names.t, labels = state_labs)

# for the susceptiblity and Rtx
states.t = c('wt.suscept','wt.Rtx')
states_names.t = c('Susceptibility','Rtx')
state_labs = c('Susceptibility (%)','Transmissibility')
loc_labs = paste('Truth', 1:n.sce)
states_all = outer(states.t, paste0('.sce',1:n.sce), FUN = paste0) %>% c
da2 = melt(da[,c('date', states_all), with=F], id.vars = 'date')
da2$loc = factor(da2$variable, levels = states_all, labels = rep(paste0('sce',1:n.sce),e=length(states.t)))
da2$state = factor(da2$variable, levels = states_all, labels = rep(states_names.t, n.sce))
da2[state == 'Susceptibility']$value = da2[state == 'Susceptibility']$value / N * 100
da2$date = da2$date %>% as.Date
setnames(da2,'value','obs')
est2 = dcast(res.train, loc + state + Week.start ~ variable, value.var = 'value')
est2$date = est2$Week.start %>% as.Date
tda2 = merge(da2, est2[state %in% states_names.t], by = c('loc','state','date'))
tda2$loc = factor(tda2$loc, levels = paste0('sce',1:n.sce), labels = loc_labs)
tda2$state  = factor(tda2$state, levels = rev(states_names.t), labels = rev(state_labs))

# overall estimates
truths.t = data.table(loc = rep(paste('Truth', 1:n.sce),e=2), 
                      state = rep(c('Tx', 'S'), n.sce) %>% as.factor(), 
                      obs = c(0, 80, 50, 0, 50, 80, 25, 40, 50, 0))
truths = truths.t 

p1 = getPlotMultiLoc(tda1, ytitle = 'Number per million population', ptitle = '(A) Mock data and model fit')
p2 = getPlotIndivLoc(tda2, ptitle = '(B) Weekly estimates', ncol.t=5)
p3 = getPlotOverall(newVstat = newVstat, eval.t=eval.t, ptitle = '(C) Overall estimates', truths = truths.t)
p1

pdf(paste0(dir_plot,'Fig1_syn.test_alpha',alpha.t,'.pdf'),width = 7, height = 5)
grid.arrange(
  grobs = list(p1, p3, p2),
  widths = c(2, 1),
  layout_matrix = rbind(c(1, 2),
                        c(1, 2),
                        c(3, 3),
                        c(3, 3),
                        c(3, 3))
)
dev.off()



