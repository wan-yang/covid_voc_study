getPlot = function(tda){
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    geom_point(mapping = aes(x = Week.start, y=obs)) + 
    # geom_line(mapping = aes(x = Week.start, y=threshold), color = 'red') + 
    facet_rep_wrap(~loc + state, scales = 'free_y', repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),
                 labels = format(seq(min(tda$Week.start), max(tda$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlotStates = function(res, var = parm.names){
  p = ggplot(res[state %in% var]) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    facet_rep_wrap(~ loc + state, scales = 'free_y', repeat.tick.labels = T, ncol=length(var)) + 
    labs(x = 'Week Start', y = 'Estimate (median, IQR, 95% CI)') +
    scale_x_date(breaks = seq(min(res$Week.start), max(res$Week.start), by = 'month'),
                 labels = format(seq(min(res$Week.start), max(res$Week.start), by = 'month'),'%Y %b')) +
    theme_minimal() + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 30))
  
  p
}


getPlotProj = function(train.t, proj.t){

  p = ggplot(train.t) +
    geom_line(aes(x = Week.start, y = median), color = 'blue') +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .2) +
    geom_line(data=proj.t, aes(x = Week.start, y = median), color = 'red') +  # no ctrl
    # geom_ribbon(data=proj.t, aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'red', alpha = .1) +
    geom_ribbon(data=proj.t, aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'red', alpha = .2) +
    geom_vline(data = train.t, aes(xintercept = max(train.t$Week.start)), linetype = 'dashed')+
    geom_point(data = train.t, mapping = aes(x = Week.start, y=obs)) + 
    facet_rep_wrap(~ loc, scales = 'free_y', repeat.tick.labels = T, ncol = 1) + 
    labs(x = 'Week Start', y = 'Estimate/Projection (median, IQR)') +
    scale_x_date(breaks = seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),
                 labels = format(seq(min(train.t$Week.start)-7, max(proj.t$Week.start)+7, by = 'week'),'%m/%d')) +
    theme_minimal() + theme(strip.text = element_text(size = 10), axis.title = element_text(size =10), axis.text.y = element_text(size=10), axis.text.x = element_text(size=10,angle = 45))
  
  p
}


getPlotMultiV = function(tda, y.lab = 'Estimate'){
  
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) +
    geom_line(aes(x = Week.start, y = v.median), color = 'blue', size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr), fill = 'blue', alpha = .3) +
    facet_rep_wrap(~ variant, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median, IQR)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 30))
  
  p
}

getPlotMultiVoverlay = function(tda, title.t, y.lab = 'Estimate',withObs = F, obs.t = NULL){
 
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) + ggtitle(title.t) +
    geom_line(aes(x = Week.start, y = v.median, color = variant), size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = variant), alpha = .3) +
    facet_rep_wrap(~ measure, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median, IQR)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 20))
  
  if(withObs & is.null(obs.t$variant)){
    p = p + geom_point(data = obs.t[date %in% dates.t], aes(x = date, y = value))
  } else if (withObs & !is.null(obs.t$variant)){
    p = p + geom_point(data = obs.t[date %in% dates.t], aes(x = date, y = value, color = variant, shape = variant))
  } 
  p
}

getPlotMultiVoverlayMedian = function(tda, title.t, y.lab = 'Estimate'){
  
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) + ggtitle(title.t) +
    geom_line(aes(x = Week.start, y = v.median, color = variant), size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    # geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = variant), alpha = .3) +
    facet_rep_wrap(~ measure, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 20))
  
  p
}


getPlotMultiVage.overlay = function(tda, title.t, y.lab = 'Estimate'){
  
  dates.t = unique(tda$Week.start) %>% as.Date
  p = ggplot(tda) + ggtitle(title.t) +
    geom_line(aes(x = Week.start, y = v.median, color = age.grp), size = 1) +  # no ctrl
    # geom_ribbon(aes(x = Week.start, ymin = ci95.lwr, ymax = ci95.upr), fill = 'blue', alpha = .1) +
    geom_ribbon(aes(x = Week.start, ymin = iqr.lwr, ymax = iqr.upr, fill = age.grp), alpha = .15) +
    facet_rep_wrap(~ variant, scales = 'free_y', 
                   repeat.tick.labels = T, ncol = 2) + 
    labs(x = 'Week Start', y = paste(y.lab, '(median, IQR)')) +
    scale_x_date(breaks = dates.t[seq(1, length(dates.t), by = 4)],
                 labels = format(dates.t[seq(1, length(dates.t), by = 4)],'%m/%d/%y')) +
    theme_minimal() + theme(strip.text = element_text(size = 12), axis.title = element_text(size =12), axis.text.y = element_text(size=12), axis.text.x = element_text(size=10,angle = 20))
  
  p
}

