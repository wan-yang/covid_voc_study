# combine results from multiple runs - diff eval method
state.names = c('S1', 'E1','I1',
                'death1', 'newIobs1','newItot1','beta',
                'Tei','Tir','Trs','Td.mean','Td.sd',
                'p.mob','alpha','ifr')
  
files = list.files(dir_res, full.names = T)
files = files[grepl('.RData',files)]
files = files[grepl('train', files)]

# compare performance across runs?


res.train <-  lapply(files, function(x) {
  
  try(load(x))
  
  id <- gsub(dir_res,'',x) %>% 
    strsplit('_') %>% unlist
  
  # get %S
  tmp = res.train$states_stats[state == 'S1', ]
  tmp[, c("mean","median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr")] = tmp[, c("mean","median","iqr.lwr","iqr.upr","ci95.lwr","ci95.upr")] / N * 100
  tmp$state = 'Susceptibility'
  
  # cumulative infection rate
  tmp2 = res.train$cumIperc_stats
  tmp2$state = 'cumItot'; tmp2$sd = NULL;
  
  d <- rbind(
    data.table(state = 'Rt',res.train$Rt_stats),
    data.table(state = 'Rtx',res.train$Rtx_stats),
    data.table(state = 'R0',res.train$R0_stats),
    data.table(tmp),
    data.table(tmp2),
    data.table(res.train$states_stats),
    fill = T
  ) %>%
    melt(id.vars = c('eval','state','Week.start'), variable.factor = F)
  
  
  d$loc = id[1] %>% strsplit('\\//') %>% unlist %>% tail(1)
  d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
  
  d
}) %>%
  rbindlist() %>%
  (function(d) d[, j = list(value = round(mean(value), 6)), by = list(eval,loc, state, variable, Week.start)]) # %>%

res.train$state = factor(res.train$state, levels=c('Rt','R0','Rtx', 'Susceptibility', 'cumItot', state.names), 
                         labels = c('Rt','R0','Rtx','Susceptibility', 'Cumulative infection rate', 'Susceptible', 'Exposed','Infectious',
                                    'death', 'case','infection','transmission rate',
                                    'latent period','infectious period','immunity period','Td.mean','Td.sd',
                                    'p.mob','infection detection rate','IFR') 
                         )

newVstat <-  lapply(files, function(x) {
  
  try(load(x))
  
  id <- gsub(dir_res,'',x) %>% 
    strsplit('_') %>% unlist
  
  d <- res.train$newVstat
  d2 <- res.train$hyp.best_all
  d = merge(d, d2, by = 'eval')
  
  d$loc = id[1] %>% strsplit('\\//') %>% unlist %>% tail(1)
  d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
  
  d
}) %>%
  rbindlist() 

