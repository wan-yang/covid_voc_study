# combine results from multiple runs - diff eval method
state.names = c('S1', 'E1','I1',
                'death1', 'newIobs1','newItot1','beta',
                'Tei','Tir','Trs','Td.mean','Td.sd',
                'p.mob','alpha','ifr')
  
files = list.files(dir_res, full.names = T)
files = files[grepl('.RData',files)]
files = files[grepl('train', files)]

# compare performance across runs?


cumIperc_ens <-  lapply(files, function(x) {
  
  try(load(x))
  
  id <- gsub(dir_res,'',x) %>% 
    strsplit('_') %>% unlist
  
  # get cumIens
  tmp = res.train$xpost_mean
  d = res.train$cumIperc_ens
  
  d$date = tmp$Week.start
  d$loc = id[1] %>% strsplit('\\//') %>% unlist %>% tail(1)
  d$run <- tail(id,1) %>% strsplit('\\.') %>% unlist %>% head(1) %>% gsub(pattern = 'r', replacement =  '') %>% as.integer()
  
  d
}) %>%
  rbindlist() 

# reshape
tmp = melt(cumIperc_ens, id.vars = c('eval','loc','run','date'))
cumIperc_ens = dcast(tmp, eval + loc + date ~ run + variable)
rm(tmp)

