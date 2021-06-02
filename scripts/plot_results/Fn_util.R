# some utility functions for analysis
# 2/5/20


fn_format= function(x,roundit=T,roundigt=0){
  x=unlist(x)
  if(roundit==T){
    x=round(x,roundigt)
  }
  paste0(x[1],' (',x[2],', ',x[3],')')
}
fn_formatCI= function(x,roundit=T,roundigt=0){
  x=unlist(x)
  if(roundit==T){
    x=round(x,roundigt)
  }
  paste0('(',x[1],', ',x[2],')')
}
