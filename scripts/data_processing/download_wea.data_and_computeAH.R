# download data from national weather sevices
# and compute the specific humidity

library(data.table)
library(magrittr)
library(stationaRy)
library(MMWRweek)
library(beepr)

dir_data = '../../data/'
dir_code = '../../scripts/'
dir_res = '../../results/'

{
  fn_getDaily = function(d){
    # note variable 'spec.hum' is the specific humidity (a measure of absolute humidity), unit: kg/kg
    parms=c("temp",'spec.hum',"rh" ,"dew_point", "atmos_pres")
    
    es_t0=611; # saturation vapor presure at 273.15K, Pa;
    t0=273.15; # K
    L=2.5e6; # latent heat of avaporation for water, J*kg^-1
    Rv=461.5; # gas constant for water vapor, J*kg^-1*K^-1
    Rd=286.9; # gas constant for dry air, J*kg^-1*K^-1
    Pa=101325; # atmospheric pressure, Pa;
    # tt=20; # rh=40;
    # es_t=es_t0*exp(L/Rv*(1/t0-1/(t0+tt))); # saturation vp
    # ea_t=es_t*rh/100; # actual vp
    # mr=Rd/Rv*ea_t/(Pa-ea_t); # mixing ratio
    # qmn=mr/(1+mr); # specific humidity 0.005822275
    es_t=es_t0*exp(L/Rv*(1/t0-1/(t0+as.numeric(d$temp)))); # saturation vp
    ea_t=es_t*as.numeric(d$rh)/100; # actual vp
    mr=Rd/Rv*ea_t/(Pa-ea_t); # mixing ratio
    d$spec.hum=mr/(1+mr); # specific humidity 0.005822275
    
    # aggregate the daily resolution
    d$date = format(d$time, '%Y/%m/%d')
    d = d %>% data.table()
    d = d[, lapply(.SD, mean, na.rm =T), by = c('id','date'), .SDcols = parms]
    d
  }
  fn_getPlot = function(res.d, res.w, loc.t, nsta, yrs.t){
    plot(res.d$day, res.d$temp, type = 'l', lwd = 2, ylab = 'Temperature (C)', xlab='', xaxt = 'n')
    plot(res.d$day, res.d$spec.hum, type = 'l', lwd = 2, ylab = 'Specific humidity', xlab='', xaxt = 'n')
    plot(res.d$day, res.d$rh, type = 'l', lwd = 2, ylab = 'Relative humidity (%)', xlab='', xaxt = 'n')
    axis(1)
    mtext('Day of year', outer = F, side = 1, line = 1.2)
    
    plot(res.w$week, res.w$temp, type = 'l', lwd = 2, ylab = 'Temperature (C)', xlab='', xaxt = 'n')
    plot(res.w$week, res.w$spec.hum, type = 'l', lwd = 2, ylab = 'Specific humidity', xlab='', xaxt = 'n')
    plot(res.w$week, res.w$rh, type = 'l', lwd = 2, ylab = 'Relative humidity (%)', xlab='', xaxt = 'n')
    mtext('Week of year', outer = F, side = 1, line = 1.2)
    axis(1)
    mtext(paste0(loc.t,' (n=',nsta,' stations; ', yrs.t[1],'-',tail(yrs.t, 1), ')'), outer = T, side = 3, line = .5)
    
  }
  
}
years.t = 2000:2021
stations = get_station_metadata() %>% data.table()
parms=c("temp",'spec.hum',"rh" ,"dew_point", "atmos_pres")

for(loc.t in c('UK','SF','BR')){ 
  
  if(loc.t == 'UK'){
    years.t = 2000:2021
  } else {
    years.t = 2000:2021
  }
  
  sta.t = stations[country == loc.t]
  
  print(paste(loc.t, ':', nrow(sta.t),'stations'), quote = F)
  
  res = NULL
  for(i in 1:nrow(sta.t)){
    print(paste(loc.t, i), quote = F)
    # i = 411
    da.t = try(get_met_data(station_id = sta.t[i]$id, 
                        make_hourly = F,
                        # too many stations and data, restrict to the most recent 10 years
                        years = years.t
                        ))
    if(class(da.t) == 'try-error'){
      print('no data')
      next
    }
    if(nrow(da.t) < 1) 
      next;
    
    # compute daily values to reduce data size
    da.t = fn_getDaily(da.t)
    
    res = rbind(res, da.t)
  }
  
  # aggregate to country level
  
  res$day = res$date %>% strftime(format = "%j") %>% as.numeric()
  tmp = res$date %>% MMWRweek 
  res$week = tmp[,'MMWRweek'] %>% as.numeric()
  res$year = res$date %>% strftime(format = '%Y') %>% as.numeric()
  # aggregate by day
  res.d = res[, lapply(.SD, mean, na.rm =T), by = c('day'), .SDcols = parms]
  res.d = res.d[order(day)]
  
  # aggregate by week
  res.w = res[, lapply(.SD, mean, na.rm =T), by = c('week'), .SDcols = parms]
  res.w = res.w[order(week)]
  
  write.csv(res.d, paste0(dir_data, 'wea.by.day_', loc.t, '.csv'), row.names = F)
  write.csv(res.w, paste0(dir_data, 'wea.by.week_', loc.t, '.csv'), row.names = F)
  

}

beep(3); beep(3)

