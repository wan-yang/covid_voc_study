# compare covid pandemic in the UK, SA, BR
# 2/23/21

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
library(tidyverse)

dir_data = '../../data/'
dir_proj = '../../data/'

if(! file.exists(dir_proj)) dir.create(dir_proj, recursive = T)

url.death = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv'
url.case ='https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv'


death = read_csv(url(url.death)) %>% data.table() 
case = read_csv(url(url.case)) %>% data.table() 

d.death = death[`Country/Region`  %in% c("United Kingdom", "South Africa", "Brazil") & is.na(`Province/State`)] %>% t
d.case = case[`Country/Region`  %in% c("United Kingdom", "South Africa", "Brazil") & is.na(`Province/State`)] %>% t
d.death = data.table(date = rownames(d.death), d.death)
colnames(d.death)[2:4] = d.death[date == 'Country/Region'] %>% unlist %>% tail(3)
d.death = d.death[!(date %in% c('Province/State','Lat','Long','Country/Region'))]
d.death$date = d.death$date %>% as.Date('%m/%d/%y')
# make it new case
tmp = d.death[,-1] %>%
  mutate_all(type.convert) %>%
  mutate_if(is.factor, as.numeric) 
tmp2 = rbind(t(rep(0,3)), d.death[-nrow(d.death),-1], use.names = F) %>%
  mutate_all(type.convert) %>%
  mutate_if(is.factor, as.numeric) 
d.death = data.table(date = d.death$date, tmp - tmp2)

# case
d.case = data.table(date = rownames(d.case), d.case)
colnames(d.case)[2:4] = d.case[date == 'Country/Region'] %>% unlist %>% tail(3)
d.case = d.case[!(date %in% c('Province/State','Lat','Long','Country/Region'))]
d.case$date = d.case$date %>% as.Date('%m/%d/%y')
# make it new case
tmp = d.case[,-1] %>%
  mutate_all(type.convert) %>%
  mutate_if(is.factor, as.numeric) 
tmp2 = rbind(t(rep(0,3)), d.case[-nrow(d.case),-1], use.names = F) %>%
  mutate_all(type.convert) %>%
  mutate_if(is.factor, as.numeric) 

d.case = data.table(date = d.case$date, tmp - tmp2)

cols = c('blue','darkgreen','orange')
matplot(d.case[,-1,with=F],type='l',lty = 1, col = cols, xaxt='n', xlab='', ylab='Daily cases')
axis(1, at = 1:nrow(d.case), labels = d.case$date)

matplot(d.death[,-1,with=F],type='l',lty = 1, col = cols, xaxt='n', xlab='', ylab='Daily deaths')
axis(1, at = 1:nrow(d.death), labels = d.death$date)

# make it weekly
tmp = MMWRweek(d.case$date) %>% setnames(paste0('MMWR',c('year','week','day')),c('year','week','day'))
d.case = cbind(tmp, d.case) %>% data.table()
d.case = d.case[,lapply(.SD, sum), by = c('year','week'), .SDcols=c("United Kingdom", "South Africa", "Brazil")]
d.case$date = MMWRweek2Date(MMWRyear = d.case$year, MMWRweek = d.case$week, MMWRday = 1)
# exclude last week with incompete data
d.case = d.case[as.Date(date) < Sys.Date()-10]

d.death = cbind(tmp, d.death) %>% data.table()
d.death = d.death[,lapply(.SD, sum), by = c('year','week'), .SDcols=c("United Kingdom", "South Africa", "Brazil")]
d.death$date = MMWRweek2Date(MMWRyear = d.death$year, MMWRweek = d.death$week, MMWRday = 1)
# exclude last week with incompete data
d.death = d.death[as.Date(date) < Sys.Date()-10]

par(mfrow=c(2,1), mar=c(2,2,1,1), oma=c(0,0,0,0), cex=.9, cex.lab=.9, cex.axis =.9, mgp=c(1.1,.2,0), tck=-.02)
cols = c('blue','darkgreen','orange')
matplot(d.case[,-1,with=F],type='l',lty = 1, col = cols, xaxt='n', xlab='', ylab='Daily cases')
axis(1, at = 1:nrow(d.case), labels = d.case$date)

matplot(d.death[,-1,with=F],type='l',lty = 1, col = cols, xaxt='n', xlab='', ylab='Daily deaths')
axis(1, at = 1:nrow(d.death), labels = d.death$date)

# convert to per population size
pop.uk = 68207116
pop.sa = 60041994
pop.br = 213993437
N = 1e6; # make it per 1 M for all

d = merge(d.case, d.death, by = c('date','year','week'), suffixes = c('.c','.d'))
d.uk = data.table(d[,c('date','year','week')], d[,paste0('United Kingdom',c('.c','.d')),with=F] / pop.uk * N)
d.sa = data.table(d[,c('date','year','week')], d[,paste0('South Africa',c('.c','.d')),with=F] / pop.sa * N)
d.br = data.table(d[,c('date','year','week')], d[,paste0('Brazil',c('.c','.d')),with=F] / pop.br * N)

d = merge(d.uk, d.sa, by = c('date','year','week'))
d = merge(d, d.br, by = c('date','year','week'))

pdf(paste0(dir_proj,'Fig_cp_case.death.uk.sa.br.pdf'),width = 6, height = 6)
par(mfrow=c(2,1), mar=c(2,2.5,1,1), oma=c(0,0,0,0), cex=1, cex.lab=1, cex.axis =.9, mgp=c(1.1,.2,0), tck=-.02)
cols = c('blue','darkgreen','orange')
matplot(d[,paste0(c("United Kingdom", "South Africa", "Brazil"),'.c'),with=F],type='l',lwd=2, lty = 1, col = cols, xaxt='n', xlab='', ylab=paste('Weekly cases per',N/1e6,"million"))
axis(1, at = 1:nrow(d), labels = d$date)
legend('topleft',c("United Kingdom", "South Africa", "Brazil"), lwd=2, lty=1, col = cols, bty='n')
matplot(d[,paste0(c("United Kingdom", "South Africa", "Brazil"),'.d'),with=F],type='l',lwd=2,lty = 1, col = cols, xaxt='n', xlab='', ylab=paste('Weekly deaths per',N/1e6,"million"))
axis(1, at = 1:nrow(d), labels = d$date)
dev.off()

setnames(d, c("United Kingdom.c", "United Kingdom.d","South Africa.c","South Africa.d","Brazil.c","Brazil.d"),
         c('case.uk','death.uk','case.sa','death.sa','case.br','death.br'))
write.csv(d, paste0(dir_proj, 'da_case.death.per',N/1e6,'M_uk.sa.br.csv'), row.names = F)

# look at mobility
mob.type.bus = c('retail_and_recreation_percent_change_from_baseline',
  'transit_stations_percent_change_from_baseline',
  'workplaces_percent_change_from_baseline')  
# ,'parks_percent_change_from_baseline',
# 'residential_percent_change_from_baseline'
mob.type.full = c('retail_and_recreation_percent_change_from_baseline',
                  'grocery_and_pharmacy_percent_change_from_baseline',
                  'parks_percent_change_from_baseline',
                  'transit_stations_percent_change_from_baseline',
                  'workplaces_percent_change_from_baseline',
                  'residential_percent_change_from_baseline')

url.mob = 'https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv'
d.mob = read_csv(url(url.mob)) %>% data.table() 
# d.mob = read.csv(paste0(dir_proj,'Global_Mobility_Report.csv'), stringsAsFactors = F) %>% data.table()

d.br = d.mob[country_region %in% c("Brazil")]
d.br = d.br[, lapply(.SD, median, na.rm=T), by = 'date', 
            .SD = mob.type.full]
d.br$mob.bus = rowMeans(d.br[,mob.type.bus,with=F])  # take the average for the 
d.br$mob.full = rowMeans(d.br[,mob.type.full,with=F]) 
matplot(d.br[,c('mob.bus','mob.full')], type='l', lty =1, col = c('blue','black'))
abline(h=0); legend('bottomleft',c('business','all'), lty=1, col = c('blue','black'), bty='n')
# SA
d.sa = d.mob[country_region %in% c("South Africa")]
d.sa = d.sa[, lapply(.SD, median, na.rm=T), by = 'date', 
            .SD = mob.type.full]
d.sa$mob.bus = rowMeans(d.sa[,mob.type.bus,with=F])  # take the average for the 
d.sa$mob.full = rowMeans(d.sa[,mob.type.full,with=F]) 
matplot(d.sa[,c('mob.bus','mob.full')], type='l', lty =1, col = c('blue','black'))
abline(h=0); legend('bottomleft',c('business','all'), lty=1, col = c('blue','black'), bty='n')


# UK
d.uk = d.mob[country_region %in% c("United Kingdom")]
d.uk = d.uk[, lapply(.SD, median, na.rm=T), by = 'date', 
            .SD = mob.type.full]
d.uk$mob.bus = rowMeans(d.uk[,mob.type.bus,with=F])  # take the average for the 
d.uk$mob.full = rowMeans(d.uk[,mob.type.full,with=F]) 
matplot(d.uk[,c('mob.bus','mob.full')], type='l', lty =1, col = c('blue','black'))
abline(h=0); legend('bottomleft',c('business','all'), lty=1, col = c('blue','black'), bty='n')


d.mob = merge(d.uk[,c('date','mob.bus','mob.full')], d.sa[,c('date','mob.bus','mob.full')], by = 'date', suffixes = c('.uk','.sa'))
d.mob = merge(d.mob, d.br[,c('date','mob.bus','mob.full')] %>% 
                setnames(c('mob.bus','mob.full'),paste0(c('mob.bus','mob.full'),'.br')), by = 'date')
write.csv(d.mob, paste0(dir_proj,'da_mob_uk.sa.br.csv'), row.names = F)

# make it weekly
tmp = MMWRweek(d.mob$date) %>% setnames(paste0('MMWR',c('year','week','day')),c('year','week','day'))
d.mob = cbind(tmp, d.mob) %>% data.table()
d.mob = d.mob[,lapply(.SD, mean, na.rm=T), by = c('year', 'week'), 
              .SDcols = c(paste0(c('mob.bus','mob.full'),'.uk'),
                          paste0(c('mob.bus','mob.full'),'.sa'),
                          paste0(c('mob.bus','mob.full'),'.br'))]
d.mob$date = MMWRweek2Date(MMWRyear = d.mob$year, MMWRweek = d.mob$week, MMWRday = 1)

d = merge(d, d.mob, all = T, by = c('date','year', 'week'))
write.csv(d, paste0(dir_proj, 'da_case.death.mob_uk.sa.br.csv'), row.names = F)

# compile vaccination data
fn_getDailyVx = function(date.t, da.vx){
  # distribute vx weekly dose to daily dose
  da.vx[which(date.start<= date.t & date.end >=date.t)]$n.daily
}
lagV1 = 14
lagV2 = 7
loc.t = 'uk' # dates are week ending
vx.day1 = as.Date('2020/12/08');
vx = read.csv(paste0(dir_data,'vx_',loc.t,'_weekly.csv')) %>% data.table() %>% 
  setnames(c("weeklyPeopleVaccinatedFirstDoseByVaccinationDate","weeklyPeopleVaccinatedSecondDoseByVaccinationDate"),
           c('n.v1','n.v2'))
vx$date = vx$date %>% as.Date
# convert to per million
vx$n.v1 = vx$n.v1 / pop.uk * N
vx$n.v2 = vx$n.v2 / pop.uk * N
vx1 = vx[,c('date','n.v1'), with=F]
vx1 = vx1[order(date)]
vx1$date.start = as.Date(vx1$date) - 6
vx1$date.end = as.Date(vx1$date)
vx1[date.start < vx.day1]$date.start = vx.day1
vx1$nday = (vx1$date.end - vx1$date.start + 1) %>% as.numeric()
vx1$n.daily = vx1$n.v1 / vx1$nday
vx1daily = data.table(date = seq(vx.day1, max(as.Date(vx1$date)),by = 'day'), n.v1 = 0)
# vx1daily = merge(vx1daily, vx1, all = T, by = 'date')
vx1daily = vx1daily[,list(n.v1 = fn_getDailyVx(date.t=date, da.vx=vx1)), by = date]
vx1daily$date = vx1daily$date + lagV1

vx2 = vx[,c('date','n.v2'), with=F]
vx2 = vx2[order(date)]
vx2$date.start = as.Date(vx2$date) - 6
vx2$date.end = as.Date(vx2$date)
vx2[date.start < vx.day1]$date.start = vx.day1
vx2$nday = (vx2$date.end - vx2$date.start + 1) %>% as.numeric()
vx2$n.daily = vx2$n.v2 / vx2$nday
vx2daily = data.table(date = seq(vx.day1, max(as.Date(vx2$date)),by = 'day'), n.v2 = 0)
# vx2daily = merge(vx2daily, vx2, all = T, by = 'date')
vx2daily = vx2daily[,list(n.v2 = fn_getDailyVx(date.t=date, da.vx=vx2)), by = date]
vx2daily$date = vx2daily$date + lagV2
vx = merge(vx1daily, vx2daily, all=T,by='date')
vx[is.na(vx)] = 0
vx.p1 = vx

# add new data
vx = read.csv(paste0(dir_data,'vx_',loc.t,'_daily_042821.csv')) %>% data.table() %>% 
  setnames(c("newPeopleVaccinatedFirstDoseByPublishDate","newPeopleVaccinatedSecondDoseByPublishDate"),
           c('n.v1','n.v2'))
vx$date = vx$date %>% as.Date
# convert to per million
vx$n.v1 = vx$n.v1 / pop.uk * N
vx$n.v2 = vx$n.v2 / pop.uk * N
vx1daily = vx[,c('date','n.v1'), with=F]
vx1daily = vx1daily[order(date)]
vx1daily$date = vx1daily$date + lagV1
vx2daily = vx[,c('date','n.v2'), with=F]
vx2daily = vx2daily[order(date)]
vx2daily$date = vx2daily$date + lagV2
vx = merge(vx1daily, vx2daily, all=T,by='date')
vx[is.na(vx)] = 0
vx.p2 = vx

idx.t = which(vx.p2$n.v1==0) %>% tail(1)
vx = rbind(vx.p1[date <= vx.p2[idx.t]$date], vx.p2[-c(1:idx.t),])

write.csv(vx, paste0(dir_data,'vx.lagged.per1Mpop_',loc.t,'.csv'), row.names = F)

# vaccination for SA
loc.t = 'sa'
url.t = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/country_data/South%20Africa.csv'
vx = read_csv(url(url.t)) %>% data.table()
vx = vx[people_vaccinated > 0]
vx$date = vx$date %>% as.Date
dates.t = data.table(date = seq(min(vx$date), max(vx$date), by = 'day'))
vx = merge(dates.t, vx, all = T, by = 'date')
vx$idx = 1:nrow(vx)
vx$daily.vx = 0
idx.all = which(is.na(vx$people_vaccinated))
consecutive =  idx.all[-1] - idx.all[-length(idx.all)]
i.div = which(consecutive >2)
w.div = idx.all[which(consecutive >2)]
grps = list()
if(length(i.div) == 1){
  grps[[1]] = idx.all[1]: w.div[1]
  grps[[2]] =idx.all[i.div[1]+1]: tail(idx.all,1)
} else {
  for(id in 1: length(i.div)){
    if(id == 1){
      grps[[id]] = idx.all[1]: w.div[id]
    } else if (id == length(i.div)){
      # both before and after
      
      grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]]
      grps[[id+1]] = idx.all[i.div[id]+1]: tail(idx.all,1)
    } else {
      grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]] # (idx.all[i.div[id-1]+1]) : w.div[id]
    }
    
  }
}
for(ig in 1: length(grps)){
  idx0 = grps[[ig]]
  
  n.add = ifelse(length(idx0) < 2, 1, 1)
  idx = seq((idx0[1] - n.add) %>% pmax(1), 
          (tail(idx0,1)+ n.add) %>% pmin(nrow(vx)),
          by = 1
          ) # add an additional week in case there is back adjustment
  vx0 = vx[idx]; vx0 = vx0[complete.cases(vx0)]
  fit = lm(people_vaccinated ~ idx, data = vx0)
  vx$daily.vx[idx0] = predict(fit, newdata = vx[idx0])
}
vx[is.na(people_vaccinated)]$people_vaccinated = vx[is.na(people_vaccinated)]$daily.vx
vx$n.v1 = c(vx$people_vaccinated[1], (vx$people_vaccinated[-1] - vx$people_vaccinated[-nrow(vx)]) %>% round(0)) / pop.sa * N
vx$n.v2 = 0  / pop.sa * N

lagV1 = 14 # for J&J
vx$date = vx$date + lagV1
vx = vx[,c('date', 'n.v1', 'n.v2'),with=F]
write.csv(vx, paste0(dir_data,'vx.lagged.per1Mpop_',loc.t,'.csv'), row.names = F)

loc.t = 'br'
url.t = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/country_data/Brazil.csv'
vx = read_csv(url(url.t)) %>% data.table()
vx = vx[people_vaccinated > 0]
vx[is.na(people_fully_vaccinated)]$people_fully_vaccinated = 0
vx$date = vx$date %>% as.Date
dates.t = data.table(date = seq(min(vx$date), max(vx$date), by = 'day'))
vx = merge(dates.t, vx, all = T, by = 'date')
vx$idx = 1:nrow(vx)
vx$daily.vx1 = 0
vx$daily.vx2 = 0
idx.all = which(is.na(vx$people_vaccinated))
consecutive =  idx.all[-1] - idx.all[-length(idx.all)]
i.div = which(consecutive >2)
w.div = idx.all[which(consecutive >2)]
grps = list()
if(length(i.div) == 1){
  grps[[1]] = idx.all[1]: w.div[1]
  grps[[2]] =idx.all[i.div[1]+1]: tail(idx.all,1)
} else {
  for(id in 1: length(i.div)){
    if(id == 1){
      grps[[id]] = idx.all[1]: w.div[id]
    } else if (id == length(i.div)){
      # both before and after
      
      grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]]
      grps[[id+1]] = idx.all[i.div[id]+1]: tail(idx.all,1)
    } else {
      grps[[id]] = idx.all[(i.div[id-1]+1): i.div[id]] # (idx.all[i.div[id-1]+1]) : w.div[id]
    }
    
  }
}
for(ig in 1: length(grps)){
  idx0 = grps[[ig]]
  
  n.add = ifelse(length(idx0) < 2, 1, 1)
  idx = seq((idx0[1] - n.add) %>% pmax(1), 
            (tail(idx0,1)+ n.add) %>% pmin(nrow(vx)),
            by = 1
  ) # add an additional week in case there is back adjustment
  vx0 = vx[idx]; vx0 = vx0[complete.cases(vx0)]
  fit1 = lm(people_vaccinated ~ idx, data = vx0)
  vx$daily.vx1[idx0] = predict(fit1, newdata = vx[idx0])
  fit2 = lm(people_fully_vaccinated ~ idx, data = vx0)
  vx$daily.vx2[idx0] = predict(fit2, newdata = vx[idx0])
}
vx[is.na(people_vaccinated)]$people_vaccinated = vx[is.na(people_vaccinated)]$daily.vx1
vx[is.na(people_fully_vaccinated)]$people_fully_vaccinated = vx[is.na(people_fully_vaccinated)]$daily.vx2

vx$n.v1 = c(vx$people_vaccinated[1], (vx$people_vaccinated[-1] - vx$people_vaccinated[-nrow(vx)]) %>% round(0)) / pop.br * N
vx$n.v2 =c(vx$people_fully_vaccinated[1], (vx$people_fully_vaccinated[-1] - vx$people_fully_vaccinated[-nrow(vx)]) %>% round(0)) / pop.br * N

lagV1 = 14 # for sinovav
lagV2 = 7

vx1daily = vx[,c('date','n.v1'), with=F]
vx1daily = vx1daily[order(date)]
vx1daily$date = vx1daily$date + lagV1
vx2daily = vx[,c('date','n.v2'), with=F]
vx2daily = vx2daily[order(date)]
vx2daily$date = vx2daily$date + lagV2
vx = merge(vx1daily, vx2daily, all=T,by='date')
vx[is.na(vx)] = 0
write.csv(vx, paste0(dir_data,'vx.lagged.per1Mpop_',loc.t,'.csv'), row.names = F)

