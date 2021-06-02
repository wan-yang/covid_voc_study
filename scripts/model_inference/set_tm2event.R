# set time to event valuables

dist_tm.to.detect.name = 'gamma'
dist_tm.to.hospital.name = 'gamma'
dist_tm.to.icu.name = 'gamma'
dist_tm.to.death.name = 'gamma'

# do not use relative numbers
tm.inf.to.onset = 1
tm.to.detect.mn = 6 + tm.inf.to.onset
tm.to.hospital.mn = 7 - 2 + tm.inf.to.onset
tm.to.icu.mn = 11 + tm.inf.to.onset  # from onset to icu: 11 days (IQR 7–14 days); plus +1 day for from viral sheding not onset



tm.to.diag = 4 # 5 days but b/c time from diag to death starts from 0, pad 4 lines
tm.to.deathFromDiag=T
if(tm.to.deathFromDiag){
  # loger time to death for other places
  t.death.adj = 5
  tm.to.death.mn = 14 
  tm.to.death.sd = 10 
  tm.to.death.max = 30 
  tm.from.inf.to.death.max = tm.to.diag + tm.to.death.max
} 

tm.to.detect.sd= 2
tm.to.hospital.sd = 3
tm.to.icu.sd = 5 # 3
tm.to.detect.max = 14
tm.to.hospital.max = 12 
tm.to.icu.max = 20 #



hospital.stay.mn.china = 11 # days
icu.stay.mn.china = 8 # days
percICU = .3
hospital.stay.mn.usICU = 6 + 21 + 14*2 # days
icu.stay.mn.usICU = 21 # 3 weeks!
hospital.stay.mn = hospital.stay.mn.china * (1-percICU) + hospital.stay.mn.usICU*percICU
icu.stay.mn = icu.stay.mn.usICU
hospital.stay.iqr = 7
icu.stay.iqr = 8
hospital.stay.sd = 5.2 # 2 # days
icu.stay.sd = 5.9 # days
hospital.stay.max = 45; # 11, IQR: 7 - 14
icu.stay.max = 15 * 2; # 8, IQR: 4 - 12


# update 11/27/20
vent.stay.mn = 8.6; # days [Mt Sinai study, overall from the five sites]
vent.stay.sd = 3; 
vent.stay.max = 21; # length of stay in icu
hospital.stay.mn.china.tot = 11 # days - including icu?
icu.stay.mn.china = 8 # days
percICU = .25
# assume that non-ICU pts spand x days, ICU pts spand icu.stay.mn.china + 2x days
# hospital.stay.mn.china = x * 75 + (8+2x) * .25 
# solve for x, x = 7.2 days for non-ICU pts
hospital.stay.mn.china = 7 # days - excluding icu, solved per above
icu.stay.mn.usICU = 13.3 # 2 weeks [Mt Sinai study, mean of the five sites]
hospital.stay.mn.usICU = 6 + icu.stay.mn.usICU + 8 # days
hospital.stay.mn = hospital.stay.mn.china * (1-percICU) + hospital.stay.mn.usICU*percICU
icu.stay.mn = icu.stay.mn.usICU
hospital.stay.iqr = 7
icu.stay.iqr = 8
hospital.stay.sd = 5.2 # 2 # days
icu.stay.sd = 5.9 # days
hospital.stay.max = 45; # 11, IQR: 7 - 14
icu.stay.max = 15 * 2; # 8, IQR: 4 - 12

# update 12/7  
vent.stay.mn = 7; 
percICU = .258
hospital.stay.mn.us.noICU = 5 # days
icu.stay.mn.usICU = 6 # possible error
icu.stay.mn.usICU = 11; # 
hospital.stay.mn.usICU = icu.stay.mn.usICU + 10 # days
hospital.stay.mn = hospital.stay.mn.us.noICU * (1-percICU) + hospital.stay.mn.usICU*percICU
icu.stay.mn = icu.stay.mn.usICU


ICU.adjNYC = 60 / (.04 /.15 * 100) # adjust for nyc
perc.seek.treat = c(.15, .25) #  + .05 # for nyc, probably higher
perc.seek.treat.mean = .2
severity.bounds = rbind(c(.015, .065),
                        c(.015, .065) * c(.2,.3), # 5/31/20
                        c(.005,.015),
                        c(.6,1) # update 11/27/20
                        ) 

# looks like % needing ventilator is higher
rownames(severity.bounds)=c('hospital','icu','death','vent')
severity.bounds
# severity

severity = t(lhs(num_ens,rect = severity.bounds))
rownames(severity)=rownames(severity.bounds)


tm.to.outcome.mn = t(lhs(num_ens,rect = rbind(tm.to.detect.mn+c(-1,1),
                                              tm.to.hospital.mn+c(-1,1),
                                              tm.to.icu.mn+c(-2,2),
                                              tm.to.death.mn+c(-1,1),
                                              hospital.stay.mn+c(-2,2),
                                              icu.stay.mn+c(-1,1),
                                              vent.stay.mn+c(-1,1))))
rownames(tm.to.outcome.mn) = c('tm.to.detect','tm.to.hospital','tm.to.icu','tm.to.death','hospital.stay','icu.stay','vent.stay')

tm.to.outcome.sd = t(lhs(num_ens,rect = rbind(tm.to.detect.sd+c(-1,1),
                                              tm.to.hospital.sd+c(-1,1),
                                              tm.to.icu.sd+c(-1,1),
                                              tm.to.death.sd+c(-1,1),
                                              hospital.stay.sd+c(-1,1),
                                              icu.stay.sd+c(-1,1),
                                              vent.stay.sd+c(-1,1))))
rownames(tm.to.outcome.sd) = c('tm.to.detect','tm.to.hospital','tm.to.icu','tm.to.death','hospital.stay','icu.stay','vent.stay')

