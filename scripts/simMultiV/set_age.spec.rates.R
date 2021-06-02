# set age specific rates based on nyc data
# 4/11/2020

# adjust time to death as for elderly 65+, time from diagnosis is ~3 days not 5
if(ia>=7){
  tm.to.diag = tm.to.diag - 2
  if(tm.to.deathFromDiag){
    tm.to.death.mn = 9.36 #  ~5 days from infectious to diagnosis
    tm.to.death.sd = 9.76 #  used before 5/15/20
    tm.to.death.max = 26 # tm.to.detect.max + diff.dd
    tm.from.inf.to.death.max = tm.to.diag + tm.to.death.max
  } else {
    tm.to.death.mn = tm.to.diag + 9.36 #  ~5 days from infectious to diagnosis
    tm.to.death.sd = 9.76 #  used before 5/15/20
    tm.to.death.max = tm.to.diag + 26 # tm.to.detect.max + diff.dd
    tm.from.inf.to.death.max = tm.to.death.max
  }
}



f.adj.hosps = c(0.0296, 0.0296, 0.0296, 0.0592, 0.2057, 0.6995, 2.5427, 3.9931)
f.adj.icus = c(0.0030, 0.0030, 0.0030, 0.0059, 0.2057, 0.6995, 2.0341, 1.5973)
f.adj.deaths = c(0.0100,  0.0100,  0.0100,  0.0090,  0.0758,  0.8077,  6.2538, 14.8883)

f.adj.deaths[which(f.adj.deaths==0)] = 1e-2


f.adj.death = f.adj.deaths[ia]

f.adj.hosp = f.adj.hosps[ia]

f.adj.icu = f.adj.icus[ia]

r.hosp.lwr = f.adj.hosps * .25 * .15
r.hosp.upr = f.adj.hosps * .25 * .3
r.icu.lwr = f.adj.icus * .25 * .15 * .25
r.icu.upr = f.adj.icus * .25 * .3 * .5
r.icu.lwr /r.hosp.lwr
r.icu.upr /r.hosp.upr

f.adj.alphas = c(0.2168, 0.1806, 0.1806, 0.2129, 0.2300, 0.2600, 0.3986, 0.4591)
f.adj.alpha = f.adj.alphas[ia]  # ifelse(ia %in% 7:8, 3, 1)

f.adj.betas = c(0.4365, 0.5153, 0.7976, 0.8390, 1.3827, 0.9965, 0.2820, 0.3952)

f.adj.beta0 = f.adj.betas[ia]
f.adj.beta = pmin(1.5, pmax(.5, f.adj.betas[ia])) # may be too high for groups with high contact rate

severity['hospital',] = severity['hospital',] * f.adj.hosp
severity['icu',] = severity['icu',] * f.adj.icu


if(ia %in% 1:4){
  ifr_bounds = c(0.00005, 0.00015)
} else if(ia %in% 5){
  ifr_bounds = c(0.0005, 0.0015)
} else if(ia %in% 6){
  ifr_bounds = c(0.005, 0.015)
} else if(ia %in% 7){
  ifr_bounds = c(0.01, 0.1)
} else if(ia %in% 8){
  ifr_bounds = c(0.02, 0.2)
}

# set ed visit rate
if(ia %in% 1:4){
  edr_bounds = c(0.001, 0.02)  # .01
} else if(ia %in% 5){
  edr_bounds = c(0.003, 0.03)
} else if(ia %in% 6){
  edr_bounds = c(0.006, 0.06)
} else if(ia %in% 7){
  edr_bounds = c(0.01, 0.15) 
} else if(ia %in% 8){
  edr_bounds = c(0.02, 0.25)
}

severity['death',] = runif(num_ens, min=ifr_bounds[1], max=ifr_bounds[2])
# add edr
edr = runif(num_ens, min=edr_bounds[1], max=edr_bounds[2])
severity = rbind(severity, edr)


if(ia ==6){
  tm.to.outcome.mn[c("hospital.stay","icu.stay","vent.stay"),] = tm.to.outcome.mn[c("hospital.stay","icu.stay","vent.stay"),] + 1
} else if (ia > 6){
  tm.to.outcome.mn[c("hospital.stay","icu.stay","vent.stay"),] = tm.to.outcome.mn[c("hospital.stay","icu.stay","vent.stay"),] + 
    matrix(c(5, 2, 5), 3, num_ens)
}
