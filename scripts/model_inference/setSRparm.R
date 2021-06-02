# set SRparm

SRparms = list(Spb.adj_lwr = .8,
               Spb.adj_upr = 1.2, 
               # minor change
               p.cum.dS.large.minor = .9, # .75,
               p.cum.dS.median.minor = .9, # .8,
               p.cum.dS.small.minor = .9, # .9,
               s.upr.t.minor = 1, s.lwr.t.minor = .5,
               Smajor.cut.lwr.minor = -.1,  # maybe too high, some need to wait till the 3rd time
               Smajor.cut.upr.minor = 2,
               Smajor.cnt.cut.minor = 2, # 1.5 
               cntSadj_large.tot.minor = 7,
               cntSadj_median.tot.minor = 7,
               cntSadj_small.tot.minor = 7,
               cntSadjtot.cut.minor = 10,
               Spb_large.minor = .2, # .15,
               Spb_median.minor = .15, # .1, # .2,
               Spb_small.minor = .1, # .1,
               cum.dS.cut.minor = .5,
               # major change - s.only
               p.cum.dS.large.major = .85, #  .85,
               p.cum.dS.median.major = .85, # .85,
               p.cum.dS.small.major = .9, # .9,
               s.upr.t.major = 1, s.lwr.t.major = .5,
               Smajor.cut.lwr.major = -.1,  # maybe too high, some need to wait till the 3rd time
               Smajor.cut.upr.major = 2,
               Smajor.cnt.cut.major = 3, # 1.5 
               cntSadj_large.tot.major = 15,
               cntSadj_median.tot.major = 15,
               cntSadj_small.tot.major = 15,
               cntSadjtot.cut.major = 20,
               Spb_large.major = .25,
               Spb_median.major = .2,
               Spb_small.major = .15,
               cum.dS.cut.major = .95,  # make it slightly higher than 1? 1.5 too much? 
               # major change - both major
               p.cum.dS.large.max = .85, # .85,
               p.cum.dS.median.max = .85, # .85,
               p.cum.dS.small.max = .85, # .9,
               s.upr.t.max = 1, s.lwr.t.max = .5,
               Smax.cut.lwr.max = -.1,  # maybe too high, some need to wait till the 3rd time
               Smax.cut.upr.max = 2,
               Smax.cnt.cut.max = 3, # 1.5 
               cntSadj_large.tot.max = 15,
               cntSadj_median.tot.max = 15,
               cntSadj_small.tot.max = 15,
               cntSadjtot.cut.max = 20,
               Spb_large.max = .25,
               Spb_median.max = .25,
               Spb_small.max = .2,
               cum.dS.cut.max = .95,  # make it slightly higher than 1? 1.5 too much? 
               # more flexible change
               s.upr.t.flex = 1, s.lwr.t.flex = .5,
               p.cum.dS.large.flex = .85, # .8,
               p.cum.dS.median.flex = .85,
               p.cum.dS.small.flex = .9,
               Smajor.cut.lwr.flex = -.1,  # maybe too high, some need to wait till the 3rd time
               Smajor.cut.upr.flex = 2,
               Smajor.cnt.cut.flex = 2.5, # 1.5 
               cntSadj_large.tot.flex = 15,
               cntSadj_median.tot.flex = 15,
               cntSadj_small.tot.flex = 15,
               cntSadjtot.cut.flex = 20, # 10
               Spb_large.flex = .25,
               Spb_median.flex = .15,
               Spb_small.flex = .1,
               cum.dS.cut.flex = .95,
               # slow change for a large spatial place like br
               p.cum.dS.large.slow = .9, # .75,
               p.cum.dS.median.slow = .9, # .8,
               p.cum.dS.small.slow = .9, # .9,
               s.upr.t.slow = 1, s.lwr.t.slow = .5,
               Smajor.cut.lwr.slow = 0,  # maybe too high, some need to wait till the 3rd time
               Smajor.cut.upr.slow = 2,
               Smajor.cnt.cut.slow = 2, # 1.5 
               cntSadj_large.tot.slow = 30,
               cntSadj_median.tot.slow = 30,
               cntSadj_small.tot.slow = 30,
               cntSadjtot.cut.slow = 30,
               Spb_large.slow = .08, # .1 seems too strong .12, # .1,
               Spb_median.slow = .07, # .1, # .08, # .2,
               Spb_small.slow = .06, # .08, # .05, # .1,
               cum.dS.cut.slow = .95,
               # restrict the level of change on S
               # slow change for a large spatial place like br
               p.cum.dS.large.slow.minor = .9, # .75,
               p.cum.dS.median.slow.minor = .9, # .8,
               p.cum.dS.small.slow.minor = .9, # .9,
               s.upr.t.slow.minor = 1, s.lwr.t.slow.minor = .5,
               Smajor.cut.lwr.slow.minor = 0,  # maybe too high, some need to wait till the 3rd time
               Smajor.cut.upr.slow.minor = 2,
               Smajor.cnt.cut.slow.minor = 2, # 1.5 
               cntSadj_large.tot.slow.minor = 30,
               cntSadj_median.tot.slow.minor = 30,
               cntSadj_small.tot.slow.minor = 30,
               cntSadjtot.cut.slow.minor = 30,
               Spb_large.slow.minor = .08, # .1 seems too strong .12, # .1,
               Spb_median.slow.minor = .07, # .1, # .08, # .2,
               Spb_small.slow.minor = .06, # .08, # .05, # .1,
               cum.dS.cut.slow.minor = .5
)

