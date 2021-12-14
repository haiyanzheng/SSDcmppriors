source("Bayesian SSD using robust commensurate priors.R")


library(sn)
#-----------------------------------------------------------------------------#
TestCvr = function(Sc, N = c(10, 10), simA, simB, true.diff, sig02,
                    l0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  smpdiff = mean(simA) - mean(simB)
  
  eta = cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/length(simA)+1/length(simB))*sig02)) + 
    smpdiff/(1+(1/length(simA)+1/length(simB))*sig02/cmsprior[,"v.cp"])
  sig2eta = 1/(1/cmsprior[,"v.cp"] + 1/((1/length(simA)+1/length(simB))*sig02))
  
  cvrInd = eta+l0/2 >= true.diff & eta-l0/2 <= true.diff
  
  return(cvrInd)
}


TestInL = function(Sc, N = c(10, 10), simA, simB, true.diff, sig02,
                   alpha0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  smpdiff = mean(simA) - mean(simB)
  
  eta = cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/length(simA)+1/length(simB))*sig02)) + 
    smpdiff/(1+(1/length(simA)+1/length(simB))*sig02/cmsprior[,"v.cp"])
  sig2eta = 1/(1/cmsprior[,"v.cp"] + 1/((1/length(simA)+1/length(simB))*sig02))
  
  IntLen = (qnorm(1-alpha0/2, mean = eta, sd = sqrt(sig2eta)) - true.diff)*2
  
  return(IntLen)
}


TestEV = function(Sc, N = c(10, 10), simA, simB, true.diff, sig02,
                   epsil0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  smpdiff = mean(simA) - mean(simB)
  
  eta = cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/length(simA)+1/length(simB))*sig02)) + 
    smpdiff/(1+(1/length(simA)+1/length(simB))*sig02/cmsprior[,"v.cp"])
  sig2eta = 1/(1/cmsprior[,"v.cp"] + 1/((1/length(simA)+1/length(simB))*sig02))
  
  prdMu = rnorm(1, mean = true.diff, sd = sqrt(sig2eta))
  
  return(prdMu)
}

#------------------------------- toy example --------------------------------#
MySc1 = cbind(c(-0.26, -0.24, -0.37, -0.34, -0.32), 
              c( 0.25,  0.23,  0.22,  0.36,  0.26))

MyWgt1 = c(0.103, 0.175, 0.081, 0.143, 0.077)

cpPriorI = cprior(p = pq(wk = MyWgt1, s0 = 0.05), mk = MySc1[,1], sk2 = MySc1[,2], 
                  wk = MyWgt1, dw = c(3, 3), br = c(18, 3))

MySSI = ACCknvar(R = 0.5, sig02 = 0.35, alpha = 0.05, l0 = 0.65, 
                 wk = MyWgt1, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), s0 = 0.05)

MySSEV = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                   wk = MyWgt1, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), s0 = 0.05)
#----------------------------------------------------------------------------#
nsim = 100000

CvrPrI1 = rep(0, nsim)
CvrPrI2 = rep(0, nsim)

IntLen1 = rep(0, nsim)
IntLen2 = rep(0, nsim)

PredMu1 = rep(0, nsim)
PredMu2 = rep(0, nsim)

## Normal data
true.sig02 = 0.35

SimMuI = rnorm(nsim, mean = cpPriorI[1], 
               sd = sqrt((2/MySSI + 2/MySSI)*0.35 + cpPriorI[2]))


set.seed(123)
for(i in 1:nsim){
  simAnorm = rnorm(MySSI/2, mean = SimMuI[i], sd = sqrt(true.sig02))
  simBnorm = rnorm(MySSI/2, mean = 0, sd = sqrt(true.sig02))
  
  CvrPrI1[i] = TestCvr(Sc = MySc1, N = c(MySSI/2, MySSI/2), 
                      simA = simAnorm, simB = simBnorm, true.diff = SimMuI[i],
                      sig02 = 0.35, l0 = 0.65, wk = MyWgt1,
                      dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  IntLen1[i] = TestInL(Sc = MySc1, N = c(MySSI/2, MySSI/2), 
                       simA = simAnorm, simB = simBnorm, true.diff = SimMuI[i],
                       sig02 = 0.35, alpha0 = 0.05, wk = MyWgt1,
                       dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  PredMu1[i] = TestEV(Sc = MySc1, N = c(MySSEV/2, MySSEV/2), 
                      simA = simAnorm, simB = simBnorm, true.diff = SimMuI[i],
                      sig02 = 0.35, epsil0 = 0.03, wk = MyWgt1,
                      dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}


## Skewed normal data
true.s = 0.98

SimMuI2 = rsn(nsim, dp = cp2dp(c(mean = cpPriorI[1], 
                                 s.d. = sqrt((2/MySSI + 2/MySSI)*0.35 + cpPriorI[2]), 
                                 gamma1 = true.s), family = "SN"))

set.seed(123)
for(i in 1:nsim){
  simAsn = rsn(MySSI/2, dp = cp2dp(c(mean = SimMuI2[i], s.d. = sqrt(true.sig02), 
                                     gamma1 = true.s), family = "SN"))
  simBsn = rsn(MySSI/2, dp = cp2dp(c(mean = 0, s.d. = sqrt(true.sig02), 
                                     gamma1 = true.s), family = "SN"))
  
  CvrPrI2[i] = TestCvr(Sc = MySc1, N = c(MySSI/2, MySSI/2), 
                       simA = simAsn, simB = simBsn, true.diff = SimMuI2[i],
                       sig02 = 0.35, l0 = 0.65, wk = MyWgt1,
                       dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  IntLen2[i] = TestInL(Sc = MySc1, N = c(MySSI/2, MySSI/2), 
                       simA = simAsn, simB = simBsn, true.diff = SimMuI2[i],
                       sig02 = 0.35, alpha0 = 0.05, wk = MyWgt1,
                       dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  PredMu2[i] = TestEV(Sc = MySc1, N = c(MySSEV/2, MySSEV/2), 
                      simA = simAsn, simB = simBsn, true.diff = SimMuI2[i],
                      sig02 = 0.35, epsil0 = 0.03, wk = MyWgt1,
                      dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}

mean(CvrPrI1); mean(CvrPrI2)
mean(IntLen1); mean(IntLen2)
var(PredMu1-SimMuI); var(PredMu2-SimMuI2)

