source("Bayesian SSD using robust commensurate priors.R")
#-----------------------------------------------------------------------------#
library(R2OpenBUGS)

# initial values
cmsinits <- function(){
  list(tld.theta = c(0, 0, 0, 0, 0)
  )
}

## Throughout, N = c(n_A, n_B)                                               ##
#----- to evaluate the average coverage probability of the posterior HPD -----#
tmixpstCP = function(Sc, N = c(10, 10), true.mu, sig02, l0, wk, dw, br, s0 = 0.05){
  
  p = pq(wk, s0)
  Nobs = nrow(Sc)
  a = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  b = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  
  Prior.lambda = Sc[,1]
  Prior.sk2 = Sc[,2]
  
  wMix = cbind(wk, 1-wk)
  
  a01 = dw[1]
  b01 = dw[2]
  a02 = br[1]
  b02 = br[2]
  
  prec = 1/sig02

  cmsdata <- list("Nobs", "a", "b", "Prior.lambda", "Prior.sk2", "wMix", "p",
                  "a01", "b01", "a02", "b02", "prec")
  
  tMixFit <- bugs(data = cmsdata, inits = cmsinits, model.file = "MyMod.txt", 
                  parameters = c("mu"), 
                  n.chains = 1, n.iter = 13000, n.burnin = 3000, bugs.seed = 1)
  
  eta = tMixFit$mean$mu
  sig2eta = tMixFit$sd$mu^2
  
  cvrInd = eta+l0/2 >= true.mu & eta-l0/2 <= true.mu
  
  return(cvrInd)
}


#------- to evaluate the average interval length of the posterior HPD -------#
tmixpstIL = function(Sc, N = c(10, 10), true.mu, sig02, alpha0, wk, dw, br, s0 = 0.05){
  
  p = pq(wk, s0)
  Nobs = nrow(Sc)
  a = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  b = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  
  Prior.lambda = Sc[,1]
  Prior.sk2 = Sc[,2]
  
  wMix = cbind(wk, 1-wk)
  
  a01 = dw[1]
  b01 = dw[2]
  a02 = br[1]
  b02 = br[2]
  
  prec = 1/sig02
  
  cmsdata <- list("Nobs", "a", "b", "Prior.lambda", "Prior.sk2", "wMix", "p",
                  "a01", "b01", "a02", "b02", "prec")
  
  tMixFit <- bugs(data = cmsdata, inits = cmsinits, model.file = "MyMod.txt", 
                  parameters = c("mu"), 
                  n.chains = 1, n.iter = 13000, n.burnin = 3000, bugs.seed = 1)
  
  eta = tMixFit$mean$mu
  sig2eta = tMixFit$sd$mu^2
  
  IntLen = (qnorm(1-alpha0/2, mean = eta, sd = sqrt(sig2eta)) - true.mu)*2
  
  return(IntLen)
}


#------- to evaluate the average interval length of the posterior HPD -------#
tmixpstEV = function(Sc, N = c(10, 10), true.mu, sig02, epsil0, wk, dw, br, s0 = 0.05){
  
  p = pq(wk, s0)
  Nobs = nrow(Sc)
  a = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  b = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  
  Prior.lambda = Sc[,1]
  Prior.sk2 = Sc[,2]
  
  wMix = cbind(wk, 1-wk)
  
  a01 = dw[1]
  b01 = dw[2]
  a02 = br[1]
  b02 = br[2]
  
  prec = 1/sig02
  
  cmsdata <- list("Nobs", "a", "b", "Prior.lambda", "Prior.sk2", "wMix", "p",
                  "a01", "b01", "a02", "b02", "prec")
  
  tMixFit <- bugs(data = cmsdata, inits = cmsinits, model.file = "MyMod.txt", 
                  parameters = c("mu"), 
                  n.chains = 1, n.iter = 13000, n.burnin = 3000, bugs.seed = 1)
  
  eta = tMixFit$mean$mu
  sig2eta = tMixFit$sd$mu^2
  
  prdMu = rnorm(1, mean = true.mu, sd = sqrt(sig2eta))
  
  return(prdMu)
}


#------------------------------- toy example --------------------------------#
MySc1 = cbind(c(-0.26, -0.24, -0.37, -0.34, -0.32), 
              c( 0.25,  0.23,  0.22,  0.36,  0.26))

MyWgt1 = c(0.103, 0.175, 0.081, 0.143, 0.077)
MyWgt2 = c(0.252, 0.319, 0.140, 0.306, 0.149)


cpPriorI = cprior(p = pq(wk = MyWgt1, s0 = 0.05), mk = MySc1[,1], sk2 = MySc1[,2], 
                  wk = MyWgt1, dw = c(3, 3), br = c(18, 3))
cpPriorII = cprior(p = pq(wk = MyWgt1, s0 = 0.05), mk = MySc1[,1], sk2 = MySc1[,2], 
                   wk = MyWgt2, dw = c(3, 3), br = c(18, 3))

MySSI = ACCknvar(R = 0.5, sig02 = 0.35, alpha = 0.05, l0 = 0.65, 
                 wk = MyWgt1, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                 s0 = 0.05)
MySSII = ACCknvar(R = 0.5, sig02 = 0.35, alpha = 0.05, l0 = 0.65, 
                  wk = MyWgt2, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                  s0 = 0.05)
#----------------------------------------------------------------------------#
nsim = 100000

tMixCvrProbI = rep(0, nsim)
tMixCvrProbII = rep(0, nsim)
tMixIntLenI = rep(0, nsim)
tMixIntLenII = rep(0, nsim)

set.seed(123)
SimMuI = rnorm(nsim, mean = cpPriorI[1], 
               sd = sqrt((2/MySSI + 2/MySSI)*0.35 + cpPriorI[2]))

set.seed(123)
SimMuII = rnorm(nsim, mean = cpPriorII[1], 
                sd = sqrt((2/MySSII + 2/MySSII)*0.35 + cpPriorII[2]))


set.seed(123)
for(i in 1:nsim){
  if(i%%5000==0) print(i)
  
  tMixCvrProbI[i] = tmixpstCP(Sc = MySc1, N = c(MySSI/2, MySSI/2), true.mu = SimMuI[i],
                              sig02 = 0.35, l0 = 0.65, wk = MyWgt1,
                              dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  tMixCvrProbII[i] = tmixpstCP(Sc = MySc1, N = c(MySSII/2, MySSII/2), true.mu = SimMuII[i],
                               sig02 = 0.35, l0 = 0.65, wk = MyWgt2,
                               dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  tMixIntLenI[i] = tmixpstIL(Sc = MySc1, N = c(MySSI/2, MySSI/2), true.mu = SimMuI[i],
                             sig02 = 0.35, alpha0 = 0.05, wk = MyWgt1,
                             dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  tMixIntLenII[i] = tmixpstIL(Sc = MySc1, N = c(MySSII/2, MySSII/2), true.mu = SimMuII[i],
                              sig02 = 0.35, alpha0 = 0.05, wk = MyWgt2,
                              dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}

mean(tMixCvrProbI); mean(tMixCvrProbII)
mean(tMixIntLenI); mean(tMixIntLenII)

#----------------------------------------------------------------------------#

MySSIII = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                    wk = MyWgt1, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                    s0 = 0.05)
MySSIV = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                   wk = MyWgt2, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                   s0 = 0.05)


tMixPredMuI = rep(0, nsim)
tMixPredMuII = rep(0, nsim)


set.seed(123)
SimMuIII = rnorm(nsim, mean = cpPriorI[1], 
               sd = sqrt((2/MySSIII + 2/MySSIII)*0.35 + cpPriorI[2]))

set.seed(123)
SimMuIV = rnorm(nsim, mean = cpPriorII[1], 
                sd = sqrt((2/MySSIV + 2/MySSIV)*0.35 + cpPriorII[2]))



set.seed(123)
for(i in 1:nsim){
  tMixPredMuI[i] = tmixpstEV(Sc = MySc1, N = c(MySSIII/2, MySSIII/2), true.mu = SimMuIII[i],
                             sig02 = 0.35, epsil0 = 0.03, wk = MyWgt1,
                             dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  tMixPredMuII[i] = tmixpstEV(Sc = MySc1, N = c(MySSIV/2, MySSIV/2), true.mu = SimMuIV[i],
                              sig02 = 0.35, epsil0 = 0.03, wk = MyWgt2,
                              dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}
  
var(tMixPredMuI-SimMuIII)
var(tMixPredMuII-SimMuIV)
