setwd("~/Dropbox/CRUK_Cam/SSE for studies using historical data/Code/Code_Nov2021")

source("Bayesian SSD using robust commensurate priors.R")

#-----------------------------------------------------------------------------#
##-------------------------- Known variance, sig02 --------------------------##
#-----------------------------------------------------------------------------#

## Throughout, N = c(n_A, n_B)                                               ##
#----- to evaluate the average coverage probability of the posterior HPD -----#
cmspstCP = function(Sc, N = c(10, 10), true.mu, sig02, l0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  respoA = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  respoB = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  smpdiff = mean(respoA) - mean(respoB)
  
  eta = cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/N[1]+1/N[2])*sig02)) + 
            smpdiff/(1+(1/N[1]+1/N[2])*sig02/cmsprior[,"v.cp"])
  sig2eta = 1/(1/cmsprior[,"v.cp"] + 1/((1/N[1]+1/N[2])*sig02))

  cvrInd = eta+l0/2 >= true.mu & eta-l0/2 <= true.mu
  
  return(cvrInd)
}


# cmspstCP(Sc = MySc1, N = c(40, 40), true.mu = -0.3,
#          sig02 = 0.35, l0 = 0.65, wk = MyWgt1,
#          dw = c(2, 2), br = c(18, 3), s0 = 0.05)

#------- to evaluate the average interval length of the posterior HPD -------#
cmspstIL = function(Sc, N = c(10, 10), true.mu, sig02, alpha0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  respoA = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  respoB = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  smpdiff = mean(respoA) - mean(respoB)
  
  eta = cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/N[1]+1/N[2])*sig02)) + 
           smpdiff/(1+(1/N[1]+1/N[2])*sig02/cmsprior[,"v.cp"])
  sig2eta = 1/(1/cmsprior[,"v.cp"] + 1/((1/N[1]+1/N[2])*sig02))
  
  IntLen = (qnorm(1-alpha0/2, mean = eta, sd = sqrt(sig2eta)) - true.mu)*2
  
  return(IntLen)
}


# cmspstIL(Sc = MySc1, N = c(40, 40), true.mu = -0.3,
#          sig02 = 0.35, alpha0 = 0.05, wk = MyWgt1,
#          dw = c(2, 2), br = c(18, 3), s0 = 0.05)

#----------------- to evaluate the average posterior variance -----------------#
cmspstEV = function(Sc, N = c(10, 10), true.mu, sig02, epsil0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  respoA = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  respoB = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  smpdiff = mean(respoA) - mean(respoB)
  
  eta = cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/N[1]+1/N[2])*sig02)) + 
    smpdiff/(1+(1/N[1]+1/N[2])*sig02/cmsprior[,"v.cp"])
  sig2eta = 1/(1/cmsprior[,"v.cp"] + 1/((1/N[1]+1/N[2])*sig02))
  
  prdMu = rnorm(1, mean = true.mu, sd = sqrt(sig2eta))
  
  return(prdMu)
}


# cmspstEV(Sc = MySc1, N = c(40, 40), true.mu = -0.3,
#          sig02 = 0.35, epsil0 = 0.03, wk = MyWgt1,
#          dw = c(2, 2), br = c(18, 3), s0 = 0.05)

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

MyCvrProbI = rep(0, nsim)
MyCvrProbII = rep(0, nsim)
MyIntLenI = rep(0, nsim)
MyIntLenII = rep(0, nsim)


SimMuI = rnorm(nsim, mean = cpPriorI[1], 
               sd = sqrt((2/MySSI + 2/MySSI)*0.35 + cpPriorI[2]))
SimMuII = rnorm(nsim, mean = cpPriorII[1], 
                sd = sqrt((2/MySSII + 2/MySSII)*0.35 + cpPriorII[2]))


set.seed(123)
for(i in 1:nsim){
  MyCvrProbI[i] = cmspstCP(Sc = MySc1, N = c(MySSI/2, MySSI/2), true.mu = SimMuI[i],
                           sig02 = 0.35, l0 = 0.65, wk = MyWgt1,
                           dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyCvrProbII[i] = cmspstCP(Sc = MySc1, N = c(MySSII/2, MySSII/2), true.mu = SimMuII[i],
                           sig02 = 0.35, l0 = 0.65, wk = MyWgt2,
                           dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyIntLenI[i] = cmspstIL(Sc = MySc1, N = c(MySSI/2, MySSI/2), true.mu = SimMuI[i],
                          sig02 = 0.35, alpha0 = 0.05, wk = MyWgt1,
                          dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyIntLenII[i] = cmspstIL(Sc = MySc1, N = c(MySSII/2, MySSII/2), true.mu = SimMuII[i],
                          sig02 = 0.35, alpha0 = 0.05, wk = MyWgt2,
                          dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}

mean(MyCvrProbI); mean(MyCvrProbII)
mean(MyIntLenI); mean(MyIntLenII)

#----------------------------------------------------------------------------#

MySSIII = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                 wk = MyWgt1, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                 s0 = 0.05)
MySSIV = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                  wk = MyWgt2, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                  s0 = 0.05)


MyPredMuI = rep(0, nsim)
MyPredMuII = rep(0, nsim)


set.seed(123)
for(i in 1:nsim){
  MyPredMuI[i] = cmspstEV(Sc = MySc1, N = c(MySSIII/2, MySSIII/2), true.mu = SimMuI[i],
                          sig02 = 0.35, epsil0 = 0.03, wk = MyWgt1,
                          dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyPredMuII[i] = cmspstEV(Sc = MySc1, N = c(MySSIV/2, MySSIV/2), true.mu = SimMuII[i],
                           sig02 = 0.35, epsil0 = 0.03, wk = MyWgt2,
                           dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}

var(MyPredMuI-SimMuI)
var(MyPredMuII-SimMuII)

#-----------------------------------------------------------------------------#
##------------------------- Unknown variance, sig02 -------------------------##
#-----------------------------------------------------------------------------#

#----- to evaluate the average coverage probability of the posterior HPD -----#
cmspstCP2 = function(Sc, N = c(10, 10), true.mu, true.sig02, 
                     df, l0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  respoA = rnorm(N[1], mean = true.mu, sd = sqrt(true.sig02))
  respoB = rnorm(N[2], mean = 0, sd = sqrt(true.sig02))
  smpdiff = mean(respoA) - mean(respoB)
  
  eta = function(S2){
    (cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/N[1]+1/N[2])*S2)) + 
        smpdiff/(1+(1/N[1]+1/N[2])*S2/cmsprior[,"v.cp"]))*finvgm(S2, a=df/2, b=df/2*cmsprior[,"v.cp"])
  }
  
  sig2eta = function(S2){
    1/(1/cmsprior[,"v.cp"] + 1/((1/N[1]+1/N[2])*S2))*finvgm(S2, a=df/2, b=df/2*cmsprior[,"v.cp"])
  }
  
  eta.pred = integrate(eta, 0, Inf)$value
  
  cvrInd = eta.pred+l0/2 >= true.mu & eta.pred-l0/2 <= true.mu
  
  return(cvrInd)
}


# cmspstCP2(Sc = MySc1, N = c(40, 40), true.mu = -0.3,
#          true.sig02 = 0.35, df = 5, l0 = 0.65, wk = MyWgt1,
#          dw = c(2, 2), br = c(18, 3), s0 = 0.05)

#------- to evaluate the average interval length of the posterior HPD -------#
cmspstIL2 = function(Sc, N = c(10, 10), true.mu, true.sig02, 
                     df, alpha0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  respoA = rnorm(N[1], mean = true.mu, sd = sqrt(sig02))
  respoB = rnorm(N[2], mean = 0, sd = sqrt(sig02))
  smpdiff = mean(respoA) - mean(respoB)
  
  eta = function(S2){
    (cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/N[1]+1/N[2])*S2)) + 
       smpdiff/(1+(1/N[1]+1/N[2])*S2/cmsprior[,"v.cp"]))*finvgm(S2, a=df/2, b=df/2*cmsprior[,"v.cp"])
  }
  
  sig2eta = function(S2){
    1/(1/cmsprior[,"v.cp"] + 1/((1/N[1]+1/N[2])*S2))*finvgm(S2, a=df/2, b=df/2*cmsprior[,"v.cp"])
  }
  
  eta.pred = integrate(eta, 0, Inf)$value
  sig2eta.pred = integrate(sig2eta, 0, Inf)$value
  
  IntLen = (qnorm(1-alpha0/2, mean = eta.pred, sd = sqrt(sig2eta.pred)) - true.mu)*2
  
  return(IntLen)
}


# cmspstIL2(Sc = MySc1, N = c(40, 40), true.mu = -0.3,
#          true.sig02 = 0.35, df = 5, alpha0 = 0.05, wk = MyWgt1,
#          dw = c(2, 2), br = c(18, 3), s0 = 0.05)

#----------------- to evaluate the average posterior variance -----------------#
cmspstEV2 = function(Sc, N = c(10, 10), true.mu, true.sig02, 
                     df, epsil0, wk, dw, br, s0 = 0.05){
  
  pk = pq(wk, s0)
  cmsprior = cprior(p = pk, mk = Sc[,1], sk2 = Sc[,2], wk, dw, br)
  
  respoA = rnorm(N[1], mean = true.mu, sd = sqrt(true.sig02))
  respoB = rnorm(N[2], mean = 0, sd = sqrt(true.sig02))
  smpdiff = mean(respoA) - mean(respoB)
  
  eta = function(S2){
    (cmsprior[,"m.cp"]/(1+cmsprior[,"v.cp"]/((1/N[1]+1/N[2])*S2)) + 
       smpdiff/(1+(1/N[1]+1/N[2])*S2/cmsprior[,"v.cp"]))*finvgm(S2, a=df/2, b=df/2*cmsprior[,"v.cp"])
  }
  
  sig2eta = function(S2){
    1/(1/cmsprior[,"v.cp"] + 1/((1/N[1]+1/N[2])*S2))*finvgm(S2, a=df/2, b=df/2*cmsprior[,"v.cp"])
  }
  
  eta.pred = integrate(eta, 0, Inf)$value
  sig2eta.pred = integrate(sig2eta, 0, Inf)$value
  
  prdMu = rnorm(1, mean = true.mu, sd = sqrt(sig2eta.pred))
  
  return(prdMu)
}


# cmspstEV2(Sc = MySc1, N = c(40, 40), true.mu = -0.3,
#          true.sig02 = 0.35, df = 5, epsil0 = 0.03, wk = MyWgt1,
#          dw = c(2, 2), br = c(18, 3), s0 = 0.05)

#----------------------------------------------------------------------------#
ACCSSI = ACCunknvar(R = 0.5, sk2 = MySc1[,2], wk = MyWgt1, 
                    df = 3, l0 = 0.65, alpha = 0.05, 
                    dw = c(2, 2), br = c(18, 3), s0 = 0.05) 

ACCSSII = ACCunknvar(R = 0.5, sk2 = MySc1[,2], wk = MyWgt2, 
                     df = 3, l0 = 0.65, alpha = 0.05, 
                     dw = c(2, 2), br = c(18, 3), s0 = 0.05) 

ALCSSI = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = MySc1[,2], wk = MyWgt1, 
                   df = 3, dw = c(2, 2), br = c(18, 3), 
                   alpha0 = 0.05, tlen = 0.65, s0 = 0.05)

ALCSSII = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = MySc1[,2], wk = MyWgt2, 
                    df = 3, dw = c(2, 2), br = c(18, 3), 
                    alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
#----------------------------------------------------------------------------#
MyCvrProbI2 = rep(0, nsim)
MyCvrProbII2 = rep(0, nsim)
MyIntLenI2 = rep(0, nsim)
MyIntLenII2 = rep(0, nsim)


set.seed(123)
for(i in 1:nsim){
  MyCvrProbI2[i] = cmspstCP2(Sc = MySc1, N = c(ACCSSI/2, ACCSSI/2), true.mu = SimMuI[i],
                             true.sig02 = 0.35, df = 5, l0 = 0.65, wk = MyWgt1,
                             dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyCvrProbII2[i] = cmspstCP2(Sc = MySc1, N = c(ACCSSII/2, ACCSSII/2), true.mu = SimMuII[i],
                              true.sig02 = 0.35, df = 5, l0 = 0.65, wk = MyWgt2,
                              dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyIntLenI2[i] = cmspstIL2(Sc = MySc1, N = c(ALCSSI/2, ALCSSI/2), true.mu = SimMuI[i],
                            true.sig02 = 0.35, df = 5, alpha0 = 0.05, wk = MyWgt1,
                            dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyIntLenII2[i] = cmspstIL2(Sc = MySc1, N = c(ALCSSII/2, ALCSSII/2), true.mu = SimMuI[i],
                             true.sig02 = 0.35, df = 5, alpha0 = 0.05, wk = MyWgt2,
                             dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
}

mean(MyCvrProbI2); mean(MyCvrProbII2)
mean(MyIntLenI2); mean(MyIntLenII2)

#----------------------------------------------------------------------------#
APVCSSI = APVCunknvar(R = 0.5, sk2 = MySc1[,2], wk = MyWgt1, epsil0 = 0.03, 
                      df = 3, dw = c(2, 2), br = c(18, 3), s0 = 0.05)

APVCSSII = APVCunknvar(R = 0.5, sk2 = MySc1[,2], wk = MyWgt2, epsil0 = 0.03, 
                       df = 3, dw = c(2, 2), br = c(18, 3), s0 = 0.05)



MyPredMuI2 = rep(0, nsim)
MyPredMuII2 = rep(0, nsim)


set.seed(123)
for(i in 1:nsim){
  MyPredMuI2[i] = cmspstEV(Sc = MySc1, N = c(APVCSSI/2, APVCSSI/2), true.mu = SimMuI[i],
                          sig02 = 0.35, epsil0 = 0.03, wk = MyWgt1,
                          dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  
  MyPredMuII2[i] = cmspstEV(Sc = MySc1, N = c(APVCSSII/2, APVCSSII/2), true.mu = SimMuII[i],
                           sig02 = 0.35, epsil0 = 0.03, wk = MyWgt2,
                           dw = c(2, 2), br = c(18, 3), s0 = 0.05)
}

var(MyPredMuI2-SimMuI)
var(MyPredMuII2-SimMuII)
