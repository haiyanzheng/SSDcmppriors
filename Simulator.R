setwd("~/Dropbox/CRUK_Cam/SSE for studies using historical data/Code/Code_Nov2021")

source("Bayesian SSD using robust commensurate priors.R")


##-------------------------- Known variance, sig02 --------------------------##
## N = c(n_A, n_B)                                                           ##
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
#----------------------------------------------------------------------------#
nsim = 100000

MyCvrProbI = rep(0, nsim)
MyCvrProbII = rep(0, nsim)
MyIntLenI = rep(0, nsim)
MyIntLenII = rep(0, nsim)

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
MyPredMuI = rep(0, nsim)
MyPredMuII = rep(0, nsim)


MySSIII = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                 wk = MyWgt1, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                 s0 = 0.05)
MySSIV = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                  wk = MyWgt2, sk2 = MySc1[,2], dw = c(2, 2), br = c(18, 3), 
                  s0 = 0.05)

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

#----------------------------------------------------------------------------#