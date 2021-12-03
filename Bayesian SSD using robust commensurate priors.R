#---------------- Pairwise distributional discrepancy ----------------#
# Squared Hellinger distance between any two normal distributions
NormHd2 = function(Norm1, Norm2){
  1 - sqrt(2*sqrt(Norm1[2]*Norm2[2])/(Norm1[2]+Norm2[2]))*exp(-1/4*(Norm1[1]-Norm2[1])^2/(Norm1[2]+Norm2[2]))
}


HdMat = function(mu.hist, var.hist){
  histdat = cbind(1:length(mu.hist), mu.hist, var.hist)
  colnames(histdat)[1] = "id"
  
  ## pairwise indicators
  prs = t(combn(1:length(mu.hist), 2))
  colnames(prs) = c("id1", "id2")
  
  myNorm1 = subset(histdat[prs[,1],], select = c("mu.hist", "var.hist"))
  myNorm2 = subset(histdat[prs[,2],], select = c("mu.hist", "var.hist"))
  
  a = numeric(10)
  for(i in 1:10) a[i] = NormHd2(Norm1 = myNorm1[i,], Norm2 = myNorm2[i,])
  
  myHdMat = cbind(prs, round(sqrt(a), 3))
  
  return(myHdMat)
}

#----------- Collective prior -----------#
pq = function(wk, s0){
  exp(-wk^2/s0)/sum(exp(-wk^2/s0))
}


cprior = function(p, mk, sk2, wk, dw, br){
  m.cp = sum(p*mk)
  v.cp = sum(p^2*(sk2 + wk*dw[2]/(dw[1] - 1) + (1 - wk)*br[2]/(br[1] - 1)))
  return(cbind(m.cp, v.cp))
}


#=====================================================================#
#     Implementation of the proposed Bayesian sample size formulae    #
#=====================================================================#
# R is the allocation ratio to A: n_A/N, where N = n_A + n_B 
# where N = n_A + n_B as the new study sample size
# gm1 = (a01, b01)
# gm2 = (a02, b02)
#------------------------- Known \sigma_0^2 -------------------------#
## ACC and ALC
ACCknvar = function(R, sig02, alpha = 0.05, l0, wk, sk2, dw, br, s0 = 0.15){
  p = pq(wk, s0)
  N = sig02/(R*(1-R))*(4*qnorm(1-alpha/2)^2/l0^2 - 1/sum(p^2*(sk2 + wk*dw[2]/(dw[1]-1)+(1-wk)*br[2]/(br[1]-1))))
  return(N)
}


## APVC
APVCknvar = function(R, sig02, epsil0, wk, sk2, dw, br, s0 = 0.15){
  p = pq(wk, s0)
  N = sig02/(R*(1-R))*(1/epsil0 - 1/sum(p^2*(sk2 + wk*dw[2]/(dw[1]-1)+(1-wk)*br[2]/(br[1]-1))))
  return(N)
}


#------------------------- Unknown \sigma_0^2 -------------------------#
finvgm = function(x, a, b){
  b^a/gamma(a)*x^(-a-1)*exp(-b/x)
}

## ACC
ACCunknvar = function(R = 0.5, sk2, wk, df, alpha = 0.05, l0 = 0.65,
                      dw = c(2, 2), br = c(18, 3), s0 = 0.05){
  p = pq(wk, s0)
  hisvar = sum(p^2*(sk2 + wk*dw[2]/(dw[1]-1)+(1-wk)*br[2]/(br[1]-1)))
  myfun = function(S2){
    ACCknvar(R, sig02 = S2, alpha, l0, 
             wk, sk2, dw, br, s0)*finvgm(S2, a = df/2, b = df/2*hisvar)
  }
  return(integrate(myfun, 0, Inf)$value)
}

## ALC
rl = function(df, N, R, alpha0 = 0.05, wk, sk2, dw, br, s0 = 0.15){
  p = pq(wk, s0)
  hisvar = sum(p^2*(sk2 + wk*dw[2]/(dw[1]-1)+(1-wk)*br[2]/(br[1]-1)))
  avgfun = function(S2){
    1/sqrt(1/hisvar + 1/S2*R*(1-R)*N)*finvgm(S2, a = df/2, b = df/2*hisvar)
  }
  myint = integrate(avgfun, 0, Inf)
  lx = 2*qnorm(1-alpha0/2)*myint$value
  return(lx)
}


## APVC
APVCunknvar = function(R = 0.5, sk2, wk, df, epsil0 = 0.05, 
                       dw = c(2, 2), br = c(18, 3), s0 = 0.05){
  p = pq(wk, s0)
  hisvar = sum(p^2*(sk2 + wk*dw[2]/(dw[1]-1)+(1-wk)*br[2]/(br[1]-1)))
  myfun = function(S2){
    APVCknvar(R, sig02 = S2, epsil0, 
              wk, sk2, dw, br, s0)*finvgm(S2, a = df/2, b = df/2*hisvar)
  }
  return(integrate(myfun, 0, Inf)$value)
}

ALCunkvar = function(MySS = 1:100, R = 0.5, wk = mywk, sk2 = MySc1[,2], df = 3, 
                     dw = c(2, 2), br = c(18, 3), alpha0 = 0.05, tlen = 0.65, s0 = 0.05){
  a = sapply(MySS, function(N){
    rl(df, N, R, alpha0, wk, sk2, dw, br, s0)
  }
  )
  return(MySS[which.max(a < tlen)])
}
