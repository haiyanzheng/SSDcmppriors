# Using proper Bayesian approach with an informative prior
#------------------------- Known \sigma_0^2 -------------------------#
## ACC/ALC
pbACCknvar = function(R, sig02, alpha = 0.05, l0 = 0.75, prvar){
  N = sig02/(R*(1-R))*(4*qnorm(1-alpha/2)^2/l0^2 - 1/prvar)
  return(N)
}


# APVC
pbAPVCknvar = function(R, sig02, epsil0, prvar){
  N = sig02/(R*(1-R))*(1/epsil0 - 1/prvar)
}


#------------------------- Unknown \sigma_0^2 -------------------------#
finvgm = function(x, a, b){
  b^a/gamma(a)*x^(-a-1)*exp(-b/x)
}

## ACC 
ACCbase2 = function(S2, R, alpha = 0.05, l0, prvar){
  N = S2/(R*(1-R))*(4*qnorm(1-alpha/2)^2/l0^2 - 1/prvar)
  return(N)
}

pbACCunknvar = function(R = 0.5, df, cl, prvar){
  myfun = function(S2){
    ACCbase2(S2, R, alpha = 0.05, l0 = cl, prvar)*finvgm(S2, a = df/2, b = df/2*prvar)
  }
  return(integrate(myfun, 0, Inf)$value)
}

pbACCunknvar(R = 0.5, df = 3, cl = 0.55, prvar = 0.10)


# ALC
rl2 = function(df, N, R, alpha0 = 0.05, prvar){
  avgfun = function(S2){
    1/sqrt(1/prvar + 1/S2*R*(1-R)*N)*finvgm(S2, a = df/2, b = df/2*prvar)
  }
  
  myint = integrate(avgfun, 0, Inf)
  lx = 2*qnorm(1-alpha0/2)*myint$value
  return(lx)
}

rl2(df = 3, N = 80, R = 0.5, alpha0 = 0.05, prvar = 0.2)


pbALCunkvar = function(MySS = 1:100, R = 0.5, df = 3, prvar, alpha0, tlen){
  a = sapply(MySS, function(N){
    rl2(df, N, R, alpha0, prvar)
  }
  )
  return(MySS[which.max(a < tlen)])
}

pbALCunkvar(MySS = 1:100, df = 3, prvar = 0.10, alpha0 = 0.05, tlen = 0.55)


# APVC
APVCbase2 = function(S2, R, epsil0, prvar){
  N = S2/(R*(1-R))*(1/epsil0 - 1/prvar)
  return(N)
}

pbAPVCunknvar = function(R = 0.5, df, prvar, pr0 = 0.07){
  myfun2 = function(S2){
    APVCbase2(S2, R, epsil0 = pr0, prvar)*finvgm(S2, a = df/2, b = df/2*prvar)
  }
  return(integrate(myfun2, 0, Inf)$value)
}

pbAPVCunknvar(R = 0.5, df = 3, prvar = 0.2, pr0 = 0.03)
