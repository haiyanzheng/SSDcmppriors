setwd("~/Dropbox/CRUK_Cam/SSE for studies using historical data/Code/Code_Nov2021")

source("Historical data.R")
source("Bayesian SSD using arbitrary priors.R")
# source("Bayesian SSD using robust commensurate priors.R")  ## already loaded in "Historical data.R"

#---------------- Figure 1 ----------------#
mats <- grep(x= ls(pos=1), pattern="MySc", value=TRUE)
HisDatV2 <- do.call(cbind, mget(mixedsort(mats)))[, c(2, 4, 6, 8)]

mats <- grep(x= ls(pos=1), pattern="MyWgt", value=TRUE)
AllWgt <- do.call(cbind, mget(mixedsort(mats)))

AllCb  = array(0, dim = c(5, 8))
AllHis = cbind(array(rep(MySc1, 2), dim = c(10, 2)), array(rep(MySc2, 2), dim = c(10, 2)),
               array(rep(MySc3, 2), dim = c(10, 2)), array(rep(MySc4, 2), dim = c(10, 2)))
HisDatM = AllHis[1:5,]
HisDatV = AllHis[6:10,]

ClPrior  = array(0, dim = c(2, 8))

for(i in 1:8){
  AllCb[,i] = pq(wk = AllWgt[,i], s0 = 0.05)
  ClPrior[,i] = round(cprior(p = AllCb[,i], mk = HisDatM[,i], sk2 = HisDatV[,i], 
                             wk = get(paste0("MyWgt", i)), 
                             dw = c(2, 2), br = c(18, 3)), 3)
}

round(AllCb, 3)
# All collective priors in Table 1
ClPrior

##----------- Known \sigma_0^2 ----------##
##--------------- ACC/ALC ---------------##
##---------------------------------------##

modWgt = array(0, dim = c(5, 4, 4))
for(j in 1:4)  modWgt[,,j] = cbind(AllWgt[, (2*j-1):(2*j)], 0, 1)

ACCf1  = array(0, dim = c(4, 7))

###############################################################################
# Configuration i, approach k
# Approach 1: Optim 1 as benchmark
# Approach 2: Optim 2 as benchmark
# Approach 3: proposed Bayesian SSD with robust weights 1 (Sc1)
# Approach 4: proposed Bayesian SSD with robust weights 2 (Sc1)
# Approach 5: proposed Bayesian SSD without robustification
# Approach 6: No borrowing at all
# Approach 7: Bayesian SSD using informative prior & arbitary guess of \sigma_0^2
###############################################################################
for(i in 1:4){
  # Optim 1
  ACCf1[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], alpha = 0.05, l0 = 0.65, 
                         wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                         s0 = 0.05)
  
  # Optim 2
  ACCf1[i, 2] = ACCknvar(R = 0.5, sig02 = ClPrior[2,2*i], alpha = 0.05, l0 = 0.65, 
                         wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                         s0 = 0.05)
}

for(i in 1:4){
  for(k in 1:4){
    ACCf1[i, (k+2)] = ACCknvar(R = 0.5, sig02 = 0.35, alpha = 0.05, l0 = 0.65, 
                               wk = modWgt[,k,i], sk2 = HisDatV2[,i], 
                               dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  }
}


for(i in 1:4){
  ACCf1[i, 7] = pbACCknvar(R = 0.5, sig02 = 0.35, alpha = 0.05, l0 = 0.65, prvar = min(HisDatV2[,i]))
}

# round(ACCf1, 1)

##--------- Known \sigma_0^2 ---------##
##-------------- APVC ----------------##
##------------------------------------##

APVCf1  = array(0, dim = c(4, 7))

for(i in 1:4){
  APVCf1[i, 1] = APVCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], epsil0 = 0.03, 
                           wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                           s0 = 0.05)
  
  APVCf1[i, 2] = APVCknvar(R = 0.5, sig02 = ClPrior[2,2*i], epsil0 = 0.03, 
                           wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                           s0 = 0.05)
}

for(i in 1:4){
  for(k in 1:4){
    APVCf1[i, (k+2)] = APVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, 
                                 wk = modWgt[,k,i], sk2 = HisDatV2[,i], 
                                 dw = c(2, 2), br = c(18, 3), s0 = 0.05)
  }
}


for(i in 1:4){
  APVCf1[i, 7] = pbAPVCknvar(R = 0.5, sig02 = 0.35, epsil0 = 0.03, prvar = min(HisDatV2[,i]))
}

# round(APVCf1, 1)

##-------- Unknown \sigma_0^2 --------##
##--------------- ACC ----------------##
##------------------------------------##

ACCf2  = array(0, dim = c(4, 7))

for(i in 1:4){
  # Optim 1
  ACCf2[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], alpha = 0.05, l0 = 0.65, 
                         wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                         s0 = 0.05)
  
  # Optim 2
  ACCf2[i, 2] = ACCknvar(R = 0.5, sig02 = ClPrior[2,2*i], alpha = 0.05, l0 = 0.65, 
                         wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                         s0 = 0.05)
}


for(i in 1:4){
  for(k in 1:4){
    ACCf2[i, (k+2)] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], df = 3, alpha = 0.05,
                                 dw = c(2, 2), br = c(18, 3), l0 = 0.65, s0 = 0.05)
  }
}


for(i in 1:4){
  ACCf2[i, 7] = pbACCunknvar(R = 0.5, df = 3, cl = 0.65, prvar = min(HisDatV2[,i]))
}

# round(ACCf2, 1)

##-------- Unknown \sigma_0^2 --------##
##-------------- ALC -----------------##
##------------------------------------##

ALCf2  = array(0, dim = c(4, 7))

for(i in 1:4){
  # Optim 1
  ALCf2[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], alpha = 0.05, l0 = 0.65, 
                         wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                         s0 = 0.05)
  
  # Optim 2
  ALCf2[i, 2] = ACCknvar(R = 0.5, sig02 = ClPrior[2,2*i], alpha = 0.05, l0 = 0.65, 
                         wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                         s0 = 0.05)
}

for(i in 1:4){
  for(k in 1:4){
    ALCf2[i, (k+2)] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], 
                                df = 3, dw = c(2, 2), br = c(18, 3), 
                                alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
  }
}


for(i in 1:4){
  ALCf2[i, 7] = pbALCunkvar(MySS = 1:300, R = 0.5, df = 3, prvar = min(HisDatV2[,i]), 
                            alpha0 = 0.05, tlen = 0.65)
}

# round(ALCf2, 1)

##-------- Unknown \sigma_0^2 --------##
##-------------- APVC ----------------##
##------------------------------------##

APVCf2  = array(0, dim = c(4, 7))

for(i in 1:4){
  # Optim 1
  APVCf2[i, 1] = APVCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], epsil0 = 0.03, 
                           wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                           s0 = 0.05)
  # Optim 2  
  APVCf2[i, 2] = APVCknvar(R = 0.5, sig02 = ClPrior[2,2*i], epsil0 = 0.03, 
                           wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                           s0 = 0.05)
}

for(i in 1:4){
  for(k in 1:4){
    APVCf2[i, (k+2)] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], df = 3,
                                   dw = c(2, 2), br = c(18, 3), epsil0 = 0.03, s0 = 0.05)
  }
}


for(i in 1:4){
  APVCf2[i, 7] = pbAPVCunknvar(R = 0.5, df = 3, prvar = min(HisDatV2[,i]), pr0 = 0.03)
}

# round(APVCf2, 1)

#-------------------- start plotting --------------------#
colnames(ACCf1) <- c("OptimA", "OptimB", "RbWgtA", "RbWgtB", "NoRb", "NoBrw", "InPr")
MyACCf1 = data.frame(Case = c("Configuration 1", "Configuration 2", 
                              "Configuration 3", "Configuration 4"), ACCf1)

colnames(APVCf1) <- c("OptimA", "OptimB", "RbWgtA", "RbWgtB", "NoRb", "NoBrw", "InPr")
MyAPVCf1 = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), APVCf1)

MyACCf1trs = melt(MyACCf1, id = "Case")
MyAPVCf1trs = melt(MyAPVCf1, id = "Case")

MyFig1Dat = rbind(data.frame(MyACCf1trs, Criter = "ACC"), 
                  data.frame(MyAPVCf1trs, Criter = "APVC"))
colnames(MyFig1Dat)[2] = "Approach"
head(MyFig1Dat)


MyFig1Dat2 = subset(MyFig1Dat, Approach != "OptimA" & Approach != "OptimB")

m1 <- ggplot(data = MyFig1Dat2, aes(x = Criter, y = value, fill = Approach)) + theme_bw() +
  geom_bar(stat = "identity", color = "black", width = 0.42, position = position_dodge()) + 
  scale_x_discrete("", labels = c("ACC/ALC", "APVC")) + ylab("Sample size") +
  theme(legend.position='none', legend.title = element_blank()) 

F1knvar <- m1 + scale_fill_manual(name = element_blank(), 
                                  labels = c("Robust weights I", "Robust weights II",
                                             "No robustification", "No borrowing", "Single source information"),
                                  values = c("#ca0020", "#f4a582", 
                                             "#0072B2", "lightgray", "#b2abd2")) + facet_wrap(~Case, ncol = 4) 

# Unknown variance
colnames(ACCf2) <- c("OptimA", "OptimB", "RbWgtA", "RbWgtB", "NoRb", "NoBrw", "InPr")
MyACCf2 = data.frame(Case = c("Configuration 1", "Configuration 2", 
                              "Configuration 3", "Configuration 4"), ACCf2)

colnames(ALCf2) <- c("OptimA", "OptimB", "RbWgtA", "RbWgtB", "NoRb", "NoBrw", "InPr")
MyALCf2 = data.frame(Case = c("Configuration 1", "Configuration 2", 
                              "Configuration 3", "Configuration 4"), ALCf2)

colnames(APVCf2) <- c("OptimA", "OptimB", "RbWgtA", "RbWgtB", "NoRb", "NoBrw", "InPr")
MyAPVCf2 = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), APVCf2)


MyACCf2trs = melt(MyACCf2, id = "Case")
MyALCf2trs = melt(MyALCf2, id = "Case")
MyAPVCf2trs = melt(MyAPVCf2, id = "Case")

MyFig2Dat = rbind(data.frame(MyACCf2trs, Criter = "ACC"), 
                  data.frame(MyALCf2trs, Criter = "ALC"),
                  data.frame(MyAPVCf2trs, Criter = "APVC"))
colnames(MyFig2Dat)[2] = "Approach"    
head(MyFig2Dat)


MyFig2Dat2 = subset(MyFig2Dat, Approach != "OptimA" & Approach != "OptimB")

m2 <- ggplot(data = MyFig2Dat2, aes(x = Criter, y = value, fill = Approach)) + theme_bw() +
  geom_bar(stat = "identity", color = "black", width = 0.6, position = position_dodge()) + 
  scale_x_discrete("", labels = c("ACC", "ALC", "APVC")) + ylab("Sample size") +
  theme(legend.position='none', legend.title = element_blank()) 


F1unknvar <- m2 + scale_fill_manual(name = element_blank(), 
                                    labels = c("Robust weights A", "Robust weights B",
                                               "No robustification", "No borrowing", "Single source information"),
                                    values = c("#ca0020", "#f4a582", 
                                               "#0072B2", "lightgray", "#b2abd2")) + facet_wrap(~Case, ncol = 4) 

prow <- cowplot::plot_grid( F1knvar,
                            F1unknvar,
                            align = 'v',
                            labels = c("(i)", "(ii)"),
                            hjust = -0.5,
                            ncol = 1
)

# Figure 1: (eps) 800*525
prow

## with legend:
# legend <- cowplot::get_legend(F1unknvar)
# cowplot::plot_grid(legend, prow, ncol = 1, rel_heights = c(0.15, 1))


#-------------------- Figure 2 --------------------#
## ACC ##

# ACCf3a to restore SS computed using Robust Weights A for each configuration based on ACC
ACCf3a  = array(0, dim = c(4, 7))
# ACCf3b to restore SS computed using Robust Weights B for each configuration based on ACC
ACCf3b = array(0, dim = c(4, 7))
# ACCf3c to restore SS computed with no robustification
ACCf3c  = array(0, dim = c(4, 7))
# ACCf3d to restore SS computed with no borrowing
ACCf3d = array(0, dim = c(4, 7))

for(i in 1:4){
  # Optim 1
  ACCf3a[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], alpha = 0.05, l0 = 0.65, 
                          wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                          s0 = 0.05)
  
  # Optim 2
  ACCf3b[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,2*i], alpha = 0.05, l0 = 0.65, 
                          wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                          s0 = 0.05)
}


mydf = c(3, 5, 10, 20, 30, 40)

for(i in 1:4){
  for(p in 1:6){
    ACCf3a[i, (p+1)] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], df = mydf[p], 
                                  dw = c(2, 2), br = c(18, 3), l0 = 0.65, s0 = 0.05)
    
    ACCf3b[i, (p+1)] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], df = mydf[p], 
                                  dw = c(2, 2), br = c(18, 3), l0 = 0.65, s0 = 0.05)
    
    ACCf3c[i, (p+1)] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,3,i], df = mydf[p], 
                                  dw = c(2, 2), br = c(18, 3), l0 = 0.65, s0 = 0.05)
    
    ACCf3d[i, (p+1)] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,4,i], df = mydf[p], 
                                  dw = c(2, 2), br = c(18, 3), l0 = 0.65, s0 = 0.05)
  }
}


## ALC ##

# ALCf3a to restore SS computed using Robust Weights A for each configuration based on ALC
ALCf3a  = array(0, dim = c(4, 7))
# ALCf3b to restore SS computed using Robust Weights B for each configuration based on ALC
ALCf3b = array(0, dim = c(4, 7))
# ALCf3c to restore SS computed with no robustification
ALCf3c  = array(0, dim = c(4, 7))
# ALCf3d to restore SS computed with no borrowing
ALCf3d = array(0, dim = c(4, 7))

for(i in 1:4){
  # Optim 1
  ALCf3a[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], alpha = 0.05, l0 = 0.65, 
                          wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                          s0 = 0.05)
  
  # Optim 2
  ALCf3b[i, 1] = ACCknvar(R = 0.5, sig02 = ClPrior[2,2*i], alpha = 0.05, l0 = 0.65, 
                          wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                          s0 = 0.05)
}


for(i in 1:4){
  for(p in 1:6){
    ALCf3a[i, (p+1)] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], 
                                 df = mydf[p], dw = c(2, 2), br = c(18, 3), 
                                 alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
    
    ALCf3b[i, (p+1)] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], 
                                 df = mydf[p], dw = c(2, 2), br = c(18, 3), 
                                 alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
    
    ALCf3c[i, (p+1)] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,3,i], 
                                 df = mydf[p], dw = c(2, 2), br = c(18, 3), 
                                 alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
    
    ALCf3d[i, (p+1)] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,4,i], 
                                 df = mydf[p], dw = c(2, 2), br = c(18, 3), 
                                 alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
  }
}


## APVC ##

# APVCf3a to restore SS computed using Robust Weights A for each configuration based on APVC
APVCf3a  = array(0, dim = c(4, 7))
# APVCf3b to restore SS computed using Robust Weights B for each configuration based on APVC
APVCf3b = array(0, dim = c(4, 7))
# APVCf3c to restore SS computed with no robustification
APVCf3c  = array(0, dim = c(4, 7))
# APVCf3d to restore SS computed with no borrowing
APVCf3d = array(0, dim = c(4, 7))

for(i in 1:4){
  # Optim 1
  APVCf3a[i, 1] = APVCknvar(R = 0.5, sig02 = ClPrior[2,(2*i-1)], epsil0 = 0.03, 
                            wk = modWgt[,1,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                            s0 = 0.05)
  # Optim 2  
  APVCf3b[i, 1] = APVCknvar(R = 0.5, sig02 = ClPrior[2,2*i], epsil0 = 0.03, 
                            wk = modWgt[,2,i], sk2 = HisDatV2[,i], dw = c(2, 2), br = c(18, 3), 
                            s0 = 0.05)
}

for(i in 1:4){
  for(p in 1:6){
    APVCf3a[i, (p+1)] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], df = mydf[p],
                                    dw = c(2, 2), br = c(18, 3), epsil0 = 0.03, s0 = 0.05)
    
    APVCf3b[i, (p+1)] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], df = mydf[p],
                                    dw = c(2, 2), br = c(18, 3), epsil0 = 0.03, s0 = 0.05)
    
    APVCf3c[i, (p+1)] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,3,i], df = mydf[p],
                                    dw = c(2, 2), br = c(18, 3), epsil0 = 0.03, s0 = 0.05)
    
    APVCf3d[i, (p+1)] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,4,i], df = mydf[p],
                                    dw = c(2, 2), br = c(18, 3), epsil0 = 0.03, s0 = 0.05)
  }
}


#-------------------- start plotting --------------------#
# ACC
colnames(ACCf3a) <- c("OptimA", paste0("c", 1:6))
MyACCf3a = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ACCf3a)

colnames(ACCf3b) <- c("OptimB", paste0("c", 1:6))
MyACCf3b = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ACCf3b)

colnames(ACCf3c) <- c("OptimA", paste0("c", 1:6))
MyACCf3c = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ACCf3c)

colnames(ACCf3d) <- c("OptimB", paste0("c", 1:6))
MyACCf3d = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ACCf3d)

# ALC
colnames(ALCf3a) <- c("OptimA", paste0("c", 1:6))
MyALCf3a = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ALCf3a)

colnames(ALCf3b) <- c("OptimB", paste0("c", 1:6))
MyALCf3b = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ALCf3b)

colnames(ALCf3c) <- c("OptimA", paste0("c", 1:6))
MyALCf3c = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ALCf3c)

colnames(ALCf3d) <- c("OptimB", paste0("c", 1:6))
MyALCf3d = data.frame(Case = c("Configuration 1", "Configuration 2", 
                               "Configuration 3", "Configuration 4"), ALCf3d)

# APVC
colnames(APVCf3a) <- c("OptimA", paste0("c", 1:6))
MyAPVCf3a = data.frame(Case = c("Configuration 1", "Configuration 2", 
                                "Configuration 3", "Configuration 4"), APVCf3a)

colnames(APVCf3b) <- c("OptimB", paste0("c", 1:6))
MyAPVCf3b = data.frame(Case = c("Configuration 1", "Configuration 2", 
                                "Configuration 3", "Configuration 4"), APVCf3b)

colnames(APVCf3c) <- c("OptimA", paste0("c", 1:6))
MyAPVCf3c = data.frame(Case = c("Configuration 1", "Configuration 2", 
                                "Configuration 3", "Configuration 4"), APVCf3c)

colnames(APVCf3d) <- c("OptimB", paste0("c", 1:6))
MyAPVCf3d = data.frame(Case = c("Configuration 1", "Configuration 2", 
                                "Configuration 3", "Configuration 4"), APVCf3d)

MyACCf3atrs = melt(MyACCf3a, id = "Case")
MyACCf3btrs = melt(MyACCf3b, id = "Case")
MyACCf3ctrs = melt(MyACCf3c, id = "Case")
MyACCf3dtrs = melt(MyACCf3d, id = "Case")

Fig3ACC = rbind(data.frame(MyACCf3atrs, type = "RbW1"),
                data.frame(MyACCf3btrs, type = "RbW2"),
                data.frame(MyACCf3ctrs, type = "NoRb"), 
                data.frame(MyACCf3dtrs, type = "NoBr"))

MyACCf3Dat = subset(Fig3ACC, variable != "OptimA" & variable != "OptimB")

m3a <- ggplot(data = MyACCf3Dat, aes(x = variable, y = value, group = type)) + theme_bw() +
  geom_line(aes(color = type), size = 0.8) + 
  geom_point(aes(shape = type, color = type), size = 2.1) + 
  scale_x_discrete("Degree of freedom", labels = c("3", "5", "10", "20", "30", "40")) + 
  scale_y_continuous(name = "ACC Sample size", 
                     breaks = seq(0, 240, by = 60), limits = c(0, 240)) + 
  theme(legend.position='none', legend.title = element_blank()) 


F3ACC <- m3a + scale_color_manual(name = element_blank(), 
                                  labels = c("Robust weights A", "Robust weights B",
                                             "No robustification", "No borrowing"), 
                                  values = c("#ca0020", "#f4a582", "#0072B2", "lightgray")) +
  scale_shape_manual(name = element_blank(), 
                     labels = c("Robust weights A", "Robust weights B", 
                                "No robustification", "No borrowing"), 
                     values = c(16, 17, 0, 1)) + facet_wrap(~Case, ncol = 4) 


## ALC ##
MyALCf3atrs = melt(MyALCf3a, id = "Case")
MyALCf3btrs = melt(MyALCf3b, id = "Case")
MyALCf3ctrs = melt(MyALCf3c, id = "Case")
MyALCf3dtrs = melt(MyALCf3d, id = "Case")

Fig3ALC = rbind(data.frame(MyALCf3atrs, type = "RbW1"),
                data.frame(MyALCf3btrs, type = "RbW2"), 
                data.frame(MyALCf3ctrs, type = "NoRb"), 
                data.frame(MyALCf3dtrs, type = "NoBr"))

MyALCf3Dat = subset(Fig3ALC, variable != "OptimA" & variable != "OptimB")

m3b <- ggplot(data = MyALCf3Dat, aes(x = variable, y = value, group = type)) + theme_bw() +
  geom_line(aes(color = type), size = 0.8) + 
  geom_point(aes(shape = type, color = type), size = 2.1) + 
  scale_x_discrete("Degree of freedom", labels = c("3", "5", "10", "20", "30", "40")) + 
  scale_y_continuous(name = "ALC Sample size", 
                     breaks = seq(0, 240, by = 60), limits = c(0, 240)) +
  theme(legend.position='none', legend.title = element_blank()) 


F3ALC <- m3b + scale_color_manual(name = element_blank(), 
                                  labels = c("Robust weights A", "Robust weights B", 
                                             "No robustification", "No borrowing"), 
                                  values = c("#ca0020", "#f4a582", "#0072B2", "lightgray")) +
  scale_shape_manual(name = element_blank(), 
                     labels = c("Robust weights A", "Robust weights B"), 
                     values = c(16, 17, 0, 1)) + facet_wrap(~Case, ncol = 4) 

## APVC ##
MyAPVCf3atrs = melt(MyAPVCf3a, id = "Case")
MyAPVCf3btrs = melt(MyAPVCf3b, id = "Case")
MyAPVCf3ctrs = melt(MyAPVCf3c, id = "Case")
MyAPVCf3dtrs = melt(MyAPVCf3d, id = "Case")

Fig3APVC = rbind(data.frame(MyAPVCf3atrs, type = "RbW1"),
                 data.frame(MyAPVCf3btrs, type = "RbW2"),
                 data.frame(MyAPVCf3ctrs, type = "NoRb"),
                 data.frame(MyAPVCf3dtrs, type = "NoBr"))

MyAPVCf3Dat = subset(Fig3APVC, variable != "OptimA" & variable != "OptimB")

m3c <- ggplot(data = MyAPVCf3Dat, aes(x = variable, y = value, group = type)) + theme_bw() +
  geom_line(aes(color = type), size = 0.8) + 
  geom_point(aes(shape = type, color = type), size = 2.1) + 
  scale_x_discrete("Degree of freedom", labels = c("3", "5", "10", "20", "30", "40")) + 
  scale_y_continuous(name = "APVC Sample size", 
                     breaks = seq(0, 240, by = 60), limits = c(0, 240)) +
  theme(legend.position='none', legend.title = element_blank()) 

F3APVC <- m3c + scale_color_manual(name = element_blank(), 
                                   labels = c("Robust weights I", "Robust weights II",
                                              bquote("All"~w[k] == 0), bquote("All"~w[k] == 1)), 
                                   values = c("#ca0020", "#f4a582", "#0072B2", "lightgray")) +
  scale_shape_manual(name = element_blank(), 
                     labels = c("Robust weights I", "Robust weights II",
                                bquote("All"~w[k] == 0), bquote("All"~w[k] == 1)), 
                     values = c(16, 17, 0, 1)) + facet_wrap(~Case, ncol = 4) 

prow2 <- cowplot::plot_grid(F3ACC,
                            F3ALC,
                            F3APVC,
                            align = 'v',
                            labels = c("(i)", "(ii)", "(iii)"),
                            hjust = -0.5,
                            ncol = 1
)

# Figure 2: (eps) 800*825
prow2

## with legend:
# legend <- cowplot::get_legend(F3APVC)
# cowplot::plot_grid(legend, prow2, ncol = 1, rel_heights = c(0.05, 1))

#---------------- Figure 3 ----------------#
##--------------------------------------##
##--------- Unknown \sigma_0^2 ---------##
##--------------------------------------##

## ACC
ACCf4  = array(0, dim = c(4, 6, 5))
mycvrg = c(0.85, 0.875, 0.90, 0.95, 0.975)


for(p in 1:5){
  for(i in 1:4){
    # Optim 1
    ACCf4[i, 1, p] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], 
                                df = 30, l0 = 0.65, alpha = 1 - mycvrg[p], 
                                dw = c(2, 2), br = c(18, 3), s0 = 0.05) 
    
    # Optim 2
    ACCf4[i, 2, p] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], 
                                df = 30, l0 = 0.65, alpha = 1 - mycvrg[p], 
                                dw = c(2, 2), br = c(18, 3), s0 = 0.05) 
  }
}

for(p in 1:5){
  for(i in 1:4){
    for(k in 1:4){
      ACCf4[i, (k+2), p] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], 
                                      df = 3, l0 = 0.65, alpha = 1 - mycvrg[p], 
                                      dw = c(2, 2), br = c(18, 3), s0 = 0.05) 
    }
  }
}

round(ACCf4, 1)

for(p in 1:5){
  for(k in 1:6){
    assign(paste0("ACC2mycvrg", p, "Appr", k), 
           data.frame(Case = c("Configuration 1", "Configuration 2",
                               "Configuration 3", "Configuration 4"),
                      SS = ACCf4[,k,p],
                      Appr = paste0("Approach ", k), 
                      Cvrg = paste0("cvrg = ", mycvrg[p])),
    )
  }
}


mats <- grep(x = ls(pos=1), pattern="ACC2mycvrg", value=TRUE)
Fig4ACC <- do.call(rbind, mget(mixedsort(mats)))

m4a <- ggplot(data = Fig4ACC, aes(x = Cvrg, y = SS, group = Appr)) + theme_bw() +
  geom_line(aes(color = Appr), size = 0.8) + 
  geom_point(aes(shape = Appr, color = Appr), size = 2.1) + 
  scale_x_discrete("Average coverage probability (%)", 
                   labels = c("85", "87.5", "90", "95", "97.5")) + 
  scale_y_continuous(name = "ACC Sample size", 
                     breaks = seq(0, 320, by = 80), limits = c(0, 320)) +
  theme(legend.position='none', legend.title = element_blank())

F4AAC <- m4a + scale_color_manual(name = element_blank(), 
                                  labels = c("Optim A", "Optim B", "Robust weights A", "Robust weights B",
                                             "No robustification", "No borrowing"),
                                  values = c("#008837", "#a6dba0", "#ca0020", "#f4a582", 
                                             "#0072B2", "lightgray")) + 
  scale_shape_manual(name = element_blank(), 
                     labels = c("Optim A", "Optim B", "Robust weights A", "Robust weights B",
                                "No robustification", "No borrowing"), 
                     values = c(13, 11, 16, 17, 0, 1)) + facet_wrap(~Case, ncol = 4)


## ALC
ALCf4  = array(0, dim = c(4, 6, 5))
mylen = seq(0.55, 0.75, by = 0.05)

for(p in 1:5){
  for(i in 1:4){
    # Optim 1
    ALCf4[i, 1, p] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], 
                               df = 30, dw = c(2, 2), br = c(18, 3), 
                               alpha0 = 0.05, tlen = mylen[p], s0 = 0.05)
    
    
    # Optim 2
    ALCf4[i, 2, p] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], 
                               df = 30, dw = c(2, 2), br = c(18, 3), 
                               alpha0 = 0.05, tlen = mylen[p], s0 = 0.05)
  }
}

for(p in 1:5){
  for(i in 1:4){
    for(k in 1:4){
      ALCf4[i, (k+2), p] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], 
                                     df = 3, dw = c(2, 2), br = c(18, 3), 
                                     alpha0 = 0.05, tlen = mylen[p], s0 = 0.05)
    }
  }
}


for(p in 1:5){
  for(k in 1:6){
    assign(paste0("ALC2mylen", p, "Appr", k), 
           data.frame(Case = c("Configuration 1", "Configuration 2",
                               "Configuration 3", "Configuration 4"),
                      SS = ALCf4[,k,p],
                      Appr = paste0("Approach ", k), 
                      Len = paste0("l0 = ", mylen[p])),
    )
  }
}


mats <- grep(x = ls(pos=1), pattern="ALC2mylen", value=TRUE)
Fig4ALC <- do.call(rbind, mget(mixedsort(mats)))


m4b <- ggplot(data = Fig4ALC, aes(x = Len, y = SS, group = Appr)) + theme_bw() +
  geom_line(aes(color = Appr), size = 0.8) + 
  geom_point(aes(shape = Appr, color = Appr), size = 2.1) + 
  scale_x_discrete("Average HPD interval length", labels = c("0.55", "0.60", "0.65", "0.70", "0.75")) + 
  scale_y_continuous(name = "ALC Sample size", 
                     breaks = seq(0, 200, by = 50), limits = c(0, 200)) +
  theme(legend.position='none', legend.title = element_blank())

F4ALC <- m4b + scale_color_manual(name = element_blank(), 
                                  labels = c("Optim A", "Optim B", "Robust weights A", "Robust weights B",
                                             "No robustification", "No borrowing"),
                                  values = c("#008837", "#a6dba0", "#ca0020", "#f4a582", 
                                             "#0072B2", "lightgray")) + 
  scale_shape_manual(name = element_blank(), 
                     labels = c("Optim A", "Optim B", "Robust weights A", "Robust weights B",
                                "No robustification", "No borrowing"), 
                     values = c(13, 11, 16, 17, 0, 1)) + facet_wrap(~Case, ncol = 4) 


## APVC
APVCf4  = array(0, dim = c(4, 6, 5))
myeps = seq(0.01, 0.05, by = 0.01)

for(p in 1:5){
  for(i in 1:4){
    # Optim 1
    APVCf4[i, 1, p] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], df = 30,
                                  dw = c(2, 2), br = c(18, 3), epsil0 = myeps[p], s0 = 0.05)
    # Optim 2  
    APVCf4[i, 2, p] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], df = 30,
                                  dw = c(2, 2), br = c(18, 3), epsil0 = myeps[p], s0 = 0.05)
    
  }
}

for(p in 1:5){
  for(i in 1:4){
    for(k in 1:4){
      APVCf4[i, (k+2), p] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], df = 3,
                                        dw = c(2, 2), br = c(18, 3), epsil0 = myeps[p], s0 = 0.05)
    }
  }
}


for(p in 1:5){
  for(k in 1:6){
    assign(paste0("APVC2myeps", p, "Appr", k), 
           data.frame(Case = c("Configuration 1", "Configuration 2",
                               "Configuration 3", "Configuration 4"),
                      SS = APVCf4[,k,p],
                      Appr = paste0("Approach ", k), 
                      Eps = paste0("eps0 = ", myeps[p])),
    )
  }
}


mats <- grep(x= ls(pos=1), pattern="APVC2myeps", value=TRUE)

Fig4APVC <- do.call(rbind, mget(mixedsort(mats)))

m4c <- ggplot(data = Fig4APVC, aes(x = Eps, y = SS, group = Appr)) + theme_bw() +
  geom_line(aes(color = Appr), size = 0.8) + 
  geom_point(aes(shape = Appr, color = Appr), size = 2.1) + 
  scale_x_discrete("Average posterior variance", labels = c("0.01", "0.02", "0.03", "0.04", "0.05")) + 
  scale_y_continuous(name = "APVC Sample size", 
                     breaks = seq(0, 650, by = 150), limits = c(0, 670)) +
  theme(legend.position='none', legend.title = element_blank())

F4APVC <- m4c + scale_color_manual(name = element_blank(), 
                                   labels = c("Optim A", "Optim B", "Robust weights A", "Robust weights B",
                                              "No robustification", "No borrowing"),
                                   values = c("#008837", "#a6dba0", "#ca0020", "#f4a582", 
                                              "#0072B2", "lightgray")) + 
  scale_shape_manual(name = element_blank(), 
                     labels = c("Optim A", "Optim B", "Robust weights A", "Robust weights B",
                                "No robustification", "No borrowing"), 
                     values = c(13, 11, 16, 17, 0, 1)) + facet_wrap(~Case, ncol = 4)

prow3 <- cowplot::plot_grid( F4AAC,
                            F4ALC,
                            F4APVC,
                            align = 'v',
                            labels = c("(i)", "(ii)", "(iii)"),
                            hjust = -0.5,
                            ncol = 1
)

# Figure 3: (eps) 800*825
prow3

## with legend:
# legend <- cowplot::get_legend(F4APVC)
# cowplot::plot_grid(legend, prow3, ncol = 1, rel_heights = c(0.15, 1))


#---------------- Figure 4 ----------------#

##-------- Unknown \sigma_0^2 --------##
##--------------- ACC ----------------##
##------------------------------------##
ACCf5  = array(0, dim = c(4, 6, 3))
myGamma = rbind(c(6, 3), c(18, 3), c(54, 3))

for(p in 1:nrow(myGamma)){
  for(i in 1:4){
    # Optim 1
    ACCf5[i, 1, p] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], 
                                df = 30, l0 = 0.65, 
                                dw = c(2, 2), br = myGamma[p,], s0 = 0.05) 
    
    # Optim 2
    ACCf5[i, 2, p] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], 
                                df = 30, l0 = 0.65, 
                                dw = c(2, 2), br = myGamma[p,], s0 = 0.05) 
  }
}

for(p in 1:nrow(myGamma)){
  for(i in 1:4){
    for(k in 1:4){
      ACCf5[i, (k+2), p] = ACCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], 
                                      df = 3, l0 = 0.65, 
                                      dw = c(2, 2), br = myGamma[p,], s0 = 0.05) 
    }
  }
}

# round(ACCf5, 1)

for(p in 1:nrow(myGamma)){
  for(k in 1:6){
    assign(paste0("ACC3myGM", p, "Appr", k), 
           data.frame(Case = c("Configuration 1", "Configuration 2",
                               "Configuration 3", "Configuration 4"),
                      SS = ACCf5[,k,p],
                      Appr = paste0("Approach ", k), 
                      MixC = paste0("Gamma(", myGamma[p, 1], ", ", myGamma[p, 2], ")"))
    )
  }
}


mats <- grep(x= ls(pos=1), pattern="ACC3myGM", value=TRUE)

Fig5ACC <- do.call(rbind, mget(mixedsort(mats)))

myFig5ACC = subset(Fig5ACC, Appr != "Approach 6")

m5a <- ggplot(data = myFig5ACC, aes(x = Appr, y = SS, group = MixC)) + theme_bw() +
  geom_point(aes(shape = MixC, color = MixC), size = 2.3) + 
  scale_x_discrete("", labels = c("OpI", "OpII", "RbI", "RbII",
                                  "NoRb")) + 
  scale_y_continuous(name = "ACC Sample size", 
                     breaks = seq(0, 210, by = 40), limits = c(0, 215)) +
  theme(legend.position='none', legend.title = element_blank()) + facet_wrap(~Case, ncol = 4) 


##-------- Unknown \sigma_0^2 --------##
##--------------- ALC ----------------##
##------------------------------------##
ALCf5  = array(0, dim = c(4, 6, 3))

for(p in 1:nrow(myGamma)){
  for(i in 1:4){
    # Optim 1
    ALCf5[i, 1, p] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], 
                               df = 30, dw = c(2, 2), myGamma[p,], 
                               alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
    
    # Optim 2
    ALCf5[i, 2, p] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], 
                               df = 30, dw = c(2, 2), myGamma[p,], 
                               alpha0 = 0.05, tlen = 0.65, s0 = 0.05) 
  }
}

for(p in 1:nrow(myGamma)){
  for(i in 1:4){
    for(k in 1:4){
      ALCf5[i, (k+2), p] = ALCunkvar(MySS = 1:300, R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], 
                                     df = 3, dw = c(2, 2), myGamma[p,], 
                                     alpha0 = 0.05, tlen = 0.65, s0 = 0.05)
    }
  }
}

# round(ALCf5, 1)


for(p in 1:nrow(myGamma)){
  for(k in 1:6){
    assign(paste0("ALC3myGM", p, "Appr", k), 
           data.frame(Case = c("Configuration 1", "Configuration 2",
                               "Configuration 3", "Configuration 4"),
                      SS = ALCf5[,k,p],
                      Appr = paste0("Approach ", k), 
                      MixC = paste0("Gamma(", myGamma[p, 1], ", ", myGamma[p, 2], ")"))
    )
  }
}


mats <- grep(x= ls(pos=1), pattern="ALC3myGM", value=TRUE)

Fig5ALC <- do.call(rbind, mget(mixedsort(mats)))

myFig5ALC = subset(Fig5ALC, Appr != "Approach 6")

m5b <- ggplot(data = myFig5ALC, aes(x = Appr, y = SS, group = MixC)) + theme_bw() +
  geom_point(aes(shape = MixC, color = MixC), size = 2.3) + 
  scale_x_discrete("", labels = c("OpI", "OpII", "RbI", "RbII",
                                  "NoRb")) + 
  scale_y_continuous(name = "ALC Sample size", 
                     breaks = seq(0, 210, by = 40), limits = c(0, 215)) +
  theme(legend.position='none', legend.title = element_blank()) + facet_wrap(~Case, ncol = 4) 


##-------- Unknown \sigma_0^2 --------##
##-------------- APVC ----------------##
##------------------------------------##
APVCf5  = array(0, dim = c(4, 6, 3))

for(p in 1:nrow(myGamma)){
  for(i in 1:4){
    # Optim 1
    APVCf5[i, 1, p] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,1,i], df = 30,
                                  dw = c(2, 2), br = myGamma[p,], epsil0 = 0.03, s0 = 0.05)
    
    # Optim 2
    APVCf5[i, 2, p] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,2,i], df = 30,
                                  dw = c(2, 2), br = myGamma[p,], epsil0 = 0.03, s0 = 0.05)
  }
}

for(p in 1:nrow(myGamma)){
  for(i in 1:4){
    for(k in 1:4){
      APVCf5[i, (k+2), p] = APVCunknvar(R = 0.5, sk2 = HisDatV2[,i], wk = modWgt[,k,i], df = 3,
                                        dw = c(2, 2), br = myGamma[p,], epsil0 = 0.03, s0 = 0.05)
    }
  }
}

# round(ALCf5, 1)


for(p in 1:nrow(myGamma)){
  for(k in 1:6){
    assign(paste0("APVC3myGM", p, "Appr", k), 
           data.frame(Case = c("Configuration 1", "Configuration 2",
                               "Configuration 3", "Configuration 4"),
                      SS = APVCf5[,k,p],
                      Appr = paste0("Approach ", k), 
                      MixC = paste0("Gamma(", myGamma[p, 1], ", ", myGamma[p, 2], ")"))
    )
  }
}


mats <- grep(x = ls(pos=1), pattern="APVC3myGM", value=TRUE)

Fig5APVC <- do.call(rbind, mget(mixedsort(mats)))

myFig5APVC = subset(Fig5APVC, Appr != "Approach 6")

m5c <- ggplot(data = myFig5APVC, aes(x = Appr, y = SS, group = MixC)) + theme_bw() +
  geom_point(aes(shape = MixC, color = MixC), size = 2.3) + 
  scale_x_discrete("", labels = c("OpI", "OpII", "RbI", "RbII",
                                  "NoRb")) + 
  scale_y_continuous(name = "APVC Sample size", 
                     breaks = seq(0, 210, by = 40), limits = c(0, 215)) +
  theme(legend.position='none', legend.title = element_blank()) + facet_wrap(~Case, ncol = 4) 



prow4 <- cowplot::plot_grid( m5a,
                            m5b,
                            m5c,
                            align = 'v',
                            labels = c("(i)", "(ii)", "(iii)"),
                            hjust = -0.5,
                            ncol = 1
)

# Figure 4: (eps) 800*825
prow4

## with legend:
# legend <- cowplot::get_legend(m5c)
# cowplot::plot_grid(legend, prow4, ncol = 1, rel_heights = c(0.05, 1))
