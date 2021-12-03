# setwd("~/Documents/SSE for studies using historical data/Code/Code_Nov2021")
setwd("~/Dropbox/CRUK_Cam/SSE for studies using historical data/Code/Code_Nov2021")
source("Bayesian SSD using robust commensurate priors.R")

library(ggplot2)
require(gtools)
require(reshape2)
#---------------- the collective prior in Section 3 ----------------#
mywk = c(0.15, 0.20, 0.17, 0.13, 0.20)

myp = pq(wk = mywk, s0 = 0.05)
round(myp, 2)
sum(round(myp, 2))

IllSc = cbind(c(-0.26, -0.24, -0.37, -0.34, -0.32), 
              c( 0.25,  0.23,  0.22,  0.36,  0.26))

round(cprior(p = myp, mk = IllSc[,1], sk2 = IllSc[,2], wk = mywk, 
             dw = c(2, 2), br = c(18, 3)), 3)


#---------------- configuration of historical data ----------------#
## Configuration 1: Consistent means and similar level of variability between historical datasets
MySc1 = cbind(c(-0.26, -0.24, -0.37, -0.34, -0.32), 
              c( 0.25,  0.23,  0.22,  0.36,  0.26))

## Configuration 2: Consistent means and similar level of variability between historical datasets, more informative
MySc2 = cbind(c(-0.26, -0.24, -0.37, -0.34, -0.32), 
              rep(0.10, 5))

## Configuration 3: Divergent means and levels of variability between historical datasets
MySc3 = cbind(c(-0.26, -0.17, -0.44, -0.15, 0.12), 
              c( 0.25,  0.64,  0.97,  1.54, 0.59))

## Configuration 4: Divergent means and levels of variability between historical datasets, more informative
MySc4 = cbind(c(-0.26, -0.17, -0.44, -0.15, 0.12), 
              c( 0.25,  0.15,  0.40,  0.89, 0.22))


## Visualise the pairwise hellinger distance between historical data ##
myHdMat = rbind(data.frame(HdMat(mu.hist = MySc1[, 1], var.hist = MySc1[, 2]), Sc = "Configuration 1"),
                data.frame(HdMat(mu.hist = MySc2[, 1], var.hist = MySc2[, 2]), Sc = "Configuration 2"),
                data.frame(HdMat(mu.hist = MySc3[, 1], var.hist = MySc3[, 2]), Sc = "Configuration 3"),
                data.frame(HdMat(mu.hist = MySc4[, 1], var.hist = MySc4[, 2]), Sc = "Configuration 4")
)

p0 <- ggplot(data = myHdMat, aes(x = id2, y = id1, fill = V3)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0, 0.5), space = "Lab", 
                       name = "Hellinger distance\n") +
  facet_wrap(~Sc, ncol = 2) + theme_bw() 

# Fig S4: eps 600*465
p0 + xlab("Historical data") + ylab("Historical data") +
  geom_text(aes(x = id2, y = id1, label = V3), color = "black", size = 3) 


#---------------- Table 1 ----------------#
MyWgt1 = c(0.103, 0.175, 0.081, 0.143, 0.077)
MyWgt2 = c(0.252, 0.319, 0.140, 0.306, 0.149)
MyWgt3 = c(0.103, 0.175, 0.081, 0.143, 0.077)
MyWgt4 = c(0.252, 0.319, 0.140, 0.306, 0.149)
MyWgt5 = c(0.101, 0.219, 0.385, 0.385, 0.304)
MyWgt6 = c(0.325, 0.203, 0.171, 0.180, 0.272)
MyWgt7 = c(0.066, 0.303, 0.459, 0.355, 0.115)
MyWgt8 = c(0.537, 0.306, 0.054, 0.220, 0.350)


for(i in 1:4){
  assign(paste0("Configuration", i),
         round(rbind(t(get(paste0("MySc", i))), 
              get(paste0("MyWgt", i*2-1)), pq(wk = get(paste0("MyWgt", i*2-1)), s0 = 0.05), 
              get(paste0("MyWgt", i*2)), pq(wk = get(paste0("MyWgt", i*2)), s0 = 0.05)), 3)
         )
  
  assign(paste0("ClPrior", i, "RWI"),
         round(cprior(p = pq(wk = get(paste0("MyWgt", i*2-1)), s0 = 0.05), 
                      mk = get(paste0("MySc", i))[,1], sk2 = get(paste0("MySc", i))[,2], 
                      wk = get(paste0("MyWgt", i*2-1)), 
                      dw = c(2, 2), br = c(18, 3)), 3)
         )
  
  assign(paste0("ClPrior", i, "RWII"),
         round(cprior(p = pq(wk = get(paste0("MyWgt", i*2)), s0 = 0.05), 
                      mk = get(paste0("MySc", i))[,1], sk2 = get(paste0("MySc", i))[,2], 
                      wk = get(paste0("MyWgt", i*2)), 
                      dw = c(2, 2), br = c(18, 3)), 3)
  )
}

# Scenarios and Robust weights I & II in Table 1
Configuration1; Configuration2; Configuration3; Configuration4

ClPrior1RWI; ClPrior1RWII
ClPrior2RWI; ClPrior2RWII
ClPrior3RWI; ClPrior3RWII
ClPrior4RWI; ClPrior4RWII
