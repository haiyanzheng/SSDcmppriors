Copyright (c) 2022 Haiyan Zheng

This file is part of the IDENT project.

The IDENT project can not be copied and/or distributed without the express permission of Haiyan Zheng <haiyan.zheng@mrc-bsu.cam.ac.uk>.

# SSDcmspriors

This repository contains R functions to implement Bayesian sample size determination for experiments comparing two groups, for which relevant pre-experimental information from multiple sources can be incoporated in a robust prior to support both the design and analysis. Readers of the following paper can use the R functions to reproduce the numerical results and/or figures reported.

**H Zheng, T Jaki, JMS Wason. (2021) Bayesian sample size determination using commensurate priors to leverage pre-experimental data. Submitted.**


For any issues or technical questions relating to the implementation, please contact Dr Haiyan Zheng at haiyan.zheng@mrc-bsu.cam.ac.uk

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Content of the folder
The file describes the infrastructure of the numerical studies for
the Bayesian sample size determination approaches written in the manuscript. This infrastructure is composed of the following files/folders:

-  [[./Bayesian SSD using arbitrary priors.R]], [[./Bayesian SSD using robust commensurate priors.R]] define functions used to compute the Bayesian sample sizes 
according to the respective criteria (i.e., ACC, ALC and APVC) for cases of known and unknown variances. 
Specifically, the former incorporates a prior for the mean difference, 
specified based upon pre-experimental information from a single source; 
whereas, the latter can accommodate pre-experimental information from multiple sources.

-  [[./Historical data.R]]: R code to stipulate the four configurations of historical data for sample size determination. 
The collective priors contained in Table 1 of the manuscript can be reproduced by the functions in this R script.

-  [[./Figures1-4.R]]: R code used to generate figures 1 - 4 of the manuscript. 
With this, the user will know how to implement the proposed Bayesian approach to sample size determination.

-  [[./MyMod.txt]] is the OpenBUGS code to be called in R for fitting the proposed Bayesian model (no Normal approximation involved) 
using Markov chain Monte Carlo for cases of known variance.  

-  [[./MyMod2.txt]] is the OpenBUGS code to be called in R for fitting the proposed Bayesian model (no Normal approximation involved) 
using Markov chain Monte Carlo for cases of unknown variance.  

-  [[./Simulator (exact Bayesian).R]] defines the functions to simulate e.g., 100,000 replicates of new experiments,
which are to be analysed to compute the average coverage, length of the highest density function,
and expected variance of the posterior distribution. 
The posterior distribution will be derived based on exact Bayesian inference, 
so the normal approximation in expression (4) of the manuscript will be used. 

-  [[./Simulator (MCMC).R]] defines the functions to simulate e.g., 100,000 replicates of new experiments, 
which are to be analysed to compute the average coverage, length of the highest density function, and expected variance of the posterior distribution. 
The posterior distribution will be derived based on Markov chain Monte Carlo samples, 
calling the OpenBUGS models, titled 'MyMod.txt' or 'MyMod2.txt'. 
The normal approximation in expression (4) of the manuscript is not used.


Note: The user would first run [[./Historical data.R]] and [[./Figures1-4.R]] for reproducible figures presented in the manuscript. 
The simulators, especially the ones using Markov chain Monte Carlo may cost time to simulate and analyse 100,000 replicates of the new experiment.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Run R code
#+BEGIN_SRC R :exports both :results output :session *R* :cache no
source("Bayesian SSD using robust commensurate priors.R")
source("Historical data.R")
#+END_SRC

* R and package versions
#+BEGIN_SRC R  :results output   :exports results  :session *R* :cache yes 
sessionInfo()
#+END_SRC

#+RESULTS[<2021-12-14 19:36:25>]:
#+begin_example
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.3 LTS

Matrix products: default
BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
 [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C      

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base    

other attached packages:
[1] R2OpenBUGS_3.2-3.2 sn_2.0.1           reshape2_1.4.3     gtools_3.8.1       ggplot2_3.3.5 

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7          pillar_1.6.3        compiler_3.4.4      plyr_1.8.6          tools_3.4.4        
 [6] boot_1.3-20         digest_0.6.28       lattice_0.20-35     lifecycle_1.0.1     tibble_3.1.5       
[11] gtable_0.3.0        pkgconfig_2.0.3     rlang_0.4.11        DBI_1.1.0           coda_0.19-3        
[16] withr_2.4.2         dplyr_1.0.7         stringr_1.4.0       generics_0.1.0      vctrs_0.3.8        
[21] grid_3.4.4          tidyselect_1.1.1    cowplot_0.9.4       glue_1.4.2          R6_2.5.1           
[26] fansi_0.5.0         purrr_0.3.4         farver_2.1.0        magrittr_2.0.1      scales_1.1.1       
[31] ellipsis_0.3.2      mnormt_2.0.2        assertthat_0.2.1    colorspace_2.0-2    numDeriv_2016.8-1.1
[36] labeling_0.4.2      utf8_1.2.2          stringi_1.7.5       munsell_0.5.0       tmvnsim_1.0-2      
[41] crayon_1.4.1  
#+end_example

