
# Model Performance. Simulation

This code replicates the simulations done to measure the model
performance as explained in the Model Performance section 5 of the
article.

It generates the data needed to construct the OC curves under 9
different operating conditions scenarions and 25 delta values. An extra
week (week no. 49) is simulated for each of the actual printers with an
out of control observation under a different delta each time.

This script takes long to run as it performs multiple simulations on
multiple scenarios. The results of the simulations are already stored in
data/results.RData so script 04OCcurves.Rmd may run without the need to
complete the current script. Futhermore, data/results.RData has the
simulation results presented in the article.

## 1. Library loading

``` r
library(ggplot2)
library(dplyr)
```

## 2. Data loading

``` r
load(file = 'model.RData')

nObs = 8113
nMaq = 357
```

## 3. Procedure

A. generate 9 printer operating conditions as combination of printer
usage levels (low, medium, high) and printing times (15, 45, 90)h

B. for each operating conditions and printer, estimate the posterior
predictive pmf for an in control and out of control situation using the
whole 30,000 samples from the model posterior

C. Estimate P(Zjt = 1)/P(Zjt=0) for each printer j = 1…357 as this value
will be used in the further loop and will remain constant

D. Set the loop indices and final result storage objects.

for k = 1 to 9 do: <br/> for delta = 0 to 6 steps = 0.25 do: <br/> E.
sample y_jt+1 OOC observations from the posterior predictive
distribution <br/> F. Estimate Zjt+1 for each yjt+1 observations and
compute odds <br/> save intermediate results <br/> end for <br/> merge
results with printer scenario description <br/> endfor <br/>

G. Save results

### 3.1 Define operating conditions scenario

``` r
# A. define printer operating conditions
#    -----------------------------------

# define printing times and usage levels
printingTime <- c(15, 45, 90)
printerTypes <- data.frame(PVCsecPerc = c(0.0, 0.5, 1.0),
                           jobsPerTime = c(0.5, 1.0, 3.0),
                           type = c('low','medium','high'))

# create levels-printingTime combinations and label them with a number
levels <- data.frame(setNames(expand.grid(printingTime, printerTypes$type), 
                              nm = c('timeImHour','type')),
                     levelNo = 1:(length(printingTime)*nrow(printerTypes)))

# create levels-printingTime combinations and label them with a number
levels <- merge( x = levels,
                 y = printerTypes,
                 by = 'type')
levels <- levels[order(levels$levelNo),]
```

### 3.2 Sample from posterior predictive

Estimate posterior predictive pmf for each printer and scenario in and
out of control

``` r
# B.1 create B matrix with fixed effects and prepare empty list to hold the
#     epmf for IC and OOC situation 

B <- cbind(model$BUGSoutput$sims.list$b0,
           model$BUGSoutput$sims.list$b1,
           model$BUGSoutput$sims.list$b2)

IClist <- vector(mode = 'list', length = 9)
OOClist <- vector(mode = 'list', length = 9)

k <- 1

# B.2 for each combination of fixed effects values (low, medium, high)

for (i in 1: nrow(printerTypes)) {
  
  # B.2.1 sample model$BUGSoutput$n.sims posterior IC and OOC rate for each 
  #       printer
  
  Bx <- B %*% c(1,as.numeric(printerTypes[i,c('PVCsecPerc','jobsPerTime')]))
  lambdaPostIC <- exp(as.numeric(Bx) + model$BUGSoutput$sims.list$alpha )
  lambdaPostOOC <- lambdaPostIC * exp(as.numeric(model$BUGSoutput$sims.list$delta))
  
  # B.2.2 for each printing time (15, 45, 90 h)
  
  for (Tjt in printingTime){
    
    print(paste('scenario:',k, '(T,pvc,jobs):',Tjt, 
                paste(printerTypes[i,c('PVCsecPerc','jobsPerTime')], collapse = ' ')))
  
    # B.2.2.1 sample for nSims from the IC posterior predictive for every printer
      
      thetaPostIC <- lambdaPostIC * Tjt
      y <- matrix(rpois(length(thetaPostIC), lambda = thetaPostIC), ncol = nMaq)
    
    # B.2.2.2 compute the IC pmf for each printer based on its posterior predictive
    #         use factor(y, levels = 0:max(y)) to ensure all counts, including
    #         those with 0 observations are included
    #         use stack to place table in data.frame format
    
      IClist[[k]] <- apply(y, MARGIN = 2, 
                           FUN = function(y) stack(prop.table(table(factor(y, levels = 0:max(y))))))
      
      
    # B.2.2.3 sample for nSims from the OOC posterior predictive for every printer
      
      thetaPostOOC <- lambdaPostOOC * Tjt
      y <- matrix(rpois(length(thetaPostOOC), lambda = thetaPostOOC), ncol = nMaq)
    
    # B.2.2.4 compute the pmf for each printer based on its posterior predictive
    
      OOClist[[k]] <- apply(y, MARGIN = 2, 
                            FUN = function(y) stack(prop.table(table(factor(y, levels = 0:max(y))))))
      
    # update k
      
      k <- k+1

  }
}
```

    ## [1] "scenario: 1 (T,pvc,jobs): 15 0 0.5"
    ## [1] "scenario: 2 (T,pvc,jobs): 45 0 0.5"
    ## [1] "scenario: 3 (T,pvc,jobs): 90 0 0.5"
    ## [1] "scenario: 4 (T,pvc,jobs): 15 0.5 1"
    ## [1] "scenario: 5 (T,pvc,jobs): 45 0.5 1"
    ## [1] "scenario: 6 (T,pvc,jobs): 90 0.5 1"
    ## [1] "scenario: 7 (T,pvc,jobs): 15 1 3"
    ## [1] "scenario: 8 (T,pvc,jobs): 45 1 3"
    ## [1] "scenario: 9 (T,pvc,jobs): 90 1 3"

``` r
# NOTE: all pmf's have all counts included, even the ones with no obsevations
#       recorded

# B.3 cast table factor labels (no. counts) into numeric

IClist <- lapply(IClist, function(lst) {
  lapply(lst, function(pmf) data.frame(cmf = cumsum(pmf$values), 
                                       pmf = pmf$values, 
                                       y = as.numeric(levels(pmf$ind))))
})

OOClist <- lapply(OOClist, function(lst) {
  lapply(lst, function(pmf) data.frame(cmf = cumsum(pmf$values),
                                       pmf = pmf$values, 
                                       y = as.numeric(levels(pmf$ind))))
})
```

### 3.4 Compute Odds

``` r
# C. determine Ppost(Z)
#    ------------------

load(file = 'data/printWeeks.RData') #printer id and week for the whole data set

pPostZ1 <- setNames(aggregate(model$BUGSoutput$mean$z ~ printWeeks$id, FUN = mean), nm = c('id','Z'))
pPostZ0 <- data.frame(id = pPostZ1$id, Z = 1-pPostZ1$Z)
ratioPpost <- pPostZ1$Z/pPostZ0$Z

# D. boiler platting
#    ---------------

S <- 5000 
deltas <- seq(from = 0, to = 6, by = 0.25)
nReplicas <- nrow(levels) * length(deltas)
results <- data.frame(scenario = integer(nReplicas), delta = numeric(nReplicas), 
                      tp = integer(nReplicas), fn = integer(nReplicas))

j <- 1 #current row in the results df

for(k in 1:9) {
  
  cat('\n')
  print(paste('scenario:', k, ' (PVC,jobs,time) = (', 
              paste(levels[levels$levelNo == k, 
                           c('PVCsecPerc','jobsPerTime','timeImHour')],
                    collapse =','),')', sep =''))
  print('delta:')
  
  for (delta in deltas){
    
    # start.time <- Sys.time()

    cat(delta)
    cat(', ')
    
    # E. sample y_jt+1 OOC observations from the posterior distribution
    
    # E.1 sample S random parameters from the posterior
    
    # sims.matrix has the following 8477 columns (and 3,000 rows for each
    # mcmc sample)
    # 1:357 the alphas_j for the 357 printers
    # 358:360 b0, b1, b2
    # 361 delta
    # 362 deviance
    # 363 p
    # 364 tau2
    # 364:8477 the 8113 printer-weeks Z_jt
    
    ids <- sort(sample(model$BUGSoutput$n.sims, S))
    B <- model$BUGSoutput$sims.matrix[ids,c('b0','b1','b2')]
    A <- model$BUGSoutput$sims.matrix[ids, 1:357]
    
    # E.2 sample theta_IC and theta_OOC
    
    Bx <- B %*% as.numeric(c(1, levels[levels$levelNo == k, c('PVCsecPerc','jobsPerTime')]))
    thetaIC <- exp(as.numeric(Bx) + A) * levels[levels$levelNo == k, 'timeImHour']
    thetaOOC <- thetaIC * exp(delta)
    
    # E.3 sample y_jt+1 from the posterior predictive
    
    ytilde <- matrix(rpois(n = length(thetaOOC), lambda = thetaOOC), 
                           nrow = nrow(thetaOOC))
    
    # F Estimate P(yjt+1|IC) and P(yjt+1|OOC)
    
    # F.1 chek ytilde is not larger/smaller than the posterior predictive pmf
    
    maxICpmf <- sapply(IClist[[k]], function(pmf) max(pmf$y))
    minOOCpmf <- sapply(OOClist[[k]], function(cmf) min(which(cmf$cmf>0))-1)
    
    truePositives <- 0
    falseNegatives <- 0
    
    for (i in 1:nMaq) {
      
      pY <- data.frame(y = ytilde[,i], ic = rep(999,S), ooc = rep(999,S), ratio = rep(999,S))
      
      # empirical pmf's may miss some values and will not have observations 
      # above or below their observed maximum, so caution needed as zero will
      # appear in the first case and NA in the second
      
      # if ytilde is larger than the largest IC value it is OOC
      indLarge <- ytilde[,i] > maxICpmf[i]
      pY[indLarge,'ratio'] = 2
      
      # if ytilde is smaller than the smallest OOC value it is IC
      indSmall <- ytilde[,i] < minOOCpmf[i]
      pY[indSmall,'ratio'] = 0.1
      
      # if ytilde is zero then it is IC
      indZero <- (ytilde[,i] == 0 & !(indSmall | indLarge))
      pY[indZero,'ratio'] <- 0.1
      
      # if ytilde is within limits
      ind <- !(indZero | indSmall | indLarge)
      pY[ind, 'ic'] <- IClist[[k]][[i]][ytilde[ind,i]+1,'pmf']
      pY[ind, 'ooc'] <- OOClist[[k]][[i]][ytilde[ind,i]+1,'pmf']
      
      # if ytilde is within limits but p(yIC) or p(yOOC) is zero
      # it may mean that there were no posterior predictive replicas in 
      # a given y, not that its probability was zero.
      pY[pY$ooc == 0,'ooc'] = 1/(2*model$BUGSoutput$n.sims)
      pY[pY$ic == 0,'ic'] = 1/(2*model$BUGSoutput$n.sims)
      
      pY$ratio[ind] <- pY$ooc[ind]/pY$ic[ind] * ratioPpost[i]
     
      # check results
      
      truePositive_i <- sum(pY$ratio > 1)
      truePositives <- truePositives + truePositive_i
      falseNegatives <- falseNegatives + (S - truePositive_i)
      
    }
    
    results[j,] <- c(k, delta, truePositives, falseNegatives)
    j <- j + 1
    
    #end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # print(time.taken)
  }
}
```

    ## 
    ## [1] "scenario:1 (PVC,jobs,time) = (0,0.5,15)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:2 (PVC,jobs,time) = (0,0.5,45)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:3 (PVC,jobs,time) = (0,0.5,90)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:4 (PVC,jobs,time) = (0.5,1,15)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:5 (PVC,jobs,time) = (0.5,1,45)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:6 (PVC,jobs,time) = (0.5,1,90)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:7 (PVC,jobs,time) = (1,3,15)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:8 (PVC,jobs,time) = (1,3,45)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 
    ## [1] "scenario:9 (PVC,jobs,time) = (1,3,90)"
    ## [1] "delta:"
    ## 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6,

``` r
results <- merge(x = results,
                 y = levels,
                 by.x = 'scenario',
                 by.y = 'levelNo')

# G. Save performance results

# Be careful, saving the current simulation results will change the simulation
# used to evaluate the performance as stated in the article.

# save(results, file = 'data/results.RData')
```
