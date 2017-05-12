rm(list = ls())
gc()

## Performance Test for mixed solver
## The working directory should still be the R directory
## Randomly generate initial estimates (but from same distributions)

setwd("/data/example-models/PKPD/torsten/R")

rm(list = ls())
gc()

modelName <- "fTwoCpt"
# modelName <- "fTwoCpt_mixed"
dataName <- "fTwoCpt"  ## use the same data for both models
testName <- modelName

N <- 100  # 100 # number of times we want to run the test

## Adjust directories to your settings.
## CHECK which directories are needed?
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path("deliv", "figure", modelName)
tabDir <- file.path("deliv", "table", modelName)
modelDir <- file.path(projectDir, modelName)
dataDir <- file.path(projectDir, dataName)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path("tools")
stanDir <- file.path("cmdstan")
tempDir <- file.path(modelDir, modelName, "temp")

# source(file.path(scriptDir, "pkgSetup.R"))
.libPaths("lib")
library(rstan)  ## make sure rstan 2.15 gets used
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(parallel)

library(qapply)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "cmdStanTools.R"))

rstan_options(auto_write = TRUE)

set.seed(11191951) ## not required but assures repeatable results

###############################################################################

## True parameters used to simulate the data
## (see *Simulation.R file)
nParameters <- 11
thetaTrue <- rep(0, nParameters)
thetaTrue[1] <- 10  # CL
thetaTrue[2] <- 15  # Q
thetaTrue[3] <- 35  # VC
thetaTrue[4] <- 105  # VP
thetaTrue[5] <- 2.0  # KA
thetaTrue[6] <- sqrt(0.001)  # sigma_PK
thetaTrue[7] <- 125  # MTT
thetaTrue[8] <- 5  # Circ0
thetaTrue[9] <- 3e-4  # alpha
thetaTrue[10] <- 0.17  # gamma
thetaTrue[11] <- sqrt(0.001)  # sigma_PD

dir.create(file.path(modelDir, modelName))

## Compile Stan Model
compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

## Randomly generate initial estimates
CLPrior = 10
QPrior = 15
VCPrior = 35
VPPrior = 105
kaPrior = 2
CLPriorCV = 0.10
QPriorCV = 0.18
VCPriorCV = 0.14
VPPriorCV = 0.17
kaPriorCV = 0.16

circ0Prior <- 5
circ0PriorCV <- 0.20
mttPrior <- 125
mttPriorCV <- 0.2
gammaPrior <- 0.17
gammaPriorCV <- 0.2
alphaPrior <- 2.0E-4
alphaPriorCV <- 0.2

init <- function(){
  list(CL = exp(rnorm(1, log(CLPrior), CLPriorCV)),
       Q = exp(rnorm(1, log(QPrior), QPriorCV)),
       VC = exp(rnorm(1, log(VCPrior), VCPriorCV)),
       VP = exp(rnorm(1, log(VPPrior), VPPriorCV)),
       ka = exp(rnorm(1, log(kaPrior), kaPriorCV)),
       sigma = runif(1, 0.0001, 2),
       alpha = exp(rnorm(1, log(alphaPrior), alphaPriorCV)),
       mtt = exp(rnorm(1, log(mttPrior), mttPriorCV)),
       circ0 = exp(rnorm(1, log(circ0Prior), circ0PriorCV)),
       gamma = exp(rnorm(1, log(gammaPrior), gammaPriorCV)),
       sigmaNeut = runif(1, 0.0001, 2))
}

dir.create(tempDir)  ## directory to store initial estimates for each chain

## 4 performance metrics (only the first two really matter):
# 1) fractional difference between estimated and real mean value of parameters
# 2) time to compute 1000 effective independent samples
# 3) Run times
# 4) flag: simulation iteration
nMetrics = 4

## Write function to run test
StanFit <- function(iSim) {
  
  ## Matrix to store results
  ## Add one row for the log posterior
  pMatrix <- matrix(nrow = nParameters + 1, ncol = nMetrics)
  
  ## CHECK - do I need parametersToPlot?
  parametersToPlot <- c("CL", "Q", "VC", "VP", "ka", "sigma", "mtt", "circ0", "alpha",
                        "gamma", "sigmaNeut")
  parametersToPlot <- c("lp__", parametersToPlot)
  otherRVs <- c("cObsPred", "neutPred")
  parameters <- c(parametersToPlot, otherRVs)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  nChains <- 4 # 4
  nPost <- 100 # 1000 ## Number of post-burn-in samples per chain after thinning
  nBurn <- 100  # 1000 ## Number of burn-in samples per chain after thinning
  nThin <- 1
  chains <- 1:nChains
  
  nIter <- nPost * nThin
  nBurnin <- nBurn * nThin
  
  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()
  
  mclapply(chains,
           function(chain, model, data, iter, warmup, thin, init) {
             dir.create(file.path(tempDir, iSim))
             tempDirSim <- file.path(tempDir, iSim, chain)
             dir.create(tempDirSim)
             inits <- init()
             with(inits, stan_rdump(ls(inits), file = file.path(tempDirSim,
                                                                "init.R")))
             runModel(model = model, data = data,
                      iter = iter, warmup = warmup, thin = thin,
                      init = file.path(tempDirSim, "init.R"), 
                      seed = sample(1:999999, 1),
                      chain = chain, refresh = 100,
                      adapt_delta = 0.95, stepsize = 0.01,
                      tag = iSim)
             },
           model = file.path(modelDir, modelName),
           data = file.path(dataDir, paste0(dataName, ".data.R")),
           init = init,
           iter = nIter, warmup = nBurnin, thin = nThin,
           mc.cores = min(nChains, detectCores()))

  ## Save results 
  fit <- read_stan_csv(file.path(modelDir, modelName, 
                                 paste0(modelName, "_", iSim, "_",
                                        chains, ".csv")))
  outputName <- paste0(modelName, iSim)
  dir.create(outDir)
  save(fit, file = file.path(outDir, paste0(outputName, "Fit.Rsave_", iSim)))

  #############################################################################
  ## compute performance metrics

  ## Get the run-time (sum of run time for all 4 chains)
  RunTimes <- get_elapsed_time(fit, parametersToPlot)
  RunTime <- sum(RunTimes[ ,1]) + sum(RunTimes[ , 2])
  
  ## create a table that contains the mean value of the parameters
  ## (obtained by combining the mean of all the chains, with index nChains + 1), 
  ## n_eff, and computation time required to generated 1000 independent samples.
  mean <- as.vector(get_posterior_mean(fit, 
                                       pars = parametersToPlot)[ , nChains + 1])
  n_eff <- as.vector(summary(fit, pars = parametersToPlot)$summary[ ,"n_eff"])
  comp_time <- RunTime/n_eff * 1000
  
  ## save parameters in the parameter matrix
  pMatrix[ , 1] <- mean
  pMatrix[ , 2] <- n_eff
  pMatrix[ , 3] <- comp_time
  pMatrix[ , 4] <- iSim
  
  pMatrix
  
} ## end of Stan fit definition

## Run Job
## Check status with qstat("-f"), qdel("-u charlesm")
qapply(1:N, StanFit, global = TRUE, nCores = N, tag = "StanFit",
       internalize = FALSE)
object <- qinternalize(file.path("out", "StanFit"))

pMatrix <- array(0, dim = c(dim(object[[1]]), length(object)))
for (i in 1:length(object)) {
  pMatrix[ , , i] <- object[[i]]
}

## write pMatrix into csv file (want to save raw data) -- FIX ME
save(pMatrix, file = file.path(outDir, paste0(modelName, "pMatrixSave")))

###############################################################################
## save results in performance table and output in a csv file
N <- dim(pMatrix)[3]  ## FIX ME - shouldn't need this
performance_table <- matrix(0, nrow = N, ncol = nParameters * 3)

## compute:
##  1) fractional difference with real parameters
##  2) effective number of independent samples
##  3) Run times
performance_table[ , 1:nParameters] <- t((pMatrix[2:(nParameters + 1), 1, 1:N] - thetaTrue)
                                         / thetaTrue * 100)
performance_table[ , (nParameters + 1):(2 * nParameters)] <- t((pMatrix[2:(nParameters + 1), 
                                                                    2, 1:N]))
performance_table[ , (2 * nParameters + 1):(3 * nParameters)] <- t((pMatrix[2:(nParameters + 1), 3, 1:N]))

## save file in test directory
testDeliv <- "test/deliv"

dir.create(testDeliv)
write.csv(performance_table, file.path(testDeliv, paste0(testName, ".summary.csv")))

###############################################################################
## other calculations for quick checks

## calculate the mean runtime
meanTime <- mean((pMatrix[1:nParameters, 3, 1:N]))

## histogram of the run time distribution for first parameter (CL)
hist(pMatrix[1, 3, 1:N])

## see other R script for analysis of results.



