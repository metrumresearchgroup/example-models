demo_twocpt <- function(model_path) {

## rm(list = ls())
## gc()

modelName <- "TwoCptModel"

scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- file.path(model_path)
toolsDir <- file.path("tools")

library(rstan)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(parallel)

source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
set.seed(11191951) ## not required but assures repeatable results

## read data
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))

## create initial estimates
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(15), 0.2)),
       V1 = exp(rnorm(1, log(35), 0.2)),
       V2 =  exp(rnorm(1, log(105), 0.2)),
       ka = exp(rnorm(1, log(2), 0.2)),
       sigma = runif(1, 0.5, 2))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

################################################################################################
## run Stan
nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste0(modelName, ".stan")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin,
            init = init,
            chains = nChains,
            cores = min(nChains, parallel::detectCores()))

return(fit)
}

fit <- demo_twocpt("../TwoCptModel")
data <- read_rdump(file.path("../TwoCptModel", paste0("TwoCptModel",".data.R")))
outputmcmc(fit, data, "deliv/TwoCptModel", c("CL", "Q", "lp__", "V1", "V2", "ka", "sigma"))
