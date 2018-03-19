demo_twocpt_cmd <- function(model_path, stan_path) {

modelName <- "TwoCptModel"

scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- file.path(model_path)
toolsDir <- file.path("tools")
stanDir <- file.path(stan_path)
tempDir <- file.path(modelDir, modelName, "temp")

# .libPaths(...)
library(rstan)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(parallel)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "cmdStanTools.R"))

rstan_options(auto_write = TRUE)

set.seed(11191951) ## not required but assures repeatable results

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

## Randomly generate initial estimates
init <- function()
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2))

# dir.create(tempDir)  ## directory to store initial estimates for each chain

################################################################################################
## run Stan
nChains <- 4 # 4
nPost <- 1000 # 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 # 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- nPost * nThin
nBurnin <- nBurn * nThin

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

dir.create(file.path(modelDir, modelName, "temp"))

chains <- 1:nChains
mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init) {
           tempDir <- file.path(tempDir, chain)
           dir.create(tempDir)
           inits <- init()
           with(inits, stan_rdump(ls(inits), file = file.path(tempDir,
                                                              "init.R")))
           runModel(model = model, data = data,
                    iter = iter, warmup = warmup, thin = thin,
                    init = file.path(tempDir, "init.R"),
                    seed = sample(1:999999, 1),
                    chain = chain, refresh = 100,
                               adapt_delta = 0.95, stepsize = 0.01)
           },
         model = file.path(modelDir, modelName),
         data = file.path(modelDir, paste0(modelName, ".data.R")),
         init = init,
         iter = nIter, warmup = nBurnin, thin = nThin,
         mc.cores = min(nChains, detectCores()))

fit <- read_stan_csv(file.path(modelDir, modelName, paste0(modelName, chains, ".csv")))

return(fit)

}

fit <- demo_twocpt_cmd("/Users/yiz/Work/Torsten/example-models/TwoCptModel", "~/Work/Torsten/cmdstan")
data <- read_rdump(file.path("../TwoCptModel", paste0("TwoCptModel",".data.R")))
outputmcmc(fit, data, "deliv/TwoCptModel", c("CL", "Q", "lp__", "V1", "V2", "ka", "sigma"))
