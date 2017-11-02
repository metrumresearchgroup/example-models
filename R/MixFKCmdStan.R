rm(list = ls())
gc()

modelName <- "FribergKarlssonMixed"

## Adjust directories to your settings.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path("deliv", "figure", modelName)
tabDir <- file.path("deliv", "table", modelName)
modelDir <- file.path(projectDir, modelName)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path("tools")
stanDir <- file.path("cmdstan")
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
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma",
                      "mtt", "circ0", "alpha", "gamma", "sigmaNeut")

## Additional variables to monitor
otherRVs <- c("cObsPred", "neutObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

## For deterministic test: real parameter values
inits_true <- list(CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0,
                   MTT = 125, circ0 = 5, alpha = 3e-4, gamma = 0.17)

## Randomly generate initial estimates
init <- function() 
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2),
       circ0 = exp(rnorm(1, log(5), 0.2)),
       mtt = exp(rnorm(1, log(125), 0.2)),
       gamma = exp(rnorm(1, log(0.17), 0.2)),
       alpha = exp(rnorm(1, log(2.0E-4), 0.2)))


################################################################################################
## run Stan
nChains <- 1 # 4
nPost <- 1 # 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 0 # 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- nPost * nThin
nBurnin <- nBurn * nThin

chains <- 1:nChains

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

dir.create(file.path(modelDir, modelName, "temp"))

## Run model with fixed parameters (for testing purposes)
# runModelFixed(model = file.path(modelDir, modelName),
#               data = file.path(modelDir, paste0(modelName, ".data.R")),
#               init = file.path(modelDir, paste0(modelName, ".init.true.R")),
#               iter = 1, warmup = 0, thin = 1,
#               refresh = 1, seed = sample(1:99999, 1), chain = 1)

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

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters
dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
mcmcDensity(fit, parametersToPlot)
pairs(fit, pars = parametersToPlot)

ptable <- parameterTable(fit, parametersToPlot)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive plot for PK observations
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))
data <- data.frame(data$cObs, data$time[data$iObsPK],
                   data$neutObs, data$time[data$iObsPD])
data <- plyr::rename(data, c("data.cObs" = "cObs", 
                             "data.time.data.iObsPK." = "time"))

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

## posterior predictive plot for PD observations
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))
data <- data.frame(data$neutObs, data$time[data$iObsPD])
data <- plyr::rename(data, c("data.neutObs" = "neutObs",
                             "data.time.data.iObsPD." = "time"))

pred <- as.data.frame(fit, pars = "neutObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data)

p1 <- ggplot(pred, aes(x = time, y = neutObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "Absolute Neutrophil Count") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
