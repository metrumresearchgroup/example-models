rm(list = ls())
gc()

# modelName <- "fTwoCpt"
modelName <- "fTwoCpt_mixed"

## Adjust directories to your settings.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path("deliv", "figure", modelName)
tabDir <- file.path("deliv", "table", modelName)
modelDir <- file.path(projectDir, modelName)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path("tools")
stanDir <- file.path("cmdstan")

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
parametersToPlot <- c("CL", "Q", "VC", "VP", "ka", "sigma", "mtt", "circ0", "alpha",
                      "gamma", "sigmaNeut")

## Additional variables to monitor
otherRVs <- c("cObsPred", "neutPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

################################################################################################
## Deterministic Tests
## Make sure the model produces the data, when it uses the true parameters and the
## random variations are set to 0.

nChains <- 1
chains <- 1:nChains

## Having issues using this on Metworx.
runModelFixed(model = file.path(modelDir, modelName),
              data = file.path(modelDir, paste0(modelName, "Det.data.R")),
              init = file.path(modelDir, paste0(modelName, "Det.init.R")),
              iter = 1, warmup = 0, thin = 1,
              refresh = 1, seed = sample(1:99999, 1))

fit <- read_stan_csv(file.path(modelDir, modelName, paste0(modelName, chains, ".csv")))
data <- read_rdump(file.path(modelDir, paste0(modelName,"Det.data.R")))

## Plot data for comparisons
## PK
dataPK <- data.frame(data$cObs, data$time[data$iObsPK])
dataPK <- plyr::rename(dataPK, c("data.cObs" = "cObs", "data.time.data.iObsPK." = "time"))

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>% bind_cols(dataPK)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = value))

maxDiffPK <- max(abs(pred$value - pred$cObs) / pred$cObs)
## For numerical model: 4.735332e-6
## For mixed model: 4.735332e-06

## PD
dataPD <- data.frame(data$neutObs, data$time[data$iObsPD])
dataPD <- plyr::rename(dataPD, c("data.neutObs" = "neutObs", 
                                 "data.time.data.iObsPD." = "time"))

pred <- as.data.frame(fit, pars = "neutPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>% bind_cols(dataPD)

p1 <- ggplot(pred, aes(x = time, y = neutObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "Absolute Neutrophil Count") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = value))

maxDiffPD <- max(abs(pred$value - pred$neutObs) / pred$neutObs)
## For numerical model: 3.412185e-06
## For mixed model: 3.412185e-06

## The two Stan models spit out the same deviation from the mrgsolve model. I'll take
## it there are in very close (exact?) agreement with one another.

################################################################################################
## run Stan
nChains <- 1 # 4
nPost <- 1 # 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1  # 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1
chains <- 1:nChains

nIter <- nPost * nThin
nBurnin <- nBurn * nThin

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

## Run diagnose
mclapply(chains,
         function(chain, model, data, init)
           runDiagnose(model = model, data = data,
                       init = init, seed = sample(1:99999, 1),
                       chain = chain,
                       refresh = 100),
         model = file.path(modelDir, modelName),
         data = file.path(modelDir, paste0(modelName, ".data.R")),
         init = file.path(modelDir, paste0(modelName, ".init.R")),
         mc.cores = min(nChains, detectCores()))

## Run model
mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init)
           runModel(model = model, data = data,
                    iter = iter, warmup = warmup, thin = thin,
                    init = init, seed = sample(1:999999, 1),
                    chain = chain, refresh = 100,
                               adapt_delta = 0.95, stepsize = 0.01),
         model = file.path(modelDir, modelName),
         data = file.path(modelDir, paste0(modelName, ".data.R")),
         init = file.path(modelDir, paste0(modelName, ".init.R")),
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
## posterior predictive plot
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))

## PK
dataPK <- data.frame(data$cObs, data$time[data$iObsPK])
dataPK <- plyr::rename(dataPK, c("data.cObs" = "cObs", "data.time.data.iObsPK." = "time"))

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(dataPK)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)


## PD
dataPD <- data.frame(data$neutObs, data$time[data$iObsPD])
dataPD <- plyr::rename(dataPD, c("data.neutObs" = "neutObs", 
                                 "data.time.data.iObsPD." = "time"))

pred <- as.data.frame(fit, pars = "neutPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(dataPD)

p1 <- ggplot(pred, aes(x = time, y = neutObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "Absolute Neutrophil Count") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
