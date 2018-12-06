rm(list = ls())
gc()

modelName <- "TwoCptModel"

useRStan <- FALSE

## Adjust directories to your settings.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- scriptDir
outDir <- file.path(modelDir, "results")
figDir <- outDir
tabDir <- outDir
stanDir <- file.path(dirname(projectDir), "cmdstan")
toolsDir <- file.path(projectDir, "tools")

dir.create(outDir)

.libPaths(file.path(projectDir, "lib"))
library(rstan)
library(ggplot2)
library(bayesplot)
library(dplyr)
library(tidyr)
library(parallel)

source(file.path(toolsDir, "stanTools.R"))
if(!useRStan) source(file.path(toolsDir, "cmdStanTools.R"))

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

################################################################################################
## run Stan
nChains <- 4
nPost <- 250 ## Number of post-burn-in samples per chain after thinning
nBurn <- 250 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

if(useRStan){
  fit <- stan(file = file.path(modelDir, paste0(modelName, ".stan")),
              data = data,
              pars = parameters,
              iter = nIter,
              warmup = nBurnin,
              thin = nThin,
              init = init,
              chains = nChains,
              cores = min(nChains, parallel::detectCores()))
}else{
  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()
  
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")))
  compileModel(model = file.path(outDir, modelName), stanDir = stanDir)
  
  chains <- 1:nChains
  mclapply(chains,
           function(chain, model, data, iter, warmup, thin, init) {
             outDir <- file.path(outDir, chain)
             dir.create(outDir)
             inits <- init()
             with(inits, stan_rdump(ls(inits), file = file.path(outDir,
                                                                "init.R")))
             runModel(model = model, data = data,
                      iter = iter, warmup = warmup, thin = thin,
                      init = file.path(outDir, "init.R"),
                      seed = sample(1:999999, 1),
                      chain = chain)
           },
           model = file.path(outDir, modelName),
           data = file.path(modelDir, paste0(modelName, ".data.R")),
           init = init,
           iter = nIter, warmup = nBurnin, thin = nThin,
           mc.cores = min(nChains, detectCores()))
  
  fit <- read_stan_csv(file.path(outDir, paste0(modelName, chains, ".csv")))
}

save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

bayesplot_theme_set(theme_gray())
bayesplot_theme_update(text = element_text(size = 12, family = "sans"),
                       axis.text = element_text(size = 12))
color_scheme_set(scheme = "brightblue")

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text()

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text()

mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
mcmcDensity(fit, parametersToPlot)
pairs(fit, pars = parametersToPlot)

ptable <- parameterTable(fit, parametersToPlot)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive plot

data <- with(data,
             data.frame(cObs = cObs, time = time[-1]))

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

dev.off()

