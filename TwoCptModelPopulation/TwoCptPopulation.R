rm(list = ls())
gc()

modelName <- "TwoCptModelPopulation"

useRStan <- TRUE

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
data$nIIV <- NULL

## create initial estimates
init <- function()
  list(CLHat = exp(rnorm(1, log(10), 0.2)),
       QHat = exp(rnorm(1, log(20), 0.2)),
       V1Hat = exp(rnorm(1, log(70), 0.2)),
       V2Hat = exp(rnorm(1, log(70), 0.2)),
       kaHat = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2),
       L = diag(5),
       eta = matrix(rep(0, 5 * data$nSubjects), nrow = 5),
       omega = runif(5, 0.5, 2))

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat", "sigma", "omega")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", "theta")

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
## posterior predictive plots

xdata <- with(data,
              data.frame(id = rep(1:length(start), end - start + 1),
                         time = time,
                         cObs = NA))
xdata$cObs[data$iObs] = data$cObs

predInd <- as.data.frame(fit, pars = "cObsCond") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lbInd = quantile(value, probs = 0.05, na.rm = TRUE),
            medianInd = quantile(value, probs = 0.5, na.rm = TRUE),
            ubInd = quantile(value, probs = 0.95, na.rm = TRUE))

predPop <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lbPop = quantile(value, probs = 0.05, na.rm = TRUE),
            medianPop = quantile(value, probs = 0.5, na.rm = TRUE),
            ubPop = quantile(value, probs = 0.95, na.rm = TRUE))

predAll <- bind_cols(xdata, predInd, predPop)

p1 <- ggplot(predAll, aes(x = time, y = cObs))
p1 <- p1 + 
  geom_line(aes(x = time, y = medianPop, 
                color = "population")) +
  geom_ribbon(aes(ymin = lbPop, ymax = ubPop, 
                  fill = "population"), 
              alpha = 0.25) +
  geom_line(aes(x = time, y = medianInd, 
                color = "individual")) +
  geom_ribbon(aes(ymin = lbInd, ymax = ubInd, 
                  fill = "individual"), 
              alpha = 0.25) +
  scale_color_brewer(name  ="prediction",
                     breaks=c("individual", "population"),
                     palette = "Set1") +
  scale_fill_brewer(name  ="prediction",
                    breaks=c("individual", "population"),
                    palette = "Set1")
p1 + geom_point() +
  labs(x = "time (h)",
       y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = c(0.8, 0.25),
        strip.text = element_text(size = 8)) +
  facet_wrap(~ id)

dev.off()
