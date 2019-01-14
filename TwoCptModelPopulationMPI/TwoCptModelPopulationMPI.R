rm(list = ls())
gc()

modelName <- "TwoCptModelPopulationMPI"
simModelName <- paste0(modelName, "Sim")

useRStan <- FALSE
nslaves = 40

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
library(Rmpi)

source(file.path(toolsDir, "stanTools.R"))
if(!useRStan) source(file.path(toolsDir, "cmdStanTools.R"))

qnorm.trunc = function(p,mean=0,sd=1,lower=-Inf,upper=Inf)
  qnorm(p*pnorm(upper,mean,sd)+(1-p)*pnorm(lower,mean,sd),mean,sd)

rnorm.trunc = function(n,mean=0,sd=1,lower=-Inf,upper=Inf)
  qnorm.trunc(runif(n),mean,sd,lower,upper)

rstan_options(auto_write = TRUE)
set.seed(11191951) ## not required but assures repeatable results
set.seed(10271998)

## Set up data structure for simulation

## Parameter values

ka = 2.0
CL = 10 # L/h
Q = 15 # L/h
V1 = 35 # L
V2 = 105 # L
sigma = 0.1
omega <- c(0.2, 0.3, 0.2, 0.3, 0.25)
rho <- diag(length(omega))

## Observation and dosing times
doseTimes <- seq(0, 168, by = 12)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
time <- sort(unique(c(xpk, doseTimes)))

nId <- 10 ## Number of individuals
weight = rnorm.trunc(nId, 70, 15, 50, 100)

## Assemble data set for simulation using Stan
obsData <- data.frame(time = time) %>%
  mutate(amt = 0,
         cmt = 1,
         evid = 0)

doseData <- data.frame(time = doseTimes) %>%
  mutate(amt = 80 * 1000, # mcg
         cmt = 1,
         evid = 1)

allData <- doseData %>%
  bind_rows(obsData) %>%
  merge(data.frame(id = 1:nId)) %>%
  arrange(id, time, desc(evid))

nt <- nrow(allData)
start <- (1:nt)[!duplicated(allData$id)]
end <- c(start[-1] - 1, nt)

dataSim <- with(allData,
                list(nId = nId,
                     nt = nt,
                     amt = amt,
                     cmt = cmt,
                     evid = evid,
                     time = time,
                     start = start,
                     end = end,
                     nti = end - start + 1,
                     weight = weight,
                     CLHat = CL,
                     QHat = Q,
                     V1Hat = V1,
                     V2Hat = V2,
                     kaHat = ka,
                     nRandom = length(omega),
                     omega = omega,
                     rho = rho,
                     sigma = sigma))

## Using Stan simulate plasma drug concentrations

sim <- stan(file = file.path(modelDir, paste(simModelName, ".stan", sep = "")),
            data = dataSim,
            algorithm = "Fixed_param",
            iter = 1,
            chains = 1)

## Assemble data set for fitting via Stan

xdata <- allData %>%
  bind_cols(as.data.frame(sim, pars = "cObs") %>%
              gather(factor_key = TRUE) %>%
              select(cObs = value))

xdata <- xdata %>%
  mutate(cObs = ifelse(time %in% xpk & time != 0 & evid == 0, cObs, NA))

head(xdata)

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$id)]
end <- c(start[-1] - 1, nt)

## Indices of records containing observed concentrations
iObs <- with(xdata, (1:nrow(xdata))[!is.na(cObs) & evid == 0])
nObs <- length(iObs)

## create Stan data set
data <- with(xdata,
             list(nId = nId,
                  nt = nt,
                  nObs = nObs,
                  iObs = iObs,
                  amt = amt,
                  cmt = cmt,
                  evid = evid,
                  time = time,
                  start = start,
                  end = end,
                  nti = end - start + 1,
                  weight = weight,
                  cObs = cObs[iObs]
             ))


## create initial estimates
init <- function()
  list(CLHat = exp(rnorm(1, log(CL), 0.2)),
       QHat = exp(rnorm(1, log(Q), 0.2)),
       V1Hat = exp(rnorm(1, log(V1), 0.2)),
       V2Hat = exp(rnorm(1, log(V2), 0.2)),
       kaHat = exp(rnorm(1, log(ka), 0.2)),
       sigma = runif(1, 0.5, 2),
       L = diag(5),
       eta = matrix(rep(0, 5 * data$nId), nrow = 5),
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
  mpi.spawn.Rslaves(nslaves = nslaves)
  fit <- stan(file = file.path(modelDir, paste0(modelName, ".stan")),
              data = data,
              pars = parameters,
              iter = nIter,
              warmup = nBurnin,
              thin = nThin,
              init = init,
              chains = nChains,
              cores = min(nChains, parallel::detectCores()))
  mpi.close.Rslaves()
}else{
  
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")))
  compileModel(model = file.path(outDir, modelName), stanDir = stanDir)
  
  mpi.spawn.Rslaves(nslaves = nslaves)
  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()
  
  chains <- 1:nChains
  startTime <- Sys.time()
  mclapply(chains,
           function(chain, model, data, iter, warmup, thin, init) {
             outDir <- file.path(outDir, chain)
             dir.create(outDir)
             with(data, stan_rdump(ls(data), file = file.path(outDir,
                                                              "data.R")))
             inits <- init()
             with(inits, stan_rdump(ls(inits), file = file.path(outDir,
                                                                "init.R")))
             runModelMPI(model = model, data = file.path(outDir, "data.R"),
                         iter = iter, warmup = warmup, thin = thin,
                         init = file.path(outDir, "init.R"),
                         seed = sample(1:999999, 1),
                         chain = chain,
                         nslaves = nslaves / nChains)
           },
           model = file.path(outDir, modelName),
           data = data,
           init = init,
           iter = nIter, warmup = nBurnin, thin = nThin,
           mc.cores = min(nChains, detectCores()))
  endTime <- Sys.time()
  elapsedTime <- endTime - startTime
  elapsedTime
  mpi.close.Rslaves()
  
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
        legend.position = c(0.8, 0.15),
        strip.text = element_text(size = 8)) +
  facet_wrap(~ id)

dev.off()
