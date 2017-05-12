## Forced ODEs (following Sebastian's scripts)

modelName <- "fTwoCpt_mixed"

scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- file.path(projectDir, modelName)
dataDir <- file.path(projectDir, dataName)
toolsDir <- file.path("tools")

library(rstan)

expose_stan_functions("ode_forcings.stan")

## True parameters
CL = 10
Q = 15
VC = 35
VP = 105
ka = 2.0
alpha = 3e-4
mtt = 125
circ0 = 5
gamma = 0.17

init <- rep(0, 8)
init[1] <- 80000
# ts <- seq(1E-3, 48, length=101)
# time <- c(0, 0, 0.083, 0.167, 0.25, 0.5, 0.75, 1)
ts <- c(0.083)
parms <- c(CL, Q, VC, VP, ka, alpha, mtt, circ0, gamma, init[1:3])
x_r <- c(0, 0)  # dummy variable

coupled <- solve_forced_ode(ts, init, parms, x_r)


