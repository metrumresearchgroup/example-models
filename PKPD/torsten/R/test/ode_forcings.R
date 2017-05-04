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
ts <- 
  c(0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083,
    12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84, 96,
    108, 120, 132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5, 168.75, 169, 169.5,
    170, 171, 172, 174, 176, 180, 186, 192, 240, 288, 336, 384, 432, 480, 528, 576, 624,
    672)

parms <- c(CL, Q, VC, VP, ka, mtt, circ0, alpha, gamma, init[1:3])
x_r <- c(0, 0)  # dummy variable

## Only save the PD data
forced <- solve_forced_ode(ts, init, parms, x_r)
coupled <- solve_coupled_ode(ts, init, parms, x_r)[ , 4:8]

diff <- max(abs(forced - coupled))
## 2.270993e-06

