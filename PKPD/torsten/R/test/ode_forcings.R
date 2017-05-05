## Forced ODEs (following Sebastian's scripts)

modelName <- "ode_forcings"
dataName <- "fTwoCpt"

scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- file.path(scriptDir, "test")
dataDir <- file.path(projectDir, dataName)

library(rstan)

expose_stan_functions(file.path(modelDir, paste0(modelName, ".stan")))

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


## TEST 1: What do the ODEs return when integrated?

## Only save the PD data
forced <- solve_forced_ode(ts, init, parms, x_r)
coupled <- solve_coupled_ode(ts, init, parms, x_r)[ , 4:8]

max(abs(forced - coupled))
## 2.270993e-06

## TEST 2: What do the evolution operators return for dosing event?
t0 <- 0
t <- 0
init[1] <- 0
amt <- 80000
cmt <- 1
evid <- 1
theta <- parms[1:9]
forced <- feedbackModel_forced(t0, t, init, amt, cmt, evid, theta, x_r)
coupled <- feedbackModel_coupled(t0, t, init, amt, cmt, evid, theta, x_r)

max(abs(forced - coupled))
## 0

## TEST 3: What do the evolution operators return for a non-dosing event?
t0 <- 0
t <- 50
init[1] <- 80000
amt <- 0
cmt <- 0
evid <- 0
forced <- feedbackModel_forced(t0, t, init, amt, cmt, evid, theta, x_r)
coupled <- feedbackModel_coupled(t0, t, init, amt, cmt, evid, theta, x_r)

max(abs(forced - coupled))
## 1.671295e-05
## PK: -2.228724e-07 -3.223982e-06 -1.671295e-05  
## PD: -1.502012e-08 -5.448061e-08  3.113174e-08 8.608153e-08 -1.124334e-07
##


## TEST 4: Do with the event Handler -- FAILS
data <- read_rdump(file.path(dataDir, paste0(dataName,".data.R")))

if (0) {
  length(data$evid)
  data$evid <- rep(0, 79)
  data$evid[1] <- 1
}

## The last argument of the Event Handler determines which evolution
## operator gets used.
forced <- feedbackModel(data$time, data$amt, data$cmt, data$evid,
                        theta, x_r, 1)

coupled <- feedbackModel(data$time, data$amt, data$cmt, data$evid,
                        theta, x_r, 0)

max(abs(forced - coupled))
## 0.00517639

## Compare PK graphs
plot(data$time, coupled[ , 2], type="l", main="PK")
lines(data$time, forced[ , 2], type="l", lty=2, col="red")

## Compare PD graphs
plot(data$time, coupled[ , 8], type="l", main="PD")
lines(data$time, forced[ , 8], type="l", lty=2, col="red")


