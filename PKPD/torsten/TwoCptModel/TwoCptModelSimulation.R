## Template to simulate PKPD data
rm(list = ls())
gc()

modelName <- "TwoCptModel"

library(rstan)
library(mrgsolve)  ## tested with version 0.7.6.9029

set.seed(11091989) ## not required but assures repeatable results

###############################################################################################
## Simulate data using mrgsolve

code <- '
$PARAM CL=5, Q=8, V2=20, V3=70, KA=1.2

$CMT GUT CENT PERI

$GLOBAL
#define CP (CENT/V2)

$PKMODEL ncmt = 2, depot = TRUE

$SIGMA 0.01

$MAIN
pred_CL = CL;
pred_Q = Q;
pred_VC = V2;
pred_VP = V3;
pred_KA = KA;
double DV = CP * exp(EPS(1));

$CAPTURE CP DV
'

mod <- mread("accum", tempdir(),code) %>% Req(DV) %>% update(end=480,delta=0.1)

e1 <- ev(amt = 1250, ii = 12, addl = 14) # Create dosing events
mod %>% ev(e1) %>% mrgsim(end = 250) %>% plot # plot data

## Create time at which data will be observed.
## More observations around the first two and the last dosing event.
td <- 0
t1 <- c(td + 0.083, td + 0.167, td + 0.25, td + 0.5, 
        td + 0.75, td + 1, td + 1.5, td + 2, td + 3)
t2 <- seq(4, 12, 2)
td <- 12
t3 <- c(td + 0.083, td + 0.167, td + 0.25, td + 0.5, 
           td + 0.75, td + 1, td + 1.5, td + 2, td + 3)
t4 <- c(seq(16, 20, 2), 24, seq(36, 168, 12))
td <- 168
t5 <- c(td + 0.083, td + 0.167, td + 0.25, td + 0.5, 
        td + 0.75, td + 1, td + 1.5, td + 2, td + 3)
t6 <- c(172, 174, 176, 180, 186, 192)
tall <- sort(c(t1, t2, t3, t4, t5, t6))

# save data in data frame 
SimData <- 
  mod %>% 
  ev(e1) %>% 
  carry.out(cmt,ii,addl,rate,amt,evid,ss) %>%
  mrgsim(Req="CP", end=-1, add=tall,recsort=3) %>%
  as.data.frame

SimData$cmt[SimData$cmt == 0] <- 2 ## adjust cmt (adopt NONMEM convention)
SimData <- SimData[!((SimData$evid == 0)&(SimData$DV == 0)),] ## remove observation with 0 drug concentration

################################################################################################
# Format Data for Stan using RStan
nt <- nrow(SimData)
iObs <- with(SimData, (1:nrow(SimData))[evid == 0])
nObs <- length(iObs)

## create Stan data set
data <- with(SimData,
             list(nt = nt,
                  nObs = nObs,
                  iObs = iObs,
                  time = time,
                  cObs = DV[iObs],
                  amt =  amt,
                  rate = rate,
                  cmt = cmt,
                  evid = evid,
                  ii = ii,
                  addl = addl,
                  ss = ss))

## create initial estimates
init <- function() 
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2))

with(data, stan_rdump(ls(data), file = paste0(modelName, ".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName, ".init.R")))
