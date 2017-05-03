## Template to simulate PKPD data
## Use to test deterministic features of model, i.e. no random variation.
## Friberg-Karlsson model, one patient
## Data is generated for version of model which does NOT use Torsten,
## i.e the event schedule is augmented in R.
rm(list = ls())
gc()

library(mrgsolve)
library(rstan)
library(dplyr)

modelName <- "fTwoCptDet"

## Simulate ME-2 plasma concentrations and ANC values
## using mrgsolve.
code <- '
$PARAM CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0,
       MTT = 125, Circ0 = 5, alpha = 3E-4, gamma = 0.17

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI PROL TRANSIT1 TRANSIT2 TRANSIT3 CIRC

$MAIN

// Reparametrization
double k10 = CL / VC;
double k12 = Q / VC;
double k21 = Q / VP;
double ktr = 4 / MTT;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;
dxdt_PROL = ktr * (PROL + Circ0) * ((1 - alpha * CENT/VC) * pow(Circ0/(CIRC + Circ0),gamma) - 1);
dxdt_TRANSIT1 = ktr * (PROL - TRANSIT1);
dxdt_TRANSIT2 = ktr * (TRANSIT1 - TRANSIT2);
dxdt_TRANSIT3 = ktr * (TRANSIT2 - TRANSIT3);
dxdt_CIRC = ktr * (TRANSIT3 - CIRC);

// $SIGMA 0.001 0.001 

$TABLE
double CP = CENT/VC;
double DV1 = CENT/VC * exp(EPS(1));
double DV2 = (CIRC + Circ0) * exp(EPS(2));

$CAPTURE CP DV1 DV2
'

mod <- mread("acum", tempdir(), code)
e1 <- ev(amt = 80 * 1000, ii = 12, addl = 14) # Create dosing events
e1 <- ev(amt = 80 * 1000)  ## TEST
# out <- mod %>% data_set(e1) %>% carry.out(dose) %>% Req(CP,DV1,DV2) %>% mrgsim(end=500)
mod %>% ev(e1) %>% mrgsim(end = 500) %>% plot # plot data

## Observation and dosing times
# doseTimes <- seq(0, 168, by = 12)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
xneut <- seq(0, 672, by = 48)
time <- sort(unique(c(xpk, xneut)))

# save data in data frame 
xdata <- 
  mod %>% 
  ev(e1) %>% 
  carry_out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "DV1, DV2", end = -1, add = time, recsort = 3) %>%
  as.data.frame

bql = 1e-9  ## set a below quantification limit
            ## mathematically convenient: avoids error with location parameter
            ## in log-normal distribution

xdata <- xdata %>%
  mutate(DV1 = ifelse((time %in% xpk & time != 0 & evid == 0) &
                        DV1 >= bql, DV1, NA),
         DV2 = ifelse(time %in% xneut & evid == 0, DV2, NA))

xdata$cmt[xdata$cmt == 0] <- 2 # switch from mrgsolve to NONMEM convention

if (TRUE) {  ## we're not using Torsten
## Augment data set with new dosing records based on the ii and addl entries
  addl <- xdata %>%
    filter(evid == 1 & addl > 0) %>%
    rowwise() %>%
    do(as.data.frame(.)[rep(1, .$addl),]) %>%
    group_by(ID) %>%
    mutate(time = time + (1:addl[1]) * ii,
           addl = 0,
           ii = 0)

  xdata <- xdata %>%
    bind_rows(addl) %>%
    arrange(ID, time, desc(evid))
}

nt <- nrow(xdata)

## Look at simulated data using plots
p1 <- ggplot(xdata %>% filter(!is.na(DV1)), aes(x = time, y = DV1))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)", y = "ME-2 plasma concentration (ng/mL)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))

p1 <- ggplot(xdata %>% filter(!is.na(DV2)), aes(x = time, y = DV2))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)",
       y = "ANC") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(DV1) & evid == 0])
nObsPK <- length(iObsPK)
## Indices of records containing observed neutrophil counts
iObsPD <- with(xdata, (1:nrow(xdata))[!is.na(DV2) & evid == 0])
nObsPD <- length(iObsPD)

## Parameters for priors (only used for informed initial estimates)
CLPrior = 10
QPrior = 15
VCPrior = 35
VPPrior = 105
kaPrior = 2
CLPriorCV = 0.10
QPriorCV = 0.18
VCPriorCV = 0.14
VPPriorCV = 0.17
kaPriorCV = 0.16

circ0Prior <- 5
circ0PriorCV <- 0.20
mttPrior <- 125
mttPriorCV <- 0.2
gammaPrior <- 0.17
gammaPriorCV <- 0.2
alphaPrior <- 2.0E-4
alphaPriorCV <- 0.2

## create data set
data <- with(xdata,
             list(
               nt = nt,
               nObsPK = nObsPK,
               iObsPK = iObsPK,
               nObsPD = nObsPD,
               iObsPD = iObsPD,
               amt = amt,
               cmt = cmt,
               evid = evid,
               time = time,
               ii = ii,
               addl = addl,
               ss = ss,
               rate = rate,
               cObs = DV1[iObsPK],
               neutObs = DV2[iObsPD],
               
               ## Priors for PK parameters
               CLPrior = CLPrior,
               CLPriorCV = CLPriorCV,
               QPrior = QPrior,
               QPriorCV = QPriorCV,
               VCPrior = VCPrior,
               VCPriorCV = VCPriorCV,
               VPPrior = VPPrior,
               VPPriorCV = VPPriorCV,
               kaPrior = kaPrior,
               kaPriorCV = kaPriorCV,
 
               ## Priors for PD parameters
               circ0Prior = circ0Prior,
               circ0PriorCV = circ0PriorCV,
               mttPrior = mttPrior,
               mttPriorCV = mttPriorCV,
               gammaPrior = gammaPrior,
               gammaPriorCV = gammaPriorCV,
               alphaPrior = alphaPrior,
               alphaPriorCV = alphaPriorCV
             ))

## create initial estimates using real values of parameters
init <- function(){
  list(CL = 10,
       Q = 15,
       VC = 35,
       VP = 105,
       ka = 2.0,
       sigma = 1e-6,
       alpha = 3e-4,
       mtt = 125,
       circ0 = 5,
       gamma = 0.17,
       sigmaNeut =1e-6)
}

with(data, stan_rdump(ls(data), file = paste0(modelName,".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName,".init.R")))

