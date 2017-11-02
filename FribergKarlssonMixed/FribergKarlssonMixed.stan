functions{
    real[] FK_Ode(real t,
			            real[] y,
			            real[] y_pk,  // forcing function
		            	real[] theta,
			            real[] x_r,
			            int[] x_i) {
	 // PK variables
   real VC = theta[3];

   // PD variables
   real MTT = theta[6];
   real circ0 = theta[7];
   real alpha = theta[8];
   real gamma = theta[9];
   real ktr = 4 / MTT;
   real prol = y[1] + circ0;
   real transit1 = y[2] + circ0;
   real transit2 = y[3] + circ0;
   real transit3 = y[4] + circ0;
   real circ = fmax(machine_precision(), y[5] + circ0);
   real conc = y_pk[2] / VC;
   real Edrug = alpha * conc;
   real dydt[5];
   
   conc = y_pk[2] / VC;
   Edrug = alpha * conc;
   
   dydt[1] = ktr * prol * ((1 - Edrug) * ((circ0 / circ)^gamma) - 1);
   dydt[2] = ktr * (prol - transit1);
   dydt[3] = ktr * (transit1 - transit2);
   dydt[4] = ktr * (transit2 - transit3);
   dydt[5] = ktr * (transit3 - circ);

   return dydt;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  real<lower = 0> rate[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
}

transformed data{
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int nTheta = 9;
  int nCpt = 8;
  int nOde = 5;
  real F[nCpt];
  real tlag[nCpt];
  
  for (i in 1:nCpt) {
    F[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  // PK parameters
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> sigma;

  // PD parameters
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigmaNeut;
}

transformed parameters{
  vector[nt] cHat;
  vector[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector[nObsPD] neutHatObs;
  matrix[nt, nCpt] x;
  real theta[nTheta] = {CL, Q, VC, VP, ka, mtt, circ0, alpha, gamma};

  x = mixOde2CptModel_rk45(FK_Ode, nOde,
                          time, amt, rate, ii, evid, cmt, addl, ss,
                          theta, F, tlag,
                          1e-6, 1e-6, 1e+6);
                             
  cHat = x[, 2] / VC;  // Divide by volume to get concentration
  neutHat = x[, 8] + circ0;  // Add baseline
  
  for(i in 1:nObsPK) cHatObs[i] = cHat[iObsPK[i]];
  for(i in 1:nObsPD) neutHatObs[i] = neutHat[iObsPD[i]];
}

model{
  // Priors for PK parameters (weakly informative)
  CL ~ normal(0, 20);
  Q ~ normal(0, 20);
  VC ~ normal(0, 100);
  VP ~ normal(0, 1000);
  ka ~ normal(0, 5);
  sigma ~ cauchy(0, 1);

  // Priors for PD parameters (weakly informative)
  mtt ~ lognormal(log(100), 20);
  circ0 ~ lognormal(log(5), 10);
  alpha ~ lognormal((2E-4), 1.5E-4);
  gamma ~ lognormal((0.10), 0.05);

  // observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities{
  real cObsPred[nObsPK];
  real neutObsPred[nObsPD];

  for (i in 1:nObsPK)
    cObsPred[i] = exp(normal_rng(log(cHatObs[i]), sigma));
  
  for (i in 1:nObsPD)
    neutObsPred[i] = exp(normal_rng(log(neutHatObs[i]), sigmaNeut));
}
