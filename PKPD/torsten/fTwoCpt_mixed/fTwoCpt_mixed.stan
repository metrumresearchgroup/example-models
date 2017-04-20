## fOneCpt.stan: example to test mixed solver
## Does NOT use Torsten

functions {
  real[] feedbackODE (real t,
			                real[] x,
			                real[] parms,
			                real[] rate,
			                int[] idummy) {
    ## PK variables
    real VC = parms[3];
    real conc;
    real Edrug;

    ## PD variables
    real MTT = parms[6];
    real circ0 = parms[7];
    real alpha = parms[8];
    real gamma = parms[9];
    real ktr = 4 / MTT;
    real prol = x[1] + circ0;
    real transit1 = x[2] + circ0;
    real transit2 = x[3] + circ0;
    real transit3 = x[4] + circ0;
    real circ = fmax(machine_precision(), x[5] + circ0);

    vector[3] initPK;
    vector[3] predPK;
    real dxdt[5];
    
    initPK[1] = parms[10];
    initPK[2] = parms[11];
    initPK[3] = parms[12];
    predPK = fTwoCpt(t, to_vector(parms), initPK, rate);

    conc = predPK[2] / VC;
    Edrug = alpha * conc;

    dxdt[1] = ktr * prol * ((1 - Edrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[2] = ktr * (prol - transit1);
    dxdt[3] = ktr * (transit1 - transit2);
    dxdt[4] = ktr * (transit2 - transit3);
    dxdt[5] = ktr * (transit3 - circ);

    return dxdt;
	}

  real[] feedbackModel1(real t0, real[]  t, real[] init,
                        real amt, int cmt, int evid,
                        real[] parms,
                        real[] rate, int[] idummy) {
    real x[8];
    real temp[1, 5];

    if (t0 == t[1]) x = init;
    else {
      x[1:3] = to_array_1d(fTwoCpt(t[1] - t0,
                           to_vector(parms[1:5]),
                           to_vector(init[1:3]),
                           rate));
      temp = integrate_ode_rk45(feedbackODE, init[4:8], t0, t, parms,
                                rate, idummy);
      x[4:8] = to_array_1d(temp);
    }

    if (evid == 1) x[cmt] = x[cmt] + amt;

    return x;
  }
  
  matrix feedbackModel(real[] time, real[] amt, int[] cmt, int[] evid,
                       real[] parms, real[] rate, int[] idummy) {
    real init[8];
    int nt = size(time);
    matrix[nt, 8] pred;
    real t0 = time[1];
    real full_parms[12];
    
    init = rep_array(0, 8);
    
    full_parms[1:9] = parms;
    
    for (i in 1:nt) {
      full_parms[10:12] = init[1:3];
      init = feedbackModel1(time[max(1, i - 1)], time[i:i], init, amt[i], cmt[i],
                            evid[i], full_parms, rate, idummy);
      for (j in 1:8) pred[i, j] = init[j];
    }
    return pred;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  real<lower = 0> rate[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  
  ## data for priors
  real<lower = 0> circ0Prior;
  real<lower = 0> circ0PriorCV;
  real<lower = 0> mttPrior;
  real<lower = 0> mttPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaPrior;
  real<lower = 0> alphaPriorCV;
  
  real<lower = 0> CLPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> QPrior;
  real<lower = 0> QPriorCV;
  real<lower = 0> VCPrior;
  real<lower = 0> VCPriorCV;
  real<lower = 0> VPPrior;
  real<lower = 0> VPPriorCV;
  real<lower = 0> kaPrior;
  real<lower = 0> kaPriorCV;
}

transformed data{
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int nParms = 9;
  int nCmt = 8;
  int idummy[0];
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
}

transformed parameters {
  vector[nt] cHat;
  vector[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector<lower = 0>[nObsPD] neutHatObs;
  real theta[nParms];  # ODE parameters
  matrix[nt, nCmt] x;
  
  theta[1] = CL;
  theta[2] = Q;
  theta[3] = VC;
  theta[4] = VP;
  theta[5] = ka;
  theta[6] = mtt;
  theta[7] = circ0;
  theta[8] = alpha;
  theta[9] = gamma;
  
  x = feedbackModel(time, amt, cmt, evid, theta, rate, idummy);
  
  cHat = x[ , 2] / VC;
  neutHat = x[ , 8] + circ0;

  cHatObs = cHat[iObsPK];
  neutHatObs = neutHat[iObsPD];
}

model{
  ## Priors
  CL ~ lognormal(log(CLPrior), CLPriorCV);
  Q ~ lognormal(log(QPrior), QPriorCV);
  VC ~ lognormal(log(VCPrior), VCPriorCV);
  VP ~ lognormal(log(VPPrior), VPPriorCV);
  ka ~ lognormal(log(kaPrior), kaPriorCV);
  sigma ~ cauchy(0, 1);

  mtt ~ lognormal(log(mttPrior), mttPriorCV);
  circ0 ~ lognormal(log(circ0Prior), circ0PriorCV);
  alpha ~ lognormal(log(alphaPrior), alphaPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  ## observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities {
  real cObsPred[nObsPK];
  real neutPred[nObsPD];

  for (i in 1:nObsPK)
    cObsPred[i] = exp(normal_rng(logCObs[i], sigma));

  for (i in 1:nObsPD)
    neutPred[i] = exp(normal_rng(logNeutObs[i], sigma));
}
