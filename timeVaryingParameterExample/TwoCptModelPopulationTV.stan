data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  row_vector<lower = 0>[nt] weight;
  real<lower = 0> time[nt];
  vector<lower = 0>[nObs] cObs;
}

transformed data{
  real<lower = 0> rate[nt] = rep_array(0.0, nt);
  real<lower = 0> ii[nt] = rep_array(0.0, nt);
  int<lower = 0> addl[nt] = rep_array(0, nt);
  int<lower = 0> ss[nt] = rep_array(0, nt);
  // Integers required to specify dimensions
  int<lower = 1> nRandom = 5; // Number of random effects
  int<lower = 1> nCmt = 3; // Number of model compartments
  int<lower = 1> nParms = 5; // Number of parameters passed to Torsten function

  // Fixed value parameters, e.g.,
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  // Population-level model parameters
  // These are the parameters for which you specify prior distributions
  // and initial estimates, e.g., 
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;

  // Inter-Individual variability
  cholesky_factor_corr[nRandom] L;
  vector<lower = 0>[nRandom] omega;
  
  // residual variability
  real<lower = 0> sigma;
  
  // Individual-level model parameters
  //  vector[nRandom] logtheta[nId];
  matrix[nRandom, nId] eta;
}

transformed parameters{
  // Vector of PK parameter typical values -- only those with IIV
  vector<lower = 0>[nRandom] thetaHat = [CLHat, QHat, V1Hat, V2Hat, kaHat]';

  // Matrix of individual-level model parameters
  matrix<lower = 0>[nId, nRandom] theta;

  // Individual-level model parameters with recognizable names, e.g.,
  row_vector<lower = 0>[nt] CL;
  row_vector<lower = 0>[nt] Q;
  row_vector<lower = 0>[nt] V1;
  row_vector<lower = 0>[nt] V2;
  row_vector<lower = 0>[nt] ka;
  
  // Covariance matrix
  //  cov_matrix[nRandom] Omega;

  // Predicted concentrations (without residual variation)
  vector<lower = 0>[nt] cHat; // All events

  // Amounts in each compartment at each event
  matrix[nt, nCmt] x;

  // Array used to pass parameters to the Torsten function
  real<lower = 0> parms[nt, nParms];

  //  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  theta = (rep_matrix(thetaHat, nId) .* 
          exp(diag_pre_multiply(omega, L * eta)))';

  for(j in 1:nId){
    
    // Calculation of individual parameter values given theta and covariates, e.g.
    CL[start[j]:end[j]] = theta[j, 1] * exp(0.75 * log(weight[start[j]:end[j]] / 70));
    Q[start[j]:end[j]] = theta[j, 2] * exp(0.75 * log(weight[start[j]:end[j]] / 70));
    V1[start[j]:end[j]] = theta[j, 3] * weight[start[j]:end[j]] / 70;
    V2[start[j]:end[j]] = theta[j, 4] * weight[start[j]:end[j]] / 70;
    ka[start[j]:end[j]] = rep_row_vector(theta[j, 5], end[j] - start[j] + 1);

    // Pack individual PK parameters into parms array, e.g.
    
    parms[start[j]:end[j],] = to_array_2d([CL[start[j]:end[j]], 
      Q[start[j]:end[j]], 
      V1[start[j]:end[j]], 
      V2[start[j]:end[j]], 
      ka[start[j]:end[j]]]');

    x[start[j]:end[j],] = PKModelTwoCpt(time[start[j]:end[j]], 
					amt[start[j]:end[j]],
					rate[start[j]:end[j]],
					ii[start[j]:end[j]],
					evid[start[j]:end[j]],
					cmt[start[j]:end[j]],
					addl[start[j]:end[j]],
					ss[start[j]:end[j]],
					parms[start[j]:end[j],], F, tLag);

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] ./ V1[j];
  }
}

model{
  // informative prior
/*
  CL ~ lognormal(log(10), 0.25);
  Q ~ lognormal(log(15), 0.5);
  V1 ~ lognormal(log(35), 0.25);
  V2 ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ cauchy(0, 1);
*/

// weakly informative priors
  CLHat ~ normal(0, 50);
  QHat ~ normal(0, 50);
  V1Hat ~ normal(0, 100);
  V2Hat ~ normal(0, 500);
  kaHat ~ normal(0, 5);
  sigma ~ cauchy(0, 1);
  omega ~ cauchy(0, 1);
  
  //  rho ~ lkj_corr(1); 
  L ~ lkj_corr_cholesky(1);
  
  // Inter-individual variability
  //  logtheta ~ multi_normal(log(thetaHat), Omega);
  to_vector(eta) ~ normal(0, 1);

  cObs ~ lognormal(log(cHat[iObs]), sigma);
}

generated quantities{
  //  vector[nRandom] logthetaPred[nId];
  matrix[nRandom, nId] etaPred;
  matrix<lower = 0>[nId, nRandom] thetaPred;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nt] cHatPred;
  vector[nt] cObsCond;
  vector[nt] cObsPred;

  // Individual-level model parameters with recognizable names, e.g.,
  row_vector<lower = 0>[nt] CLPred;
  row_vector<lower = 0>[nt] QPred;
  row_vector<lower = 0>[nt] V1Pred;
  row_vector<lower = 0>[nt] V2Pred;
  row_vector<lower = 0>[nt] kaPred;

  matrix[nt, nCmt] xPred;
  real<lower = 0> parmsPred[nt, nParms];

  rho = L * L';
  for(j in 1:nId) 
    for(i in 1:nRandom)
      etaPred[i, j] = normal_rng(0, 1);

  thetaPred = (rep_matrix(thetaHat, nId) .* 
              exp(diag_pre_multiply(omega, L * etaPred)))';

  for(j in 1:nId){

    // Population predictions

    //    logthetaPred[j] = multi_normal_rng(log(thetaHat), Omega);

    // Calculation of individual parameter values given logtheta and covariates, e.g.
    CLPred[start[j]:end[j]] = thetaPred[j, 1] * exp(0.75 * log(weight[start[j]:end[j]] / 70));
    QPred[start[j]:end[j]] = thetaPred[j, 2] * exp(0.75 * log(weight[start[j]:end[j]] / 70));
    V1Pred[start[j]:end[j]] = thetaPred[j, 3] * weight[start[j]:end[j]] / 70;
    V2Pred[start[j]:end[j]] = thetaPred[j, 4] * weight[start[j]:end[j]] / 70;
    kaPred[start[j]:end[j]] = rep_row_vector(thetaPred[j, 5], end[j] - start[j] + 1);

    // Pack individual PK parameters into parms array, e.g.
    
    parmsPred[start[j]:end[j],] = to_array_2d([CLPred[start[j]:end[j]], 
      QPred[start[j]:end[j]], 
      V1Pred[start[j]:end[j]], 
      V2Pred[start[j]:end[j]], 
      kaPred[start[j]:end[j]]]');

    xPred[start[j]:end[j],] = PKModelTwoCpt(time[start[j]:end[j]], 
					    amt[start[j]:end[j]],
					    rate[start[j]:end[j]],
					    ii[start[j]:end[j]],
					    evid[start[j]:end[j]],
					    cmt[start[j]:end[j]],
					    addl[start[j]:end[j]],
					    ss[start[j]:end[j]],
					    parmsPred[start[j]:end[j],], F, tLag);

    // Calculate target concentration for specified compartment.
    // Change compartment number and distribution volume as appropriate.

    cHatPred[start[j]:end[j]] = xPred[start[j]:end[j], 2] ./ V1Pred[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObsCond[i] = 0;
      cObsPred[i] = 0;
    }else{
      cObsCond[i] = lognormal_rng(log(cHat[i]), sigma);
      cObsPred[i] = lognormal_rng(log(cHatPred[i]), sigma);
    }
  }

}
