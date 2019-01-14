functions{
  
    vector twoCptAll(vector phi, 
		   vector theta,
		   data real[] x_r, 
		   int[] x_i){
    int nt = x_i[1];
    int evid[nt] = x_i[(1+1):(1+nt)];
    int cmt[nt] = x_i[(1+nt+1):(1+nt+nt)];
    int addl[nt] = x_i[(1+2*nt+1):(1+2*nt+nt)];
    int ss[nt] = x_i[(1+3*nt+1):(1+3*nt+nt)];
    real time[nt] = x_r[1:nt];
    real amt[nt] = x_r[(nt+1):(nt+nt)];
    real rate[nt] = x_r[(2*nt+1):(2*nt+nt)];
    real ii[nt] = x_r[(3*nt+1):(3*nt+nt)];
    matrix[nt, 3] x;
    real F[3] = rep_array(1.0, 3);
    real tLag[3] = rep_array(0.0, 3);

    x = PKModelTwoCpt(time, 
					amt,
					rate,
					ii,
					evid,
					cmt,
					addl,
					ss,
					to_array_1d(theta), F, tLag);

/*
    x = PKModelTwoCpt(x_r[1:nt], 
					x_r[(nt+1):(nt+nt)],
					x_r[(2*nt+1):(2*nt+nt)],
					x_r[(3*nt+1):(3*nt+nt)],
					evid,
					cmt,
					addl,
					ss,
					to_array_1d(theta), 
					x_r[(4*nt + 1):(4*nt + 3)], 
					x_r[(4*nt + 4):(4*nt + 6)]);
*/

    return to_vector(x[,2] ./ theta[3]);
  }

}

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
  int<lower = 1> nti[nId];
  real<lower = 0> weight[nId];
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
  
  int x_i[nId, 4 * max(nti) + 1];
  real x_r[nId, 4 * max(nti) + 6];
  
  for(i in 1:nId){
    x_i[i,1] = nti[i];
    x_i[i, (1+1):(1+nti[i])] = evid[start[i]:end[i]];
    x_i[i, (1+nti[i]+1):(1+nti[i]+nti[i])] = cmt[start[i]:end[i]];
    x_i[i, (1+2*nti[i]+1):(1+2*nti[i]+nti[i])] = addl[start[i]:end[i]];
    x_i[i, (1+3*nti[i]+1):(1+3*nti[i]+nti[i])] = ss[start[i]:end[i]];
    x_r[i, 1:nti[i]] = time[start[i]:end[i]];
    x_r[i, (nti[i]+1):(nti[i]+nti[i])] = amt[start[i]:end[i]];
    x_r[i, (2*nti[i]+1):(2*nti[i]+nti[i])] = rate[start[i]:end[i]];
    x_r[i, (3*nti[i]+1):(3*nti[i]+nti[i])] = ii[start[i]:end[i]];
    x_r[i, (4*nti[i] + 1):(4*nti[i] + 3)] = F;
    x_r[i, (4*nti[i] + 4):(4*nti[i] + 6)] = tLag;
  }
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

  // Individual-level model parameters adjusted for covariates
  vector<lower = 0>[nRandom] thetaAdj[nId];
  
  // Individual-level model parameters with recognizable names, e.g.,
  real<lower = 0> CL[nId];
  real<lower = 0> Q[nId];
  real<lower = 0> V1[nId];
  real<lower = 0> V2[nId];
  real<lower = 0> ka[nId];
  
  // Covariance matrix
  //  cov_matrix[nRandom] Omega;

  // Predicted concentrations (without residual variation)
  vector<lower = 0>[nt] cHat; // All events

  // Amounts in each compartment at each event
  matrix[nt, nCmt] x;

  // Array used to pass parameters to the Torsten function
  //real<lower = 0> parms[nParms];

  //  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  theta = (rep_matrix(thetaHat, nId) .* 
          exp(diag_pre_multiply(omega, L * eta)))';

  for(j in 1:nId){
    
    // Calculation of individual parameter values given logtheta and covariates, e.g.
    CL[j] = theta[j, 1] * (weight[j] / 70)^0.75;
    Q[j] = theta[j, 2] * (weight[j] / 70)^0.75;
    V1[j] = theta[j, 3] * weight[j] / 70;
    V2[j] = theta[j, 4] * weight[j] / 70;
    ka[j] = theta[j, 5];

    thetaAdj[j] = [CL[j], Q[j], V1[j], V2[j], ka[j]]';
  }

  // Calculate target concentration for specified compartment.
  // Change compartment number and distribution volume as appropriate.

  cHat = map_rect(twoCptAll, rep_vector(0.0, 1), 
    thetaAdj, x_r, x_i);

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
  real<lower = 0> CLPred[nId];
  real<lower = 0> QPred[nId];
  real<lower = 0> V1Pred[nId];
  real<lower = 0> V2Pred[nId];
  real<lower = 0> kaPred[nId];

  matrix[nt, nCmt] xPred;
  real<lower = 0> parmsPred[nParms];

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
    CLPred[j] = thetaPred[j, 1] * (weight[j] / 70)^0.75;
    QPred[j] = thetaPred[j, 2] * (weight[j] / 70)^0.75;
    V1Pred[j] = thetaPred[j, 3] * weight[j] / 70;
    V2Pred[j] = thetaPred[j, 4] * weight[j] / 70;
    kaPred[j] = thetaPred[j, 5];

    // Pack individual PK parameters into parms array, e.g.
    parmsPred = {CLPred[j], QPred[j], V1Pred[j], V2Pred[j], kaPred[j]};

    xPred[start[j]:end[j],] = PKModelTwoCpt(time[start[j]:end[j]], 
					    amt[start[j]:end[j]],
					    rate[start[j]:end[j]],
					    ii[start[j]:end[j]],
					    evid[start[j]:end[j]],
					    cmt[start[j]:end[j]],
					    addl[start[j]:end[j]],
					    ss[start[j]:end[j]],
					    parmsPred, F, tLag);

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
