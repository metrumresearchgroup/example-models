data{
  int<lower = 1> nId;
  int<lower = 1> nt;
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  row_vector<lower = 0>[nt] weight;
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  int<lower = 0> nRandom;
  corr_matrix[nRandom] rho;
  vector<lower = 0>[nRandom] omega;
  real<lower = 0> sigma;
}

transformed data{
  real<lower = 0> rate[nt] = rep_array(0.0, nt);
  real<lower = 0> ii[nt] = rep_array(0.0, nt);
  int<lower = 0> addl[nt] = rep_array(0, nt);
  int<lower = 0> ss[nt] = rep_array(0, nt);
  // Integers required to specify dimensions
  int<lower = 1> nCmt = 3; // Number of model compartments
  int<lower = 1> nParms = 5; // Number of parameters passed to Torsten function
  vector<lower = 0>[nRandom] thetaHat =
    [CLHat, QHat, V1Hat, V2Hat, kaHat]';
  cov_matrix[nRandom] Omega;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
  
  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)

}

parameters{
}

transformed parameters{
}

model{
}

generated quantities{
  vector[nRandom] logtheta[nId];
  row_vector<lower = 0>[nt] CL;
  row_vector<lower = 0>[nt] Q;
  row_vector<lower = 0>[nt] V1;
  row_vector<lower = 0>[nt] V2;
  row_vector<lower = 0>[nt] ka;
  vector[nt] cHat;
  real cObs[nt];
  matrix[nt, nCmt] x;
  real<lower = 0> parms[nt, nParms];

  for(j in 1:nId){
    logtheta[j] = multi_normal_rng(log(thetaHat), Omega);
    CL[j] = exp(logtheta[j, 1]) * (weight[j] / 70)^0.75;
    Q[j] = exp(logtheta[j, 2]) * (weight[j] / 70)^0.75;
    V1[j] = exp(logtheta[j, 3]) * weight[j] / 70;
    V2[j] = exp(logtheta[j, 4]) * weight[j] / 70;
    ka[j] = exp(logtheta[j, 5]);
    CL[start[j]:end[j]] = exp(logtheta[j, 1] + 0.75 * log(weight[start[j]:end[j]] / 70));
    Q[start[j]:end[j]] = exp(logtheta[j, 2] + 0.75 * log(weight[start[j]:end[j]] / 70));
    V1[start[j]:end[j]] = exp(logtheta[j, 3]) * weight[start[j]:end[j]] / 70;
    V2[start[j]:end[j]] = exp(logtheta[j, 4]) * weight[start[j]:end[j]] / 70;
    ka[start[j]:end[j]] = rep_row_vector(exp(logtheta[j, 5]), end[j] - start[j] + 1);

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

    cHat[start[j]:end[j]] = x[start[j]:end[j], 2] / V1[j];
  }

  for(i in 1:nt){
    if(time[i] == 0){
      cObs[i] = 0;
    }else{
      cObs[i] = lognormal_rng(log(fmax(machine_precision(), cHat[i])),
			       sigma);
    }
  }

}

