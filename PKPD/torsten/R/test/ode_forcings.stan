## ODE Forcing

functions {
  ## Analytical solution to the two compartment model
  vector twoCptModel1(real dt, vector init, vector parms) {
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    real ksum = k10 + k12 + k21;
    vector[3] alpha;
    vector[3] a;
    vector[3] x = rep_vector(0.0, 3);

    alpha[1] = (ksum + sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[2] = (ksum - sqrt(ksum * ksum - 4.0 * k10 * k21))/2.0;
    alpha[3] = ka;

    if(init[1] != 0.0){
      x[1] = init[1] * exp(-alpha[3] * dt);
      a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1]) * (alpha[2] - alpha[1]));
      a[2] = ka * (k21 - alpha[2]) / ((ka - alpha[2]) * (alpha[1] - alpha[2]));
      a[3] = -(a[1] + a[2]);
      x[2] = init[1] * sum(a .* exp(-alpha * dt));
      a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[2] - alpha[1]));
      a[2] = ka * k12 / ((ka - alpha[2]) * (alpha[1] - alpha[2]));
      a[3] = -(a[1] + a[2]);
      x[3] = init[1] * sum(a .* exp(-alpha * dt));
    }
    
    if(init[2] != 0){
      a[1] = (k21 - alpha[1]) / (alpha[2] - alpha[1]);
      a[2] = (k21 - alpha[2]) / (alpha[1] - alpha[2]);
      x[2] = x[2] + init[2] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
      a[1] = k12 / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[3] = x[3] + init[2] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
    }

    if(init[3] != 0){
      a[1] = k21 / (alpha[2] - alpha[1]);
      a[2] = -a[1];
      x[2] = x[2] + init[3] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
      a[1] = (k10 + k12 - alpha[1]) / (alpha[2] - alpha[1]);
      a[2] = (k10 + k12 - alpha[2]) / (alpha[1] - alpha[2]);
      x[3] = x[3] + init[3] * sum(segment(a, 1, 2) .* exp(-segment(alpha, 1, 2) * dt));
    }
    // if(evid == 1) x[cmt] = x[cmt] + amt;
    return x;
  }


  ## ODE sytem for mixed solver
  real[] feedbackODE_forced(real t, 
                            real[] x,
                            real[] parms,
                            real[] rate, 
                            int[] idummy) {
    ## PK variables
    real VC = parms[3];

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

    vector[3] initPK = to_vector(parms[10:12]);
    vector[3] predPK;
    real conc;
    real Edrug;
    real dxdt[5];

    // predPK = fTwoCpt(t, to_vector(parms[1:5]), initPK, rate);
    predPK = twoCptModel1(t, initPK, to_vector(parms[1:5]));
    // predPK = twoCptModel2(t, initPK, to_vector(parms[1:5]));

    conc = predPK[2] / VC;
    Edrug = alpha * conc;


    dxdt[1] = ktr * prol * ((1 - Edrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[2] = ktr * (prol - transit1);
    dxdt[3] = ktr * (transit1 - transit2);
    dxdt[4] = ktr * (transit2 - transit3);
    dxdt[5] = ktr * (transit3 - circ);

    return dxdt;
  }
  
    ## ODE system for numerical solver
  real[] feedbackODE (real t,
			                real[] x,
			                real[] parms,
			                real[] rdummy,
			                int[] idummy) {
    ## PK variables
    real CL = parms[1];
    real Q = parms[2];
    real VC = parms[3];
    real VP = parms[4];
    real ka = parms[5];
    real k10 = CL / VC;
    real k12 = Q / VC;
    real k21 = Q / VP;
    real conc;
    real Edrug;

    ## PD variables
    real MTT = parms[6];
    real circ0 = parms[7];
    real alpha = parms[8];
    real gamma = parms[9];
    real ktr = 4 / MTT;
    real prol = x[4] + circ0;
    real transit1 = x[5] + circ0;
    real transit2 = x[6] + circ0;
    real transit3 = x[7] + circ0;
    real circ = fmax(machine_precision(), x[8] + circ0);

    real dxdt[8];

    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2] / VC;
    Edrug = alpha * conc;

    dxdt[4] = ktr * prol * ((1 - Edrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;
	}  

	## ODE system for numerical solver
  real[] feedbackODE_coupled (real t,
			                        real[] x,
			                        real[] parms,
			                        real[] rdummy,
			                        int[] idummy) {
    ## PK variables
    real CL = parms[1];
    real Q = parms[2];
    real VC = parms[3];
    real VP = parms[4];
    real ka = parms[5];
    real k10 = CL / VC;
    real k12 = Q / VC;
    real k21 = Q / VP;
    real conc;
    real Edrug;

    ## PD variables
    real MTT = parms[6];
    real circ0 = parms[7];
    real alpha = parms[8];
    real gamma = parms[9];
    real ktr = 4 / MTT;
    real prol = x[4] + circ0;
    real transit1 = x[5] + circ0;
    real transit2 = x[6] + circ0;
    real transit3 = x[7] + circ0;
    real circ = fmax(machine_precision(), x[8] + circ0);

    real dxdt[8];

    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2] / VC;
    Edrug = alpha * conc;

    dxdt[4] = ktr * prol * ((1 - Edrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;
	}
	
	matrix solve_coupled_ode(real[] ts, real[] init, real[] theta, real[] x_r) {
	  int x_i[0];
	  
	  return(to_matrix(integrate_ode_rk45(feedbackODE_coupled, init, 0, ts, theta[1:9], x_r, x_i)));
	}
	
	matrix solve_forced_ode(real[] ts, real[] init, real[] theta, real[] x_r) {
	  int x_i[0];

	  return(to_matrix(integrate_ode_rk45(feedbackODE_forced, init[4:8], 0, ts, theta,
	                                       x_r, x_i)));
	}
	
}
  
model { }
