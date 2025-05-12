
void thermal_conduction(scalar T, bool VAP, double rhocp, double lambda) {
  scalar f2[];
  scalar thetacorr[];
  face vector f2s[];
  face vector lambdaf[];

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    if (VAP){
      f2[] = 1. - f[];
    }else{
      f2[] = f[]; 
    }
  }

  face_fraction (f2, f2s);

  foreach_face() {
    if (VAP){
      lambdaf.x[] = lambda_V/rho_V/cp_V*f2s.x[]*fm.x[];
    }else{
      lambdaf.x[] = lambda_L/rho_L/cp_L*f2s.x[]*fm.x[];
    }
  }

  foreach() {
    thetacorr[] = cm[]*max(f2[], F_ERR);
  }
  
  scalar r_temp[];
  foreach(){
    if (VAP){
      r_temp[] = sT_V[];
    }else{
      r_temp[] = sT_L[];
    }
  }
  diffusion (T, dt, D=lambdaf, r=r_temp, theta=thetacorr);
}

void conduction_step(int i){
  bool VAP = false;
  thermal_conduction(T_L, VAP, rhocp_L, lambda_L);
  VAP = true;
  thermal_conduction(T_V, VAP, rhocp_V, lambda_V);
}


