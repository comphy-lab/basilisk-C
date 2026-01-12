#if dimension == 2
#define dA Delta
#define dV sq(Delta)
#endif

#if dimension == 3
#define dA sq(Delta)
#define dV pow(Delta,3)
#endif

double nusselt_top (scalar T){
  double flux = 0.;

  foreach_boundary (top, reduction(+:flux))
    flux += dA*(T[] - T[0,1])/Delta;

  return flux;
}

double nusselt_bot (scalar T){
  double flux=0.;

  foreach_boundary (bottom, reduction(+:flux))
    flux += dA*(T[0,-1] - T[])/Delta;

  return flux;
}

double nusselt_right (scalar T){
  double flux=0.;
  
  foreach_boundary (right, reduction(+:flux))
    flux += dA*(T[] - T[1,0])/Delta;
  
  return flux;
}

double nusselt_left (scalar T){
  double flux=0.;
  
  foreach_boundary (left, reduction(+:flux))
    flux += dA*(T[-1,0] - T[])/Delta;

  return flux;
}

double nusselt_vol (scalar T, vector u){
  double flux=0.;
  
  foreach(reduction(+:flux))
    flux += dV*(sqrt(Ra)*u.y[]*T[] - (T[0,1] - T[0,-1])/(2*Delta));
  
  return flux;
}

double nusselt_viscous (vector u){
  double flux=0.;  

  // Viscous dissipation  
  foreach(reduction(+:flux)){
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);    
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double Sxx = dudx;
    double Sxy = 0.5*(dudy + dvdx);
    double Syx = Sxy;
    double Syy = dvdy;
    double S2 = sq(Sxx) + sq(Sxy) + sq(Syx) + sq(Syy);
    #if dimension == 3
      double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
      double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
      double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      double Szz = dwdz;
      double Sxz = 0.5*(dwdx + dudz);
      double Syz = 0.5*(dwdy + dvdz);
      double Szx = Sxz;
      double Szy = Syz;
      S2 += sq(Szz) + sq(Sxz) + sq(Syz) + sq(Szx) + sq(Szy);
    #endif
    flux += dV * S2;
  }
  
  return flux;
}

double nusselt_tmp (scalar T){
  double flux=0.;
  
  foreach(reduction(+:flux)){
    foreach_dimension(){
      flux += dV*sq((T[1] - T[-1])/(2.*Delta));
    }
  }

  return flux;
}

double nusselt_mix (scalar T, scalar Yref){
  double flux=0.;
  
  foreach(reduction(+:flux)){
    foreach_dimension(){
      flux += dV*(T[1] - T[-1])*(Yref[1] - Yref[-1])/sq(2.*Delta);
    }
  }

  return flux;
}