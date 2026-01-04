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
