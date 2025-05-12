/**
#Finite difference (non-conservative) WENO5.

Here is a finite difference implementation of a WENO5. Why? Because in my
level-set method, my velocity is non-divergence-free and the level-set function
and the velocity are colocated.

Taken from [F. Couderc PhD Thesis](#https://tel.archives-ouvertes.fr/tel-00143709), pp56-58.


*/
foreach_dimension()
static double WENO5_x(Point point, scalar q, int i){
  static const double coeff[3][3]= {
    {1./3., -7./6., 11./6.},
    {-1/6., 5./6., 1./3.},
    {1./3., 5./6., -1./6.}
  };
  static const double wIS = 13./12;
  static const double weights[3] = {0.1,0.6,0.3};
  static const double eps = 1.e-12;

  double ENO3[3], alpha[3], IS[3];
  ENO3[0] =  q[i*2]*coeff[0][0] + q[i]*coeff[0][1]  + q[]*coeff[0][2]     ;
  ENO3[1] =  q[i]*coeff[1][0]   + q[]*coeff[2][1]   + q[-i]*coeff[1][2]  ;
  ENO3[2] =  q[]*coeff[2][0]    + q[-i]*coeff[2][1] + q[-2*i]*coeff[2][2];

  IS[0] = wIS*sq(q[2*i]- 2*q[i]  + q[]    ) + sq(q[2*i]-4*q[i]  + 3*q[]  )*0.25;
  IS[1] = wIS*sq(q[i]  - 2*q[]   + q[-i]  ) + sq(q[i]           - q[-i]  )*0.25;
  IS[2] = wIS*sq(q[]   - 2*q[-i] + q[-2*i]) + sq(3*q[] -4*q[-i] + q[-2*i])*0.25;

  alpha[0] = weights[0]/sq(eps + IS[0]);
  alpha[1] = weights[1]/sq(eps + IS[1]);
  alpha[2] = weights[2]/sq(eps + IS[2]);

  double sum = alpha[0]+ alpha[1] + alpha[2];

  return (alpha[0]*ENO3[0] + alpha[1]*ENO3[1] + alpha[2]*ENO3[2])/sum;
}

foreach_dimension()
static double WENO3_x (Point point, scalar f, int dir)
{
  static double epsilon = 1.e-12;
  static double coeff[2][2] = {
    {-1./2., 3./2.},
    {1./2., 1./2.}
  };    
  static double weights[2] = {1./3.,2./3.};

  double p[2], beta[2], ar[2];

  beta[0] = sq(f[] - f[-dir]);
  beta[1] = sq(f[dir]  - f[]);
  
  p[0] = coeff[0][0]*f[-dir] + coeff[0][1]*f[];
  p[1] = coeff[1][0]*f[]     + coeff[1][1]*f[dir];
  
  double arsum = 0.;
  for (int i = 0; i < 2; i++) {
    ar[i] = weights[i]/sq(epsilon + beta[i]);
    arsum += ar[i];
  }
  
  double fweno = 0.;
  for (int i = 0; i < 2; i++)
    fweno += ar[i]/arsum*p[i];    

  return fweno;
}

static void FE_WENO5(scalar phi, scalar fi, vector u, double dt, double maxd){
  vector upfluxp[], upfluxm[];
  foreach(){
    foreach_dimension(){
      upfluxp.x[] = (fi[1] - fi[])/Delta;
      upfluxm.x[] = (fi[] - fi[-1])/Delta;
    }
  }
  boundary((scalar *){upfluxm,upfluxp});
  restriction((scalar *){upfluxm,upfluxp});

  foreach(){
    if(fabs(phi[])<maxd){
      foreach_dimension(){
        phi[] -= dt*u.x[]*
        (u.x[] > 0 ? WENO5_x(point, upfluxm.x, -1) : WENO5_x(point,upfluxp.x, 1));
      }
    }
  }
  boundary({phi});
  restriction({phi});
}

static void FE_WENO3(scalar phi, scalar fi, vector u, double dt, double maxd){
  vector upfluxp[], upfluxm[];
  foreach(){
    foreach_dimension(){
      upfluxp.x[] = (fi[1] - fi[])/Delta;
      upfluxm.x[] = (fi[] - fi[-1])/Delta;
    }
  }
  boundary((scalar *){upfluxm,upfluxp});
  restriction((scalar *){upfluxm,upfluxp});

  foreach(){
    if(fabs(phi[])<maxd){
      foreach_dimension(){
        phi[] -= dt*u.x[]*
        (u.x[] > 0 ? WENO3_x(point, upfluxm.x, -1) : WENO3_x(point,upfluxp.x, 1));
      }
    }
  }
  boundary({phi});
  restriction({phi});
}

void get_weno_derivative (scalar f, vector df)
{
  foreach()
    foreach_dimension()
      df.x[] = (WENO5_x (point, f, -1) - WENO5_x (point, f, 1));
  boundary ((scalar* ){df});
}