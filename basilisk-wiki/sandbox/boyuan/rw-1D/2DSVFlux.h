/** A copy from **[WENO scheme in rajarshi's sandbox](http://basilisk.fr/sandbox/rajarshi/WENO_CODES/Saint_Venant_Test/)**.
*/

/**
#WENO 2D Flux computation for Saint venant equation
*/

#define epsilon 1e-30
#define epsilon_WENO 1E-06

void kurganovM (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double cp = sqrt(G*hp), cm = sqrt(G*hm);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);
  if (a > epsilon) {
    *fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fq = (ap*(qm*um + G*sq(hm)/2.) - am*(qp*up + G*sq(hp)/2.) + 
	    ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/(a*20);
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = *fq = 0.;
}

double WENO5L(Point point, scalar *Xl, double gradL){

  scalar X = Xl[0];
  double T_Stencil[3], Beta[3], Weights[3];
  int i;
  double sum;

  double gammaL[3] = {1./10.,3./5.,3./10.};
  T_Stencil[0]  =   1.*(X[-1] - 2.*Delta*gradL)/3. - 7.*X[-2]/6. + 11.*X[-1]/6.;
  Beta[0]       =  13.*sq(X[-1] - 2.*Delta*gradL - 2.*X[-2] + X[-1])/12. + sq(X[-1] - 2.*Delta*gradL - 4.*X[-2] + 3.*X[-1])/4.;
  T_Stencil[1]  =  -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
  Beta[1]       =  13.*sq(X[-2] - 2.*X[-1] + X[])/12. + sq(X[-2]-X[])/4.;
  T_Stencil[2]  =   1.*X[-1]/3. + 5.*X[]/6. - 1.*X[1]/6.;
  Beta[2]       =  13.*sq(X[-1] - 2.*X[] + X[1])/12. + sq(3.*X[-1] - 4.*X[] + X[1])/4.;
     
  for(i=0; i<=2; i++)
    Weights[i] = gammaL[i]/sq(epsilon_WENO + Beta[i]);
  
  sum = Weights[0] + Weights[1] + Weights[2];
  for(i=0; i<=2; i++)
    Weights[i] /= sum;

  sum = 0.;
  for(i=0;i<=2;i++)
   sum += T_Stencil[i]*Weights[i]; 
 
  return(sum);
}

double WENO5R(Point point, scalar *Xl, double gradR){

 scalar X = Xl[0];
 double T_Stencil[3], Beta[3], Weights[3];
 int i;
 double sum;

 double gammaR[3] = {3./10.,3./5.,1./10.};
 T_Stencil[0] = -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
 Beta[0]      = 13.*sq(X[-2] - 2.*X[-1] + X[])/12. + sq(X[-2] - 4.*X[-1] + 3.*X[])/4.;
 T_Stencil[1] =  1.*X[-1]/3. + 5.*X[]/6. - 1.*X[1]/6.;
 Beta[1]      = 13.*sq(X[-1] - 2.*X[] + X[1])/12. + sq(X[-1]-X[1])/4.;
 T_Stencil[2] = 11.*X[]/6. - 7.*X[1]/6. + (X[0]+2.*Delta*gradR)/3.;
 Beta[2]      = 13.*sq(2.*X[] - 2.*X[1] + 2.*Delta*gradR)/12. + sq(4.*X[] - 4.*X[1] + 2.*Delta*gradR)/4.;
 for(i=0; i<=2; i++)
   Weights[i] = gammaR[i]/sq(epsilon_WENO + Beta[i]);
  
 sum = Weights[0] + Weights[1] + Weights[2];
 for(i=0; i<=2; i++)
   Weights[i] /= sum;

 sum = 0;
 for(i=0;i<=2;i++)
   sum += T_Stencil[i]*Weights[i]; 
 
 return(sum);
}


double WENO3L(Point point, scalar *Xl, double gradL){

  scalar X = Xl[0];
  double T_Stencil[2], Beta[2], Weights[2];
  int i;
  double sum;
  
  double gammaL[2] = {1./3.,2./3.};
  T_Stencil[0] = -1.*X[-2]/2. + 3.*X[-1]/2.;
  Beta[0]      = sq(X[-1]-X[-2]);
  T_Stencil[1] = 1.*X[-1]/2. + 1.*X[]/2.;
  Beta[1]      = sq(X[]-X[-1]);
  for(i=0; i<=1; i++)
   Weights[i] = gammaL[i]/sq(epsilon_WENO + Beta[i]);
  
  sum = Weights[0] + Weights[1];
  for(i=0; i<=1; i++)
   Weights[i] /= sum;

  sum = 0;
  for(i=0;i<=1;i++)
    sum += T_Stencil[i]*Weights[i]; 
 
  return(sum);
}


double WENO3R(Point point, scalar *Xl, double gradR){
 
  scalar X = Xl[0];
  double T_Stencil[2], Beta[2], Weights[2];
  int i;
  double sum;

  double gammaR[2] = {2./3.,1./3.};
  T_Stencil[0] = 1.*X[-1]/2. + 1.*X[]/2.;
  Beta[0]      = sq(X[]-X[-1]);
  T_Stencil[1] = 3.*X[]/2. - 1.*X[1]/2.;
  Beta[1]      = sq(X[1]-X[]);
  for(i=0; i<=1; i++)
   Weights[i] = gammaR[i]/sq(epsilon_WENO + Beta[i]);
  
  sum = Weights[0] + Weights[1];
  for(i=0; i<=1; i++)
   Weights[i] /= sum;

  sum = 0;
  for(i=0;i<=1;i++)
    sum += T_Stencil[i]*Weights[i]; 
 
  return(sum);
  
}


double LIMITER2L(Point point, scalar *Xl, double gradL){
  scalar X = Xl[0];
  double grad = gradient (X[-2], X[-1], X[])/Delta;
  return (X[-1] + Delta*grad/2.);
}

double LIMITER2R(Point point, scalar *Xl, double gradR){
  scalar X = Xl[0];
  double grad = gradient (X[-1], X[0], X[1])/Delta;
  return (X[] - Delta*grad/2.);
} 
