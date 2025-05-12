/**
# Saramito elasto-viscoplastic model

Drop coalescence. Propeties adapted for AXIS, VOF and Bingham model. This code is for a 3 phase flow.

*/
#include "axi.h"
#include "vof.h"
#include "tension.h"
#include "log-conform.h"

double tau_y;      // Yield stress
double K;          // Consistency index
double epsilon;    // Regularization parameter
double lamb_c;     // Characteristic relaxation time
double mupp;       // Characteristic polymeric viscosity
double rho1, rho2, mu1, mu2;
double fi;         // Flow index for the Herschel–Bulkley model

scalar rhov[];
face vector alphav[], muv[];
scalar f1[], f2[], * interfaces = {f1, f2};
scalar tau_norm[], deviatoric[], yielded[], mupv[], lambv[];



event defaults (i = 0)
{
  rho = rhov;
  alpha = alphav;
  mu = muv;

  lambda = lambv;   // Relaxation time
  mup = mupv;      // Polymeric viscosity
}




 /* Calculating $\mu_p$ and $\lambda$ in the event properties */
event properties (i++)
{
#if 1
  double ff = 0.01;
  foreach()
  {
    tau_norm[] = sqrt(0.5)*sqrt( 2*sq(tau_p.x.y[]) + sq((2*tau_p.x.x[]-tau_p.y.y[]-tau_qq[])/3) + sq((2*tau_p.y.y[]-tau_p.x.x[]-tau_qq[])/3)+ sq((2*tau_qq[]-tau_p.y.y[]-tau_p.x.x[])/3) );

    if ((f1[] + f2[]) > (1. - ff))
      mupv[] = 0.;
    else if (tau_norm[] <= tau_y && (f1[] + f2[]) < ff)
      mupv[] = tau_y / epsilon;  // Bingham fluid
    else if (tau_norm[] > tau_y && (f1[] + f2[]) < ff)
      mupv[] = K * tau_norm[]/(tau_norm[] - tau_y + K * epsilon);     // Bingham model
//    	mupv[] = pow( (K * pow(tau_norm[], fi)/(tau_norm[] - tau_y + K* pow(epsilon, fi))), 1/fi);  // Herschel-Bulkley model
    else if ((f1[] + f2[] < 1. - ff) && ((f1[] + f2[]) > ff)) 
    {
      int c1, c2, c3, c4;
      c1 = c2 = c3 = c4 = 1;
      if (f1[-1,0] + f2[-1,0] > 1. - ff && f1[-1,0] + f2[-1,0] < ff)
       c1 = 0;
      if (f1[1,0] + f2[1,0] > 1. - ff && f1[1,0] + f2[1,0] < ff)
       c2 = 0;
      if (f1[0,-1] + f2[0,-1] > 1. - ff && f1[0,-1] + f2[0,-1] < ff)
       c3 = 0;
      if (f1[0,1] + f2[0,1] > 1. - ff && f1[0,1] + f2[0,1] < ff)
       c4 = 0;

    if (c1 == 0 && c2 == 0 && c3 == 0 && c4 ==0)
    {
      mupv [] = 0.;
    }
    else
      mupv[] = max((mupv[-1,0]*c1 + mupv[1,0]*c2 + mupv[0,-1]*c3 + mupv[0,1]*c4) * (f1[]+f2[]) / (c1 + c2 + c3 + c4), 0.);
    }

  //  lambv[] = max(lamb_c*mupv[]/mupp, 0.);

  //  rhov[] = cm[]*((1 - f1[] - f2[])*rho1 + f1[]*rho2 + f2[]*rho2);
  }
  boundary ({tau_norm, mupv, lambv, rhov});
  #endif


  scalar fa1[], fa2[];
  foreach()
  {
    fa1[] =  (4.*f1[] + 2.*(f1[-1,0] + f1[1,0] + f1[0,-1] + f1[0,1]) + f1[1,1] + f1[-1,1] + f1[1,-1] + f1[-1,-1])/16.;
    fa2[] =  (4.*f2[] + 2.*(f2[-1,0] + f2[1,0] + f2[0,-1] + f2[0,1]) + f2[1,1] + f2[-1,1] + f2[1,-1] + f2[-1,-1])/16.;
  }
  boundary ({fa1, fa2});
  
  foreach_face()
  {
    double fm1 = (fa1[] + fa1[-1])/2.;
    double fm2 = (fa2[] + fa2[-1])/2.;
    
//    muv.x[] = fm.x[] / ((1. - fm1 - fm2)/mu1 + fm1/mu2 + fm2/mu2);
    muv.x[] = fm.x[] * ((1. - fm1 - fm2)*mu1 + fm1*mu2 + fm2*mu2);
    alphav.x[] = fm.x[]/ ((1 - fm1 - fm2)*rho1 + fm1*rho2 + fm2*rho2);
  }
  boundary ({muv,alphav});

  foreach()
  {
    rhov[] = cm[]*((1 - f1[] - f2[])*rho1 + f1[]*rho2 + f2[]*rho2);
    lambv[] = max(lamb_c*mupv[]/mupp, 0.)*(1 - fa1[] - fa2[]);
    mupv[] = (1 - fa1[] - fa2[])*mupv[];
  //  lambv[] = lamb_c*(1 - fa1[] - fa2[]);
  }
  boundary ({mupv, lambv, rhov});


 foreach()
    yielded[] = (f1[] + f2[] < 0.01) ? (tau_norm[] < tau_y ? 0 : 1) : 1;
}
