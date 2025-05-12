#include "grid/multigrid1D.h"
#include "spherisym.h"
#include "compressible/thermal.h"
#include "compressible/NASG.h"

#define LEVEL 16

#define PR 50

#define pdim 1e4
#define rhodim 998.21
#define Tdim 293.15
#define Rdim 1e-4
#define cvv1 3610*rhodim*Tdim/pdim
#define cpp1 4285*rhodim*Tdim/pdim
#define cvv2 719.18*rhodim*Tdim/pdim
#define cpp2 1006.85*rhodim*Tdim/pdim
#define PII1 7028e5/pdim
#define bb1 6.61e-4*rhodim
#define qqq1 -1177788*rhodim/pdim
#define cst1 (1. - bb1)
#define cst2 cst1*PII1
#define rcp (cpp1 - cvv1)
#define pin 1
#define Tin (cst1*pin + cst2)/rcp

double rhoL = 1., rhoR = pin/(Tin*(cpp2 - cvv2));
double p0L = PR;
double p0 = pin;
double tend = 1.5;
double Rbub = 1.;
double lambda = 64.;
double tr;

q.n[right]  = neumann(0.);

int main() 
{
  L0 = lambda;
  CFLac = 0.5;
  DT = HUGE [0];
  
  tr = 0.91468*sqrt(1./(p0L - p0));
  tend *= tr;

  f.gradient = zero;

  gamma1 = 1.187;
  PI1 = 7028e5/pdim;
  b1 = 6.61e-4*rhodim;
  q1 = -1177788*rhodim/pdim;

  cv1 = 3610*rhodim*Tdim/pdim; cv2 = 719.18*rhodim*Tdim/pdim;
  cp1 = 4285*rhodim*Tdim/pdim; cp2 = 1006.85*rhodim*Tdim/pdim;

  /**
For adiabatic simulations, comment out the following definitions of the thermal
conductivities.
  */

  /* kappaT1 = 0.59846028987077/(Rdim/Tdim*sqrt(cube(pdim)/rhodim)); */
  /* kappaT2 = 25.685e-3/(Rdim/Tdim*sqrt(cube(pdim)/rhodim)); */

  init_grid(1 << LEVEL);

  TOLERANCE = 1e-6;
  
  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) {
  
    foreach() {
      f[] = (Rbub - x) > Delta/2. ? 0. :
	fabs(Rbub - x) < Delta/2. ? 1. - 0.5 - (Rbub - x)/Delta :
	(Rbub - x) == Delta/2. ? 0. :
	(Rbub - x) == -Delta/2. ? 1. :
	1.;
      
      frho1[]  = f[]*rhoL;
      frho2[]  = (1. - f[])*rhoR;
    
      double pL = p0L*(1. - Rbub/x) + p0*Rbub/x;
    
      p[] = pL*f[] + p0*(1. - f[]);

      double fc = clamp (f[],0.,1.);
      double rhocpmcvavg = (cp1 - cv1)*frho1[] + (cp2 - cv2)*frho2[];
      double const1 = (fc - frho1[]*b1) + (1. - fc - frho2[]*b2);
      double const2 = (fc - frho1[]*b1)*PI1 + (1. - fc - frho2[]*b2)*PI2;
      T[] = (const1*p[] + const2)/rhocpmcvavg;
    
      fE1[]   = (pL + gamma1*PI1)/(gamma1 - 1.)*(f[] - frho1[]*b1) + frho1[]*q1;
      fE2[]   = (1. - f[])*(p0/(gamma2 - 1.));
    }
  }
}

event centroid (t += 0.0001*0.91468*sqrt(1./(PR - 1.))) {
  scalar ff[];
  double radius = 0.;
  double pressure = 0., area = 0.;
  
  foreach(reduction(+:radius) reduction(+:pressure) reduction(+:area)) {
    ff[] = 1. - f[];
    radius += ff[]*Delta;
    pressure += frho2[]*(cp2 - cv2)*T[]*Delta;
    area += ff[]*Delta;
  }

  if (pid() == 0) {
    char name[80];
    sprintf(name,"radius.txt");
    FILE * fp = fopen(name,"a");
    char str[80];
    sprintf(str,"%g %0.16f %g\n",t/tr,radius,pressure/area);
    fputs(str,fp);
    fclose(fp);
  }
}

event logfile (t += 0.01*0.91468*sqrt(1./(PR - 1.))) {
  stats sp = statsf (p);
  stats su = statsf (q.x);
  stats sT = statsf (T);
  fprintf (stderr,"t = %g, i = %d, dt = %g, min(p) = %g, max(p) = %g, min(T) = %g, max(T) = %g, min(u) = %g, max(u) = %g\n", t/tr, i, dt/tr, sp.min, sp.max, sT.min, sT.max, su.min, su.max);
}

event ending (t = tend) {
  return 1.;
}
