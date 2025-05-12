#include "grid/multigrid1D.h"
#include "spherisym.h"
#include "compressible/thermal.h"
#include "compressible/NASG.h"

#define LEVEL 12

#define Req 1e-5
#define kappaa 1.050878134564148
#define omega 1.7888921295827648e6

/* Physical properties of water */
#define pinf 101325.
#define rhoo 998.21
#define Tinf 293.15
#define mul 0.
#define st 0.

/* Physical properties of air */
#define peq  (pinf + 2*st/Req)
#define cp 1006.85
#define cv 719.18
#define rhoeq (peq/(354.8799982*(cp - cv)))

#define deltaR 0.001
#define Rbub (1. + deltaR)

#define p0 peq*pow(1./Rbub,3.*kappaa)/pinf

double rhoL = 1., rhoR = rhoeq*cube(1./Rbub)/rhoo;
double tend = 1.5;
double lambda = 256;
double tr;

uf.n[right] = neumann(0.);
q.n[right]  = neumann(0.);

uf.n[left] = 0.;

int main() {
  L0 = lambda;
  CFLac = 0.5;
  DT = HUGE [0];
  
  f.gradient = zero;

  gamma1 = 1.187;
  PI1 = 7028e5/pinf;
  b1 = 6.61e-4*rhoo;
  q1 = -1177788*rhoo/pinf;

  cv1 = 3610*rhoo*Tinf/pinf; cv2 = cv*rhoo*Tinf/pinf;
  cp1 = 4285*rhoo*Tinf/pinf; cp2 = cp*rhoo*Tinf/pinf;

  kappa1 = 0.59846028987077/(Req/Tinf*sqrt(cube(pinf)/rhoo));
  kappa2 = 25.685e-3/(Req/Tinf*sqrt(cube(pinf)/rhoo));

  tr = 2.*M_PI/(omega*Req*sqrt(rhoo/pinf));
  
  tend *= tr;

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
    
      double pL = (1. - Rbub/x) + p0*Rbub/x;
    
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

event centroid (t += 0.001*2.*M_PI/(omega*Req*sqrt(rhoo/pinf))) {
  scalar ff[];
  double radius = 0.;
  foreach(reduction(+:radius)) {
    ff[] = 1. - f[];
    radius += ff[]*Delta;
  }

  if (pid() == 0) {
    char name[80];
    sprintf(name,"radius.txt");
    FILE * fp = fopen(name,"a");
    char str[80];
    sprintf(str,"%g %0.16f\n",t/tr,radius);
    fputs(str,fp);
    fclose(fp);
  }
}

event logfile (t += 0.01*2.*M_PI/(omega*Req*sqrt(rhoo/pinf))) {
  stats sp = statsf (p);
  stats su = statsf (q.x);
  stats sT = statsf (T);
  fprintf (stderr,"t = %g, i = %d, dt = %g, min(p) = %g, max(p) = %g, min(T) = %g, max(T) = %g, min(u) = %g, max(u) = %g\n", t/tr, i, dt/tr, sp.min, sp.max, sT.min, sT.max, su.min, su.max);
}

event ending (t = tend) {
  return 1.;
}

/**
~~~gnuplot Bubble radius as a function of time
set xlabel 't'
set ylabel 'R'
plot "radius.txt" u 1:2 w p
~~~
*/
