#include "grid/quadtree.h"
#include "axi.h"
#include "src/compressible-thermal.h"
#include "src/compressible-tension.h"

#define LEVEL 14

double rhoL = 1., rhoR = 4.972194208811462e-5;
double p0L = 1. + 1./19.;
double p0 = 1./19.;
double tend = 1.5;
double Rbub = 1.;
double lambda = 64.;
double tr;

double pdim = 101325.*19./20.;
double rhodim = 998.21;
double Tdim = 293.15;
double Rdim = 1e-4;

face vector cs[];

scalar centroid_x[], centroid_y[], frac_area[];

scalar pdump[];

uf.t[left] = dirichlet(0.);
uf.n[left] = dirichlet(0.);
q.n[left]  = dirichlet(0.);
q.t[left]  = dirichlet(0.);

uf.n[top]  = neumann(0.);
q.n[top]   = neumann(0.);

uf.n[right] = neumann(0.);
q.n[right]  = neumann(0.);

uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0.);

double CFLac = 5.;
double dtmin = HUGE;

event stability (i++) {
  foreach() {
    double cl = f[] ? sqrt(gamma1*(PI1 + p[])/rhoL/(1. - rhoL*b1)) : 0.;
    double cg = (1. - f[]) ? sqrt(gamma2*p[]/rhoR) : 0.;
    double cs = fmax(cl,cg);
    dtmin = min(Delta*CFLac/cs, dtmin);
  }
  DT = dtmin;
  dtmax = dtmin;
}

int main() {
  L0 = lambda;
  X0 = -4.*Rbub;

  tr = 0.915;
  tend *= tr;

  double Weber  = 1000.;
  f.sigma = 1./Weber;
  f.gradient = zero;

  gamma1 = 1.187;
  PI1 = 7028e5/pdim;
  b1 = 6.61e-4*rhodim;
  qq1 = -1177788*rhodim/pdim;

  cv1 = 3610*rhodim*Tdim/pdim; cv2 = 719.18*rhodim*Tdim/pdim;
  cp1 = 4285*rhodim*Tdim/pdim; cp2 = 1006.85*rhodim*Tdim/pdim;

  kappaT1 = 0.59846028987077/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));
  kappaT2 = 25.685e-3/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));

  init_grid(1 << 5);

  TOLERANCE = 1e-6;
  
  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) {
    int maxlevel = LEVEL;
    refine ( level <= (maxlevel - sqrt(sq(x) + sq(y))/4./Rbub));

    fraction (f, - (sq(Rbub) - sq(x) - sq(y)));
  
    foreach() {      
      frho1[]  = f[]*rhoL;
      frho2[]  = (1. - f[])*rhoR;

      double pL = p0L*(1. - Rbub/sqrt(sq(x) + sq(y))) + (p0 - 2.*f.sigma)*Rbub/sqrt(sq(x) + sq(y));
    
      p[] = pL*f[] + p0*(1. - f[]);

      double fc = clamp (f[],0.,1.);
      double rhocpmcvavg = (cp1 - cv1)*frho1[] + (cp2 - cv2)*frho2[];
      double const1 = (fc - frho1[]*b1) + (1. - fc - frho2[]*b2);
      double const2 = (fc - frho1[]*b1)*PI1 + (1. - fc - frho2[]*b2)*PI2;
      T[] = (const1*p[] + const2)/rhocpmcvavg;
    
      fE1[]   = (pL + gamma1*PI1)/(gamma1 - 1.)*(f[] - frho1[]*b1) + frho1[]*qq1;
      fE2[]   = (1. - f[])*(p0/(gamma2 - 1.));
      q.x[] = 0.;
      q.y[] = 0.;
    }
    boundary ((scalar *){q,frho1,frho2,p,fE1,fE2});
  }
}

event centroid (t += 0.0001*0.915) {
  scalar ff[];
  
  foreach() {
    ff[] = 1. - f[];
    
    double xc = x, yc = y;
    if (ff[] > 0. && ff[] < 1.) {
      coord n = facet_normal (point, ff, cs), p;
      double alpha = plane_alpha (ff[], n);
      line_center (n, alpha, ff[], &p);
      xc += p.x*Delta, yc += p.y*Delta;
    }
    centroid_x[] = xc; centroid_y[] = yc;
    frac_area[] = ff[]*Delta*Delta;

    double Ek = 0.;
    foreach_dimension()
      Ek += sq(q.x[]);

    double fc = clamp (f[],0.,1.);
    double invgammaavg = (fc - frho1[]*b1)/(gamma1 - 1.) +
      (1. - fc - frho2[]*b2)/(gamma2 - 1.);
    double PIGAMMAavg = PI1*gamma1*(fc - frho1[]*b1)/(gamma1 - 1.) + frho1[]*qq1 +
      PI2*gamma2*(1. - fc - frho2[]*b2)/(gamma2 - 1.) + frho2[]*qq2;

    pdump[] = (fE1[] + fE2[] - Ek/(frho1[] + frho2[])/2. - PIGAMMAavg)/invgammaavg;
  }
  
  double Volume = 0., area = 0., pressure = 0.;
  foreach(reduction(+:area) reduction(+:Volume) reduction(+:pressure)) {
    area += frac_area[];
    Volume += 2*M_PI*centroid_y[]*frac_area[];
    pressure += pdump[]*frac_area[];
  }

  if (pid() == 0) {
    FILE * fp = fopen("volume.txt","a");
    char str[80];
    sprintf(str,"%g %g %g\n",t/tr,Volume,pressure/area);
    fputs(str,fp);
    fclose(fp);
  }
}

event logfile (t += 0.01*0.915) {
  stats sp = statsf (p);
  stats su = statsf (q.x);
  stats sT = statsf (T);
  fprintf (stderr,"t = %g, i = %d, dt = %g, min(p) = %g, max(p) = %g, min(T) = %g, max(T) = %g, min(u) = %g, max(u) = %g\n", t/tr, i, dt/tr, sp.min, sp.max, sT.min, sT.max, su.min, su.max);
}

event output (t += 0.01*0.915) {
  char name[80];
  sprintf (name,"dump-%g",t/0.915);
  dump (name, list = (scalar *){f,pdump,T});
}

event ending (t = tend) {
  return 1.;
}
