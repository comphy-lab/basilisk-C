/**
# Axisymmetric expansion of a bubble
 */

#include "grid/quadtree.h"
#include "axi.h"
#include "compressible/thermal.h"
#include "compressible/NASG.h"
#include "view.h"

#define LEVEL 8

/**
The initial density of the gas is chosen in way such as the initial temperature
inside the bubble is half the far-field temperature $T_\infty$.
 */

double rhoL = 1., rhoR = 0.080470865772519;
double p0L = 1.;
double p0 = 1.;
double tend = 1.;
double Rbub = 1.;
double lambda = 8.;
double tr;

/**
The problem is rendered dimensionless using the ambient pressure, the liquid
density, the far-field temperature and the bubble initial radius. The values
employed for this simulation are respectively listed.
 */

double pdim = 5e6;
double rhodim = 975.91;
double Tdim = 350;
double Rdim = 1e-4;

scalar pdata[], tdata[];

uf.n[right] = neumann(0.);
p[right]    = dirichlet(p0L);
q.n[right]  = neumann(0.);

uf.n[left] = neumann(0.);
p[left]    = dirichlet(p0L);
q.n[left]  = neumann(0.);

uf.n[top] = neumann(0.);
p[top]    = dirichlet(p0L);
q.n[top]  = neumann(0.);

uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0.);

/**
Although the thermal solver is implicit and unconditionally stable, a diffusive
CFL condition is employed for better accuracy.
 */

event stability (i++) {
  dtmax = rhoR*cp2*sq(L0/pow(2,LEVEL))/kappa2/2.;
}

int main() {
  L0 = lambda;
  X0 = -L0/2.;
  
  f.gradient = zero;

  /**
Liquid water parameters in the Noble-Abel Stiffened Gas (NASG) equation of
state.
   */

  gamma1 = 1.187;
  PI1 = 7028e5/pdim;
  b1 = 6.61e-4*rhodim;
  q1 = -1177788*rhodim/pdim;

  /**
Specific heats and thermal conductivity of the fluids.
   */
  
  cv1 = 3610*rhodim*Tdim/pdim; cv2 = 729.1*rhodim*Tdim/pdim;
  cp1 = 4285*rhodim*Tdim/pdim; cp2 = 1063*rhodim*Tdim/pdim;

  kappa1 = 0.6705/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));
  kappa2 = 0.03153/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));

  mu1 = 3.7e-4/(Rdim*sqrt(rhodim*pdim));
  mu2 = mu1*1e-2;

  /**
The diffusive time scale $\tau$ based on the gas properties.
   */

  tr = rhoR*cp2/kappa2;
  
  tend *= tr;

  init_grid(1 << LEVEL);

  TOLERANCE = 1e-6;

  system ("mkdir dumps");
  
  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) {
    fraction (f, - (sq(Rbub) - sq(x) - sq(y)));
  
    foreach() {
      frho1[]  = f[]*rhoL;
      frho2[]  = (1. - f[])*rhoR;
    
      double pL = p0L;
    
      p[] = pL*f[] + p0*(1.-f[]);

      /**
Initialization of the initial temperature using the the NASG EoS in
compatibility with the initialization of both the density and pressure fields.
       */

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

/**
Outputting the bubble radius in an $0.01*\tau$ increment.
 */

event centroid (t += 0.01*18951.2) {
  scalar ff[];
  foreach()
    ff[] = 1. - f[];
  
  double Volume = statsf(ff).sum;

  if (pid() == 0) {
    FILE * fp = fopen ("radius.txt","a");
    char str[80];
    sprintf(str,"%g %g\n",t/tr,pow(3.*Volume/2.,1./3.));
    fputs(str,fp);
    fclose(fp);
  }
}

/**
Outputting some statistics about the fields.
 */

event logfile (i++) {
  stats sp = statsf (p);
  stats su = statsf (q.x), sv = statsf (q.y);
  stats sT = statsf (T);
  fprintf (stderr,"t = %g, i = %d, dt = %g, max(p) = %g, max(T) = %g, max(u) = %g, max(v) = %g\n",
	   t/tr, i, dt/tr, sp.max, sT.max, su.max, sv.max);
}

/**
Saving dump files for post-processing. (Uncomment)
 */

/* event outputs (t += 0.01*18951.2) { */
/*   char name[80]; */
/*   sprintf (name,"dumps/dump-%g",t/18951.2); */
/*   dump (name, list = (scalar *){f,p,T}); */
/* } */

/**
On the fly movie generation.
 */

event movie (t += 0.01*18951.2) {
  char name[80], s[80];
  sprintf (name,"T0.5.mp4");
  sprintf (s, "t = %.2f", t/18951.2);

  clear();
  view(fov = 12.5, quat = {0,0,-cos(M_PI/4.),cos(M_PI/4.)}, width = 1280, height = 1980);
  draw_vof("f");
  squares("p", min = 1., map = cool_warm);
  mirror({0,1}) {
    draw_vof("f");
    squares("T", min = 0.5548, max = 1.0896, map = cool_warm);
    draw_string (s, pos = 2, size = 16, lc = {255,255,255}, lw = 4);
    draw_string ("Temperature", pos = 3, size = 25, lc = {255,255,255}, lw = 4);
    draw_string ("Pressure", size = 25, lc = {255,255,255}, lw = 4);
  }
  save(name);
}

/**

<center>
<video width="426" height="660" controls>
<source src="AXI_expansion/T0.5.mp4" type="video/mp4">
</video>
</center>

 */

event ending (t = tend) {
  double ending = 0.;
  ending++;
  return ending;
}

/**
~~~gnuplot Bubble radius as a function of time
set xlabel 't'
set ylabel 'R'
plot "radius.txt" u 1:2 w p
~~~
*/
