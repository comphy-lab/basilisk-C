/**
# Propagation of an acoustic disturbance in a tube

In this test, We quantify the numerical dissipation of the all-Mach solver in the acoustic limit
by simulating the propagation of a gaussian pulse of small amplitude. The results show the conclusions of the previous test(http://www.basilisk.fr/sandbox/fuster/Allmach3.0/gaussianaxi.c) quantitatively.
*/

#include "grid/multigrid.h"
/** We use the two-phase flow formulation */

#include "two-phase-compressible.h"

/** Parameters of the problem */
double tend = 4.;
double cflac = 0.01;
double sigmaP = 1.4;
double cson;
double deltaP = 1.e-3;

event stability (i++) {
  double dt = 100.;
  foreach ()
    dt = min(dt,Delta/sqrt(gamma1));
  dtmax = dt*cflac;
  DT = dt*cflac;
}


int main()
{  

  /** Size of the domain: */
  size (20.);
  X0 = -L0/2.;


  /** The EOS for an adiabatic perfect gas is defined by its polytropic coefficient $\Gamma = \gamma = 1.4$ */
  gamma1 = 1.4;

  /** We use an upwind method for the tracer advection associated to the VOF tracer f*/
  f.gradient = zero; 

  /** We perform a convergence study */
  N = 256;
  //for (sigmaP = 0.4; sigmaP <= 1.4; sigmaP += 0.05) {
  //for (cflac = 0.1; cflac <= 100.; cflac *= 1.4) {
  run();
  //}
  //}

}

event init (i = 0)
{   
  cson = sqrt(gamma1);

  foreach() {
    f[] = 1.;
    p[] = (1.+ deltaP*exp(-x*x/sq(sigmaP)));
    frho1[] = (1. + (p[] - 1.)/sq(cson));
    q.x[] = 0.;
    q.y[] = 0.;
    fE1[] = p[]/(gamma1 - 1.) + 0.5*pow(q.x[],2)/frho1[];
  }
  boundary ((scalar *){q, frho1, p, fE1});
}

event endprint (t = tend) {

  char name[80];
  sprintf (name, "snapshot-%g-%3.2f", sigmaP, cflac);
  FILE * fp = fopen(name, "w");  

  foreach () 
    if ( y < Delta && x > 0.) 
    fprintf(fp, "%g %g \n", x, p[]);

    fclose(fp);


  double pmax = 0;
  double xmax = 0;
  foreach () {
    if (p[]>pmax && x > 0.){
      pmax = p[];
      xmax = x;
    }
  }
  /** For a plane wave propagating in x direction, we calculate the CFL using the relation $ P = /rho c_{son} u_{x} $ in the acoustic limit */
  printf("%g %g %g %g %g \n", sigmaP, cflac, deltaP/2./cson*cflac, fabs(xmax/(cson*tend)-1.), fabs((pmax - 1.)*2./deltaP - 1.0));
}

/**

~~~gnuplot Numerical error of attenuation
set dgrid3d 30,30
set pm3d map
set logscale y
unset key
set title "Numerical error of attenuation"
set ylabel "CFL"
set autoscale xfix
set xlabel "sigmaP"
set autoscale yfix
sp "out" u 1:3:5 palette
~~~ 

*/

/**

~~~gnuplot Numerical error of sound speed
set dgrid3d 30,30
set pm3d map
set logscale y
unset key
set title "Numerical error of sound speed"
set ylabel "CFL"
set autoscale xfix
set xlabel "sigmaP"
set autoscale yfix
sp "out" u 1:3:4 palette
~~~ 

*/