/**
# 2D Kelvin-Helmholtz 

This is a shear flow instability simulation created by flows with different velocities and densities : the shearflow simulation in 2D.
*/
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tag.h"
#include "two-phase.h"

/**
Uncomment for conserving method, which implements the momentum-conserving VOF advection of the velocity components.
*/

//#include "navier-stokes/conserving.h"

/**
We define air and water densities in both two-phase.h and this simulation.
*/

double Pi = 3.141592653589793;
double rho_1 = 1000.;
double rho_2 = 1.;

/**
Initialization of the parameters.
*/

int maxlevel;
double vis, Re; 
FILE * fpn;
char fnamev[99];
char fnameg[99];
scalar omega[], lev[], k[];

/**
We setup the shape function which will allow the fraction function display.
*/

double shape (double x, double y)
  {
    double p1 = 0.1 - y ;
    double p2 = 0.1 + y ;
    return min ( p1 , p2 ) ;
  }

/**
The domain is periodic in $x$ and resolved using 128$^2$ points. */
int main(){
  periodic(left); 

  X0 = Y0 = -0.5*L0; // flow position
  size (1.); // domain size 

  maxlevel = 7 ; 
  init_grid (1<< 4); //grid
  display_control (maxlevel, 4, 12);
  rho1 = rho_1 ;
  rho2 = rho_2 ;
  N = 128;
  run();
}

/**
We don't involve viscosity effects.
*/
event init(t = 0){
  vis = 0.;
  const face vector muc[] = {0., 0.};
  mu = muc ;

  /**
  Initial flows' velocities : 
  $u_1 = 15$ if $|y|>1/10$ ; $u_1 = 1$ if $|y|<1/10$.
  */
  
  foreach()
    {
      u.x[] = ( fabs(y) < 0.1 ? 1. : 15. ); 
      u.y[] = 0.01*sin(2.*Pi*x)*exp((-20.)*pow((y),2.));// perturbation 
    }
  
 boundary({u.x});
  boundary({u.y});
  fraction (f,shape(x,y));
}

/**
We compute the vorticity.
*/
event output (t = 0.; t += 0.001){

  int n = 0;
  double m = 0;
  foreach(reduction(+:n) reduction(max:m)){
    n++;
    lev[] = level;
    omega[] = (u.x[0,1]-u.x[0,-1] - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
    if (fabs(omega[]) > m)
      m = fabs(omega[]);
  }
  boundary({omega});
}

/**
We now compute the kinetic energy and stop the simulation when it is exploding.
*/

static double energy()
{
  double ke = 0.;
    foreach(reduction(+:ke))
      ke += (sq(u.x[]) + sq(u.y[]))/2.*sq(Delta);
      if (ke > 500.) exit(1); 
  return ke;
}


event kinetic_energy (i += 25)
{
  fprintf (stderr, "%g %g\n", t, energy());
}

event end (t = 1.)
  fclose(fpn);

/**
# Kinetic energy comparison
We log the evolution of the kinetic energy.

*/


/*![[Kinetic energy comparison](sandbox/Benoit/Kelvin-Helmholtz/kh.c)(/sandbox/nlemoine/README)](sandbox/Benoit/Kelvin-Helmholtz/ke.png)*/
  
/**
~~~gnuplot Evolution of the kinetic energy
set xlabel 't'
set ylabel 'Kinetic energy'
set title 'Kinetic energy as a function of time'
set autoscale
plot 'log128t=1g' u 1:2 w point pointtype 6 ps 1 lc rgb "blue" title 'Conserving N = 128', \
     'log128d' u 1:2 w point pointtype 6 ps 1 lc rgb "red" title 'Centered N = 128', \
     'log256t=1g' u 1:2 w point pointtype 6 ps 1 lc rgb "green" title 'Conserving N = 256', \
     'log256d' u 1:2 w point pointtype 6 ps 1 lc rgb "orange" title 'Centered N = 256'
~~~
*/

/**/
   