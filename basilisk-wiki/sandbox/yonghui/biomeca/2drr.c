/**
# An poiseuille flow  With output resistance 
Here, we study the effect of the output resistance in the axisymmetric tube and compare our results with the theoretical values. Mesh adaptation is used by default.

Here $ P_ {out} = R_1 * Q $ instead of 0, $ R_1 $ is a constante coefficient resistance, and $ Q $ is the value of the mass flow. A more complex parallel RC conditional model will be applied in another code later.

This method is commonly used to simulate the effects of external vascular models. Details in [KVmodel](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Voigt_material) 
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
#define ADAPT 1 // 0: no adaption
#define alpha_w 10.
#define ORR 0.15
#define MIN_LEVEL 3
#define LEVEL 5
#define MAX_LEVEL 7
#define tmax 1*2.*M_PI
FILE * fp1;
FILE * fp2;
FILE * fp3;
double	Pout; //$P_{out}$
double	DEBIT; // $Q$ 

/**
## main fonction
*/
int main(){
  N = 2 << LEVEL;
  fp1 = fopen ("testpoint", "w");
  fp2 = fopen ("pressure.dat", "w");
  run();
}

/**
we set the no slip Boundary conditions for lateral wall*/
u.n[top]  = dirichlet(0.);
u.t[top]  = dirichlet(0.);
//We apply a poiseuille flow in the x direction
u.n[left] = dirichlet( 0.25 * (1. - sq(y)));
u.t[left] = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = dirichlet(0.);
//pf is the facet pressure
p[left]  = neumann(0.);
pf[left]  = neumann(0.);
/** 
##initial event 
here we calculate the $\mu$ based on womersley number $\alpha$*/
event init (t = 0) {
  double viscosity = 1. / sq(alpha_w);
  const face vector muc[] = {viscosity , viscosity};
  mu = muc;
  fprintf(stderr,"Corresponding viscosity: %g\n", viscosity);  
}

/**
## Time integration
In each iteration, we calculate the value of the mass flow and apply a new value to the output pressure based on it.
We record the pressure value of a choosen point in order to check if the result is converged after the simulation (periodic form).*/
event midpressure(t <= tmax; i += 1) {
  DEBIT = 0.;
  // we calculate the mass flow 
  foreach_boundary(right)
    DEBIT += u.x[]*Delta;
  // we calculate the $P_{out}$
  Pout = ORR * DEBIT;
  double px = L0/2.;	
  double py = 0.;
  fprintf(fp1,"%g %g %g %g %g\n" , t, interpolate(u.x, px, py), DEBIT, Pout, ORR);
  // We set the new value of outlet pressure to BC
  p[right]  = dirichlet(Pout);
  pf[right]  = dirichlet(Pout);
}

/**
We output the value of pressure for a stright line $y = 0.5$, which we can compare it with the results of p_{out} = 0.*/
event tracer (t = end) {
  for (double xx = 0.; xx < 1.; xx += L0/100. ){
    double pyy = L0/2.;
    fprintf(fp2,"%d %g %g %g %g\n", N, t, xx, interpolate(p , xx, pyy), Pout);
  }
}

/**
## Mesh adaptation
We adapt the mesh according to the error on the volume fraction field and the velocity. */
#if ADAPT
event adapt (i++) {
  adapt_wavelet((scalar *){u}, (double[]){5e-4,1e-3}, MAX_LEVEL, MIN_LEVEL) ;
}
#endif

/**
# Results
~~~gnuplot testpoint
plot 'testpoint' u 1:2 w l t'u.x' ,\
'testpoint' u 1:3 w l t'debit' ,\
'testpoint' u 1:4 w l t'Pressure' ,\
'testpoint' u 1:5 w l t'ORR'
~~~

~~~gnuplot pressure at y = L0/2
plot 'pressure.dat' u 3:4 w l t'pressure' ,\
'pressure.dat' u 3:5 w l t'pressure out'
~~~

*/