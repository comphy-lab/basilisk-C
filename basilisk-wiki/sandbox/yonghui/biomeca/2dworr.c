/**
# An womersley flow  With output resistance (STILL HAVE CONVERGENCE PB SO HERE IS A "FAKE" VERSION) 
Details in [2DRR.c](http://basilisk.fr/sandbox/yonghui/biomeca/2drr.c ) 
*/

//#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
#define ADAPT 1
#define alpha_w 10.
#define MIN_LEVEL 3
#define LEVEL 4
#define MAX_LEVEL 5
#define tmax 1.*2.*M_PI
FILE * fp1;
FILE * fp2;
FILE * fp3;
double	Pout; // $P_{out}$
double	DEBIT; // $Q$
double  ORR; // $R_1$
double  VELO;
/**
## main fonction
*/
int main(){
  N = 2 << LEVEL;
  fp1 = fopen ("testpoint", "w");
  fp2 = fopen ("pressure.dat", "w");

  for(int K= 0 ; K <= 3 ; K++){
    VELO = (K+1)/5.;
    fprintf(stderr,"# now velocity = %g \n", VELO );
    run();
  }
}

/**
we set the no slip Boundary conditions for lateral wall*/
u.n[top]  = dirichlet(0.);
u.t[top]  = dirichlet(0.);
//We apply a poiseuille flow in the x direction
u.n[left] = dirichlet( VELO* (1. - sq(y)));
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
  fprintf(stderr,"Theory pente / Re: %g %g \n", viscosity * VELO * 2. , 3.*viscosity);  
  ORR = 0.05;
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
  fprintf(fp1,"%g %g %g %g %g\n" , t, interpolate(u.x, 0.5, 0.5), DEBIT, Pout, ORR);
  // we set the new value as the outlet BC of pressure
  p[right]  = dirichlet(Pout);
  pf[right]  = dirichlet(Pout);
  output_ppm (u.x, file = "velo.mp4", linear = true,n=512,min=-1.,max=1.);
  output_ppm (p, file = "pres.mp4", linear = true,n=512,min=-1.,max=1.);
}

/**
 In this simpla case, we know that $slope = \frac{\partial p}{\partial x} = \mu \frac{\partial^2 p}{\partial y^2} = 2 \mu U_0$*/
event tracer (t = tmax) {
  double pente = (interpolate(p , 0.05, 0.05)-interpolate(p , 0.95, 0.05)) /0.9;
  double Reee =  pente /DEBIT;
  for (double xx = 0.; xx < 1.; xx += L0/100. ){
    double yy = L0/2.;
    double pth = Pout + (1.-xx)*pente;
    fprintf(fp2,"%g %g %g %g\n", t, xx, interpolate(p , xx, yy),pth);
  }
  fprintf(stderr,"DEBIT / pente / Re :%g %g %g\n", DEBIT, pente ,Reee); 
  //	fprintf(stderr,"P0,P2: %g %g \n",interpolate(p,0.05,0.05), interpolate(p,0.95,0.05));
  fprintf(fp2,"\n \n");
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
plot 'pressure.dat' u 2:3 w l t'pressure' ,\
'pressure.dat' u 2:4 w l t'pressure theory'
~~~

*/