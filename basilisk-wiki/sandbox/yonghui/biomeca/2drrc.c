/**
# An poiseuille flow  With output resistance 
Here, we study the effect of the output resistance in the axisymmetric tube and compare our results with the theoretical values. Mesh adaptation is used by default.

Here P_ {out} is similar to a paraller RC circuit (Details in [KVmodel](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Voigt_material) ) instead of 0,  ($R_1, R_2, C$) are resistors and capacitors ,and Q is the value of the mass flow. This method is commonly used to simulate the effects of external vascular models.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
#define ADAPT 1 
#define alpha_w 10. 
#define MIN_LEVEL 3
#define LEVEL 5
#define MAX_LEVEL 6
#define tmax 1*2.*M_PI
#define  ORR 0.1 //R1
#define  ORC 0.1 //R2
#define  OCC 10  //C

FILE * fp1; FILE * fp2; FILE * fp3;
double	Pout,Pold;
//double	Pold;
double	DEBIT, Qold;
//double	Qold;
double	PPPP;

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
u.n[left] = dirichlet(0.25 * (1. - sq(y)));
//u.n[left] = dirichlet(0.25*cos(t) * (1. - sq(y)));
u.t[left] = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = dirichlet(0.);
p[left]  = neumann(0.);
pf[left]  = neumann(0.);

/** 
##initial event 
here we calculate the $\mu$ based on womersley number $\alpha$.*/
event init (t = 0) {
  double viscosity = 1. / sq(alpha_w);
  const face vector muc[] = {viscosity , viscosity};
  mu = muc;
  fprintf(stderr,"Corresponding viscosity: %g\n", viscosity);   
  Pout =0.;
  Pold =0.;
  Qold =0.;
}

/**
## Time integration
In each iteration, we calculate the value of the mass flow and apply a new value to the output pressure based on it.
We record the pressure value of a choosen point in order to check if the result is converged after the simulation (periodic form).*/
event midpressure(t <= tmax; i += 1) {
  DEBIT = 0.; 
  // we calculate the Q
  foreach_boundary(right){
    DEBIT += u.x[]*Delta;
  }

  /**
  DETAILS OF THE THEORY CALCULATION

Here $R_e$ corresponding the resistance equivelent of our model,$R_1$ is the first resistance outside, $R_2 \& C $ are the paraller RC after,  We can define 2 variable:

$\alpha = \frac{R_e+R_1}{R_2}+1$  and $\beta = 1 - \frac{R_e}{R_{tot}}$

We can see in the video that  

$\mathbf{i}_{tot} = \frac{V}{R_{tot}}(1-exp(\frac{-\alpha t}{(R_e+R_1)C})) 
                    +\frac{V}{R_e+R_1}exp(\frac{-\alpha t}{(R_e+R_1)C})$

With simple caculation we can get our theory outlet pressure

$V_{out} = \beta V (1.-exp(\frac{-\alpha t}{(R_e+R_1)C}))$
*/
  double RRe =0.8;
  double alpha = (RRe + ORR) / ORC + 1.;
  double beta =  1. - RRe/(RRe + ORR +ORC);
  double QQQ = beta * (1.-exp( -alpha*t/((RRe+ORR)*OCC)))* DEBIT;
  /**
  DETAILS OF THE OUTPUT PRESSURE CALCULATION

The equation od the KV model is 

$P + R_2 C \partial_t P = (R_1 + R_2) Q + R_1 R_2 C \partial_t Q $ 

which can be write in the form

$P (1 + \frac{R_2 C}{dt} )= \frac{R_2 C}{dt} P_{old} + (R_1 +R_2)Q + R_1 R_2 C  \frac{Q - Q_{old}}{dt}$
  */

  double term1 = ORC * OCC * Pold / dt;
  double term2 = (ORR + ORC) * DEBIT ;
  double term3 = ORR * ORC * OCC * (DEBIT - Qold) / dt;
  double deno = 1. + ORC * OCC / dt ;
  //Pout = t < tmax / 10.? (ORR * DEBIT) : ((term1 + term2 + term3) / deno);
  Pout =  (term1 + term2 + term3) / deno;
  fprintf(fp1,"%g %g %g %g %g %g\n" , t, DEBIT, Qold, Pout, Pold,QQQ);
  PPPP = (ORR +ORC) * DEBIT; 

  /** 
  We set the new value of outlet pressure to BC, and give the P Q to the Pold Qold*/
  p[right]  = dirichlet(Pout);
  pf[right]  = dirichlet(Pout); 
  Pold = Pout;
  Qold = DEBIT;
}

/**
We output the pressure in position $y=0.5$.*/
event tracer (t = end) {
  for (double xx = 0.; xx < 1.; xx += L0/100. ){
    double pyy = L0/2.;
    fprintf(fp2,"%d %g %g %g %g %g\n", N, t, xx, interpolate(p , xx, pyy), Pout,PPPP);
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
plot 'testpoint' u 1:2 w l t'debit' ,\
'testpoint' u 1:3 w l t'Q' ,\
'testpoint' u 1:4 w l t'Pressure' ,\
'testpoint' u 1:5 w l t'Pold'
~~~

~~~gnuplot convergence
plot 'testpoint' u 1:5 w l t'Pressure-time' ,\
'testpoint' u 1:6 w l t'theory'
~~~


~~~gnuplot pressure at y = L0/2
plot 'pressure.dat' u 3:4 w l t'pressure' ,\
'pressure.dat' u 3:5 w l t'pressure out'
~~~
*/
