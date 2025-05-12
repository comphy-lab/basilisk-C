/**
# Behavior near the front of the water train

## Problem

 The question is : what is the behavior of a perturbation of the free surface at large time and distance?
 
The response may be found using asymptotic behavior of a wave packet.
 It may be obtained from the superposition of all the waves :
 $$\eta = \int_{-\infty}^{\infty}F(k) e^{i(kx-\omega t) }dk$$
 considering that we deal with long wave approximation of the dispersion relation of water wave
 $$\omega =   \sqrt{gh_0}(k - \frac{1}{6} k^3h_0^2+...)$$
 and using the steepest descent method one obtains the asymptotic behavior
 (see Whitham)
 $$ \eta \sim 2 \sqrt{\pi } F(0)  \left(   \frac{2}{( h_0^2 c_0 t ) (c_0t  - x)}\right)^{1/4} \sin \left( \frac{2}{3}
 \frac{(c_0t - x)^{3/2}}{  ( h_0^2c_0 t /2)^{1/2}} + \frac{\pi}{4}
 \right).$$
 A simpler way is to start by the linearized KdV (corresponding to $\omega =   c_0 k - \frac{c_0 h_0^2}{6} k^3)$):
 $$\frac{\partial \eta}{\partial t}+c_0 \frac{\partial \eta}{\partial x} + \frac{c_0h_0^2}6 \frac{\partial^3 \eta}{\partial x^3}=0.$$
 it is possible to obtain a selfsimilar solution from this equation
 $$\eta(x,t)=(\frac{2}{c_0h_0^2 t})^{1/3}Ai[(\frac{2}{c_0h_0^2})^{1/3} \frac{(x-c_0t)}{t^{1/3}} ].$$
 see  Whitham  13.6 page 441 for details (of course the asymptotic behavior with a sine is an asymptotic approximation of the Airy function)
 
 
 
 # Code
 
 This can be done by the Multilyer Euler Lagrange `layered/hydro.h`
  with `layered/nh.h`the non hydrostatic correction, or w
 
 
 */
#define ML 1
#include "grid/multigrid1D.h"
#if ML
  #include "layered/hydro.h"
  #include "layered/nh.h"
 // #include "layered/remap.h"
#else
  #include "green-naghdi.h"
#endif

/**
The bassin needs to be long enough so as to minimise the influence of
wave reflection at the outlet. Relatively high resolution is needed to
capture the dynamics properly. */

int main() {
  N = 2048*8*2;
  L0 = 400;
  X0= -20;
  G = 1;
  DT = HUGE;
#if ML
  N=N/8;
  nl = 4;
  //breaking = 0.1;
#endif
  run();
}

/**
We use ["radiation"
conditions](/src/elevation.h#radiation-boundary-conditions) at the
inlet and outlet. */

event init (i = 0) {

#if ML
    u.n[left]  = - radiation (0);
    u.n[right] = + radiation (0);
 
#else
  u.n[left]  = - radiation (0);
  u.n[right] = + radiation (0);  
#endif
  
  /**
  Here we define the flat bathymetry and the initial perturbation, depth is $h_0=1$  */
  foreach() {
    zb[] = -1.;
#if ML
    foreach_layer()
      h[] =  (- zb[] + 0.05*exp(-pi*x*x))/nl ;
#else
    h[] = - zb[] + 0.05*exp(-pi*x*x);
    u.x[] = 0.*exp(-pi*x*x);
#endif
  }
}

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=350$.

![Snapshot of the water front wave, comparison with the analytical Airy solution.](airy_watertrainfront/snapshot.png)
*/
#if 1
void plot_profile (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
             "p [-10:5][-1:1]'-' u 1:(2*$2/.05) w l lc 3 t 'num' , airy(x) t'Airy(x)' \n", t);
  foreach()
    fprintf (fp, "%g %g \n", (x-t)/pow((t>0?t/2:1e-4),1./3),eta[]*pow(t/2,1./3) );
  fprintf (fp, "e\n\n");
  fflush (fp);
}
#else
void plot_profile (double t, FILE * fp)
{
    fprintf (fp,
             "set title 't = %.2f'\n"
               "p [-20:380][-0.02:0.01]'-' u 1:3:2 w filledcu lc 3 t ''\n", t);
    foreach()
    fprintf (fp, "%g %g %g\n", x,eta[],zb[]);
    fprintf (fp, "e\n\n");
    fflush (fp);
}

#endif

event profiles (t += 1) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  plot_profile (t, fp);
  fprintf (stderr, "%g %f\n", t, interpolate (eta, 17.3, 0.));
}

event gnuplot (t = end) {
  FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
           "set term png enhanced size 640,200 font \",8\"\n"
           "set output 'snapshot.png'\n");
  plot_profile (t, fp);
   // getchar();
}

/**
 Water gauges  in $x=50,...300$  */

Gauge gauges[] = {
  {"WG50", 50},
  {"WG100", 100},
  {"WG150", 150},
  {"WG200", 200},
  {"WG250", 250},
  {"WG300", 300},
  {NULL}
};
/**

There is now a bug 

`warning: could not conserve barotropic flux at -14.4336,0,3

[Detaching after vfork from child process 630414]
64 /home/basilisk-src/basilisk/src/hessenberg.h: No such file or directory.
src/hessenberg.h:64:error: Program received signal SIGFPE, Arithmetic exception``

*/
event output (t += 1; t <= 180) //350
  output_gauges (gauges, {eta});


/**

## Compilation

~~~bash
 qcc airy_watertrainfront.c -o airy_watertrainfront
 ./airy_watertrainfront
 
 
 make airy_watertrainfront.tst
 make airy_watertrainfront/plots;make airy_watertrainfront.c.html
 

 sed -i - 's/\\)//g' airy_watertrainfront.c.html
 sed -i - 's/\\(//g' airy_watertrainfront.c.html
 sed -i - 's/\\\[//g' airy_watertrainfront.c.html
 sed -i - 's/\\\]//g' airy_watertrainfront.c.html
 open airy_watertrainfront.c.html
 
~~~
 

## Results

Comparison of theoretical Airy solution and numerical timeseries in $x=50,...300$ for $t$ rescaled.

~~~gnuplot Comparison of theoretical and numerical timeseries
 p[-5:10]\
 'WG50' u (($1-50)/(50**(1./3))):($2*(50**(1./3))) not ,\
 'WG100' u (($1-100)/(100**(1./3))):($2*(100**(1./3))) w l not,\
 'WG150' u (($1-150)/(150**(1./3))):($2*(150**(1./3))) w l not,\
 'WG200' u (($1-200)/(200**(1./3))):($2*(200**(1./3))) w l not,\
 'WG250' u (($1-250)/(250**(1./3))):($2*(250**(1./3))) w l not,\
 'WG300' u (($1-300)/(300**(1./3))):($2*(300**(1./3))) w l not,\
 0.5*0.05*1.2599*airy(-1.2599*x) t'Analytical Airy solution'

 
~~~

## Links
 * [Popinet (2019)](/Bibliography#popinet2019)  from [bar.c](http://basilisk.fr/src/test/bar.c) test case
 
 * [PYL](http://basilisk.fr/sandbox/M1EMN/Exemples/boussinesqc.c) linearised example in C of Boussinesq, boussinesqc.c
 
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/airy_watertrainfront.c]() this example

 * [http://basilisk.fr/sandbox/M1EMN/Exemples/ressaut_mascaret.c]()
 Poor’s man dispersive model





 
## References

 * Hinch "Perturbation methods" page 30
 
 * Whitham "Linear and non linear waves" page 441
 
 * Mei  "Ocean Surface Waves" 2005
 
 * Noda JGR 76 (30) 1971,Water Waves Generated by a Local Surface Disturbance.

 * [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEhoule.pdf) details for front computation.
 
 * [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/M2MHP/kdv.pdf) about KdV equation.
 
 
 
*/
