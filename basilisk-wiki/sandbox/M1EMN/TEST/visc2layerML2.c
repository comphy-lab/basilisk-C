/**
# collapse of a rectangular viscous column on a slope,

 
From the paper: Huppert 82 "Flow and instability of a viscous current along a slope"

the 1D theory is presented  [with shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c), here we solve the full profiles

Here we solve it using the Multilayer Shallow Water (Saint Venant Multi Couches) strategy of Audusse Sainte-Marie  et al 2011. See De Vita 2020 for details, this example is presented there as a test case.
 
 */
#define ML 1
#define IMPLICIT 1

#include "grid/multigrid1D.h"
#if IMPLICIT

#if ML
  #include "layered/hydro.h"
  #include "layered/implicit.h"
  //#include "layered/remap.h"
#else 
  #include "saint-venant-implicit.h"
#endif //ML
#else //IMPLICIT

#if ML
  #include "layered/hydro.h"
  //#include "layered/remap.h"
#else 
  #include "saint-venant-implicit.h"
#endif //ML

#endif //IMPLICIT

const double HR = 0.5;
double Cf;
double alpha;  

int main() {
  X0 = -5.;
  L0 = 30;
  G  = 1.;
  alpha=0.5;
  N  = 128;//512;
  nl = 2;
  nu = 1.;
  Cf = 3.;

  #if ML //fixme: need to refine dt when using hydro.h 
   CFL_H = 0.125;
  #endif  
  run();
}

/**
We impose boundary condition for $h$ and $\eta$. */
// h[left] = neumann (0);
// eta[left] = neumann (0);
// u.n[left] = dirichlet(0);

// h[right] = neumann (0);
// eta[right] = neumann (0);

/**
## Initialization  */

event init (i = 0) {
  /**
  We set a zero velocity at the inlet and a free outlet. */

  #if !ML
    for (vector u in ul) { 
      u.n[left] = dirichlet(0); //0;
      u.n[right] = neumann(0.);
    }
  #else
    u.n[left] = dirichlet(0); //0.;
    u.n[right] = neumann(0.);
  #endif

/**
We impose boundary condition for $h$ and $\eta$. */  

 // h[left] = neumann (0);
  //h[right] = neumann (0);
  //eta[left] = neumann (0);
  //eta[right] = neumann (0);
  //u.n[left] = dirichlet(0);

/** we initialize *h* and *z*. */
  foreach(){
    zb[] = -(x-X0)*alpha; 
  #if !ML  
    h[] = (fabs(x)<HR);//(fabs(x)<.5);
  #else
    foreach_layer()
      //h[] =  (fabs(x)<HR)/nl;
      h[] =  (fabs(x)<HR) ? 1./nl : dry; //x<.5  ? 1./nl : dry; //
  #endif
  }
}


// event friction (i++) {
//   // Poiseuille friction (implicit scheme)
//   foreach(){
//     double H=0.;
//     foreach_layer()  
//       H+=h[];
//     //double ff = h[] < dry ? HUGE : (1. + Cf*dt/h[]/h[]);
//     double ff = h[] < dry ? HUGE : (1. + Cf*dt/H/H);
//     u.x[] /= ff;
//   }
// }


/**
## Output
We print the elevation and the stress. */

event output (t += 10;  t <= 1000) {
	#if !ML  
    vector u0 = ul[0];
    foreach()
      fprintf (stderr, "%g %g %g %g %g\n", x, h[], 2.*u0.x[]/((h[]+dry)/nl),t, dt );
	#else
    foreach(){  
      double H=0.; 
      foreach_layer()
        H += h[];
      fprintf (stderr, "%g %g %g %g %g %g\n", x, H, 2.*u.x[]/(h[]+dry),t, dt, CFL_H); 
      //fprintf (stderr, "%g %g %g %g %g %g\n", x, H, 2.*u.x[]/h[],t, dt, CFL_H); 
    }   	
  #endif 
  fprintf (stderr, "\n");
}

/**

## Run

To run 

~~~bash  
qcc -O2 -o viscous_collapsesqrt_ML viscous_collapsesqrt_ML.c 
 ./viscous_collapsesqrt_ML  2>log
~~~

To make
 
~~~bash
 make viscous_collapsesqrt_ML.tst
 make viscous_collapsesqrt_ML/plots
 make viscous_collapsesqrt_ML.c.html
 
 source c2html.sh viscous_collapsesqrt_ML
~~~

## Results


~~~gnuplot 
 set xlabel "x"
 set ylabel "h(x,t)"
p[-2:]'log' w l
~~~

the selfsimilar collapse over the selfsimilar solution: in $\sqrt(\eta)$: (the overshoot desappears with *N*)

~~~gnuplot 
 set xlabel "x/t^{1/3}"
 set ylabel "h(x,t) t^{1/3}"
 n=1./3
 p [-0.5:1.5]'log' u ($1/($4**n)):($4>2000? $2*($4**n): NaN) t'comp.' w l,sqrt(x/(.5)+.2) t'anal''
~~~



Montpellier 07/17


## Links

* related example   [with shallow watzer](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c)

* see other multilayer examples

# Bibliography
 
 
* [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper49.pdf)
 "Flow and instability of a viscous current along a slope"
 Nature volume 30 1982 p 427
 
* Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
 [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
 Volume 79, January–February 2020, Pages 233-246, European Journal of Mechanics - B/Fluids
 https://doi.org/10.1016/j.euromechflu.2019.09.010

*/
