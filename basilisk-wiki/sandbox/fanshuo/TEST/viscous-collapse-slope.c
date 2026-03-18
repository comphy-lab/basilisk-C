/**
# collapse of a rectangular viscous column on a slope,

 
From the paper: Huppert 82 "Flow and instability of a viscous current along a slope"

the 1D theory is presented  [with shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c), here we solve the full profiles

Here we solve it using the Multilayer Shallow Water (Saint Venant Multi Couches) strategy of Audusse Sainte-Marie  et al 2011. See De Vita 2020 for details, this example is presented there as a test case.
 
 */
#define ML 0

# include "grid/multigrid1D.h"
# if ML
  //#include "layered/hydroMT.h"
#  include "layered/hydro.h"
//#  include "layered/remap.h"
  //#include "layered/implicit.h"
# else 
#  include "saint-venant.h"
#endif //ML

const double HR = 0.5, T0 = 100;
double alpha, slope;  

int main() {
  X0 = -5.;
  L0 = 30;
  G  = 1.;
  alpha = 0.5;
  slope = atan(-alpha);
  N  = 128;
  nl = 15;
  nu = 1.;

  #if ML //fixme: need to refine dt  
     CFL_H = 0.0625; // fixme, CFL has to be very small otherwise the tail slips until it reaches the boundary 
  #endif  
  run();
}

/**
We impose boundary condition for $h$ and $\eta$. */
h[left] = neumann (0);
eta[left] = neumann (0);
u.n[left] = dirichlet(0);

h[right] = neumann (0);
eta[right] = neumann (0);

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

/** we initialize *h* and *z*. */
  foreach(){
    zb[] = 0.;
    //zb[] = -(x-X0)*alpha; 
  #if !ML  
    h[] = (fabs(x)<HR );
  #else
    foreach_layer()
      //h[] =  (fabs(x)<HR)/nl;
      h[] = (fabs(x)<HR) ? 1./nl : dry; 
  #endif
  }
}

#if ML  
event acceleration (i++, t<=T0)
{
  foreach_face()
    foreach_layer()
      ha.x[] += - G*sin(slope)*hf.x[]; 
}
#else
event acc (i++, t<=T0){
  foreach () 
    for (vector u in ul) {
      u.x[] = u.x[] - G*sin(slope)*dt;
  }
}
#endif

/**
## Output
We print the elevation and the stress. */

event output (t += 10;  t <= T0) {
#if !ML  
  vector u0 = ul[0];
  foreach()
    fprintf (stdout, "%g %g %g %g %g %g\n", x, h[], 2.*u0.x[]/((h[]+dry)/nl),t, dt, CFL );
#else
  foreach(){  
    double H=0.;//zb[]; 
    foreach_layer()
      H += h[];
      fprintf (stdout, "%g %g %g %g %g %g\n", x, H, 2.*u.x[]/(h[]+dry),t, dt, CFL_H); 
  }   	
#endif 
  fprintf (stderr, "\n");
}

// fixme, it does not work as well as saint-venant.h
/**

## Results


~~~gnuplot 
 set xlabel "x"
 set ylabel "h(x,t)"
p[-2:]'out' w l
~~~

The selfsimilar collapse over the selfsimilar solution : in $\sqrt(\eta)$: (the overshoot desappears with *N* but it works only with saint venant for now) 

~~~gnuplot 
 set xlabel "x/t^{1/3}"
 set ylabel "h(x,t) t^{1/3}"
 n=1./3
 p [-0.5:1.5]'out' u ($1/($4**n)):($4>2000? $2*($4**n): NaN) t'comp.' w l,sqrt(x/(.5)+.2) t'anal''
~~~

## Links

* related example   [with shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c)

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
