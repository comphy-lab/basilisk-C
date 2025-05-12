/**
# collapse of a rectangular viscous column,

or collapse of a viscous fluid (double viscous dam break)
From the paper: Huppert 82 “The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface”
it is done with [only mass equation and lubrication ](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and 
  [with shallow watzer](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c)
 
Here we solve it using the Multilayer Shallow Water (Saint Venant Multi Couches) strategy of Audusse Sainte-Marie  et al 2011. See De Vita 2020 for details, this example is presented there as a test case.
 
 
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

int main() {
  X0 = 0.;
  L0 = 5;
  G  = 1.;
 // N  = 512;
/**
 one needs details in $y$ */
 // nl = 256;
  N  = 128;
  nl = 15;
  nu = 1.;
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

  for (vector u in ul) {
    u.n[left] = 0;
    u.n[right] = neumann(0.);
  }
  
  /**
  We initialize *h*. */
  foreach()
    h[] =   (x<1);
}

/**
## Output
We print the elevation and the stress. */
 
event output (t += 5; t <= 100) {
  vector u0 = ul[0];
  foreach()
    fprintf (stderr, "%g %g %g %g\n", x, h[], 2.*u0.x[]/((h[]+dry)/nl),t );
     fprintf (stderr, "\n");
}
/**

## Run

To run 

~~~bash  
qcc -O2 -o viscous_collapse_ML viscous_collapse_ML.c 
 ./viscous_collapse_ML  2>log
~~~

 
 To run
 
~~~bash
make viscous_collapse_ML.tst
make viscous_collapse_ML/plots
make viscous_collapse_ML.c.html

source c2html.sh viscous_collapse_ML
 ~~~

## Results

the selfsimilar collapse over the selfsimilar solution:
 $h(x,t)t^{1/5}$ plotted as a function of $(xt^{-1/5})$

~~~gnuplot 
 set xlabel "x/t^{1/5}"
 set ylabel "h(x,t) t^{1/5}"
 p [0:1.5]'log' u ($1/($4**.2)):($2*($4**.2)) t'comp.' w l,(9./10*(1.28338-x*x))**(1/3.) t'analytic'
~~~



Montpellier 07/17

## Links
related examples [only mass equation and lubrication ](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and 
  [with shallow watzer](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c)
see as well Navier Stokes solution.



## Bibliography
 
* Huppert
"The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface" 
[JFM 82](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 
* Lagrée  [M2EMN
 Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)

* Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
 [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
 Volume 79, January–February 2020, Pages 233-246
 European Journal of Mechanics - B/Fluids
 https://doi.org/10.1016/j.euromechflu.2019.09.010

*/
