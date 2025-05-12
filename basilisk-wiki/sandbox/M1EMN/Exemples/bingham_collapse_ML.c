/**
# collapse of a rectangular Bingham column on a slope,
 
## Problem
 
 A heap of fluid following Bingham rheology is released along a constant slope.
 We see the front moving to the right, and the left front going slowly up hill to the left. A real Bingham flow should stop. We do not solve NS but we solve RNSP equations. 
 
 
## RNSP 
We solve the boundary layer (RNSP) equations corresponding to a thin layer approximation of the Navier Stokes equations:
 $$\left\{ 
\begin{aligned}
   \dfrac{\partial  u}{\partial   x}  + \dfrac{\partial  v}{\partial   y} &= 0  \\
 \rho (\dfrac{\partial   u}{\partial   t} + \dfrac{\partial   u^2}{\partial   x}  + \dfrac{\partial  u  v}{\partial   y} )&=  -   \rho g  Z'_b -
  \dfrac{\partial   p}{\partial   x} + 
   \dfrac{\partial   \tau_{xy}}{\partial   y} \\
0 &=   -\rho g -\dfrac{\partial  p}{\partial  y}.
\end{aligned}
\right.$$
with  $\tau_{xy}= (\tau_y + \mu \dfrac{\partial u}{\partial z})$
with the Multilayer Solver technique from Audusse Sainte Marie.
We tune the solver  to consider a non newtonian viscosity function of the shear
$$\nu_{eq}=\nu + \frac{\tau_y/\rho}{|\frac{\partial u}{\partial z}|}$$
the `multilayer.h` is changed in `saint-venantNN.h` (all is in `saint-venantB.h`).

## Link with 1D model

If we integrate over the depth of flow $h$ the incompressiblility equation, we obtain 
$\frac{\partial h}{\partial t} + \frac{\partial Q}{\partial x}=0$,
where $Q=\int_0^h udy$. If we neglect inertia in the momentum equation, we solve for $u$ in 
$$0 =  -   \rho g  Z'_b -
  \dfrac{\partial   p}{\partial   x} + 
   \dfrac{\partial   \tau_{xy}}{\partial   y} $$
 were the pressure is hydrostatic $p=\rho g (h-y)$, and with the Bingham rheology. This gives then $Q=\int_0^h udy$ that we put in the mass conservation, and hence we obtain the [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c) from the paper of Balmforth. 
Here we use the Multilayer Solver to see the influence of the inertia neglected in the 1D model.

# Code
*/
#include "grid/cartesian1D.h"
double BBingham;
#include "saint-venantB.h"

/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear:
 */

double nu_eq(double shear,double pipi){
    double nu_eq;
    nu_eq = nu*(1 + BBingham/sqrt(sq(shear) + sq(.1e-8)));
    return nu_eq;}
/**



*/

double alpha;  
char s[80];

int main() {
  X0 = -.5;
  L0 = 1.5;
  G  = 1.;
  alpha = 1;
  N  = 512*4;
  nl = 128/2;
  nu = 1.;
 /**
 a loop for the three values of Bingham parameter
 */ 
 for (B = 0.5 ; B <= 2 ; B += 0.75){ //  B = 1.25;  //.5 1.25 2
  BBingham = B;
  sprintf (s, "xML-%.2f.txt", B);
  FILE * fp = fopen (s, "w"); 
  fclose(fp); 
  sprintf (s, "shapeML-%.2f.txt", B);
  FILE * fs = fopen (s, "w");  
  fclose(fs); 
  run();
 }
}
/**
## Initialization  */
event init (i = 0) {
/**
We impose boundary condition for $h$ and $\eta$. */
  h[left] = neumann (0);
  eta[left] = neumann (0);
  u.n[left] = dirichlet(0);
  h[right] = neumann (0);
  eta[right] = neumann (0);
  /**
  We set a zero velocity at the inlet and a free outlet. */

  for (vector u in ul) {
    u.n[left] = 0;
    u.n[right] = neumann(0.);
  }
  /**
  We initialize *h*. */
  foreach(){
    h[] =   (fabs(x)<.25);
    zb[] = -(x-X0)*alpha;}
}
  /**
## Output
  We print the elevation  */
event outputfront (t += .1 , t< 5 ) {   
    double  xf=0,xe=0; 
/**
tracking the front and the end of the heap
*/   
  foreach(){
   xf = h[] > 1e-4 ?  max(xf,x) :  xf ;
   xe = h[] > 1e-4 ?  min(xe,x) :  xe ;
  } 
 /**
save them
 */  
  sprintf (s, "frontML-%.2f.txt", B);
  FILE * f = fopen (s, "w");  
  foreach()
    fprintf (f, "%g %g %g \n", fmin((x-xf),0), h[], xe-xf);   
  fclose(f);
  fprintf (stdout, "%g %g %g \n", t, xf, xe);
  sprintf (s, "xML-%.2f.txt", B);
  FILE * fp = fopen (s, "a");   
   fprintf (fp, "%g %g  \n", t, xf); 
  fclose(fp);
}
/**
save the hight the flux and the yield surface as a function of time
*/ 
  event output  (t = {0, 0.0625, 0.25, 1, 4}){
    sprintf (s, "shapeML-%.2f.txt", B);
    FILE * fp = fopen (s, "a"); 
    foreach()
      fprintf (fp, "%g %g %g\n", x, h[],t );
      fprintf (fp, "\n");
    fclose(fp);
  }
  /**

# Results
## Run

  To run 

~~~bash  
  qcc -O2 -o bingham_collapse_ML bingham_collapse_ML.c 
    ./bingham_collapse_ML  2>log
~~~

or with `make`

~~~bash
 make bingham_collapse_ML.tst; make bingham_collapse_ML/plots
 make bingham_collapse_ML.c.html ;
 source ../c2html.sh bingham_collapse_ML

~~~

## Comparisons RNSP - 1D

  Compare with Balmforth results $B=0.5$ (last $t=100$ not plotted in 2D)


~~~gnuplot B=0.5
   set xlabel "x"
   set ylabel "h(x,t)"
   p[-0.4:1.][:1.6]'../bingham_collapse_noSV/shape-0.50.txt' t '1D' w l,'shapeML-0.50.txt' t'2D' w lp
~~~

 Compare with Balmforth results $B=1.25$ (last $t=100$ not plotted in 2D)

~~~gnuplot B=1.25
   set xlabel "x"
   set ylabel "h(x,t)"
   p[-0.4:1.][:1.6]'../bingham_collapse_noSV/shape-1.25.txt' t '1D' w l,'shapeML-1.25.txt' t'2D' w lp
~~~

 Compare with Balmforth results $B=2.0$ (last $t=100$ not plotted in 2D)

~~~gnuplot  B=2
   set xlabel "x"
   set ylabel "h(x,t)"
   p[-0.4:1.][:1.6]'../bingham_collapse_noSV/shape-2.00.txt' t '1D' w l,'shapeML-2.00.txt' t '2D'w lp
~~~

 Compare with Balmforth results for position of front as function time
  
~~~gnuplot compare front position
p[:10]'xML-0.50.txt't'B=0.5 2D' w l,'../bingham_collapse_noSV/x-0.50.txt't'B=0.5 1D' w l,'xML-1.25.txt't'B=1.25 2D' w l,'../bingham_collapse_noSV/x-1.25.txt't'B=1.25 1D' w l,'xML-2.00.txt't'B=2.0 2D' w l,'../bingham_collapse_noSV/x-2.00.txt't'B=1.25 2D' w l 
~~~


Montpellier 07/17


# Links  

* see the related example in 1D  
   [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)
* see the related example in 1D  for Herschel-Bulkley 
   [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/herschel-column-noSV.c)
* see the related example in 2D for comparison
   [with hydro.h](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_NH.c)
   
# Bibliography   

* Audusse Sainte-Marie 2011

* Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet (2020)
   [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
   Volume 79, January–February 2020, Pages 233-246
   European Journal of Mechanics - B/Fluids
   https://doi.org/10.1016/j.euromechflu.2019.09.010

      
*/
