/**
#  Saint Venant / Shallow Water "transcritical flow on a bump", 

 
## Problem
 
For a given bump at the bottom of a river, find the free surface deformation.
 
 We test the analytical linearized solution :
 $$h = 1 + \frac{z_b(x)}{Fr^2-1}$$
 
 
 
 
## Shallow-Water Saint-Venant model
The classical Shallow Water Equations in 1D with a topography and no friction. 
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x Q=0\\
  \partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= - gh \partial_x z_b
        \end{array}\right. 
$$
with  a given topography $z_b(x)$, and we vary the Froude number $Fr$. As we are without dimension, $g=1$, the depth is measured with unit, and the initial velocity will be $Fr$,  in practice $z_b=\alpha exp(-x^2)$ with $\alpha=0.1$

Depending on the Froude number, the free surface will present various deformations.

## Code
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"
double tmax,Fr;
/** 
free output
*/
u.n[left]   = neumann(0);
u.n[right]  = neumann(0); 
h[left]   = neumann(0);
h[right]  = neumann(0);  
/** 
start by a constant flow
*/
event init (i = 0)
{
  foreach(){
    zb[] =  .1*exp(-x*x);
    h[] =  1 - zb[];
    u.x[] = Fr;
   }
   boundary({zb,h,u});

#ifdef gnuX
  printf("\nset grid\n");
#endif    
}
/**
position of domain `X0`, length `L0`, no dimension `G=1` 
run with 512 points (less is not enough). Do a loop on Froude
*/
int main() {
  X0 = -10.;
  L0 = 20.;
  G = 1;
  N = 512;
  tmax=50;
  
  double FF[4]={ .2 , .4 , 0.65  , 2};  // loop on Froude
  
  for (int iF=0;iF<=3;iF++){ 
    Fr = FF[iF];
  run();}
}
/**
Output in gnuplot if the flag `gnuX` is defined
*/
#ifdef gnuX
event plot (t<tmax;t+=0.2) {
    printf("set title 'Ressaut en 1D --- t= %.2lf  Fr= %g '\n"
      "Fr= %g ;Z(x)=0.1*exp(-x*x); h(x)=1+Z(x)/(Fr*Fr-1);   ; "
      "p[%g:%g][-.25:2.25]  '-' u 1:($2+$4) t'free surface' w l lt 3,"
      "'' u 1:($2*$3) t'Q' w l lt 4,\\\n"
      "'' u 1:4 t'topo' w l lt -1,\\\n"
      "h(x) t 'theo h(x) w p'\n",
           t,Fr,Fr,X0,X0+L0);
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}
/**
Output at the end if not defined
*/
#else
event end(t=tmax ) {
    foreach()
    printf ("%g %g %g %g %g %g \n", x, h[], u.x[], zb[], t, Fr);
}
#endif
/**
 
## Run
To compile and run with gnuplot:

~~~bash
  qcc -DgnuX=1   -O2   -o bump_trans bump_trans.c  -lm;  ./bump_trans | gnuplot
~~~

To compile and  plot gnuplot at the end :

~~~bash
  qcc  -O2   -o bump_trans bump_trans.c  -lm 
~~~


## Plots
 
Plot of free surface and comparison with the nalaytical linearized solution : 
$$h = 1 + \frac{z_b(x)}{Fr^2-1}$$
 
 
### sub critical case $Fr=0.2$
 
In this case, the water depth presents a depression, flowx is accelerated over the bump  

~~~gnuplot  result, free surface (blue) and bottom (black)  Fr=0.2 
 set xlabel "x"
  h(x) = (1+Z(x)/(Fr*Fr-1))
 Z(x)=0.1*exp(-x*x); 
 Fr=0.2
 p [:][-.5:1.5] 'out' u 1:($6==Fr? $2: NaN) t'free surface' w lp lc 3, '' u 1:4 t'zb' w l lc -1, '' u 1:($6==Fr? (1+$4/($6*$6-1)) : NaN)  w l lc 1 t 'analytic'
 
~~~

### sub critical case $Fr=0.4$

In this case, the water depth presents a depression, flowx is accelerated over the bump

~~~gnuplot  result, free surface (blue) and bottom (black)  Fr=0.4 
 set xlabel "x" 
 Fr=0.4
 p [:][-.5:1.5] 'out' u 1:($6==Fr? $2: NaN) t'free surface' w lp lc 3, '' u 1:4 t'zb' w l lc -1, '' u 1:($6==Fr? (1+$4/($6*$6-1)) : NaN)  w l lc 1 t 'analytic'
~~~

 
### Trans critical case $Fr=0.65$,

in this case $Fr=1$ at the top of the bump, then the flow in supercritical in the lee side. note the recompression hydrolic jump.

~~~gnuplot  result, free surface (blue) and bottom (black)  Fr=0.65 
 set xlabel "x" 
 Fr=0.65
 p [:][-.5:1.5] 'out' u 1:($6==Fr? $2: NaN) t'free surface' w lp lc 3, '' u 1:4 t'zb' w l lc -1, '' u 1:($6==Fr? (1+$4/($6*$6-1)) : NaN)  w l lc 1 t 'analytic'
~~~


### Super critical case $Fr=2$,  

In this case, the water depth presents a hump, flowx is decelerated over the bump
  
~~~gnuplot  result, free surface (blue) and bottom (black)  Fr=0.65 
 set xlabel "x" 
 Fr=2
 p [:][-.5:1.5] 'out' u 1:($6==Fr? $2: NaN) t'free surface' w lp lc 3, '' u 1:4 t'zb' w l lc -1, '' u 1:($6==Fr? (1+$4/($6*$6-1)) : NaN)  w l lc 1 t 'analytic'
~~~



##Bibliography
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC



Version 1: Montpellier 2017

*/
