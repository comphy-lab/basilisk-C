/**
# collapse of a rectangular  column of dry granular (with cohesion) on a slope,
 
## Problem
 
 A heap of dry granular (with cohesion) fluid is released along a constant slope.
 
 
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
 following $\mu(I)$ rheology, with $I=\dfrac{d}{\sqrt{p/\rho}}\dfrac{\partial u}{\partial y}$ this gives
 $\tau_{xy}= \tau_Y + (\mu( \dfrac{d}{\sqrt{p/\rho}}\dfrac{\partial u}{\partial y})p),$
 there is cohesion with additional $\tau_Y$  (Yield stress) to the classical $(\mu(I)p)$.
 
 
 
We solve this with the Multilayer Solver technique from Audusse & Sainte Marie.
We tune the solver  to consider a non newtonian viscosity function of the shear
$$\nu_{eq}=  \frac{(\tau_Y+ \mu(I)p  )/\rho}{|\frac{\partial u}{\partial y}|}$$
the `multilayer.h` is changed in `saint-venantNN.h` (all is in `saint-venantB.h`).

## Link with 1D model

If we integrate over the depth of flow $h$,  the pressure is hydrostatic $p=\rho g (h-y)$.
Integrating over the depth of flow $h$ the incompressiblility equation, we obtain
$$\frac{\partial h}{\partial t} + \frac{\partial Q}{\partial x}=0$$
where $Q=\int_0^h udy$.  For momentum it gives:
 $$\frac{\partial Q}{\partial t} + \frac{\partial }{\partial x}(\frac{Q^2}{h} +g \frac{h^2}{2})=
 - g h Z'_b - \mu g h - \tau_Y$$
 
 
 this is not simple to put $\tau_Y$,
 Note that we cannot  neglect inertia in the momentum equation. Indeed to solve for $Q$ in
 $$0+g \frac{\partial h }{\partial x}=
 - g Z'_b - \mu g  $$
with $\mu$   constant is impossible.
 
 
 
But,  if $\mu$ is function of a mean $I_m=(2/5)(d/h)Q/\sqrt{g h }$
obtained for a  Bagnold profile,  then we can solve for $Q$
in
 $$0+g \frac{h^2}{2}=
 - g h Z'_b - \mu(I_m) g h $$
 and find a kind of 1D kinetic wave with
 $Q = \frac{2 }{5} (\frac{-Z_b' - \mu_s -\partial_x h}{\Delta \mu}) \sqrt{g h} \frac{h^2}{d}$. This flux is  substituted  in mass conservation;
 $$\frac{\partial h}{\partial t} +  \frac{\partial Q}{\partial x}=0$$
 
 
 
 
 
 as in [Huppert  problem](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) or [Balmforth problem](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c). This is done
 for this Bagnold case [here](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_ondesimple_noSV.c)
 
In this file we use the Multilayer Solver to see the influence of the inertia neglected in the 1D model.

# Code
*/
#define pourtoutes(item, array) \
for(int keep = 1, \
count = 0,\
size = sizeof (array) / sizeof *(array); \
keep && count != size; \
keep = !keep, count++) \
for(item = (array) + count; keep; keep = !keep)


#include "grid/cartesian1D.h"
double tauY;
#include "saint-venantB.h"



/**
 
### Non Newtonian viscosity
 The definition of viscosity as a function of shear ($z$ is transverse if we use $x,y$ in the plane):
 $I=  \dfrac{d}{\sqrt{p/\rho}}\dfrac{\partial u}{\partial z},$
 there is cohesion with additional $\tau_Y$ to the classical $(\mu(I)p)$.
the  non newtonian viscosity function of the shear and pressure is
 $$\nu_{eq}=  \frac{(\tau_Y+ \mu(I)p  )/\rho}{|\frac{\partial u}{\partial z}|}$$
 
 $\tau_Y$ without dimension is refered as a "Bingham", this is a Yield stress
 */

double nu_eq(double shear,double press){
    double nu_eq,In;
    In=(press>0?  1./32*sqrt(sq(shear))/sqrt((press)) : 1e30);
    nu_eq = min( ( tauY + (.4 + .28*In/(In+.4)) *  press)/sqrt(sq(shear)+dry),1e30 ) ;
    return nu_eq;}
/**
 
 
### main
 
 */

double alpha,xfront;
char s[80];

int main() {
  X0 = 0;
  L0 = 8;
  G  = 1.;
  alpha = 0.;
  N  = 512*2;
  nl = 32*2;
  nu = 1;
  DT =0.01;
  FILE * f = fopen ("tauYxmax.txt", "w");
    
 /**
 a loop for the three values of Bingham parameter
 */
   double values[] = { 0, 0.01, 0.02, 0.05,  0.1 , 0.2 , 0.3, 0.5 , 1 };
   pourtoutes(double *v, values) {
   tauY=*v;
   printf("# tauY %lf\n",tauY);
  
  sprintf (s, "xML-%.3f.txt", tauY);
  FILE * fp = fopen (s, "w"); 
  fclose(fp);
  sprintf (s, "shapeML-%.3f.txt", tauY);
  FILE * fs = fopen (s, "w");  
  fclose(fs);
        
  run();
   fprintf (f, "%g %g  \n", tauY, xfront );
   }
    fclose(f);
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
    zb[] =  (x<2?-(x-2)*alpha:0);
    h[] =   (fabs(x)<2)-zb[];
    }
}
  /**
## Output
  We print the elevation  */
event outputfront (t += .1 ; t<=5) {
    double  xf=0,xe=0; 
/**
tracking the front and the end of the heap
*/   
  foreach(){
   xf = h[] > 1e-10 ?  max(xf,x) :  xf ;
   xe = h[] > 1e-10 ?  min(xe,x) :  xe ;
  } 
 /**
save them
 */
  xfront = xf;
  sprintf (s, "frontML-%.3f.txt", tauY);
  FILE * f = fopen (s, "w");  
  foreach()
    fprintf (f, "%g %g %g \n", fmin((x-xf),0), h[], xe-xf);   
  fclose(f);
  fprintf (stdout, "%g %g %g \n", t, xf, xe);
  sprintf (s, "xML-%.3f.txt", tauY);
  FILE * fp = fopen (s, "a");   
   fprintf (fp, "%g %g  \n", t, xf); 
  fclose(fp);
}
/**
save the hight the flux and the surface as a function of time
*/ 
  event output  (t+=1){
    sprintf (s, "shapeML-%.3f.txt", tauY);
    FILE * fp = fopen (s, "a"); 
    foreach()
      fprintf (fp, "%g %g %g %g \n",x,h[],zb[],t);
      fprintf (fp, "\n");
    fclose(fp);
  }
  /**

# Results
## Run

  To run 

~~~bash  
  qcc -O2 -o bingham_muI_collapse_ML bingham_muI_collapse_ML.c
    ./bingham_muI_collapse_ML  2>log
~~~

or with `make`

~~~bash
 make cohesive_muI_collapse_ML.tst; make cohesive_muI_collapse_ML/plots
 make cohesive_muI_collapse_ML.c.html ;
   
 source ../c2html.sh cohesive_muI_collapse_ML

~~~

## Comparisons RNSP - 1D

  Compare with 1D result and variable $\mu$


~~~gnuplot
   set xlabel "x"
   set ylabel "h(x,t)"
   p[:8][:1.6]'../savagestaron/log' t '1D' w l,'shapeML-0.000.txt' u 1:($2+$3) t'2D' w l,'' u 1:3 not w l linec -1
~~~

Compare    Kzerswell Balmforth results $\tau_Y=0$ $\mu=cst$.
   
~~~gnuplot cmp B&K 05
   set term png ;   set output 'bk.png';
   
   Xc=520
   L1=907-Xc
   H1=640-185
   Yc=185
   unset tics
   set key right top
   plot [0:1400]  '../Img/Balmforth_Kerswell05.png' binary filetype=png with rgbimage not,\
   '../savagestaron/log' u (Xc+($1-2)*L1*.4):(Yc+$2*(H1)) w l not,\
   'shapeML-0.00.txt' u (Xc+($1-2)*L1*.4):(Yc+($2)*(H1)) t'2D' w l
~~~
 
 
# with cohesion  

Compare with and without 

~~~gnuplot
   reset
   set multiplot layout 2,2
   set xlabel "x"
   set ylabel "h(x,t)"
   p[:8][:1.6]'../savagestaron/log' t '1D' w l,\
   'shapeML-0.000.txt' u 1:($2+$3) t'2D 0.0' w l,\
   'shapeML-0.100.txt' u 1:($2+$3) t'2D 0.1' w l,\
   '' u 1:3 not w l linec -1
   
   p[:8][:1.6]'../savagestaron/log' t '1D' w l,\
   'shapeML-0.000.txt' u 1:($2+$3) t'2D 0.0' w l,\
   'shapeML-0.200.txt' u 1:($2+$3) t'2D 0.2' w l,\
   '' u 1:3 not w l linec -1
   
   p[:8][:1.6]'../savagestaron/log' t '1D' w l,\
   'shapeML-0.000.txt' u 1:($2+$3) t'2D 0.0' w l,\
   'shapeML-1.000.txt' u 1:($2+$3) t'2D 1.0' w l,\
   '' u 1:3 not w l linec -1
   
   p[:8][:1.6]'../savagestaron/log' t '1D' w l,\
   'shapeML-0.000.txt' u 1:($2+$3) t'2D 0.0' w l,\
   'shapeML-0.500.txt' u 1:($2+$3) t'2D 0.5' w l,\
   '' u 1:3 not w l linec -1
     unset multiplot
~~~



  Compare with Balmforth results for position of front as function time
  
   
~~~gnuplot compare front position
   reset
   mu=.4
   p[:5]'xML-0.000.txt' t'B=0.000 2D' w l,\
        'xML-0.001.txt' t'B=0.001 2D' w l,\
        'xML-0.010.txt' t'B=0.010 2D' w l,\
        'xML-0.020.txt' t'B=0.020 2D' w l,\
        'xML-0.050.txt' t'B=0.050 2D' w l,\
        'xML-0.100.txt' t'B=0.100 2D' w l,\
        'xML-0.200.txt' t'B=0.200 2D' w l,\
        'xML-0.300.txt' t'B=0.300 2D' w l,\
        'xML-1.000.txt' t'B=1.000 2D' w l,\
         (2+x*(2-x*mu/2))*(x<2/mu?1:NaN) t'K','../savagestaron/x.txt' t'num SH'
~~~


 Value of the final position ($x_{front}-2$) as function of cohesion.
 Plain line: comparison wit full 1D constant $\mu$ 
 [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_savagehutter.c]()
 
   
~~~gnuplot
   set xlabel "tauY"
   set ylabel "extend"
   p'tauYxmax.txt' t '2D','../cohesive_savagehutter/tauYxmax.txt' t'1D' w l 
~~~
   

 
 
 
 


## Links

* see the related example in 1D  
   [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)
* [Balmforth problem](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c).
* for a Bagnold case [here]((http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_ondesimple.c))
* [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapse_ML.c]() 2D viscous `multilayer.h`
* [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapse_NH.c]() 2D viscous with `hydro.h`
* [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]() 1D with `saintvenant.h`
   
   
# Bibliography   

* Audusse Sainte-Marie 2011

* Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet (2020)
   [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
   Volume 79, January–February 2020, Pages 233-246
   European Journal of Mechanics - B/Fluids
   https://doi.org/10.1016/j.euromechflu.2019.09.010

      
*/
