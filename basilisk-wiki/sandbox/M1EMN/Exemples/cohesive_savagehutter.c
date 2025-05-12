/**
# Savage Hutter Collapse of granular heap with friction and cohesion

## Problem

This is the problem of the collapse of a initial rectangular heap of cohesive grains.
It is a kind of Dam Break problem.

 
 ![animation of the collapse](cohesive_savagehutter/animate.gif)
 

## Equations

We ccompute the  "Shallow Water" or "Savage Hutter equations" or "depth averaged equations" or "Saint Venant" simplified equations with a  constant basal friction with cohesion and friction $\tau_Y+\mu P$.
Here $P=\rho gh$  .
The system is
 $$
 \left\{\begin{array}{l}
\frac{\partial }{\partial t}  h \; +\; \frac{\partial }{\partial x} uh=0\\
 \frac{\partial }{\partial t} hu +
 \frac{\partial }{\partial x}  \dfrac{(hu)^2}{h} +
 \frac{\partial }{\partial x}g\dfrac{h^2}{2}
 = - gh \frac{\partial }{\partial x} Z-(\frac{\tau_Y}{\rho} +  \mu g h)\frac{u}{|u|}
 \end{array}\right.
 $$
 
 It is may be convinient to define $\ell_c = \frac{\tau_Y}{\mu \rho g }$, the friction is in its definition, so that the equation rescale with $\mu$ (see Balmforth & Kerswell 05  and Kerswell 05).
 
 
 
 We solve the problem by splitting, first (note that we use $Q=hu$)
 $$\frac{h^*-h^n}{\Delta t} + \frac{\partial Q^n}{\partial x} =0, \text{ and }
 \frac{Q^*-Q^n}{\Delta t}+ \frac{\partial }{\partial x}\frac{Q^n}{h^n}
 +
 \frac{\partial }{\partial x}g\dfrac{h^{n2}}{2}
 = - g \frac{\partial }{\partial x} Z$$

is solved with [http://basilisk.fr/src/saint-venant.h]().
 
Second, friction gives
 $$ \frac{u^{n+1}-u^*}{\Delta t} =- \mu g (1+ \frac{\ell_c}{h}) \frac{u}{|u|}$$
is solved like in [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]()
 
 
 

## Code

*/
 
#include "grid/multigrid1D.h"
#include "saint-venant.h"
#define pourtoutes(item, array) \
for(int keep = 1, \
count = 0,\
size = sizeof (array) / sizeof *(array); \
keep && count != size; \
keep = !keep, count++) \
for(item = (array) + count; keep; keep = !keep)
/**
The domain is 7 long, the height is the unit of length.
The problem is without dimensions. The problem is
solved in one dimension. 
*/

/** 
Wall symmetry at the left
Neumann conditions at the exit
*/
u.n[left] = dirichlet(0);
h[left] = neumann(0);
u.n[right] = neumann(0);
h[right] = neumann(0);
/**
Main and parameters
*/
double Ltas,Htas,tmax,xfront,tauY;
char s[80];
int main()
{
  X0 = 0.;
  L0 = 7;
  G = 1.; 
  N = 1024;
//  source = 0;
  Htas = 1;
  Ltas = 2.; 
  tmax=6;
  DT=0.001;
  FILE * f = fopen ("tauYxmax.txt", "w");
  sprintf (s, "shape-%.3f.txt", tauY);
  FILE * fs = fopen (s, "w");
  fclose(fs);
  sprintf (s, "x-%.3f.txt", tauY);
  FILE * fp = fopen (s, "w");
  fclose(fp);
  double values[] = { 0, 0.01, 0.02, 0.05,  0.1 , 0.2 , 0.3, 0.5 , 1 };
  pourtoutes(double *v, values) {
  tauY=*v;
  printf("# tauY %lf\n",tauY);
    
  run();
  fprintf (f, "%g %g  \n", tauY, xfront );
  }
}
/**
 The initial conditions are  a given heap of length $L_{tas}$ (the double by symetry) of height $H_{tas}$.
to left. */
event init (i = 0){
  double h0;
/**
initial heap
*/ 
   h0=Htas;
   foreach(){
    zb[] = 0;	
    h[] = (x < Ltas) ? h0 : 0;
    u.x[]= 0;}
}

/**
## Friction

We use a simple implicit scheme to implement coulomb bottom friction i.e.
$$\frac{d\mathbf{u}}{dt}=-(\frac{\tau_Y}{\rho h}+\mu g) \frac{\mathbf{u}}{|\mathbf{u}|}$$
with $\mu$ is constant (it may be funtion  of a mean $I$, written with $Q/h^{3/2}$).
 
 We have here a $\mu_{eff}=(\frac{\tau_Y}{\rho h}+\mu g)=\mu g (\frac{\ell_c}{h}+1)$
 Of course we have splitted, first the Rieman Problem to solve the waves, second the friction problem.

 
For friction, we noticed that there is an exact solution
  (see granular sand glass with friction)
 of this splited problem,
 defining the norm of velocity
 $U=|\mathbf{u}|$ and $\overrightarrow{T}=\frac{\mathbf{u}}{|\mathbf{u}|}$
 the equation
 $\frac{d (U \overrightarrow{T})}{dt}$=$-\mu_{eff} g \overrightarrow{T}$
 is solved explicitely for $t'>t$,
 obviously, as $\frac{d \overrightarrow{T}}{ds}=\frac{\overrightarrow{N}}{R}$
 we have just to solve:
 $$\frac{d U}{dt}=-\mu_{eff} g$$
 it is linear and  the solution is $U(t')=U(t)- \mu_{eff}  g (t'-t)$ down to stop, so:
 $$U(t+\Delta t) = max(U(t) - \mu_{eff}  \Delta t g,0)$$
 This gives a real stop. But finally it is not so acurate.
 
*/
event coulomb_friction (i++) {
  double mu,U;
  foreach() {
    U=norm(u);
    mu=.4;
      if(U>0){
          foreach_dimension()
          u.x[] = max(U -dt *( mu * G + tauY/h[]),0)*u.x[]/U;}
  }
  boundary ({u.x});
}
/**
## output
 
*/

event outputfront (t += .1 ) {
    double  xf=0,xe=0;
    static FILE * ff = fopen("x.txt", "w");
    /**
     tracking the front and the end of the heap
     */
    foreach(){
        xf = h[] > 1e-20 ?  max(xf,x) :  xf ;
        xe = h[] > 1e-20 ?  min(xe,x) :  xe ;
    }
     fprintf (ff, "%g %g %g \n", t, xf , xe);
      xfront = xf;
    sprintf (s, "x-%.3f.txt", tauY);
    FILE * fp = fopen (s, "a");
    fprintf (fp, "%g %g  \n", t, xf);
    fclose(fp);
}

/**
 generate the gif
 */

event output (t += .02; t < tmax) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
#ifdef gnuX
    
#else
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .5\n");
#endif
    fprintf (fp,"set title ' collapse tauY=%.3lf--- t= %.2lf '\n"
             "p[0:5][-.5:2]  '-' u 1:($2) t'free surface' w lp lt 3,0 not w l linec -1\n",tauY,t);
  foreach()
      fprintf (fp, "%g %g \n", x, h[]);
  fprintf (fp,"e\n\n");
}
  event outputlog (t += 1; t < tmax) {
   foreach()
     fprintf (stderr, "%g %g %g \n", x, h[], t);
   fprintf (stderr, "\n");
   char s[80];
   sprintf (s, "shape-%.3f.txt", tauY);
   FILE * fp = fopen (s, "a");
    foreach()
      fprintf (fp, "%g %g %g %g \n",x,h[],zb[],t);
      fprintf (fp, "\n");
    fclose(fp);
}
/**
# Results

to run and plot with gnuplot during the run (note the `-DgnuX`)

~~~bash
qcc -g -O2 -DTRASH=1 -Wall -DgnuX=1 -o cohesive_savagehutter cohesive_savagehutter -lm
~~~
 
 with `make`
 
~~~bash
 make cohesive_savagehutter.tst;make cohesive_savagehutter/plots; make cohesive_savagehutter.c.html
 
 
 source ../c2html.sh cohesive_savagehutter
~~~

plot with gnuplot of the results, dynamics is not very good, but deposit is not so bad for the dry case.

~~~gnuplot Fluid depth profile, Contact Dynamics vs Savage Hutter
set xlabel 'x'
set ylabel 'h'
set xlabel "x"
set ylabel "h(x,t)"
d=0.005
h0=0.149
p[0:6][0:1.5]'../../granular_column/ShapeTime.A-01.dat' u (($1*d)/h0):(($2*d)/h0) t'DCM'w l,\
   'shape-0.000.txt' w l t 'h'
~~~

 Comparison with Balmforth & Kerswell 05
 who rescaled to have a canonical problem in the case of Constant friction:
 
 ~~~gnuplot cmp B&K 05
 reset
 set term png ;   set output 'bk.png';
 
 Xc=520
 L1=907-Xc
 H1=640-185
 Yc=185
 mu=0.4
 unset tics
 set key center top
 plot [0:1300]  '../Img/Balmforth_Kerswell05.png' binary filetype=png with rgbimage not,\
 'shape-0.000.txt' u (Xc+($1-2)*L1*mu):(Yc+$2*(H1)) w l not
 ~~~
 
 
 ## effect of $\tau_Y$
 
 
 ~~~gnuplot
 reset
 set multiplot layout 2,2
 set xlabel "x"
 set ylabel "h(x,t)"
 p[:8][:1.6]'shape-0.000.txt' t 'tauY=0.000' w l,\
 'shape-0.010.txt' u 1:($2+$3) t'tauY=0.010' w l
 
 p[:8][:1.6]'shape-0.000.txt' t 'tauY=0.000' w l,\
 'shape-0.100.txt' u 1:($2+$3) t'tauY=0.100' w l
 
 p[:8][:1.6]'shape-0.000.txt' t 'tauY=0.000' w l,\
 'shape-0.200.txt' u 1:($2+$3) t'tauY=0.200' w l
 
 p[:8][:1.6]'shape-0.000.txt' t 'tauY=0.000' w l,\
 'shape-0.500.txt' u 1:($2+$3) t'tauY=0.500' w l

 unset multiplot
 ~~~
 
 
 
Comparison with Balmforth & Kerswell 05 for the runout proposed by Kerswell 05
 
 $u =  \frac{dx}{dt}=2 - \mu t$ so that $x  = 2 (t - \mu t^2/4)$ and
  $t_{max} =\frac{2}{\mu}$ and the variation of positio of the front is
 $\Delta x_{max} =\frac{2}{\mu}$

 position of the front as a funtion of time for various $\tau_Y$
~~~gnuplot
 reset
 set xlabel "t"
 set ylabel "xmax-xinit"
 mu =.4
 set key bottom
 p [][0:]'x-0.000.txt'u 1:($2-2) w l,'x-0.010.txt'u 1:($2-2) w l ,\
 'x-0.020.txt'u 1:($2-2) w l ,'x-0.050.txt'u 1:($2-2) w l ,'x-0.100.txt'u 1:($2-2) w l ,\
 (x*(2-x*mu/2))*(x<2/mu?1:NaN) t'caract pred.'
~~~
 
 We note that even for $\tau_Y$ the result is not exactly the parabola predicted by Kerswell.
 We plot now the runaout as function of $\tau_Y$

 
~~~gnuplot
 set xlabel "tauY"
 set ylabel "(x max- x init)"
 p'tauYxmax.txt' u ($1):($2-2)  w lp not
~~~
 
 
 ~~~gnuplot
 set xlabel "tauY"
 set ylabel "x max"
 set logscale x
 p'tauYxmax.txt' u ($1):($2-2)  w lp not
 ~~~
 
 
 We plot here the extent measured with $\ell_c$ function of $h_0/\ell_c$
 
 ~~~gnuplot
 reset
 set xlabel "h_0/l_c"
 set ylabel "(x max-x init)/l_c"
 val = `awk 'END{print $2}' ../cohesive_savagehutter/x-0.000.txt`
 p[0:200]'tauYxmax.txt' u (1/$1):(($2-2)/$1)  w lp not, val*x t' dry'
 ~~~
 
 
 
*/ 


/**
## Links
 
* Ideal fluid dam break [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c)
 
* Savage Hutter  [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]()
 
* [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass_muw.c]() granular sand glass with friction
 
* a version in [python](https://colab.research.google.com/drive/19oLUumzwTzTiRxnF3VgHvp0j8OjE4MlD) of this file with no cohesion
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapse_ML.c]() the RNSP 
 
 
## Bibliography
 
* Savage Hutter
 
* Larrieu Staron Hinch
"Raining into shallow water as a description of the collapse of a column of grains"
J. Fluid Mech. (2006), vol. 554, pp. 259–270.
 
* Balmforth & Kerswell 05 "Granular collapse in two dimensions"
 J. Fluid Mech. (2005), vol. 538, pp. 399–428
 
 *  R. R. Kerswell
 Dam break with Coulomb friction: "A model for granular slumping?"
PHYSICS OF FLUIDS 17, 057101 2005


v1: Montpellier 11 juillet 15
*/
