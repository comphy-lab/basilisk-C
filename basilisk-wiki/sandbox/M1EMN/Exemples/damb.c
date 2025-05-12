/**
#    Ideal fluid "dam break problem", 

## Problem
 
What happens if the dam which retains a lake is suddenly removed? That is the dam break problem ("Problème de la rupture de Barrage").
Here an animation of the final result:

 ![animation of the dam break](damb/animate.gif)
 

## Equations
 
 
The classical Shallow Water (Saint-Venant) Equations in 1D with no friction:
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x Q=0\\
	\partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= 0
        \end{array}\right. 
$$
with at  $t=0$, a lake on the left $h(|x|<0,t=0)=1$ and zero water  at the righjt$h(|x|>0,t=0)=0$ for $x>0$.
As the flow is not viscous, the convective term is important. The solution is the Ritter 1862 solution with simple waves.
 
This is solved with [http://basilisk.fr/src/saint-venant.h]()

## Code*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax;

u.n[left]   = neumann(0);
u.n[right]  = neumann(0); 
 
/** 
start by a lake on the left and dry on the right, the soil is flat $z_b=0$
*/
event init (i = 0)
{
  foreach(){
    zb[] =  0;
    h[] =  (x<0)*1  +.0;
    u.x[] =0;
   }
}

/**
 position of domain `X0`, length `L0`, no dimension `G=1`,  
run with 512 points (less is not enough)
*/
int main() {
  X0 = -5.;
  L0 = 10.;
  G = 1;
  N = 512/4;
  tmax=3;
  DT = HUGE;
  run();
}
/**
The following event is just to do a nice animation with gnuplot.
Output in gnuplot if the flag `gnuX` is defined, put in a `gif` if not
*/
event plot (t<tmax;t+=0.01) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
#ifdef gnuX
#else
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
#endif
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Ressaut en 1D --- t= %.2lf '\n"
      "h(x)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t));t= %.2lf  ; "
      "p[%g:%g][-.5:2]  '-' u 1:($2+$4) t'free surface' w lp lt 3,"
      "'' u 1:($2*$3) t'Q' w l lt 4,\\\n"
      "'' u 1:4 t'topo' w l lt -1,\\\n"
      "'+' u (0):(0.296296296296296) t'theo Qmax ' , h(x) t 'theo h(x)'\n",
           t,t,X0,X0+L0);
    foreach()
    fprintf (fp,"%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    fprintf (fp,"e\n\n");

    if(t==tmax){
        fprintf (fp,"! cp animate.gif a2.gif \n");}
}
/**
Output at the end  
*/
event end(t=tmax ) {
    foreach()
    fprintf (stderr,"%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    
}

/**

## Run
To compile and run with gnuplot (and X11 window):

~~~bash
  qcc -DgnuX=1   -O2   -o damb damb.c -lm;  ./damb  
~~~

## Plots
 
Plot of the surface and flux at time $t=t_{max}$:
Indeed the analytical solution 
$$ x<-t, \;\; h(x,t)=1$$
$$ -t<x<2t,  \;\; h(x,t)= (\frac{2}{3}(1-\frac{x}{2t}))^2$$
$$ 2t < x , \;\;  h(x,t)=0$$
 and the numerical one are superposed, here at time $t=3$

 
~~~gnuplot  result, free surface (blue) and bottom (black)
set xlabel "x"
h(x)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t))
t=3
 p [:][-1:3] 'log' t'free surface' w p lc 3,'' u 1:($3*$2) w l t'Q', '' u 1:4 w l lc -1,h(x) w l lc 1
 
~~~

## Links
 
* see the same
  [with turbulent friction](http://basilisk.fr/sandbox/M1EMN/Exemples/damb_dressler.c)

* see non viscous dam break with [standard C](http://basilisk.fr/sandbox/M1EMN/Exemples/svdb.c)
and with [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c) (this file)
 
* see all the viscous collapse examples ([with laminar friction](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c))
 
 * a version in [python](https://colab.research.google.com/drive/1960Q9Cgu9anAv6wfB9MblsBFo0BvWdrO) of this file

## Bibliographie
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC





Version 1: Montpellier 2017

*/

   

 
 
 
