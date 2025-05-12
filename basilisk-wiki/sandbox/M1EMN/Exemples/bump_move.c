/**
#    Saint Venant "fluvial flow over an erodible  bump", 

 
 
 
## Problem

How moves a dune in a subcritical flow?
 
 
## Model Equations : Saint-Venant + Exner
 
 The classical Shallow Water Equations in 1D, we first neglect the viscous dissipation
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x Q=0\\
  \partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= - gh \partial_x z_b
        \end{array}\right. 
$$
coupled with Exner equation
$$ \partial_t h+\partial_x q_s=0 $$
where the flux of sediments is $q_s= Q_0( \frac{u}{h} - \tau_{threshold})_+$

with  a given in initial topography $z_b(x)$. 
As we are without dimension, $g=1$, the depth is measured with unit, and the initial velocity will be $Fr$,  
in practice $z_b=\alpha exp(-x^2)$ with $\alpha$ small in order to linearize.

##Code
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax,Fr,Cf,Q0,tauth,alpha;
scalar tau[],qs[];
/**
 Impose a characteristic BC, outcoming waves leave the domain $u-2\sqrt{gh}$ is constant
*/

u.n[left]   = u.x[] - 2*(sqrt (G*h[]) - sqrt(G*1));
u.n[right]  = neumann(0);

h[left]   = dirichlet(1) ;  
h[right]  = neumann(0);

zb[left]   = dirichlet(0);
zb[right]  = neumann(0);

qs[left]   = dirichlet(0);
qs[right]  = neumann(0);

/**
position of domain `X0`, length `L0`, no dimension `G=1` 
run with 512 points (less is not enough)
*/

int main() {
  X0 = -2.5;
  L0 = 10.;
  G = 1;
  N = 128*8;
  tmax=1000;
  Cf=0.01*0;
  Q0=0.01;
  alpha=.0125;
  tauth=0.05*0;
  Fr=.2;
  run();

}
/** 
Initialisation, 
start  by a bump $z_b$, initial linearized velocity and linearized height
$$
u = Fr (1 + \frac{z_b}{1-Fr^2}),\;\;\;
h = (1 - \frac{z_b}{1-Fr^2})\;\;\;
$$
this is the $O(\alpha)$ solution
*/
event init (i = 0)
{
  foreach(){
    zb[] =  alpha*exp(-x*x);
    h[] =  1 - zb[]/(1-Fr*Fr);
    u.x[] = Fr*(1+zb[]/(1-Fr*Fr));
   }
   boundary({zb,h,u});

#ifdef gnuX
  printf("\nset grid\n");
#endif    
}
/**
 Compute friction:
*/ 
event friction (i++) {
/**
  NOTE that in fact we put no friction here in the flow which remains "ideal"
 Nevertheless, we
 computate   an equivalent "friction" at the wall  $\tau = \frac{u}{h}$
 which will drive the sediments
*/
  foreach(){
    tau[]=u.x[]/h[];
    }
}
/**
Exner equation, first computation of the sediment flux $q_s$ as a function of the previous  shear stress.
We take a simple law with a threshold :
$$q_s=Q_0(\frac{u}{h} - \tau_{threshold})_+$$
 we even take tauth =0
*/
event exner (i++){
 foreach()
   { 
    qs[]= (tau[] > tauth ? tau[]-tauth : 0 );
  }
  boundary ({qs}); 
/** 
discretisation is just the simple downstream first order derivative $\partial_x q_s = \frac{q_{si}-q_{si-1}}{\Delta x} + O(\Delta x)$
*/
  foreach(){
    zb[] = (x> X0+.1 ? zb[] - Q0*(dt/Delta)*(qs[0]-qs[-1])  :0);
   }
 boundary ({zb}); 
}
/**
Output in gnuplot if the flag `gnuX` is defined
*/
#ifdef gnuX
event plot (t<tmax;t+=1) {
    printf("set title 'Moving dune 1D --- t= %.2lf  Fr= %g '\n"
      "Fr= %g ;Z(x)=%g*exp(-x*x); h(x)=1+Z(x)/(Fr*Fr-1);   ; "
      "p[%g:%g][-.05:.05]  '-' u 1:($2+$4-1) t'free surface' w l lt 3,"
      "'' u 1:($3) t'Qs' w l lt 4,\\\n"
      "'' u 1:4 t'actual topo' w l lt -1,\\\n"
      "'' u 1:5 t'u' w l lt -1,\\\n"
      "Z(x-%g) t 'theo Z(x-ct) 'w l linec 1,\\\n"
      " Z(x)  t 'init Z(x) 'w l linec 2\n",
           t,Fr,Fr,alpha,X0,X0+L0,2*t*Q0*Fr/(1-Fr*Fr));
    foreach()
    printf ("%g %g %g %g %g %g\n", x, h[], qs[] , zb[], u.x[], t);
    printf ("e\n\n");
}
/**
Output at the end if not defined
*/
#else
event end(t=tmax ) {
    foreach()
    printf ("%g %g %g %g %g %g  %g\n", x, h[], u.x[], zb[], t, Fr,qs[]);
}
#endif
/**
 
## Run
To compile and run with gnuplot:

~~~bash
  qcc -DgnuX=1   -O2   -o bump_trans bump_move.c  -lm;  ./bump_move | gnuplot -persist
~~~

To compile and  plot with gnuplot at the end :

~~~bash
  qcc  -O2   -o bump_move bump_move.c  -lm ;  ./bump_move > out
~~~


## Plots
 
Plot of free surface and comparison with the analytical linearized solution : 
$$h = 1 + \frac{z_b(x)}{Fr^2-1}$$
 
 
Sub critical case $Fr=0.2$, 
 note that if we are over threashold  
 $$\partial_x q_s= \partial_x  (\frac{u}{h} )$$
gives by linearisation
as $\frac{1}{h} = 1 - \frac{z_b(x)}{Fr^2-1}$ and  $u = Fr (1 + \frac{z_b}{1-Fr^2})$:
 $$\partial_x q_s=   (\frac{2 Q_0 Fr}{1-Fr^2} ) \partial_xz_b $$
so that the topography moves at the constant velocity $c_0$,
$\partial_t z_b +  c_0 \partial_xz_b=0$
with $c_0=(\frac{2 Q_0Fr}{1-Fr^2})$
 
~~~gnuplot  result, the computed bump and teh alaytical linearized solution,
set xlabel "x" 
Q0=0.01;
alpha=.0125 
Fr=.2 
tmax=1000
Z(x)=alpha*exp(-x*x)
p'out' u 1:4 t'comp. pos. at time tmax',Z(x-2*tmax*Q0*Fr/(1-Fr*Fr)),Z(x) t 'initial position'
 
~~~

## Conclusion
 
 In a fluvial flow,
the bump moves slowly without changing its shape (up to order two effects and numerical errors).


## Bibliography
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC


*/
