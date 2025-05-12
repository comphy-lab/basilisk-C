/** 
#   Stewartson's impulsive problem   


At time $t=0$ a infinite layer of fluid is impulsively set into motion. On the left, a constant velocity is imposed.
The Stewartson's problem consists to solve the unsteady Prandtl equations:
$$
\left\{
\begin{aligned}
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}&=0,
\\ \frac{\partial u}{\partial t}+u \frac{\partial u}{\partial x}+v \frac{\partial u}{\partial y}&=
  \frac{\partial^2 u}{\partial y^2},\\
u(x,0,t)=v(x,0,t)&=0, \\
u(x,y>0,t=0)&=1\\
v(x,y>0,t=0)&=0\\
\; \mbox{and}\; u(x,\infty,t>0)&=1.
\end{aligned}
\right.
$$




To solve this with multilayer, we reduce to a 
layer of hight $h$ of fluid which is impulsively set into motion. 
The Reynolds is so large that it is negligible compared to the depth of fluid.
On the left, a constant veloicty is imposed.


Behaviours of the solution are different depending on two asymptotic regimes: 
for large $t$ or small $x$ we recover Blasius steady solution;
$$ \delta_1 = 1.73 \sqrt{\frac{x}{Re}},\qquad  \qquad  \qquad  \tau = \frac{.3321}{\sqrt{Re x}},$$   
  
inversely, for large $x$ or small $t$, convective terms in momentum equation  are negligible (the flow is not aware of the beginning of the flow)
so Prandtl system reduces to Stokes' first problem (also called Rayleigh problem by Stewartson).
 The solution for the velocity profile can be expressed with the error function erf:
$\frac{\bar u}{u_e}=\textrm{erf}\left(\frac{\bar y}{2 \sqrt{\nu t}}\right),$
this gives      
$$  \delta_1 = 2 \sqrt{\frac{\nu t}{\pi}},\qquad  \qquad  \qquad  \tau = \frac{1}{\sqrt{\pi \nu t}},$$   
  
Transition between these two solutions occurs at $x=O(t)$, this is "Stewartson's problem" 
see Stewartson (1951) and Stewartson (1973).

## Code

We use the Multilayer solver
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h" 

double tmax,Fr,hout;
scalar tau[],delta1[];
 
/**
position of domain `X0=0`, length `L0`, no dimension `G=1` 
run with 64 points. We define `hout` from the analytical solution, 
it will be usefull for the boundary condition.
Note that the velocity is $Fr$ the Froude number. 
The Reynolds is $Re=500$ ($\nu=0.001$ and $Fr=0.5$)
*/

int main() {
  system("rm xproff.txt");
  X0 = 0;
  L0 = 1.-X0;
  G = 1;
  N = 64;
  tmax=5; 
  Fr=.5;
  nu=0.001;
  nl = 100;
  hout = 1+1.7*sqrt(L0*nu/Fr)*(Fr*Fr)/(Fr*Fr-1);
  run();
}
/** 
Initialisation,  we initialize *h* and *hc*.  
We define a scalar field *hc* to check for convergence on $h$. */

scalar hc[];

event init (i = 0)
{
  foreach(){
    zb[] = 0;
    h[] = hc[] = 1; 
   }
   boundary({zb,h,hc});

 for (vector u in ul) {
  foreach()
  u.x[] =Fr;
 }

/**
Boundary condtion, neumann for $h$ and impose  flux  and a Rieman invariant 
$u - 2 * \sqrt{gh} $ is constant  around  $h=1$
*/
  h[left] = neumann(0);
  eta[left] = neumann(0);
 
  for (vector u in ul) {
    u.n[left] = dirichlet(h[left] ? Fr /h[left] : 0.)  - 2*sqrt(G*h[]) +  2*sqrt(G*1.);  
    u.n[right] = neumann(0.);
  }
/**
Boundary condition, neumann for $u$ and impose  the height from the analytical linearized solution
*/ 
 eta[right] = dirichlet(hout); 
}
/**
 Compute friction:
*/ 
event friction (i++) {
/**
 Poiseuille friction: implicit scheme (time splitting, here commented)
 $$\frac{\partial u}{\partial t} = - C_f \frac{u}{h^2}$$
 hence, if implicited for $u$, explicited for $h$:
 $$ u^{n+1} = u^n -C_f \Delta t \frac{u^{n+1}}{(h^{n})^2}$$
*/ 
 if(nl==1){ 
  foreach()
  {    double ff = h[] < dry ? HUGE : (1. +   3*nu*dt/h[]/h[]);
    u.x[] /= ff;
    }

   foreach()
    tau[]=3*u.x[]/(h[]); 
  }
/**
Multilayer case, the derivative is twice the value in the first layer, 
*/   
else{
  vector u0 = ul[0];
  foreach()
    tau[]=2*u0.x[]/(h[]/nl);
   }
}

/**
We check for convergence. */
event logfile (t += 0.1; i <= 50000) {
  double dh = change (h, hc);
  fprintf (stderr,"coucou t= %g  -- var=  %g, Re = %lf \n", t, dh,  Fr/nu);
  if (i > 0 && dh < 1e-5)
    return 1;

/** }
And save fields

event sauv (t<=tmax;t+=.1) { */


 foreach() { 
      double z = zb[],psi=0;
      vector ue;
      ue=ul[nl-1];
      int l = 0;
      z -= (layer[0]/2)*h[];
      for (vector u in ul) {
  z += (layer[l++])*h[];
  psi += (1-u.x[]/ue.x[])*(layer[l])*h[];
      }
  delta1[]=  psi; 
}


 FILE * fq = fopen("xproff.txt","a");
 foreach() { 
      fprintf(fq,"%g %g %g %g %g %g %g %g\n", x, h[], tau[], u.x[],t,Fr,nu,delta1[]);}
fprintf(fq," \n");
}


/**
save for comaprison with analytical solution
*/
event sauf1 (t = 0.25) {
 
FILE * fp = fopen("uprof1.txt","w");
 foreach() { 
      double z = zb[];
      int l = 0;
      z -= (layer[0]/2)*h[];
     // fprintf (fp,"%g %g %g %g\n", x, z, u.x[],tau[]*z);
      for (vector u in ul) {
  z += (layer[l++])*h[];
  fprintf (fp,"%g %g %g %g \n", x, z, u.x[],tau[]*z);
      }
      fprintf (fp,"\n");
}
}

event sauf2 (t = end) {
/** 

numerical solution of 
$$2 f'''(\eta) + f(\eta) f''(\eta)=0$$
with 
$f'(0)= 0.332096$, $f'(0)=0$, and $f'(\infty)=1$

table of $\eta \qquad f'(\eta)$
*/

 FILE * fb = fopen("ublas.txt","w");
  fprintf (fb,"0.0  0.0 \n");
  fprintf (fb,"0.5  0.165905\n");
  fprintf (fb,"1.0  0.329819\n");
  fprintf (fb,"1.5  0.486845\n");
  fprintf (fb,"2.0  0.629836\n");
  fprintf (fb,"2.5  0.75134 \n");
  fprintf (fb,"3.0  0.84613\n");
  fprintf (fb,"3.5 0.913127\n");
  fprintf (fb,"4.0  0.955603\n");
  fprintf (fb,"4.5 0.979597\n");
  fprintf (fb,"5.0 0.991623\n");
  fprintf (fb,"6.0 0.9990523\n");
  fprintf (fb,"7.0 1.\n");


FILE * fp = fopen("uprof2.txt","w");
 foreach() { 
      double z = zb[];
      int l = 0;
      z -= (layer[0]/2)*h[];
     // fprintf (fp,"%g %g %g %g\n", x, z, u.x[],tau[]*z);
      for (vector u in ul) {
  z += (layer[l++])*h[];
  fprintf (fp,"%g %g %g %g \n", x, z, u.x[],tau[]*z);
      }
      fprintf (fp,"\n");
}
}

/**
 
## Run
To compile and run  :

~~~bash
   qcc -O2  -o stewartson_imp  stewartson_imp.c  -lm;  ./stewartson_imp 
~~~
 

## Plots
 
Plot of free surface and comparison with the analytical linearized solution : 
$$h = 1 + Fr^2 \frac{\delta_1^0}{Fr^2-1}$$
with $Re=500$
 
 ~~~gnuplot steady h
 set xlabel "x"
 set ylabel "h"
 p'xproff.txt' u ($5>6?$1:NaN):2 w l t'num.',1.0+1.7*sqrt(x/500)*(.5*.5)/(.5*.5-1) t'analytical',1
 ~~~

Plot of velocity profile compared to erf at time 0.25 and $x=.75$ with $\eta=y/\sqrt{\nu t}$

~~~gnuplot erf solution
set xlabel "eta"
 set ylabel "u"
 p[0:5][0:1.5] 'uprof1.txt' u ($1>0.7&&$1<.8? $2/sqrt(0.001*.25) :NaN):($3/.5) w lp,erf(x/2)
~~~


Plot of velocity profile compared to Blasius at final time  and $x<.25$ with $\eta=y\sqrt{Re}/sqrt{x}$

~~~gnuplot Blasius
set xlabel "eta"
 set ylabel "u"
p[0:8][0:1.5] 'uprof2.txt' u ($1>0.1&&$1<.25? $2/sqrt($1/500) :NaN):($3/.5) w lp,'ublas.txt' w l

~~~


Plot of $\delta_1$ as a function of $x$, for $x>u t/H$ Stokes $\delta_1=2*\sqrt{(\nu t/\pi)}$ for
 $x<u t/H$ Blasius  $\delta_1=1.732 \sqrt{x/Re}$ : 

 
 ~~~gnuplot
set xlabel "x"
set ylabel "delta1"
p'xproff.txt'  u 1:($8) w lp,'' u 1:(1.732*$1/sqrt($1*500)),'' u 1:($1>(1/2.55*.5*$5)? 2*sqrt($5*.001/pi): NaN) w l

~~~



Plot of transition between Blasius and Stokes-Rayleigh:, plot of $\sqrt{x}\tau$ as function of $u t/x$:

~~~gnuplot
 set xlabel "t/x"
 set ylabel "tau*sqrt(x)"  
 set logscale ;p[0.01:100][.1:10]'xproff.txt' u ($5*$6/$1):($1>4./256?(sqrt($7/$6)*($3/$6)*sqrt($1)):NaN) t'num.',(x<3?1/sqrt(pi*x):NaN)t'Stokes Rayleigh',0.332*(x>3? 1:NaN) t'Blasius'
~~~

## Bibliography


*  K. Stewartson. On the impulsive motion of a flat plate in a viscous fluid. ii.
 The Quarterly Journal ofMechanics and Applied Mathematics, 26(2):143–152, 1973.

* K. Stewartson.  On the impulsive motion of a flat plate in a viscous fluid.
 The Quarterly Journal ofMechanics and Applied Mathematics, 4(2):182–198, 1951
 
* [Lagree P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Boundary Layer" Cours MSF12, M2 SU

* [Lagree P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/CISM/blasius_CISM.pdf)
"A rem(a)inder on Ideal Fluid - Boundary Layer decomposition" Cours MF2A, M2 SU

* [James Lagree Legrand ](https://hal-univ-tlse3.archives-ouvertes.fr/IJLRDA/hal-01341563v1) A viscous layer model for a shallow water free surface flow, HAL

* James Lagree Le Legrand





dans le Marseille-Paris Mars 2018 
*/