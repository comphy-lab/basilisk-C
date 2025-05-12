/**
# Resolution of KPZ equation in 1D
A classic model for the evolution of the profile of a growing interface 
is the Kardar–Parisi–Zhang (KPZ) equation :
$$\frac{\partial h}{\partial t} =
-v_0 \frac{\partial h}{\partial x} +
 \frac{\lambda}{2} (\frac{\partial h}{\partial x})^2 + \nu \frac{\partial^2 h}{\partial x^2}
 +v_c$$
without here the stochastic source.
$v_c$ is a constant feeding, $\lambda$ an effect of curvature, $\nu$ a diffusion effect, $v_0$ an extra slope effect.
We propose here an explicit 1D resolution of KPZ equation compared with an analytical solution.

##Code
mandatory declarations:
*/
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the height of interface
its O(Delta) derivative and it O(Delta^2) derivative, time step
*/
scalar h[];
scalar hp[];
scalar hp2[]; 
double dt;
/**
Main with definition of parameters
*/
int main() {
  L0 = 10.;
  X0 = -L0/2;
  N = 200;
#define nu 1.    
  DT = (L0/N)*(L0/N)/nu/10 ;
#define lambda 1. 
#define EPS 0.1  
#define v0 1.
#define vc 2.
  run();
}
/** 
initial elevation: a "triangle"
*/
event init (t = 0) {
  foreach()
    h[] =   (1-fabs(x))*(fabs(x)<1);
  boundary ({h});
  }
/** 
print data
*/
event printdata (t += 0.1; t < 1) {
  foreach()
    fprintf (stdout, "%g %g %g\n", x, h[], t);
  fprintf (stdout, "\n\n");
}
/** integration 
*/
event integration (i++) {
  double dt = DT;
  scalar dh[];
/**
finding the good next time step
*/
  dt = dtnext (dt);
/**
  $O(\Delta)$ derivative
*/
  foreach()
    hp[] =  ( h[1,0] - h[0,0] )/Delta;
  boundary ({hp});
/**
  centered derivative $O(\Delta^2)$
*/ 
  foreach()
    hp2[] =  ( h[1,0] - h[-1,0] )/2/Delta;
  boundary ({hp2});
/** 
explicit step
$$\frac{ (h(x+\Delta x) - h(x))/\Delta  -(h(x) - h(x-\Delta x))/\Delta }{\Delta  } \simeq \frac{\partial^2 h}{\partial x^2}$$
*/  
  foreach()
    dh[] =  ( hp[] - hp[-1,0] )/Delta;
  boundary ({dh});
/** 
update viscous term
*/
  foreach()
    h[] += dt*nu*dh[];
  boundary ({h});
/** 
update uniform growth
*/  
  foreach()
    h[] += dt*vc;
  boundary ({h});
/** 
update normal effect 
*/
  foreach()
    h[] += -dt*v0*hp2[];
  boundary ({h});
/** 
update curvature effect
*/
  foreach()
    h[] += dt*(lambda/2.*hp2[]*hp2[]);
  boundary ({h});
  
}
/**
## Run
Then compile and run:

~~~bash
qcc  -g -O2 -DTRASH=1 -Wall  kpz.c -o kpz
./kpz > v.out
~~~

or with `make` 

~~~bash
 make kpz.tst;make kpz/plots    
 make kpz.c.html ; open kpz.c.html 
~~~


## Results
The analytical (for $v_0=0$) solution is 
$$h(x,t) = 2 \log \left(-\frac{1}{2} e^{-x/2} \left(e^{\frac{t}{4}+\frac{1}{2}}
   \text{erf}\left(\frac{t-x}{2 \sqrt{t}}\right)-e^{\frac{t}{4}+\frac{1}{2}}
   \text{erf}\left(\frac{t-x+1}{2 \sqrt{t}}\right)-e^{x/2}
   \text{erf}\left(\frac{x-1}{2 \sqrt{t}}\right)+e^{\frac{t}{4}+x+\frac{1}{2}}
   \text{erf}\left(\frac{t+x}{2 \sqrt{t}}\right)-e^{\frac{t}{4}+x+\frac{1}{2}}
   \text{erf}\left(\frac{t+x+1}{2 \sqrt{t}}\right)-e^{x/2}
   \text{erfc}\left(\frac{x+1}{2 \sqrt{t}}\right)-e^{x/2}\right)\right)$$
in gnuplot type

~~~bash
 h(x,t)=2*log(-(-exp(x/2.) + exp(0.5 + t/4.)*erf((t - x)/(2.*sqrt(t))) -  exp(0.5 + t/4.)*erf((1 + t - x)/(2.*sqrt(t))) -  exp(x/2.)*erf((-1 + x)/(2.*sqrt(t))) +   exp(0.5 + t/4. + x)*erf((t + x)/(2.*sqrt(t))) -   exp(0.5 + t/4. + x)*erf((1 + t + x)/(2.*sqrt(t))) - exp(x/2.)*erfc((1 + x)/(2.*sqrt(t))))/(2.*exp(x/2.)))
  p'out' t'num'w lp,'' u 1:(h($1-$3,$3)+2*$3) t'exact' w l
~~~
which gives $h(x,t)$ plotted here for t=0 .1 .2 .3 ... .9 1.0 and $-5<x<5$ 
~~~gnuplot
 h(x,t)=2*log(-(-exp(x/2.) + exp(0.5 + t/4.)*erf((t - x)/(2.*sqrt(t))) -  exp(0.5 + t/4.)*erf((1 + t - x)/(2.*sqrt(t))) -  exp(x/2.)*erf((-1 + x)/(2.*sqrt(t))) +   exp(0.5 + t/4. + x)*erf((t + x)/(2.*sqrt(t))) -   exp(0.5 + t/4. + x)*erf((1 + t + x)/(2.*sqrt(t))) - exp(x/2.)*erfc((1 + x)/(2.*sqrt(t))))/(2.*exp(x/2.)))
 p'out' t'num'w lp,'' u 1:(h($1-$3,$3)+2*$3) t'exact' w l
~~~



# Bibliography 
M. Kardar, G. Parisi, Yi-C. Zhang, 
"Dynamic Scaling of Growing Interfaces"
Phys. Rev. Lett. 56 (1986) 889.

M.T. Batchelor, R.V. Burne, B.I. Henryc; S.D. Watt
"Deterministic KPZ model for stromatolite laminae"
Physica A 282 (2000) 123–136




ready for new site 09/05/19
*/

