/**
# Resolution of simple transport equation, centered and decentered
Explicit resolution of
$$\frac{\partial h}{\partial t}+ v_0 \frac{\partial h}{\partial x}=0$$

with various schemes for comparison: centered derivative $O(\Delta^2)$:

 $$\frac{h(t+\Delta t,x) - h(t ,x)}{\Delta t}+ v_0 \frac{h(t,x+\Delta x) - h(t ,x-\Delta x)}{2 \Delta x}=0$$
  it is not a good idea, it is unstable, 
down stream decentered derivative $O(\Delta)$, is  stable, but damped:
 
  $$\frac{h(t+\Delta t,x) - h(t ,x)}{\Delta t}+ v_0 \frac{h(t,x) - h(t ,x-\Delta x)}{\Delta x}=0$$

 

## Code
mandatory declarations:
*/
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the height of interface
its `O(Delta)` derivative and its `O(Delta^2)` derivative, time step
*/
scalar h0[];
scalar hc[];
scalar hd[];
scalar hp[];
scalar hdp[];
scalar hp2[]; 
double dt;
/**
Main with definition of parameters
*/
int main() {
  L0 = 8.;
  X0 = -L0/2;
  N = 128/2;
  DT = (L0/N)/10 ;
#define v0 .5
  run();
}
/** 
initial elevation: two "triangles"
*/
event init (t = 0) {
  foreach(){
    h0[] =   -(1-fabs(x-2))*(fabs(x-2)<1) + (1-fabs(x+2))*(fabs(x+2)<1);
    hc[]=h0[];
    hd[]=h0[];
    }
  boundary ({hc,hd});
  }
/** 
print data
*/
event printdata (t += .5; t <= 2) {
  foreach()
    fprintf (stdout, "%g %g %g %g \n", x, hc[],hd[], t);
  fprintf (stdout, "\n");
}
/** 
integration 
*/
event integration (i++) {
  double dt = DT;
/**
finding the good next time step
*/
  dt = dtnext (dt);
/**
  $O(\Delta)$ derivative
*/
  foreach()
    hp[] =  ( hc[1,0] - hc[0,0] )/Delta;
  boundary ({hp});
/**
  down stream decentered derivative $O(\Delta)$, stable, but damped
*/   
  foreach()
    hdp[] =  ( hd[0,0] - hd[-1,0] )/Delta;
  boundary ({hp});
/**
  centered derivative $O(\Delta^2)$, it is not a good idea, it is unstable
*/ 
  foreach()
    hp2[] =  ( hc[1,0] - hc[-1,0] )/2/Delta;
  boundary ({hp2});
/** 
update of the centered and decentered fields  
$$ h(t+\Delta t,x)=  h(t ,x)- \Delta t v_0 \frac{h(t,x+\Delta x) - h(t ,x-\Delta x)}{2 \Delta x}$$
  and
  $$h(t+\Delta t,x)=  h(t ,x)- \Delta t v_0 \frac{h(t,x) - h(t ,x-\Delta x)}{\Delta x}$$
*/
  foreach(){
    hc[] += -dt*v0*hp2[];
    hd[] += -dt*v0*hdp[]; 
  }
  boundary ({hc,hd});  
}
/**
## Run
Then compile and run:

~~~bash
qcc  -g -O2 -DTRASH=1 -Wall  advecte.c -o advecte
./advecte > out
~~~

or better 

~~~bash
 ln -s ../../Makefile Makefile
 make advecte.tst;make advecte/plots    
 make advecte.c.html ; open advecte.c.html
~~~
 
~~~bash
 source ../c2html.sh advecte
~~~



## Results
The analytical (for $v_0=0$) solution is 
$$h(x,t) =  h_0(x-v_0t)$$
in gnuplot type

~~~bash
  v0=.5
  h0(x,t)=-(1-abs(x-v0*t-2))*(abs(x-v0*t-2)<1)  + (1-abs(x-v0*t+2))*(abs(x-v0*t+2)<1) 
  p'out' t'cent.'w lp,'' u 1:3 t'dec.','' u 1:(h0($1,$4)) t'exact' w l
~~~
which gives $h(x,t)$ plotted here for t=0 .5 ... 2.0 and $-8<x<8$ 

~~~gnuplot
 v0=.5
 h0(x,t)=-(1-abs(x-v0*t-2))*(abs(x-v0*t-2)<1)  + (1-abs(x-v0*t+2))*(abs(x-v0*t+2)<1) 
 p'out' t'cent.'w lp,'' u 1:3 t'dec.','' u 1:(h0($1,$4)) t'exact' w l
~~~

 


# Links
  
* [advecte1.c]() explains the notions of advection, testing the flux, coded with Basilisk 

* [advecte1c.c]() the same than the previous one but with standard C  

* [advecte.c]() advection with Basilisk, compares centered and decentered

# Bibliography 
 
* [LeVeque](https://www.cambridge.org/core/books/finite-volume-methods-for-hyperbolic-problems/finite-volume-methods/CB7B0A27A6D37AE3B906D4AE7C60A87E)  chapitre 4, Finite Volume methods



ready for new site 09/05/19
*/

