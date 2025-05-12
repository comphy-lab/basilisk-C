/**
# Resolution of Advection equation in 1D
 
All the problem consists to solve
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
 The equation is written as well 
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial U} 
\frac{\partial U}{\partial x} = 0\;\text{  or } \frac{\partial U}{\partial t}+c_F
\frac{\partial U}{\partial x} = 0,$$
the advection equation 
with $c_F= \frac{\partial F(U)}{\partial U}.$

 We consider here the simple case $F=U$,
hence $c_F=1$.

 
 

## Finite volumes

Space is decomposed in small sgments of  a priori different length, but here we suppose that the length is a constant $\Delta x$.
The same for time increment: $\Delta t$ is constant.  
Those small segments are kind of  "volumes", as we are in  dimension one.


Consider the system written in conservative form:
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
and integrate in $x$ on a small interval between $x_i-\Delta x /2$ and $x_{i}+\Delta x /2$  and consider two times  $t^n$ and $t^{n+1}$.
 

A mean value of $U$ around $x_i$ between $x_i-\Delta x /2$ and 
 $x_i+\Delta x/2$ may be defined as 
  $U_i^n$ the integral between $x_i-\Delta x /2$ and
 $x_i+\Delta x/2$ is so by definition
$$
 U_i^n =\dfrac{1}{\Delta x} \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U(x,t_n)dx
$$
index $i$ is for the segment
$C_i=(x_{i}-\Delta x /2,x_{i}+\Delta x /2)$,  centered in $x_{i}$
 
index $n$ corresponds to time  $t_n$ with $t_{n+1}-t_{n}=\Delta t$. So that with the definition:
$$
 U_i^{n+1} =\dfrac{1}{\Delta x} \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U(x,t_{n+1})dx
$$


Hence by the integration on the "segment/ Volume" of 
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
we have if we integrate (on the fixed volume centered in $x_i$):

 $\int_{x_i-\Delta x /2}^{x_{i}+\Delta x/2} \frac{\partial U}{\partial t}  dx  =\frac{d}{dt}  \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U dx$
 and for the flux 
 $\int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} \partial_x F(U) dx = F^n_{i+1}-F^n_{i}.$
 This is  $\frac{d}{dt}  \int_{x_{i}-\Delta x/2}^{x_{i}+\Delta x/2} U dx = F^n_{i+1}-F^n_{i}.$ Note that sometimes we prefer the notation 
  $\int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} \partial_x F(U) dx = F^n_{i+1/2}-F^n_{i-1/2}$ to insist on the values taken on the faces $i-1/2$. But as the index of tables is an integer we write it as $F^n_{i+1}-F^n_{i}.$ 
 
 The "exact" integration from the flux on the  "volume" for the mean value gives:
  $$\frac{d}{dt}(\frac{1}{\Delta x} \int_{x_{i}-\Delta x/2}^{x_{i}+\Delta x/2} U dx )
 + \frac{F^n_{i+1}-F^n_{i}}{\Delta x}=0.$$
 
  


 Taylor expansion:
$$ U(t + \Delta t) = U(t )+\Delta t \partial_t U  + (1/2)(\Delta t)^2 \partial^2_t U  +O(\Delta t)^3, $$
allows us to write the approximation (at first order in time) whis is   the finite volume method. 
$$
 \dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1}-F^n_{i}}{\Delta x}=0,
$$
Numerical flux $F_{i}$ is an approximation of  $F(U)$ at left interface (in $x_{i-1/2}$) of segment  $C_i$ centered in $x_i$.
Numerical flux $F_{i+1}$ is an approximation of  $F(U)$ at right interface
(in $x_{i+1/2}$) of segment  $C_i$ centered in $x_i$.

$F_{i}$ is a function of the value of  $U_i$ in the considered segment, which begins in  $i-1/2$ 
 and of the   value $U_{i-1}$ from the previuous one which ends in $i+1/2$.
 $F_i$  enters the volume $i$, trough the face in $x_{i-1/2}$ 
 
 
$F_{i+1}$ is a function of the value of  $U_i$ in the considered segment, which begins in  $i-1/2$ 
 and of the   value $U_{i+1}$ from the next one which begins in $i+1/2$.
 $F_{i+1}$ is the flux that leaves the cell $i$ through the face$x_{i+1/2}$.


The numerical flux across face i+1/2 is denoted $F_{i+1}$ (or $F^n_{i+1}$ at time $n$), it is function (say $f$) of values 
before and after the face (i+1) which are  $U_i$ and $U_{i+1}$   
$$
 F_{i+1}=f(U_i,U_{i+1}).
$$
The position of the center of the cell is $x_{i}$. 



~~~gnuplot
set samples 9 
set label "U i-1" at 1.5,3.1
set label "U i" at 2.5,3.15
set label "U i+1" at 3.5,2.5
set xtics ("i-2" 0.5, "i-1" 1.5, "i" 2.5,"i+1" 3.5,"i+2" 4.5,"i+3" 5.5)
set arrow from 2,1 to 2.5,1
set arrow from 3,1 to 3.5,1
set label "F i" at 2.1,1.25
set label "F i+1" at 3.1,1.25

set label "x i-1/2" at 1.5,0.25
set label "x i" at 2.4,0.25
set label "x i+1/2" at 3.,0.25

set label "x"  at 0.5,2+sin(0) 
set label "x"  at 1.5,2+sin(1)
set label "x"  at 2.5,2+sin(2) 
set label "x"  at 3.5,2+sin(3) 
set label "x"  at 4.5,2+sin(4) 
set label "x"  at 5.5,2+sin(5) 
p[-1:7][0:4] 2+sin(x) w steps not,2+sin(x) w impulse not linec 1
~~~


So, if the length of domain is $L$, and if the domain starts in $x_0$, and if we take $N$ points, hence $\Delta x=L/N$
faces are $x_0 + (i-1) \Delta x$.

A the center of the cell, $x_{i}=x_0 + (i-1/2) (\Delta x)$, we have the mean value $U^n_i$. 


The finite volume method is 
$$
 \dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1}-F^n_{i}}{\Delta x}=0,
$$
It is explicit:  compute new values  $U_i^{n+1}$ as a function of the old ones $U_i^{n}$ .
 
Nota 1:

Do not confuse the $f(,)$ function, numerical flux across face $F_i$ and actual flux function $F$ comming from the physics.

Nota 2:

This is a simplified version point of view, Basilisk uses 
"foreach_face()"  and  "face vector"

 

### Flux BCG 
 
 
 In the previous examples we wrote 
 $$F_i=\dfrac{F(U_{i-1})+F(U_i)}{2} 
  -  c_\Delta \frac{({(U_{i})-(U_{i-1})})}{2}$$
 with an $ad$ $hoc$ $c_\Delta$. Let us find BCG's one.
 
 
The value of $U$ at intermediate $n+1/2$ time on face ($i-1/2$) is say  $U_{i-1/2}^{n+1/2}$. We write it as a double Taylor expansion starting from $U_{i-1}$. The derivative $\frac{\partial U_i^n}{\partial x}$ is approximated by $\frac{  U_i^n -U_{i-1}^n}{\Delta x}$ (values after  and before $i-1/2$) hence :
$$U_{i-1/2}^{n+1/2}=U_{i-1}^n +\frac{\Delta x}{2}\frac{\partial U_i^n}{\partial x} +\frac{\Delta t}{2}\frac{\partial  U_i^n}{\partial t}+ O(\Delta t^2,\Delta t^2)$$

Using  $\frac{\partial U}{\partial t}= - c_F\frac{\partial U}{\partial x}$ 
 temporal derivative can be replaced by spatial derivatives yielding 

$$U_{i-1/2}^{n+1/2}=U_{i-1}^n +(\frac{\Delta x}{2}- \frac{c_F\Delta t}{2}) \frac{\partial U_i^n}{\partial x}  $$
this is a simplified upwind scheme, 
 which gives the flux  (as $F_i=c_F U_i^{n+1/2}$)
 
 
 $$F_i=c_F  (U_{i-1}^n +\frac{\Delta x}{2} (1 - \frac{c_F\Delta t}{\Delta x}) (\frac{  U_i^n -U_{i-1}^n}{\Delta x}))$$
 
 we rearrange : 
  $$F_i=c_F (\frac{ U_i^n +U_{i-1}^n}{2}) -  c_F ( \frac{c_F\Delta t}{\Delta x}) (\frac{  U_i^n -U_{i-1}^n}{2})$$

as $F=c_FU$  this is the usual formaulation of the mean flux plus a correction
 $$F_i=\dfrac{F(U_{i-1})+F(U_i)}{2} 
  -  c_\Delta \frac{({(U_{i})-(U_{i-1})})}{2}$$
 with $c_\Delta = c_F ( \frac{c_F\Delta t}{\Delta x})$
 
 
 
## Code
mandatory declarations:
*/
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the field h, the flux, its derivative, time step and 
*/
scalar U[];
scalar F[];
double dt;  
double cF=1;
/**

Boundary conditions

*/
U[left]  = neumann(0);
U[right] = neumann(0);
/**
Main with definition of parameters
*/
int main() {
  L0 = 12.;
  X0 = -L0/4;
  N = 64;    
  DT = (L0/N)/2;
  run();
}
/** 
initial elevation: an exponential "bump"
*/
event init (t = 0) {
  foreach()
    U[] = exp(-x*x);
  boundary ({U});
  }
/** 
print data

first point is in `X0+1/2*(L0/N))`, 
ith point is in `X0+(i-1/2)*(L0/N))`, 
last point is `X0+(N-1/2)*(L0/N))` which is `X0+L0-1/2*(L0/N))`
*/
event printdata (t += 1; t <=3) {
  foreach()
    fprintf (stdout, "%g %g %g  \n", x, U[], t);
  fprintf (stdout, "\n\n");
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

 We use the BCG flux 
 $$F_i=U_{i-1}^n +\frac{\Delta x}{2} (1 - \frac{c_F\Delta t}{\Delta x}) \frac{  U_i^n -U_{i-1}}{\Delta x} $$
  the commented line corresponds to this notation.
 
But as it written as well 
  $$F_i=c_F (\frac{ U_i^n +U_{i-1}^n}{2}) -  c_F ( \frac{c_F\Delta t}{\Delta x}) (\frac{  U_i^n -U_{i-1}^n}{2})$$

as $F=c_FU$  this is  
 $$F_i=\dfrac{F(U_{i-1})+F(U_i)}{2} 
  -  (c_F ( \frac{c_F\Delta t}{\Delta x})) \frac{({U_{i}-U_{i-1}})}{2}$$
 which is the code line


*/ 
  foreach() {
   // F[] = U[-1] + (Delta/2)*(1-dt*cF/Delta)*(U[ ]-U[-1])/Delta;
    F[] = (cF*U[]+cF*U[-1])/2  - dt*cF/Delta*cF*(U[ ]-U[-1])/2;
   }
  boundary ({F});
/** 
explicit step
update 
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i})}{\Delta x}
$$*/
  foreach()
    U[] +=  - dt* ( F[1] - F[0] )/Delta;
  boundary ({U});
}
/**
## Run
Then compile and run:

~~~bash
 qcc  advecteB-C-G.c -o advecteB-C-G ;./advecteB-C-G > out
~~~

or better 

~~~bash
 ln -s ../../Makefile Makefile
 make advecteB-C-G.tst;make advecteB-C-G/plots    
 make advecteB-C-G.c.html ; open advecteB-C-G.c.html 
~~~

 



## Results
The analytical solution is 
$$U(x,t) =  exp(-(x-t)^2)$$
in gnuplot type

~~~bash
 U(x,t)= exp(-(x-t)*(x-t))
 p'out' u ($1):($2)t'num'w l,'' u 1:(U($1,$3)) t'exact' w l
~~~
which gives $U(x,t)$ plotted here for t=0 1 2 3   and $-3<x<9$  
$\Delta=L0/N=12/64=0.1875$
first point is $-3 + \Delta/2= -2.90625$ second point is in $-3 + \Delta/2+\Delta=-2.71875$ next is in -2.53125
up to the previous last one $-3 +1/2 \Delta +(N-2)\Delta = 8.71875$  and finally the 64th point ($N$) is  $-3 +1/2 \Delta + (N-1) \Delta =8.90625.$

(`x1=X0+(1-1/2) L0/N, x2=X0 +(2-1/2) L0/N, ...,xi= X0 +(i-1/2) L0/N, .. xN=X0 +(N-1/2) L0/N= X0 +L0-1/2 L0/N`)  

  




~~~gnuplot
 set output 'dessin.svg'; 
 reset
 set xlabel "x"
 U(x,t)= exp(-(x-t)*(x-t))
 p'out' u ($1):($2)t'num.'w p,'' u 1:(U($1,$3)) t'exact' w l
 
~~~


 

# Links
  
* [advecte1.c]() explains the notions of advection, testing the flux, coded with Basilisk 

* [advecte1c.c]() the same than the previous one but with standard C  

* [advecte.c]() advection with Basilisk

* [advecteB-C-G.c]() advection with B-C-G flux in Basilisk

# Bibliography 
 
* [Popinet](https://hal.science/hal-01445436) 
Gerris: a tree-based adaptive solver for the incompressible Euler equations in complex geometries 
 
 * J. B. Bell, P. Colella, and H. M. Glaz. A second-order projection method for the incompressible navier-stokes equations. J. Comput. Phys., 85:257–283, 1989.
 
 
*/

