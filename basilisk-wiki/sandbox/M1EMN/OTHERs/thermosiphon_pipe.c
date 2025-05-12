/**
# Resolution of Advection equation for heat in 1D
 
All the problem consists to solve heat equation in 1D
$$\frac{\partial T}{\partial t}+U\frac{\partial T}{\partial x} = - K(T-T_0)$$
the advection equation plus a source therm. We consider here the simple case $F=U T$
splitted in two parts, the first part is the convection
$$\frac{\partial T}{\partial t}+\frac{\partial F(U)}{\partial x} = 0
\text{ and  }\;\frac{\partial T}{\partial t}  = - K(T-T_w)$$
the second part is heat transfer from the wall (say $q_w=  K(T-T_w)$). 
We start by an adiabtic wall $K=0$ for $-1<x<0$, then the hot wall has 
$T_w=1$ for $0<x<1$,
 then an adiabatic wall has $K=0$
 for $1<x<2$, finaly the cold one has $T_w=0$ for $2<x<3$.
 
 
## Code
mandatory declarations:
*/
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the field $T,U$, the flux, its derivative, time step and 
*/
scalar T[];
scalar U[];
scalar K[];
scalar Tw[];
scalar qw[];
scalar F[];
scalar dU[];
double dt;
double cDelta;
double Kn;
/**

Boundary conditions

*/
U[left]   =  neumann(0);
U[right]  = neumann(0);
T[left]   =  neumann(0);
T[right]  = neumann(0);
qw[left]  = neumann(0);
qw[right] = neumann(0);
dU[left]  = dirichlet(0);
dU[right] = neumann(0);
/**
Main with definition of parameters
*/
int main() {
  L0 = 4.;
  X0 = -1;
  N = 256;
  DT = (L0/N)/2;
  Kn=5;
  run();
}
/** 
initial temperature 0 

Wall temperature is 1 for $0<x<1$
$K=5$ for $0<x<1$ and $2<x<3$, else $K=0$
*/
event init (t = 0) {
    foreach(){
    U[] = 1;
    T[] = 0;
    Tw[] = (x<0)?0:(x<1)?1:0;
    K[]=(x<0)?0:(x<1)?Kn:(x<2)?0:(x<3)?Kn:0;
    }
  boundary ({U,T});
  }
/** 
print data

first point is in `X0+1/2*(L0/N))`, 
ith point is in `X0+(i-1/2)*(L0/N))`, 
last point is `X0+(N-1/2)*(L0/N))` which is `X0+L0-1/2*(L0/N))`
*/
event printdata (t += .1; t <=8) {
  foreach()
    fprintf (stdout, "%g %g %g %g %g \n", x, T[],U[],dU[], t);
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

 
Pour fixer les idées, dans le cas de l'équation d'advection simple
$$\partial_t T+\partial_x (U T)=0$$
on a $F(T)=U T$, la valeur propre est $dF/dU=U_0$ tout simplement, si on utilise le flux
$$ 
 F_{i}=  \dfrac{U_{i-1}T_{i-1}+U_{i}T_i}{2} -c_\Delta (\dfrac{T_{i}-T_{i-1}}{2})
$$
  ou  si on prend $c_{\Delta}=U_0$, c'est le flux  upwind:
et le nouveau $T$
$$
 T_i^{n+1}=T_i^{n} -{\Delta t} \dfrac{(F_{i+1}-F_{i})}{ \Delta x}
$$ 



  Expression of the flux  :
$$
F_i=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2} 
% - c_\Delta ({F(U_{i+1})-F(U_{i-1})})
  - c_\Delta \frac{({(T_{i})-(T_{i-1})})}{2}
$$ 
avec   
dans le cas Lax
 $c_\Delta=\Delta/\Delta t$
 dans le cas centré 
 $c=0$ 
 dans le cas upwind  $c=1$ 
*/ 
  foreach() {
    cDelta = (U[]+U[-1])/2.;
    F[] =(T[]*U[]+T[-1]*U[-1])/2. - cDelta *(T[]-T[-1])/2;}
  boundary ({F});
/** 
explicit step
update 
$$
 T_i^{*}=T_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i})}{\Delta x}
$$*/
  foreach()
    T[] +=  - dt* ( F[1] - F[] )/Delta;
  boundary ({T});
/** 
source term explicit 
$$
 T_i^{n+1}=T_i^{*} -{\Delta t} K (T^{*}-T_w)   
$$
but 
implicit here
$$
 T_i^{n+1}=T_i^{*} -{\Delta t} K (T^{n+1}-T_w)  
$$


*/  
  foreach()
    T[] = (T[]  - K[]*dt*( - Tw[]))/(1 + K[]*dt);
  boundary ({T});

  foreach()
    qw[]=K[]*(T[]-Tw[]);
  boundary ({qw});

/**

solves variation of velocity 
$$
\frac{\partial u'}{\partial x} = q_w
$$
it is solved as teh steady solution of 
$$
\frac{\partial u'}{\partial t}= - ( \frac{\partial u'}{\partial x} - q_w)
$$
*ie* a fixed point method (implicit form) 
$$
 u'_{i,n+1} -u'_{i,n} = -\alpha ( \frac{u'_{i,n+1} -u'_{i-1,n}}{\Delta} - q_w)
 \text{ so we iterate }\; u'_{i,n+1}  = 
 \frac{ u'_{i,n}-\alpha (\frac{-u'_{i-1,n}}{\Delta} - q_w)}{1+ \frac{\alpha}{\Delta}}
$$
(or explicit form) 
*/
#if 1
  foreach()
    dU[]=0; 
  double  e=1,alphat=1000;
  while(fabs(e)>.0000001){
    foreach(){
      dU[]=(dU[]- (alphat)* (-dU[-1]/Delta-qw[]))/(1+(alphat)*1/Delta);
      e=(dU[]-dU[-1])/Delta-qw[]; }
    boundary ({dU});}
#else
  foreach()
    dU[] = -x;
  double ddU=1;
   while(fabs(ddU)>.00001){
    foreach(){
      ddU =- (.01)* ((dU[]-dU[-1])/Delta-qw[]);
      dU[] += ddU;}
    boundary ({dU});}
#endif
/**
then one couples  $U=1+u'$ check the signs of $q_w$ 
*/    
  foreach()
    U[] = 1 - 0.1*dU[];

}
/**
## Run
Then compile and run:

~~~bash
 qcc  -g -O2 -DTRASH=1 -Wall  advecte1.c -o thermosiphon ;./thermosiphon > out
~~~

or better 

~~~bash
 ln -s ../../Makefile Makefile
 make thermosiphon.tst;make thermosiphon/plots    
 make thermosiphon.c.html ; open thermosiphon.c.html 
~~~

 

## Results: plot of temperature
The analytical solution for large $t$ is 
$$T(x,t) =  1 - exp(-K (x-x_0))\text{   or  } T(x,t) =    exp(-K (x-x_2))$$
 
  on the graph, we plot the temperature every $\Delta t=.1$

The difference between the analytical curve and the steady computed is due to the fact that velocity changes. With no coupling, the curves are superposed.


~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "T" 
 T1(x)=(x<1? (1-exp(-5*x)):NaN)
 T2(x)=(x>2?exp(-5* (x-2)):NaN)
 p[-1:][0:1.4] 'out' t'computed temp' w l,T1(x) w l,T2(x) w l linec -1
 
~~~


at fixed $x$, not to close to 0 

The analytical solution  is 
$$T(x,t) =  1 - exp(-K t) $$
small difference is due to the coupling 

~~~gnuplot
 reset
 set xlabel "t"
 set ylabel "T" 
p[:1][0:1.4] 'out' u ($1>.85&&$1<1?($5):NaN):2 t'computed temp' w p,1-exp(-5*x) t'analytic'

~~~


## Results: plot of velocity

The velocity changes due to the surece term (thermal flux)

~~~gnuplot
 reset
 set xlabel "x" 
 set ylabel "U" 
 p[-1:][:] 'out' u 1:3 not
 
~~~

 

# Links
  
* [advecte1.c]() explains the notions of advection, testing the flux, coded with Basilisk 
* [http://basilisk.fr/sandbox/Antoonvh/integrator2.h]() and [http://basilisk.fr/sandbox/Antoonvh/ti2.c]() solve 
$\frac{\partial^2 u'}{\partial x^2} = \frac{\partial }{\partial x} q_w$
 

# Bibliography 
 
* [LeVeque](https://www.cambridge.org/core/books/finite-volume-methods-for-hyperbolic-problems/finite-volume-methods/CB7B0A27A6D37AE3B906D4AE7C60A87E) 
 chapitre 4, Finite Volume methods
* [Toro](https://www.springer.com/gp/book/9783540252023) "Riemann Solvers and Numerical Methods"  Chap 5 Springer
* [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/ENSTA/coursENSTA.html)   cours ENSTA thermique
*  [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MECAVENIR/index.html)   cours EPU thermique


Version juin 2020
*/
