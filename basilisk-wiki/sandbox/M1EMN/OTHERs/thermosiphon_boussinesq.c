/**
# Resolution of thermosiphon 1D with Boussinesq approximation
 
 Configuration : heated pipe and cooled pipe, vertical, joined by two adiabatic pipes, horizontal.
 
 ~~~gnuplot
 set arrow 1 from 0.01,0.01 to 0.01,.99 front
 set arrow 2 from .99,.99 to .99,0.01 front
 set arrow 3 from .5,.8 to 0.5,0.2 front
 set arrow 4 from .99,0.01 to 0.01,0.01 front
 set arrow 5 from 0.01,.99 to .99,0.99 front
 set label 1 "hot" at 0.1,.5 front
 set label 2 "cold" at .9,.5 front
 set label 3 "adiabatic " at .5,.95 front
 set label 4 "adiabatic " at .5,.05 front
 set label 5 "g" at 0.5,.7 front
 p [0:1]0 not,1 not w l
 unset arrow 1; unset arrow 2; unset arrow 3
 unset arrow 4; unset arrow 5;
 unset label 1; unset label 2; unset label 3
 unset label 4; unset label 5
 
 ~~~
 
All the problem consists to solve incompressibility (Boussinesq approximation),  momentum and  heat equation in 1D.

 
Incompressibility in the Boussinesq framework
$$ \frac{\partial  U}{\partial x} =0$$
 
 
The momentum equation is
 
$$0  =  (-\frac{\partial \Pi_{dyn}}{\partial x} -  (T-T_0) H_\theta  - U)$$

 
The reference temperature is $T_0$, $T_w$ is the wall temperature, the advection equation plus a source therm:
 
$$U\frac{\partial T}{\partial x} = - K(T-T_w)$$
 
$q_w=  K(T-T_w)$
 
 The domain is a bit strange! it is from $x=-1$ to $x=3$ and corresponds to the perimeter of the square. We start by:
 
 an adiabtic wall $K=0$ for $-1<x<0$ (bottom from right corner to left bottom corner: on the sketch $0<x<1$ $y=0$),
 
 then the hot wall has
 $T_w=1+\Delta T/2$ for $0<x<1$ (it correspond to the left pipe, from lower left corner to upper left corner:  on the sketch  $x=0$ $0<y<1$ ),
 
 then an adiabatic wall has $K=0$ for $1<x<2$, (top pipe:  on the sketch  $0<x<1$ $y=1$)
 
 finaly the cold one has $T_w=1-\Delta T/2$ for $2<x<3$
 (descending pipe:  on the sketch  $x=1$ $1>y>0$).
 
 
 We define $H_\theta$ on each pipe, it corresponds to the relative angle of the gravity:
 $H_\theta(-1<x<0)=0$, $H_\theta(0<x<1)=1$ $H_\theta(1<x<2)=0$ $H_\theta(2<x<3)=-1$,
 
 


 
## Code
mandatory declarations:
*/

#include "grid/multigrid.h"
#include "run.h"
#include "poisson.h"
/** definition of the field $T,U$, the flux, its derivative, time step and 
*/
scalar T[];
scalar U[];
scalar K[];
scalar Tw[];
scalar qw[];
scalar F[];
scalar Pid[];
double dt;
double cDelta;
double Kn,Tc,Tf,tmax,DeltaT;
mgstats mgpoi;
/**



 
Main with definition of parameters
*/
int main() {
  periodic (right);
  
    
  DeltaT = .5;
  Kn=1;
    
  L0 = 4.;
  X0 = -1;
  N = 128;
  DT = (L0/N)/2*(DeltaT/2);
  tmax = 100;
  
  run();
}
/** 
initial temperature 0 

Wall temperature is 1 for $0<x<1$
$K=5$ for $0<x<1$ and $2<x<3$, else $K=0$
*/
event init (t = 0) {
    foreach(){
    U[] = .05;
    T[] = 1;
    Tc = 1.+DeltaT/2;
    Tf = 1.-DeltaT/2;
    Tw[] = (x<=0)?0:(x<=1)?Tc:(x<=2)?0:(x<=3)?Tf:0;
    K[]  = (x<=0)?0:(x<=1)?Kn:(x<=2)?0:(x<=3)?Kn:0;
    }
  boundary ({U,T});
  }
/** 
print data

first point is in `X0+1/2*(L0/N))`, 
ith point is in `X0+(i-1/2)*(L0/N))`, 
last point is `X0+(N-1/2)*(L0/N))` which is `X0+L0-1/2*(L0/N))`
*/
event printdata (t += 5; t <=tmax) {
  foreach()
    fprintf (stdout, "%g %g %g %g %g %g  \n", x, T[],Pid[],U[],qw[], t);
    fprintf (stdout, "\n\n");
    
     fprintf (stderr, "%g %g \n",t,interpolate(U,1.5));
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
 
 
 We solve the unstaeady heat equation
 
 $$\frac{\partial T}{\partial t}+U\frac{\partial T}{\partial x} = - K(T-T_w)$$

 We consider here the simple case $F=U T$
 splitted in two parts, the first part is the convection
 $$\frac{\partial T}{\partial t}+\frac{\partial F(U)}{\partial x} = 0
 \text{ and  }\;\frac{\partial T}{\partial t}  = - K(T-T_w)$$
 the second part is heat transfer from the wall (say $q_w=  K(T-T_w)$).
 
 
 
 
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
    qw[]=-K[]*(T[]-Tw[]);
  boundary ({qw});

/**

 solves variation of velocity with a $\tau$ for stability,  $(T-T_0)$ is the Boussinesq source;
 
 $H_\theta( -1<x<0)=0$, $H_\theta(0<x<1)=1$ $H_\theta(1<x<2)=0$ $H_\theta(2<x<3)=-1$
 
 derivation of momentum
 $$
 0  =  (-\frac{\partial \Pi_{dyn}}{\partial x} -  H_\theta (T-T_0) - U)$$
 
 gives
 $$
   \frac{\partial^2 \Pi_{dyn}}{\partial x^2} =  \frac{\partial  }{\partial x}  ( H_\theta (T-T_0))$$
 
 
*/
    Pid[right]  = dirichlet(0.);
    Pid[left]  = dirichlet(0.);
    scalar bouss[],source[];

    double s1,s2;
    foreach() {
        s1= 1*(T[]-1);
        s2=-s1;
        bouss[] = ((x<0)?0:(x<1)?s1:(x<2)?0:(x<3)?s2:0);
    }
    boundary ({source});
    
    foreach() {
        source[] = 0*(qw[]) + (bouss[1]-bouss[0])/Delta;
    }
    boundary ({source});
    
    mgpoi = poisson (Pid,source);
    boundary ({Pid});

    
   /**
   Once we have $\Pi_{dyn}$, we do a relaxation (again an unsteady loop)
    $$
    \frac{\partial u }{\partial t} = \frac{1}{\tau} (-\frac{\partial \Pi_{dyn}}{\partial x} -  H_\theta (T-T_0) - U)
    $$
    
    Solved in implicit way
   */
    

   double tau = 20.;
    foreach() {
      U[] =  (U[]  +  (1./tau)*dt*( - (Pid[0]-Pid[-1])/Delta + bouss[] ))/(1+ (1./tau)*dt);
    }
      boundary ({U});
}
/**
## Run
Then compile and run:

~~~bash
 qcc  -g -O2 -DTRASH=1 -Wall  thermosiphon.c    ;./a.out > out
~~~

or better 

~~~bash
 ln -s ../../Makefile Makefile
 make thermosiphon.tst;make thermosiphon/plots    
 make thermosiphon.c.html ; open thermosiphon.c.html 
~~~

 

## Results: plot of temperature
The analytical solution for large $t$ is 
$$T(x,t) =  1- \Delta T/2 +(\Delta T )(1-exp(-K (x-x_0)))\text{   or  } T(x,t) = 1 + \Delta T/2 -(\Delta T )(1-exp(-K (x-x_2)))$$
 
  on the graph, we plot the temperature every $\Delta t=5$

The difference between the analytical curve and the steady computed is due to the fact that velocity changes. With no coupling, the curves are superposed.


~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "T" 
 T1(x)=(x<1? (1-exp(-5*x)):NaN)
 T2(x)=(x>2?exp(-5* (x-2)):NaN)
 TT=0
 p[-1:][0:1.4] 'out' u 1:($6>TT?$2:NaN)t'computed temp' w l
 
~~~


 Pressure $\Pi_{dyn}$ and assocaited velocity

~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "Pidyn, U"
  TT=70
 p[-1:][:] 'out' u 1:($6>TT?$3:NaN) t'Pi_dyn' w l,'' u 1:($6>TT?$4:NaN) t 'U' w l
 
~~~
 
 Velocity as function  of time

 
 
 ~~~gnuplot
 reset
 set xlabel "t"
 set ylabel "U(0.5)"
 p[:][:] 'log' u 1:2 w l,0,0.111113
 
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

* [https://math.unice.fr/~massonr/MAM5/Cours-VF-Boyer.pdf]()

Version juin 2021
*/
