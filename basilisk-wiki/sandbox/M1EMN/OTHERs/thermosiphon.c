/**
# Resolution of thermosiphon 1D (Low Mach developpment of 1D NS equations)
 
 Configuration : heated pipe and cooled pipe, vertical, joined by two adiabatic pipes, horizontal. By natural convection, it creates a flow.
 
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
 
All the problem consists to solve mass,  momentum and   energy conservation  in 1D (with low Mach approximation, constant pressure, 0 mean flux here).
 This is not boussinesq approximation.

 
 Mass conservation (should add source term $s=A P' - B q_w$  ($B>0$; $A<0$), thermodynamic pressure $P$  is supposed constant by symmetry of the problem)
$$ \frac{\partial  U}{\partial x} =s$$

The steady momentum equation is
 $$0  =  (-\frac{\partial \Pi_{dyn}}{\partial x} - \frac{H_\theta}{T}  - c_f U)$$
the term $1/T$ is density, $H_\theta=0,\pm 1$ corresponds to direction of gravity, the friction is laminar and there is a strong approximation  on the dependance of viscosity with temperature so that we obtain a constant coeffcient $c_f$.
 
The reference temperature is $T_0=1$, $T_w$ is the wall temperature, imposed hott and cold $T_w=1 \pm \Delta T/2$,
the heat equation:
 
$$U\frac{\partial T}{\partial x} = - K(T-T_w)$$
 with wall flux
$q_w=  K(T-T_w)$
 
 The domain is a bit strange because we follow the perimeter of the square. It is from $x=-1$ to $x=3$ and corresponds to the perimeter of the square. We start by:
 
 an adiabtic wall $K=0$ for $-1<x<0$ (bottom from right corner to left bottom corner: on the sketch $0<x<1$ $y=0$),
 
 then the hot wall has an imposed temperature
 $T_w=1+\Delta T/2$ for $0<x<1$ (it corresponds to the left pipe, from lower left corner to upper left corner:  on the sketch  $x=0$ $0<y<1$ ),
 
 then an adiabatic wall has $K=0$ for $1<x<2$, (top pipe:  on the sketch  $0<x<1$ $y=1$)
 
 finaly the cold one has $T_w=1-\Delta T/2$ for $2<x<3$
 (descending pipe:  on the sketch  $x=1$ $1>y>0$).
 
 
 We define $H_\theta$ on each pipe, it corresponds to the relative angle of the gravity:
 $H_\theta(-1<x<0)=0$, $H_\theta(0<x<1)=1$, $H_\theta(1<x<2)=0$, $H_\theta(2<x<3)=-1$.
 
 


 
## Code
mandatory declarations:
*/

#include "grid/multigrid.h"
#include "run.h"
#include "poisson.h"
/** definition of the field $T,U$, the flux, its derivative, time step and 
*/
scalar T[];
scalar U[],Uold[],Ustar[];
scalar K[];
scalar Tw[];
scalar qw[];
scalar F[];
scalar Pid[];
double dt;
double cDelta;
double Kn,Tc,Tf,tmax,DeltaT,cf,B,err;
mgstats mgpoi;
/**
 
Main with definition of parameters ($\Delta T, K, c_f, -B<0$), domain size, the CFL is heuristical, $B=0$ to test a simple analytical solution,
*/
int main() {
  periodic (right);
  
    
  DeltaT = .5;
  Kn = 1;
  cf = 1./3;
  B = .2*0;
    
  L0 = 4.;
  X0 = -1;
  N = 128;
  DT = (L0/N)*(DeltaT/2)/50;
  tmax = 20;
  
  run();
}
/** 
initial temperature 1

Wall temperature is 1 for $0<x<1$
$K=1$ for $0<x<1$ and $2<x<3$, else $K=0$
*/
event init (t = 0) {
    foreach(){
    U[] = .5;
    Uold[] = U[];
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
event printdata (t += 1; t <=tmax) {
  foreach()
    fprintf (stdout, "%g %g %g %g %g %g  \n", x, T[],Pid[],U[],qw[],t);
    fprintf (stdout, "\n\n");
    
     fprintf (stderr, "%g %g %g \n",t,interpolate(U,1.5),sqrt(err)/N);
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

 solves variation of velocity with a $\tau$ for stability,  $1/(T)$ is the source;
 
 $H_\theta( -1<x<0)=0$, $H_\theta(0<x<1)=1$ $H_\theta(1<x<2)=0$ $H_\theta(2<x<3)=-1$
 
 
 
  derivation of momentum
 $$
 \frac{\partial U }{\partial t} =  (-\frac{\partial \Pi_{dyn}}{\partial x} -  \frac{H_\theta}{ T} - c_f U)
 $$
 we do   a split in two steps.
 First a friction step:
 $$
 \frac{U^{*}-U^n }{\Delta t} =   -  c_f  U^{*}
 $$
*/
    foreach() {
        Uold[] =U[];
        Ustar[] =  U[]/(1+ (cf)*dt);
    }      boundary ({Ustar});
/**
 Second, a pressure step.
 $$
 \frac{U^{n+1}-U^* }{\Delta t} =  (-\frac{\partial \Pi_{dyn}}{\partial x} -  \frac{H_\theta}{ T})$$
For the pressure part,
we want $\frac{\partial U^{n+1}}{\partial x}=s$, (we implement here $s=0$).   Hence, the "projection step" gives
$$\frac{\partial^2 \Pi_{dyn}}{\partial x^2} = - \frac{\partial  }{\partial x} (\frac{H_\theta}{ T}) +  \frac{1}{\Delta t}(\frac{\partial U^*}{\partial x}-s)$$
                    */

    
    Pid[right]  = dirichlet(0.);
    Pid[left]  = dirichlet(0.);
    scalar bouss[],source[],div[];

    double s1,s2;
    foreach() {
        s1= -1/T[];
        s2=-s1;
        bouss[] = ((x<0)?0:(x<1)?s1:(x<2)?0:(x<3)?s2:0) ;
    }
    boundary ({source});
 
    
    double s;
    foreach() {
        s = -B*qw[];
        source[] = (bouss[1]-bouss[0])/Delta  +  ((Ustar[1]-Ustar[0])/Delta - s)/dt;
    }
    boundary ({source,div});
    
    mgpoi = poisson (Pid,source);
    boundary ({Pid});

    
   /**
   Once we have $\Pi_{dyn}$, we update $U$
 $$ U^{n+1}= U^*  - {\Delta t}(-\frac{\partial \Pi_{dyn}}{\partial x} -  \frac{H_\theta}{ T})$$
   */
    
    err=0;
    foreach() {
        Uold[] =U[];
        U[] =  (Ustar[]  + dt*( - (Pid[0]-Pid[-1])/Delta + bouss[]));
        err=+sq( Uold[]-U[]);
    }  boundary ({U});
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

 
 
## Analytical solution
 
 the analytical solution is obtained in writing that along a loop  (the reason for the 4 in the equations) the pressure is periodic.
 
 $$\frac{4 c_f U}{\lambda }-\frac{2 \varepsilon}{\left(1-\varepsilon^2\right) \lambda }+\frac{2 \log
 \left(\frac{1+\varepsilon \tanh \left(\frac{1}{2 \lambda }\right)}{1-\varepsilon \tanh \left(\frac{1}{2
 \lambda }\right)}\right)}{1-\varepsilon^2}=0$$
 with $\varepsilon = \Delta T/2$ and $\lambda=K/U$
 
 Examples of values with $C_f=1/3$,  $K=1$
 
 
 $\Delta T=.5$, $U=0.222404$
 
 $\Delta T= 1$ $U=0.34048$
 
 $\Delta T= .25$ $U=0.137786$
 
 
## Results: plot of temperature
The analytical solution for large $t$ is 
$$T(x,t) =  1- \Delta T/2 +(\Delta T )(1-exp(-K/U (x-x_0)))\text{   or  } T(x,t) = 1 + \Delta T/2 -(\Delta T )(1-exp(-K/U (x-x_2)))$$
 
  on the graph, we plot the temperature every $\Delta t=5$

The difference between the analytical curve and the steady computed is due to the fact that velocity changes. With no coupling, the curves are superposed.


~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "T" 
 T1(x)=(x>0 && x<2 ? .75+ .5*(1-exp(-x/.2224)):NaN)
 TT=10
 p[-1:][0:1.4] 'out' u 1:($6>TT?$2:NaN)t'computed temp' w l,T1(x) t'exacte'
 
~~~

## Results: plot of pressure

 Pressure $\Pi_{dyn}$, this pressure is the deviation of pressure from thermodynamical pressure.
 
~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "Pidyn"
 TT=10
 p[-1:][:] 'out' u 1:($6>TT?$3:NaN) t'Pi_dyn' w l
 
~~~
 
## Results: plots of velocity
 
 velocity, note that $\partial_x U = -B q_w$, here $B=0$
 
 velocity at a  point  as function  of time compared with analytical solution
 $C_f=1/3$, $\Delta T=.5$, $K=1$ we have $U=0.222404$
 
 
 
~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "U"
 TT=15
 p[-1:][0:] 'out'  u 1:($6>TT?$4:NaN) t 'U' w l,0.222404 t'B=0'
 
~~~

 
 
 
 
 Error on stationarisation of velocity as function  of time
 
 
 
 ~~~gnuplot
 set logscale
 set xlabel "t"
 set ylabel "erreur"
 p[10:]'log' u 1:3 w l,1e-9,1e-12
 ~~~


 
 ~~~gnuplot
 reset
 set xlabel "t"
 set ylabel "U(1.5)"
 p[:][:] 'log' u 1:2 w l,0,.22240
 
 ~~~
 
 ~~~gnuplot
 reset
  set logscale
 set xlabel "t"
 set ylabel "U(1.5)-0.222404"
 p[1:][:] 'log' u 1:($2-0.222404) w l,0
 
 ~~~
 
# Links
  
* [../BASIC/advecte1.c]() explains the notions of advection, testing the flux, coded with Basilisk 
* [http://basilisk.fr/sandbox/Antoonvh/integrator2.h]() and [http://basilisk.fr/sandbox/Antoonvh/ti2.c]() solve 
$\frac{\partial^2 u'}{\partial x^2} = \frac{\partial }{\partial x} q_w$
* [thermosiphon_pipe.c]() explains the notions of advection, testing the flux, coded with Basilisk
* [http://basilisk.fr/sandbox/M1EMN/OTHERs/thermosiphon_boussinesq.c]() same with Boussinesq approximation
 
 
 
# Bibliography 
 
* [LeVeque](https://www.cambridge.org/core/books/finite-volume-methods-for-hyperbolic-problems/finite-volume-methods/CB7B0A27A6D37AE3B906D4AE7C60A87E) 
 chapitre 4, Finite Volume methods
* [Toro](https://www.springer.com/gp/book/9783540252023) "Riemann Solvers and Numerical Methods"  Chap 5 Springer
* [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/ENSTA/coursENSTA.html)   cours ENSTA thermique
*  [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MECAVENIR/index.html)   cours EPU thermique
* [https://math.unice.fr/~massonr/MAM5/Cours-VF-Boyer.pdf]() (about dirac source term, not used here)

Version juin 2021
*/
