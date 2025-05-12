/**
# Resolution of flood wave numericaly and with characteristics
 
 
## Problem
 
how a shock forms with the transport equation for a flood wave
 
 
# Model equations
 
The flood wave or kinetic wave, see  Whitham p 82,  Fowler  p 76,  or Chanson, is a long wave approximation
 of shallow water equations were inertial and pressure gradient are neglected.
 THe waves go only downstream.
  We have only mass conservation, aband momentum equationb reduces to a balance between friction and mild slope only:
 $$\frac{\partial h}{\partial t} +  \frac{\partial  Q}{\partial x}  =0\;\;\;\ \;\;\;0 =  - g h Z'_b   - c_f \frac{Q^2}{h^2}$$
 the slope $-Z'_b>0$ is constant, hence the flux is  $Q\sim h^{3/2}$.

This is the turbulent case, in the laminar one
 $Q= -gZ_b' h^3/(3 \nu)$ and mass equation is
$$\frac{\partial h} {\partial t} +   \frac{\partial }{\partial x}((\frac{-gZ_b'}{3\nu}) h_{\;}^3)   =0,$$
See example from [viscous collapse on a slope (Huppert's problem)](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c).
 
 
Then, for the turbulent case, without dimension we have to solve:
$$\frac{\partial}{\partial \bar t} \bar h  + \frac{\partial}{\partial \bar x}(\bar h^{3/2})= 0 $$
 
the method is similar to [advection](http://basilisk.fr/sandbox/M1EMN/BASIC/advecte1.c) which explains the notions of advection, testing the flux, coded with Basilisk.
 The flux is here $\bar Q=\bar h^{3/2}$, and we solve:
$$ \frac{\partial}{\partial \bar t} \bar h  + \frac{\partial}{\partial \bar x} \bar Q = 0
 \text{ solved as } \frac{\partial}{\partial \bar t} \bar h  + \bar c \frac{\partial}{\partial \bar x} \bar h = 0 \text{ with }  \bar c = \partial \bar Q/\partial \bar h$$
 
 
## Code
 mandatory declarations:
 */
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the field $h$, the flux $Q$, Boundary conditions
 */
scalar h[];
scalar Q[];

h[left] = neumann(0);
h[right] = neumann(0);

/**
 the flux is $\bar Q =\bar h^{3/2}$.
 Let wfrite it $\bar Q =\bar h^{m}$ so that we can test $m=1$ (pure advection) and $m=2$ (Burgers).
*/
double m=3./2;  // change for tests, m=1 advection, m=2 Burgers
double flux(double z)
{
    return pow(fabs(z),m);
}
/**
 the velocity $\bar c = \partial \bar Q/\partial \bar h$:
*/
double celerity(double z)
{
    return fmin(m*pow(fabs(z),m-1.),10000);
}
/**
 Main with definition of parameters, note that time step is small
 */
int main() {
    L0 = 20.;
    X0 = -4;
    N = 256*2;
    DT = (L0/N)/16;

    run(); 
}
/**
 initial elevation: a constant level plus a gaussian perturbation
 */
event init (t = 0) {
    foreach()
      h[] = .5+exp(-x*x) ;//fabs(x)<1;
}
/**
 begin the time loop, print data, in practice a max time of 5 is enough.
 */
event printdata (t += 1; t <=100)
{
    foreach()
    fprintf (stdout, "%g %g %g  \n", x, h[], t);
    fprintf (stdout, "\n\n");
}
/**
 integration step, at each time step
 */
event integration (i++) {
    double dt = DT;
    double cDelta = 1;
    /**
     finding the good next time step
     */
    dt = dtnext(dt);
    /**
     the algorithm is based on the flux.
     Approximation of the numerical flux taking into account
     $$Q_i = \frac{h_i^{3/2}+h_{i-1}^{3/2}}2 - c \frac{(h_i-h_{i-1})}2 $$
     */
    foreach()
    {
        cDelta = (celerity(h[])+celerity(h[-1]))/2.;
        Q[] = (flux(h[])+flux(h[-1]))/2.  - cDelta *(h[]-h[-1])/2;
    }

    /**
    explicit step update
     $$
     h_i^{n+1}=h_i^{n} -{\Delta t} \dfrac{Q_{i+1}-Q_{i}}{\Delta x}
     $$*/
    foreach()
      h[] +=  - dt* ( Q[1] - Q[] )/Delta;
}
/**
## Run
 Then compile and run, either by `qcc` either with `make`:
 
~~~bash
 qcc  -g -O2 -DTRASH=1 -Wall  floodwave.c -o floodwave ;./floodwave > out
 
 make floodwave.tst
 make floodwave/plots
 make floodwave.c.html
 
 source ../c2html.sh floodwave
~~~
 
## Results
 
`Time evolution with splot:
~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "t"
 set zlabel "h"
 sp [-5:10][0:5][0:2]'out' u 1:3:2 w l
~~~

 The equation is of advection type, with $\bar c = \partial \bar Q/\partial \bar h$ we have
 $$ \frac{\partial}{\partial \bar t} \bar h  + \frac{\partial}{\partial \bar x} \bar Q = 0
 \text{ is } \frac{\partial}{\partial \bar t} \bar h  + \bar c \frac{\partial}{\partial \bar x} \bar h = 0$$
In the $x,t$ plane, along the lines $\bar c = (d \bar x/d \bar t)$ the value of the  solution is constant.
 This constant is the value in  $t=0$, at the position $x=\xi$, given by function $h_0$.
 We write  $x=ct+\xi$, and  the value of the  solution is constant along this line.
 This constant value   along the line   $x-ct$, it is $h(x,t=0)=h(\xi,0)$.
Given at time  $t=0$ an initial  value $h(x,0)=h_0(x)$, then we construct the solution.
 Every  $x$ as an initial  $\xi$, so that the solution is
$$x=\xi +t c(\xi,0)\;\;\;\text{and}\;\; h(x,t)=h_0(x-t c(h_0(\xi))).$$

 On the graph we plot the numerical result and the analitycal one.
 
~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "h"
 h(x)=.5+exp(-x*x)
 c(x)=3./2*sqrt(h(x))
 set parametric
 set dummy x

 p[:][:]'out' u 1:($3<=5? $2:NaN) not w l ,x+1*c(x),h(x) not w l lc -1,\
 x+2*c(x),h(x) not w l lc -1,\
 x+3*c(x),h(x) not w l lc -1,\
 x+4*c(x),h(x) not w l lc -1,\
 x+5*c(x),h(x) not w l lc -1, x ,0 not w l lc -1
 ~~~
 
We see the formation of the shock (red computed values) when the solution is no more a function (analytical solution in black).

## As exercice

change tmax, and start by a heap  fabs(x)<1

show that solution is $H=4 \eta^2/9$ and $h=t^{-2/3}H$ and $\eta=x/t^{2/3}$


 

## Links
 
 * [http://basilisk.fr/sandbox/M1EMN/BASIC/advecte1.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c]().
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/floodwaveC.c]()
 
## Bibliography
 
 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 "Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC
 
 * G. B. Whitham "Linear and Nonlinear Waves" Wiley-Interscience, p 82
 
 * Fowler ["Mathematcis of the environment"](https://people.maths.ox.ac.uk/fowler/courses/mathenvo/envonotes.pdf)
 
 */
