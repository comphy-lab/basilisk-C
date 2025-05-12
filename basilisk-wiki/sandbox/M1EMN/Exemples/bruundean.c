/**
# Resolution of the shape of a Beach by Bruun 64 & Dean 91 model

It is believed that the shore cross profile is a power law $Z=Ax^n$ where $x$ is distance offshore ($x<100m$), $0<Z<10m$ is water depth, 
averaged measurements give $n=0.66$ (in fact $0.03<x<1.4$ see Pilkey 94, we did himself a mistake in the equation). Bruun obtained that from fits of measurements.

To establish this,  Dean supposes that there is a uniform wave energy dissipation constant per unit volume $D_*$, and that within the breaking zone the flux of energy is steady and dissipated at a rate $D_*$:
$$D_*= \frac{1}{h}\frac{d}{dx}(E C_g)$$
where $E$ and $C_g$ are the energy density and group velocity.
The waves come from the right to the left, we suppose that the mean sea level defines level 0, 
as $E \propto Z^2$ and $C_g \propto \sqrt{Z}$, hence 
the equation  without dimension is :
$$ \frac{d}{dx}(Z^{5/2})= Z $$
we solve it with a relaxation as:
$$ \frac{\partial}{\partial x}Z  + \frac{\partial}{\partial x}(Z^{5/2})= Z $$
the method is similar to  [http://basilisk.fr/sandbox/M1EMN/BASIC/advecte1.c]() which explains the notions of advection, testing the flux, coded with Basilisk (see the end, maybe $\frac{\partial}{\partial x}Z^2 $ is better).

Starting by a slope and a flat bottom, we wil compute evolution to the steady power law 
~~~gnuplot 
set xlabel "x"
set ylabel "depth"
set arrow from 3,-(3.*3./5)**(2./3) to 3,0
set label "depth -Z(x)" at  3.5,-1
set key left bottom
p[0:10][]  -0.25*(x<1? x:1) t'initial', -(3*x/5)**(2./3) t'expected' w l  linec 3 
~~~
 Note that here $Z$ is positive, in regular Shallow Water configuration we use $z_b$ which is $-Z$, as here $\eta=0$.
 
##Code
mandatory declarations:
*/
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the field Z, the flux F, Boundary conditions 
*/
scalar Z[];
scalar F[];
 
Z[left] = dirichlet(0);
Z[right] = dirichlet(5);;  

double m=5./2;  // change for tests, m=1 advection, m=2 Burgers
double flux(double z)
{ 
  return pow(fabs(z),m);
}
double celerity(double z)
{ 
  return fmin(m*pow(fabs(z),m-1.),10000);
}
/**
Main with definition of parameters, note that time step is small
*/
int main() {
  L0 = 10.;
  X0 = 0;
  N = 256;    
  DT = (L0/N)/16/4; 
  run();
 }
/** 
initial elevation: a "slope and a flat bottom"
*/
event init (t = 0) {
  foreach(){
    Z[] = .25*fmin(x,1);
  }
  boundary ({Z});
  }
/**
 begin the time loop, print data, in practice a max time of 5 is enough.
*/
event printdata (t += .5; t <=5)  
{
  foreach()
    fprintf (stdout, "%g %g %g  \n", x, Z[], t);
  fprintf (stdout, "\n\n");
}
/** 
integration step, at each time step
*/
event integration (i++) {
  double dt = DT;
  double cDelta;
/**
finding the good next time step
*/
  dt = dtnext(dt);
/** 
the algorithm is based on a split in time, the flux and the source term.
*/    
  foreach()
   {
    cDelta = (celerity(Z[])+celerity(Z[-1]))/2.;
    F[] = (flux(Z[])+flux(Z[-1]))/2.  - cDelta *(Z[]-Z[-1])/2;
    }
    boundary ({F});
/** 
split: explicit step update 
$$
 Z_i^{n+1}=Z_i^{n} -{\Delta t} \dfrac{F(Z_{i+1})-F(Z_{i})}{\Delta x}
$$*/ 
  foreach()
    Z[] +=  - dt* ( F[1] - F[] )/Delta;
  boundary ({Z});  
/**
  split: 
  implicit resolution of 
  $\dfrac{\partial Z}{\partial t} = Z$ :
 $$
 Z_i^{n+1}=Z_i^{n} +{\Delta t} Z_i^{n+1}
$$ 
 */  
  foreach()
    Z[] = Z[]/(1 - 1.0 * dt); //Z[]*(1 + 1.0 * dt);
/**
 Boundary condition
*/
  boundary ({Z});
 } 
/**
##Run
Then compile and run:

~~~bash
qcc  -g -O2 -DTRASH=1 -Wall  bruundean.c -o bruundean ;./bruundean > out
~~~
##Results
The analytical solution solves the equilibrium flux/dissipation
$$dZ(x)^{m}/dx =  Z  $$
so that 
$Z=(\frac{(m-1)}{m} x)^{\frac{1}{m-1}}$
with m=5/3, so that we obtain
$Z(x)= (3 x/5)^{2/3}$
the Bruun's solution.

Note that there is obviously a problem of B.C. at the right:


~~~gnuplot
 set xlabel "x"
 set ylabel "-Z"
 set key left bottom
 p[:][-3.5:0]'out' u ($1):(-$2)t'num'w l,-(3*x/5)**(2./3) t'exact' w l linec 3 
~~~

which is  $-Z(x,t)$ plotted here for t=0 .5 1 1.5  2 ...


time evolution with splot:
~~~gnuplot
 set xlabel "x"
 set ylabel "t"
 set zlabel "Z"
 sp [:][:][0:5]'out' u 1:3:2 w l,(3*(x)/5)**(2./3) t'exact' w l    
~~~

## Exercise
  
Solve now with $\partial_t Z^2$ instead of $\partial_t Z$, to have an energy like equation in $Z^2$
  
$$ \frac{\partial}{\partial x}Z^2  + \frac{\partial}{\partial x}(Z^{5/2})= Z $$

which is as well, $m=5/2$
  
  $$ \frac{\partial}{\partial x}Z + \frac{\partial}{\partial x}(\frac{m}{2 (m-1)}Z^{m-1})= \frac{1}{2} $$


##Bibliography

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC

 
* [Robert G. Dean](http://journals.fcla.edu/jcr/article/download/78405/75816)
"Equilibrium Beach Profiles: Characteristics and Applications" 
Journal of Coastal Research, Vol. 7, No. 1, 1991

* [O. Pilkey, Jr.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.527.8205&rep=rep1&type=pdf)
"Mathematical Modeling of Beach Behavior Doesn't Work"
 
*/
