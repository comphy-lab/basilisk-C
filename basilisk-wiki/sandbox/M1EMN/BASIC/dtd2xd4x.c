/**
# Resolution of growth equation

Explicit resolution of  equation which models the growth of an interface
$$\frac{\partial T}{\partial t}=  \phi_0 + 
\alpha  \frac{\partial^2 T}{\partial x^2}
 - \beta \frac{\partial^4 T}{\partial x^4}$$

the flux of material $\phi_0$ is zero for $x>0$, unit for $x<0$, there is diffusion and curvature effects

*/
#include "grid/cartesian1D.h"
#include "run.h"

scalar T[];
scalar dT4[];
scalar dT2[];
double dt;

int main() {
  L0 = 10.;
  X0 = -L0/2;
  N = 256;
/* very small space step due to the explicit scheme */  
  DT = (L0/N)*(L0/N)/2*(L0/N)*(L0/N)/2;
#define EPS 0.1  
  run();
}
/** 
initial surface : flat
*/
event init (t = 0) {
  foreach()
   T[] =  0;
  boundary ({T});
}
/** 
print data
*/
event printdata (t += 0.1; t <= 2) {
  foreach()
    fprintf (stdout, "%g %g %g %g\n", x, T[], t , dT4[]);
  fprintf (stdout, "\n\n");
}
/** integration 
*/
event integration (i++) {
  double dt = DT;
  

/**
finding the good next time step
*/
  dt = dtnext ( dt);
/** 
explicit step
$$\frac{ T(x+\Delta x) - 2 T(x) +T(x-\Delta x)}{\Delta x^2 } \simeq \frac{\partial^2 T}{\partial x^2}$$
 
 and the same for the fourth order derivative
*/
  foreach()
    dT2[] = ( T[-1,0] - 2 * T[0,0] + T[1,0] )/Delta/Delta;
  boundary ({dT2});
    
 foreach()
    dT4[] = ( dT2[-1,0] - 2 * dT2[0,0] + dT2[1,0] )/Delta/Delta;
  boundary ({dT4});
    /**
update
*/
  foreach()
    T[] += dt*( (x<0) + 0.1 * dT2[] - .25 * dT4[]);
  boundary ({T});
}
/**
Then compile and run:

~~~bash
qcc  -g -O2 -Wall  dtd2xd4x.c -o dtd2xd4x
./dtd2xd4x > out
~~~
or better 

~~~bash
 make dtd2xd4x.tst;make dtd2xd4x/plots    
 make dtd2xd4x.c.html ; open dtd2xd4x.c.html 
~~~



Note the overshoot in the growth due to the curvature term

~~~gnuplot surface
 p[-5:5][:2]'out' u 1:($3*($1<0)) t'pure growth' w l ,''  u  1:2 t'with dx^2 and dx^4 ' w l linec -1
~~~
 


ready for new site 09/05/19
*/
