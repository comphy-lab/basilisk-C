/**
# Non-trivial in and outflow. 

The usage of a [*Robin* boundary condition](robin.c) is explored for
non-trivial in and outflow.
 */
#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta))) + ((neumann (0))*((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))
/**
   The Navier-Stokes equations as solved in Two-Dimensions,
 */
#include "navier-stokes/centered.h"
#include "utils.h"
/**
The test aims to loosly advect periodic vortices into the domain. 
 */
#define U_LEFT  (cos(y)*sin(t)   + 3*sin(t))
#define U_RIGHT (cos(2*y)*sin(t) + 3*sin(t))
/**
Boundary conditions for the normal-flow component are set. For the
pressure(s), a condition is *guessed*.
 */
u.n[left]  = robin (1, 0.1, U_LEFT);
u.n[right] = robin (1, 0.1, U_RIGHT);
p[right]   = robin (1, 1, 0);
p[left]    = robin (1, 1, 0);
pf[right]  = robin (1, 1, 0);
pf[left]   = robin (1, 1, 0);

int main() {
  L0 = 2*pi;
  DT = 0.1;
  run();
}

event mov (t += 0.1) {
  stats f = statsf (u.x);
  printf ("%g %g\n", t, f.sum/f.volume);
  scalar omg[];
  vorticity (u, omg); 
  output_ppm (omg, file = "omg.mp4", n = 300,
	      min = -1, max = 1, map = cool_warm);
}

event stop (t = 5*2*pi) ;

/**
## Results

![Vorticity](rob_in_and_out/omg.mp4)

~~~gnuplot Flow rate
set grid
set xlabel 'time [-]'
set ylabel 'Mean horizontal velocity [-]'
set key right outside
plot 'out' t 'Result' , 3*sin(x) t 'Goal function'
~~~
*/
  
