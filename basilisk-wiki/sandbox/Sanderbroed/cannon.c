/**
# A simple vortex cannon

On this page, we discuss and setup a vortex-ring cannon. This may help
to answer prominent [questions](README). 

## The setup and its dimensionless numbers.

A essence of a vortex-ring cannon consists of a round orifice
(introducing radius $R_o$), through which a fluid (e.g. air) is pushed
at a velocity $U_o$, for a time period $t_o$. The resulting vortex
structure owes it vorticity to the detachment of the thin viscous
boundary layer at the orifice' edge. As such, the fluid's viscousity
($\nu$) is an important system parameter. From the system parameters
$R_o,\ U_o,\ t_o \&\ \nu$, we can construct two dimensionless
numbers. First,

$$Re_o = \frac{U_oR_o}{\nu},$$

that we call the "orifice Reynolds number", and second,

$$\Pi_1 = \frac{t_oU_o}{R},$$ 

which is closely related to the socalled "stroke-bore ratio" or just
"stroke ratio". It is well known that $\Pi_1$ should not exceed 8 in
order omit a turbulent wake. It is important to realize that all
infinite combinations in the four-dimensional parameter space, span up
by $R_o, U_o, t_o \&\ \nu$, for which the values of $Re_o$ and $\Pi_1$
are the same, behave dyanmically identical.

If we were to introduce a cross-flow wind-speed magnitude, $U_c$, a new
dimensionless parameter could be constructed. An example is,

$$\Pi_2 = \frac{U_o}{U_c}.$$

## Numerical Setup 

The Navier-Stokes flow solver is used in three dimensions. Further, a
switch is placed for 3D visualization (see below). 
*/
#include "grid/octree.h" 
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h" //Trace performance
		      
#define BVIEW 1 //0 is off.
#if BVIEW
#include "view.h"
#endif
/**
Some free-parameter values are set, along side the two dimensionless
number values.
  */
double Uo = 1, Ro = 1;      // Normalized values
double Pi1 = 6, Re = 5000;  // Dimensionless numbers  
double to, nu;              // Dependend varuables

/**
the cannon injects in the left-to-right direction through a round
orifice in the $y-z$ plane
*/
#define RADIUS (sqrt(sq(y) + sq(z)))

u.n[left] = dirichlet ((t <= to)*(RADIUS <= Ro)*Uo); // Injection

u.n[right]  = neumann (0);   // free Outflow
p[right]    = dirichlet (0); 
/**
The maximum level or refinement is set and the domain size is chosen
large, such that the exact location of the boundaries does not affect
the dynamics.
 */
int maxlevel = 9;
int main() {
  L0 = 15*Ro;
  X0 = Y0 = Z0 = -L0/2;
  nu = Uo*Ro/Re;
  to = Pi1*Ro/Uo;
  const face vector muc[] = {nu, nu, nu};
  mu = muc;
  run();
}
/**
## Initialization

The grid is refined close to the opening.
*/
event init (t = 0) {
  refine (x < X0 + Ro    && RADIUS < 2*Ro   && level < maxlevel - 1);
  refine (x < X0 + Ro/2. && RADIUS < 1.1*Ro && level < maxlevel);
}
/**
## Output

We output a simple movie, depicting a slice of the simulation
 */
event movie (t += 0.1) {
  scalar omg[], lev[];
  vorticity (u, omg);
  boundary ({omg});
  foreach()
    lev[] = level;
  output_ppm (omg, file = "o.mp4", n = 300, linear = true, min = -1, max = 1);
  output_ppm (lev, file = "l.mp4", n = 300, min = 1, max = maxlevel);
}
/**
   If we use bview, we can visualize the vortex structure in 3D using
   the $\lambda_2$ vortex-detection algorithm.

![The chosen isosurface-value seems reasonable](cannon/l2.mp4)
*/
#include "lambda2.h"
double val = -1;

#if BVIEW
event bviewer (t += 0.1) {
  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});
  view (theta = 0.4, phi = 0.3);
  isosurface ("l2", val);
  translate (z = Z0)
    cells();
  save ("l2.mp4");
}
#endif

/**
The vortex structure's trajectory is also tracked. The coordinates are
printed to a file called `data`. 

~~~gnuplot A straight trajectory is diagnosed
  set xr [-7.5:0.5]
  set yr [-4:4]
  set zr [-4:4]
  set grid
  set view equal xyz
  splot 'data' u 2:3:4:1 palette t 'Time'
~~~

Distance from the opening as a function of time:

~~~gnuplot The velocity is quite constant.
  reset
  set xlabel ' tR_o/U_o'
  set ylabel 'distance from opening d/R_o'
  set key top left
  plot 'data' u 1:(($2+7.5)**2+$3**2+$4**2)**0.5 t 'data', x/2.6 t 'U_o t / 2.6'
~~~
 */
event track_and_trace (t += 0.5) {
  scalar l2[];
  lambda2 (u, l2);
  double xp = 0, yp = 0, zp = 0, tot = 0;
  foreach (reduction (+:xp) reduction (+:yp) reduction (+:zp) reduction (+:tot)) {
    if (l2[] < val) {
      tot +=   l2[]*dv();
      xp  += x*l2[]*dv();
      yp  += y*l2[]*dv();
      zp  += z*l2[]*dv();
    }
  }
  if (tot != 0) { 
    static FILE * fp = fopen ("data", "w");
    fprintf (fp, "%g %g %g %g\n", t, xp/tot, yp/tot, zp/tot);
  }
}
/**
## Grid adaptation

These adaptation values were readily derived from cheap experiments in
cylindrical coordinates.
 */
event adapt (i++) {
  double uc = Uo*0.04;
  adapt_wavelet ((scalar*){u}, (double[]){uc, uc, uc}, maxlevel);
}
/**
The simulation is stopped at some point.
 */
event stop (t = 4*to);

/**
## Running this example

The simulation' speed performance may benefit from multithreading. To
enable this feature, compile with the `-fopenmp` option. For using
three threads do something like;

~~~literatec
$qcc -O2 -Wall -fopenmp cannon.c -o cannon -lm
$OMP_NUM_THREADS=3 ./cannon
~~~

Or if you have bview installed, 

~~~literatec
$qcc -O2 -Wall -fopenmp cannon.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
$OMP_NUM_THREADS=3 ./cannon
~~~

## To do

Introduce a cross flow in a jet-perpendicular direction and charterize
the behaviour as a function of $\Pi_2$.
 */
