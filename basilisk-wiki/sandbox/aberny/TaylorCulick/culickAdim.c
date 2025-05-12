/**
# Taylor--Culick retraction of a liquid sheet

This is an example of 2D instability first introduced by Culick in [Culick,
1960] and numerically described by [Agbaglah, Josserand & Zaleski 2009].  This
instability is about the evolution of a viscous liquid sheet, in a 2D plane. 
Before launching the simulation, we decide the type of what we are simulating.
axi will be for an axisymmetric simulation, imposeVelocity will define the
boundary condition if we want to impose a veocity in the code. Finally,
gfsOutput is for an output into a gfs file.*/

#define axi 0
#define imposeVelocity 0
#define gfsOutput 0


#if axi
  #include "axi.h"
#endif

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

/**
We define the geometrical parameter, and the dimensionless time of the
simulation. */

#define LEVEL 10
#define L0 20
#define R1 0.5

double tEnd = 6;

#define Oh 1.

/**
We define the value of the fluid parameter, we will use them to compute the
Taylor-Culick velocity. */

#define rhoInt 1.
double muInt = Oh;
#define rhoRatio 100.
#define muRatio 55.
#if axi
  #define uCulick (4*muInt/(6*R1*rhoInt*Oh)*sqrt((2+sq(2*R1)/2)*1/(2*R1*M_PI)))
#else
  #define uCulick (muInt/(rhoInt*R1*Oh*sqrt(2)))
  // #define uCulick 0
#endif

/**
We setup boundary condition. Indeed, we want to simulate the retracting flow
from a moving observation point, moving at the Taylor Culick velocity. In our
case, this velocity is: $u_c = \mu_{Int}/(\rho_{Int}*R1*Oh*\sqrt{2})$ */

#if imposeVelocity
  u.n[left] = dirichlet(uCulick);
  p[left] = neumann(0.);
  pf[left] = neumann(0.);

  u.n[right] = neumann(0.);
  p[right] = dirichlet(0.);
  pf[right] = dirichlet(0.);

  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.);
  pf[top] = dirichlet(0.);
#endif

int main() {

  /**
  The domain will be $[-10, 10]\times[0:20]$ and will be resolved with
  $512\times 512$ grid points at the beginning. */

  size (L0);
  origin (-L0/2., 0.);
  init_grid (1 << LEVEL);

  /**
  Definition of the parameter of the fluids.  In this case, internal
  fluid for 1, external fluid for 2.  We define the surface tension from
  the Ohnesorge number. */

  rho1 = rhoInt, mu1 = muInt;
  rho2 = rhoInt/rhoRatio, mu2 = muInt/muRatio;
  f.sigma = 1.;
    
  run();
}

/**
We define the geometry with boolean operations. We define a
rectangular area (the intersection of Line and Line_vert). We also
define a half circle.  Then, our geometry will be the union of the
rectangle with the half circle. */

double geometry (double x, double y) {
  double Line = y-R1;
  double Line_vert = x-2.5;
  double Circle = sq(x-2.5)+sq(y)-sq(R1);
  double right_part = max(-Circle,-Line_vert);
  return min(-Line, right_part);
}

double uemax = 0.1;
double intemax = 0.01;

event init (t = 0) {

  /**
  We refine the original grid, to have a good initialise mesh. */
  double iteration = 0;
  do {
    iteration++;
    fraction(f, geometry(x,y));
  } while( adapt_wavelet({f,u}, (double []){intemax,uemax,uemax,uemax},
     maxlevel = LEVEL, 5).nf != 0 && iteration <= 10);

  #if imposeVelocity
    foreach()
      u.x[] = uCulick;
    fprintf(stderr, "Impose velocity: %f\nSigma: %f\n", uCulick, sq(mu1/Oh)*1./(rho1*R1*2.));
  #endif

  /**
  We output the initial situation*/

  static FILE * fp = fopen("initial.png", "w");
  output_ppm(f, fp, 256, min = 0, max = 1, box = {{-L0/2,0},{5,3}});
  static FILE * fp2 = fopen("initial.gfs", "w");
  output_gfs(fp2);
}

/**
We use an adaptative mesh with a maximum level of refinement of 9*/


event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){intemax,uemax,uemax,uemax}, LEVEL, 5);
}

/**
We output the tracer to have an idea of the evolution of the instability. */

event interface (i += 20) {
  static FILE * fp = popen ("ppm2mpeg > taylorCulick.mpeg", "w");
  static FILE * fp2 = popen ("ppm2mpeg > grid.mpeg", "w");
  output_ppm(f, fp, min=0, max=1);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm(l,fp2, min=5, max=LEVEL);
}

/**
We output the interface of the fluid to track the evolution of the
liquid bulge. */

event plotInterface (t += tEnd/10; t <= tEnd) {

  /**
  We open a new file for each selected time step. */
  
  char name[80];
  sprintf (name, "interface-%f.txt", t);
  FILE* fp = fopen (name, "w");
  
  /**
  We only output the interface in this function. */
  
  output_facets (f, fp);
}

/**
We output the maximum x position of the liquid finger. Indeed this
position should linearly evolve in time, with a few variations coming
from the capillary wave. We choose to have an output every 10 steps which 
should be enough to capture the mean velocity of the interface with a good
precision. */

double xPrev = -1, tPrev = -1;

event extractPosition (i ++) {
  static FILE * fp = fopen("interfacePosition.txt", "w");
  /**
  We define a new vector field, h. */
  
  vector h[];
  
  /**
  We reconstruct the height function field and take the corresponding 
  maximum along the x axis. */
  
  heights (f, h);
  double xMax = -HUGE;;
  foreach()
    if (h.x[] != nodata) {
      double xi = x + height(h.x[])*Delta;
      if (xi > xMax)
  xMax = xi;
  }
  double veloTip = tPrev >= 0 ? (xMax - xPrev)/(t - tPrev) : 0.;
  fprintf(fp, "%g %g %g\n", t, xMax, veloTip);
  fflush (fp);
  tPrev = t, xPrev = xMax;
}

/**
We output, in the standard output file, the step with the corresponding
time. */

event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
}

event end (t = tEnd) {}

/**
If gfsview is installed on the system, we can output the initial and
the final state of the simulation. */

#if gfsOutput
event gfsview (i += 100) {
  char name[80];
  sprintf (name, "taylorCulick-%05ld.gfs", i);
  FILE* fp = fopen (name, "w");
  output_gfs (fp);
}
#endif

/**
## Results

~~~gnuplot Evolution of the $x$ position of the end of the interface
set xlabel 'time'
set ylabel 'x position'
f(x) = a*x + b
fit [0.5:]f(x) 'interfacePosition.txt' via a,b
plot [0:] 'interfacePosition.txt' every 10 title "Basilisk Data",\
 f(x) t sprintf("%.2f t + %.2f", a, b)
~~~

According to Taylor and Culick, the retraction velocity should be $u_c
= \sqrt{2}$. With the linear fit, we obtain $u_c =
1.19$ (for $Oh=0.1$), which is close to the theoretical value. Note that this 
is only true for the 2D case.

If we display the evolution allong the time we obtain: 

![The sheet retraction](culickAdim/taylorCulick.gif)

We can also observe the mesh evolution:

![Evolution of the mesh](culickAdim/grid.gif)

## Comments on the results:

The retraction occurs at constant velocity. However the retraction velocity 
mesured for the 2D case is different than the theoretical velocity. More: if 
we simulate several case by changing the Oh number, we observe a strange
behaviour.

![Velocity of the retraction for the 2D case](data/retractionVelocity.jpg) 

The velocity seems to drop for high Oh number while we would expect that she 
goes up to the theoretical value from Taylor and Culick. This behaviour was 
expected. Indeed, the end time remain the same for each simulation, while the 
transiant regime from the beginning of the retraction to the equilibrium 
between inertia and capillarity increase in time. This transient regime is 
define by the viscous characteristic time $\tau_{vis} = \mu*H/2\gamma$. If the 
simulation time is high enough (idealy $30*\tau_{vis}$), the transient regime 
is pratically over.

We plot the evolution of the velocity.

~~~gnuplot Evolution of the velocity of the end of the interface
set xlabel 't^*'
set ylabel 'u^*'
set yrange [0:1]
plot 'interfacePosition.txt' u ($1/(0.1/2)):(-$3/sqrt(2)) title "retraction velocity"
~~~

Here, $t^* = t/\tau_{vis}$ and $u^*=u/u_c$. We can clearly observe the 
difference between the retraction velocity and the culick retraction velocity, 
for $Oh=0.1$. This difference between the 2 velocity was already observe in a 
previous article from Song et Al in 1999 (see [here](http://aip.scitation.org/doi/abs/10.1063/1.870113)).*/
