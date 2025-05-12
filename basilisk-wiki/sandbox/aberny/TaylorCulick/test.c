/**
# Taylor-Culick Retraction of a Liquid Sheet

This is an example of 2D instability first introduced by Culick in [Culick,
1960](/sandbox/aberny/ref.bib#culick1960) and by Taylor in [Taylor, 1959](/sandbox/aberny/ref.bib#taylor1959).
This instability is about the evolution of a viscous liquid sheet, in a 2D plane.

The theoretical velocity of retraction should be $\sqrt(2\sigma/\rho H)$ where 
H is the total height of the liquid sheet.

We non-dimensionlise the parameters. Since we are working with a 2-phase 
flow, the liquid sheet will be the reference for the dimensionless value. The 
liquid density $\rho_1$ is 1, the surface tension $\sigma$ is 1 and the height 
of the liquid sheet. The viscosity $\mu_1$ is then equal to the Ohnesorge 
number ($Oh = \mu/\sqrt(\rho\sigma H)$).

The retraction velocity is not reaching instantaneously. There is indeed a 
transient regime. The characteristic time scale is called the viscous time 
scale. Its definition is: $\tau_{vis} = \mu H/2\sigma$.

With the dimensionless parameters, the values are: 
\begin{itemize}
  \item $u_c=\sqrt(2)$
  \item $\tau_{vis} = Oh/2$
\end{itemize}

Bu default, we will choose $Oh=0.1$ which lead to $\tau_{vis}=0.05$. To 
overpass practically all the transient regime, we need at least 
$t_{end}=30\tau_{vis} = 1.5$. We choose $t_{end} = 6$


In this simulation, there is the possibility to run the retraction of an 
axisymmetric liquid finger or a semi-infinite liquid sheet, with or without an 
impose velocity. This velocity corresponds to the theoretical velocity from 
Taylor and Culick (in 2D or in 3D axisymmetric).*/

#define axi 0
#define imposeVelocity 1
#define gfsOutput 0


#if axi
  #include "axi.h"
#endif

#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

/**
The level of refinement is 10. The size of the computation domain is 20.*/

#define LEVEL 10
#define L0 20

/**
We are simulating only one half of the liquid sheet.*/

#define R1 0.5

double tEnd = 6;


/**
We define the parameters of the 2 fluids in a dimensionless way.*/

#define rhoInt 1.

/**
The ratios are the one for a water sheet which is retracting in the air. */

#define rhoRatio 100.
#define muRatio 55.

int main(int argc, char *argv[]) {

  /**
  The domain will be $[-10, 10]\times[0:20]$ and will be resolved with
  $1024\times 1024$ grid points at the beginning. */

  size (L0);
  origin (-L0/2., 0.);
  init_grid (1 << 9);

  /**
  By default the Ohnesorge number is equal to 0.1, but we can change that
  with an input argument.*/

  double Oh = 0.1;
  if (argc >=2)
    Oh = atof(argv[1]);

  /**
  We define the parameters of the fluids. The internal fluid (the water) is 
  the fluid 1. The external fluid (the air) is the fluid 2.*/

  double muInt = Oh;

  rho1 = rhoInt, mu1 = muInt;
  rho2 = rhoInt/rhoRatio, mu2 = muInt/muRatio;
  f.sigma = 1.;
    
  run();
}

/** 
Depending on the case (an axisymmetric simulation or not), the Culick velocity 
is different. */

#if axi
  #define uCulick ((4./3.)*sqrt(5./(2.*M_PI)))
#else
  #define uCulick (sqrt(2.))
#endif

/**
If we impose a velocity to the fluid, then we have to set the boundary 
conditions. We allow the fluid to exit the domain on all the border, except 
the bottom one.*/

#if imposeVelocity
  u.n[left] = dirichlet(uCulick);
  p[left] = neumann(0.);

  u.n[right] = neumann(0.);
  p[right] = dirichlet(0.);

  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.);
#endif


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

/**
We are simulating the following geometry:

~~~gnuplot Initial profile of the retracting liquid
set xlabel 'x'
set ylabel 'y'
set size ratio -1
plot 'interface-0.000000.txt' w l title "initial profil"
~~~
*/

event init (t = 0) {

  /**
  We initialise the geometry of the interface.*/

  fraction(f, geometry(x,y));

  /**
  If we also want to impose an incoming velocity to the fluid, then we also 
  impose this velocity to all the domain.*/

  #if imposeVelocity
    foreach()
      u.x[] = uCulick;
  #endif
}

/**
We use an adaptive mesh with a maximum level of refinement of 10. We adapt 
the mesh with respect to the interface the velocity.*/


event adapt (i++) {
  double uemax = 0.05;
  adapt_wavelet ({f,u}, (double[]){0.005,uemax,uemax,uemax}, LEVEL, 5);
}

/**
We output the tracer to have an idea of the evolution of the instability. */

event interface (t += 0.01) {
  static FILE * fp1 = popen ("ppm2mp4 interface.mp4", "w");
  view (fov = 5.2, tx = 0.0111485, ty = -0.027554, width = 612, height = 144);
  clear();
  box();
  draw_vof ("f");
  save (fp = fp1);
  dump ("dump");
  
  static FILE * fp2 = popen ("ppm2mp4 grid.mp4", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp2, 512, min=5, max=LEVEL, box = {{X0,Y0},{X0+L0,Y0+L0/4}});
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
from the capillary wave. */

double xPrev = -1, tPrev = -1;

event extractPosition (i ++) {

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

  /**
  We also output the velocity of the end of the bulge.*/

  double veloTip = tPrev >= 0 ? (xMax - xPrev)/(t - tPrev) : 0.;
  fprintf(stderr, "%g %g %g\n", t, xMax, veloTip);
  fflush (stderr);
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
## Results

~~~gnuplot Evolution of the $x$ position of the end of the interface
set xlabel 'time'
set ylabel 'x position'
f(x) = a*x + b
fit [4:]f(x) 'log' u 1:2 via a,b
plot [0:] 'log' u 1:2 every 10 title "Basilisk Data",\
 f(x) t sprintf("%.2f t + %.2f", a, b)
~~~

TTo avoid the most part of the transient regime, we fit a linear curve for $t>
4$. We obtain a retraction velocity of $1.22$. The difference with the 
theoretical results is about $15\%$. This difference was observed by 
[Savva, 2009](/sandbox/aberny/ref.bib#savva2009) but also by [Song, 1999](/sandbox/aberny/ref.bib#song1999).

If we observe the evolution allong the time we obtain: 

<p><center>
<video width="612" height="144" controls>
  <source src="culickTest/taylorCulick.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
<br>
The sheet retraction.
</center></p>

We can also observe the mesh evolution:

<p><center>
<video width="512" height="128" controls>
  <source src="culickTest/grid.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
<br>
Evolution of the mesh.
</center></p>

![Evolution of the mesh](culickTest/grid.mp4)

The mesh is not saturate, our refinement criteria are good.

## Comments on the Results

The retraction seems to occur at a constant velocity. But this is not really 
the case. We are still in the transient regime, even if we spend more than 30 
times the viscous time scale.

We plot the evolution of the velocity, by rescaling it with the Taylor-Culick 
velocity.

~~~gnuplot Evolution of the velocity at the end of the interface
set xlabel 't^*'
set ylabel 'u^*'
set yrange [0:1]
plot 'log' u ($1/(0.1/2)):((sqrt(2)-$3)/sqrt(2)) title "retraction velocity"
~~~

Here, $t^* = t/\tau_{vis}$ and $u^*=u/u_c$. We can clearly observe the 
difference between the retraction velocity and the culick retraction velocity, 
for $Oh=0.1$. */
