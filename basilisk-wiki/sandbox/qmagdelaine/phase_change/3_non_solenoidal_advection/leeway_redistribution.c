/**
# Non-solenoidal advection minimal working exemple 

This code tests the functions defined in
[mixtures.h](../mixtures.h) to advect a tracer
with respect to a velocity defined just at the interface.
We define a planar interface thanks to VOF tracer and a patch of a non-diffusive
tracer, thereafter called solute. The planar interface receeds and grabs the
solute.

![Advection of patch of solute by an interfacial velocity, concentration
field](leeway_redistribution/video_solute.mp4)
*/

#define VIDEO 1

/**
This code is called leeway.c because, with the advection scheme defined in
[vof.h](/src/vof.h), the patch of tracer moved a little bit along the interface
during the advection, instead of being just advected perpendicularly to it, see
[leeway_classic.c](leeway_classic.c).
The lateral displacement was around 5% of the advection velocity and depended
on the tilt of the interface. The variable *mytheta* defines this tilt. */

#define mytheta (17.5)
#define theta_rad (mytheta*pi/180.)

/**
We define the geometrical, temporal and resolution parameters: */

#define LEVEL 8
#define L 3.
#define T_END 1000.
#define DT_MAX 1.
#define SPEED 1e-3

/** 
[advection_Q.h](../advection_Q.h) essentially sets the order of the events,
[vof.h](/src/vof.h) is used to advect the interface and
[mixtures.h](../mixtures.h) to advect the solute.
*/

#include "../advection_Q.h"
#include "vof.h"
#include "curvature.h"
#include "../mixtures.h"
#include "view.h"

/**
We allocate several scalars and vector fields to describe both the
interface and the concentration field. */

scalar f[], solute[];
scalar * interfaces = {f}, * tracers = NULL;

/**
We add the diffusion coefficient as an attribute to the scalar fields. It is
defined in
[elementary_body.h](../elementary_body.h)
but this file is not included in this minimal working example. We
do not make diffuse any tracer here, but
[mixtures.h](../mixtures.h) needs it to compile
(even if the function which needs it is not used in this code). */

attribute {
  double D;
}

/**
In the main function, we set the domain geometry. */

int main() {
  N = 1 << LEVEL;
  size(2.*L);
  origin (- L, - L);
  DT = DT_MAX;
  init_grid (N);
  run();
}

/**
We define functions to set the initial position of both the interface and
the patch of solute. */

#define shape(x, y) (- x - tan(theta_rad)*y + 1.123e-2)
#define patch(x, y, R) (sq(R) - sq(x) - sq(y))

double initial_solute_mass, initial_tangential_position;

event init (i = 0) {

  /**
  We initialise the interface and to patch of solute.*/

  fraction (f, shape(x, y));
  fraction (solute, patch(x+0.2, y+0.2, 0.1));

  foreach_face()
    uf.x[] = 0.;
  boundary({f, solute, uf});
  
  /**
  We compute the initial mass of solute and the position of its center. */
  
  initial_solute_mass = statsf(solute).sum;
  coord center;
  foreach_dimension()
    center.x = 0.;
  foreach() {
    coord p = {x, y};
    foreach_dimension()
      center.x += dv()*solute[]*p.x;
  }
  foreach_dimension()
    center.x /= initial_solute_mass;
  
  /**
  We project the center of the patch on the direction parallel to the interface.
  */  
  
  initial_tangential_position = center.y*cos(theta_rad) - center.x*sin(theta_rad);
  
  CFL = 0.2;
}

/**
We define a uniform velocity normal to the interface. */

event stability (i++) {
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});
  normal_velocity (f, uf, speed = SPEED);
  boundary((scalar*){uf});
}

static scalar * interfaces_save = NULL;

event vof (i++) {

#if 1
  scalar f_previous[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    f_previous[] = f[];
  }
  boundary ({f, solute, f_previous, uf});
#else
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, solute, uf});
#endif
  
  /**
  We first advect just the vof tracer alone. */
   
  f.tracers = NULL;
  vof_advection({f}, i);
  boundary ({f});
  
  /**
  We redistribute now the solutes with the functions defined in
  [mixtures](../mixtures.h). */
  
  distribution (f, solute, uf);
  #if 1
  distribution_over (f, f_previous, solute, uf);
  #else
  distribution_over_2 (f, solute);
  #endif

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f, solute});

  /**
  We set the list of interfaces to NULL to prevent vof.h to do a second
  advection.*/

  interfaces_save = interfaces;
  interfaces = NULL;
}

event tracer_advection (i++) {
    interfaces = interfaces_save;
}

/**
## Post-processings and videos

We write A post-processing event to measure the loss of solute and the leeway of 
the center of the patch, and to save a video.
*/


event outputs (t = 0.; t += max(T_END/100., DT_MAX); t <= T_END) {
  double liquid_mass = statsf(f).sum;
  double total_solute = statsf(solute).sum;

  coord center;
  foreach_dimension()
    center.x = 0.;
  foreach() {
    coord p = {x, y};
    foreach_dimension()
      center.x += dv()*solute[]*p.x;
  }
  foreach_dimension()
    center.x /= total_solute;
    
  /**
  To measure the leeway of the patch of solute we project its center on the
  interface. */ 
  
  double leeway = center.y*cos(theta_rad) - center.x*sin(theta_rad);

  static FILE * fp = fopen("output_data", "w");
  fprintf (fp, "%.17g %.17g %.17g %.17g %.17g\n", t, initial_solute_mass,
           total_solute, liquid_mass, leeway-initial_tangential_position);
  fflush(fp);

  #if VIDEO
    /**
    solute */
    view (fov = 18., width = 640, height = 640, samples = 1,
          bg = {0.7, 0.7, 0.7});
    clear();
    squares ("solute", min = -1., max = 1., linear = false, map = cool_warm);
    draw_vof("f", lw = 1.5, lc = {0.05, 0.05, 0.05});
    char legend[100];
    sprintf(legend, "t = %0.2g", t);
    draw_string(legend, 1, size = 30., lw = 2.);
    save ("video_solute.mp4");  
  #endif
}

/**
# Results

We plot the gain (or loss if it is negative) of solute and the leeway of the
patch center.

~~~gnuplot Gain of solute
set style line 1 pt 7 ps 0.7
darkgray="#666666"
blue="#5082DC"
raspberry="#FA0F50"
set format y "%.1e"
plot 'output_data' u 1:($3 - $2) ls 1 lc rgb raspberry t "gain or loss"
~~~

~~~gnuplot leeway
plot 'output_data' u 1:5 ls 1 lc rgb blue t "leeway"
~~~

Both seem ok to me. These results have to be compared with the ones of
[leeway_classic.c](leeway_classic.c) where the scheme defined in [vof.h](/src/vof.h) has been used.
*/
