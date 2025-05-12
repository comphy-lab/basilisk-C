/**
# Level_set advection of a circle

We simply advect a level set function

We will need the advection solver combined with the VOF advection
scheme and the reinitialization function of the LS function. */

#include "utils.h"
#include "advection.h"
#include "vof.h"
#include "../level_set.h"
#include "basic_geom.h"




/**
The volume fraction is stored in scalar field `f` which is listed as
an *interface* for the VOF solver. The level set function is a *tracer* `dist`.


We do not advect any *level set* with
the default (diffusive) advection scheme of the advection solver.  
 */

scalar f[], dist[];
scalar * interfaces = {f}, * tracers = {dist};
scalar * level_set = NULL;
/**
Here are the parameters for the simulation. We use a narrow band (NB) approach 
meaning that the level set function has meaning only in the direct vicinity of 
the 0 value of the level set function. For this test case, the NB is made of 
only 4 cells. 
 */
int     MAXLEVEL;
int     nb_cell_NB =  8 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
#define T  2*10*(1 << 8)

                                // will display the results
double mass_ls_init, mass_vof_init;

/**
We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
  origin (-L0/2.,-L0/2.);
  
  /**
  We then run the simulation for different levels of refinement. */

  
  int j;
  for (j=0;j<=1;j++){
    MAXLEVEL = 6 + j;
    periodic(right);
    periodic(top);
    NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);
    init_grid (1 << MAXLEVEL);
    run();
  }
}




coord center_circle ={0.,0.};
double Radius   =  0.25;

/**
We define the auxiliary levelset function $\phi$ on each vertex of the grid and
compute the corresponding volume fraction field. 

The level set function `dist` is taken positive out of the circle and we 
clamp the distance due to our NB approach. We take a 2\% overshoot that prevents
 NB cells from appearing due to spurious oscillations.
*/

event init (i = 0) {
  fraction (f, circle(x,y,center_circle,Radius));

  foreach(){
    dist[] = -clamp(circle(x,y,center_circle,Radius),
     -1.02*NB_width, 1.02*NB_width);
  }

  if (MAXLEVEL  == 7) {
    scalar l[];
    // foreach()
    //   l[] = level;
    // output_ppm (l, file = "levels_init.png", n = 400, min = 0, max = 7);

    foreach()
      l[] = f[];
    output_ppm (l, file = "f_init.png", n = 400, min = 0, max = 1);

    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_init.png", n = 400, min = -NB_width,
      max = NB_width);

  }

  foreach_face(){
    u.x[] = 0.25/(1 << MAXLEVEL);
  }

}

/**
At the start and end of the simulation we check the sum, min and max
values of the volume fraction field. The sum must be constant to
within machine precision and the volume fraction should be bounded by
zero and one. */

event logfile (t = {0,T}) {
  stats s = statsf (f);

  fprintf (stderr, "# %f %f %f %g\n", t, s.sum, s.min, s.max);
  if(N == (1<<MAXLEVEL)) {
    mass_vof_init = s.sum;
  }
  scalar f2[];
  fractions (dist, f2);
  s = statsf (f2);
  fprintf (stderr, "# %f %.12f %.9f %g\n", t, s.sum, s.min, s.max);
  if(N == (1<<MAXLEVEL)) {
    mass_ls_init = s.sum;
  }
  if (N == (1<<MAXLEVEL))
    output_facets (f);

  if(MAXLEVEL == 7)exit(1);
}

/**
To compute the errors, we reinitialise field `e` at the end of the
simulation with the initial shape and compute the differences with the
final shape. We output the norms as functions of the maximum
resolution `N`.  

Note that the error of the level set function is only studied in half of the NB width.
*/

event field (t = T) {
  scalar e[], e2[];
  fraction (e, circle(x,y,center_circle,Radius));
  foreach()
    e2[] = -circle(x,y,center_circle,Radius);
  foreach(){
    e[]  -= f[];
    e2[] -= fabs(e2[])< 0.5*NB_width ? dist[] : e2[];
  }
  norm n  = normf (e);
  norm n2 = normf (e2);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, 
    n2.avg, n2.rms, n2.max);
  if(MAXLEVEL == 8){
    scalar l[];
    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_final.png", n = 400,min = -NB_width,
      max = NB_width);
    foreach()
      l[] = f[];
    output_ppm (l, file = "f_reversed.png", n = 400, min = 0, max = 1);
  }
}

/**
Level set reinitialization event. The number of iteration of the LS_reinit2 
function is set to  $1.4 \times 2 \times nb_{cell NB}$ which is a bit more than 
the $NB_{width}$
*/
event LS_reinitialization(i++,last){
  if(i>0 && i%2==1){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      4);
  }
}


/**
## Results

~~~gnuplot Phase change velocity
set logscale
plot 'log' u 1:7
~~~



|                          |                          |
:-------------------------:|:-------------------------:
| Initial Level Set        |  Final Level set         |
|![Initial](test_advection/dist_init.png)(width="400" height="300") |  ![Final](test_advection/dist_final.png)(width="400" height="300")  |
| Initial VOF              |  Final VOF               |
| ![Initial](test_advection/f_init.png)(width="400" height="300") |  ![Final](test_advection/f_reversed.png)(width="400" height="300")  |
*/
