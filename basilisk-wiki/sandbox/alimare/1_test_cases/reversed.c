
/**
# Time-reversed VOF and Level_set advection in a vortex

This classical test advects and stretches an initially circular
interface in a non-divergent vortical flow. The flow reverts in time
and the interface should come back to its original position. The
difference between the initial and final shapes is a measure of the
errors accumulated during advection.

We will need the advection solver combined with the VOF advection
scheme and the reinitialization function of the LS function. */

#define QUADRATIC 1
#define BGHOSTS 2
#include "utils.h"
#include "vof.h"
#include "advection.h"
#include "../LS_reinit.h"
#include "../basic_geom.h"
#include "view.h"
#define Pi 3.14159265358979323846


/**
The volume fraction is stored in scalar field `f` which is listed as
an *interface* for the VOF solver. The level set function is a *tracer* `dist`.


We do not advect any *level set* with
the default (diffusive) advection scheme of the advection solver.  
 */

scalar f[], dist[];
scalar * interfaces = {f}, * tracers = NULL;
scalar * level_set = {dist};

// vector velo[];
/**
Here are the parameters for the simulation. We use a narrow band (NB) approach 
meaning that the level set function has meaning only in the direct vicinity of 
the 0 value of the level set function. For this test case, the NB is made of 
only 4 cells. 
 */

int     MAXLEVEL;
int     nb_cell_NB =  1 << 4 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

int     N_display = 8;          // maximum level of refinement for which we 
                                // will display the results
double mass_ls_init, mass_vof_init;

double smin = 1.e-10;


/**
We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
/**
We then run the simulation for different levels of refinement. 
*/

  for (MAXLEVEL = 6; MAXLEVEL <= 8; MAXLEVEL++) {
    NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);
    init_grid (1 << MAXLEVEL);
    DT = 0.5*L0/(1<< MAXLEVEL);
    run();
  }
}

/**
The initial interface is a circle of radius 0.15 centered on
(0.5,0.75) . We use the levelset
function `circle()` to define this interface. 

*/

#define T 4.

coord center_circle ={0.5,0.75};
double Radius   =  0.15;

/**
We define the auxiliary levelset function $\phi$ on each vertex of the grid and
compute the corresponding volume fraction field. 

The level set function `dist` is taken positive out of the circle and we 
clamp the distance due to our NB approach. We take a 2\% overshoot that prevents
 NB cells from appearing due to spurious oscillations.
*/

#define circle2(x,y) (sq(0.15) - (sq(x - 0.5) + sq(y - .75)))

event init (i = 0) {
  fraction (f, circle2(x,y));

  foreach(){
    dist[] = -clamp(circle(x,y,center_circle,Radius),
     -1.2*NB_width, 1.2*NB_width);
  }
}

event velocity (i++) {

  /**
  This event defines the velocity field.
  
  On trees we first adapt the grid so that the estimated error on
  the volume fraction is smaller than $5\times 10^{-3}$. We limit the
  resolution at `MAXLEVEL` and we refine the volume fraction field
  `f` and on the level set function `dist`. */

#if 1
  adapt_wavelet ({f,dist}, (double[]){1.e-3, 1.e-4}, MAXLEVEL);
#endif
  
  trash ({u});

  double sign = 1.;
  if(t > T/2.)sign = -1;
  foreach_face(x)
    u.x[] = -sign*pow(sin(Pi*x),2.)*sin(2.*Pi*y);
  foreach_face(y)
    u.y[] = sign*pow(sin(Pi*y),2.)*sin(2.*Pi*x);    
  boundary ((scalar *){u});
}

/**
We use the advection solver for co-located variables from my sandbox.
*/
event LS_advection(i++,last){
  advection (level_set, u, dt);
}

/**
At the start and end of the simulation we check the sum, min and max
values of the volume fraction field. The sum must be constant to
within machine precision and the volume fraction should be bounded by
zero and one. */

event logfile (t = {0,T}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %.9f %g\n", t, s.sum, s.min, s.max);
  if(N == (1<<N_display) && t == 0) {
    mass_vof_init = s.sum;
  }
  scalar f2[];
  vertex scalar distn[];
  cell2node(dist,distn);
  fractions (distn, f2);
  s = statsf (f2);
  fprintf (stderr, "# %f %.12f %.9f %g\n", t, s.sum, s.min, s.max);
  if(N == (1<<N_display) && t == 0) {
    mass_ls_init = s.sum;
  }
}


/**
To compute the errors, we reinitialise field `e` and `e2` at the end of the
simulation with the initial shape and compute the differences with the
final shape. We output the norms as functions of the maximum
resolution `N`.  

Note that the error of the level set function is only studied in half of the NB width.
*/

event field (t = T) {
  scalar e[], e2[], loge[];
  fraction (e, circle2(x,y));
  foreach_vertex()
    e2[] = -circle(x,y,center_circle,Radius);
  foreach(){
    e[]  -= f[];
    if(e[]!=0.)
      loge[] = log(fabs(e[]));
    else
      loge[] = nodata;
    e2[] -= fabs(e2[])< 0.5*NB_width ? dist[] : e2[];
  }
  norm n  = normf (e);
  norm n2 = normf (e2);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, 
            n2.avg, n2.rms, n2.max);
  stats s = statsf(loge);
  fprintf(stderr, "#####%g %g\n", s.min, s.max);
  view(tx = -0.5, ty = -0.5);
  squares("loge");
  save("loge.png");
  dump();

  scalar l[];
    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_final.png", n = 400,min = -NB_width,
      max = NB_width);
}

/**
We also output the shape of the reconstructed interface at regular
intervals (but only on the finest grid considered). */

event shape (t += T/4.) {
  if (N == (1<<N_display)){
    output_facets (f);

    scalar cs[];
    face vector fs[];
    vertex scalar distn[];
    cell2node(dist,distn);
    fractions (distn, cs, fs);
    static FILE * fp2 = fopen ("out2","w");
    output_facets(cs, fp2);
  }
}

event mass (t += T/100.,last) {
  if (N == (1<<N_display)){
    static FILE * fp = fopen ("mass","w");
    scalar f2[];
    vertex scalar distn[];
    cell2node(dist,distn);
    fractions (distn, f2);
    stats s  = statsf(f2);
    stats s2 = statsf (f);
  /**
  We create a temporary VOF variable to calculate loss of mass during advection
  of the level set function */
  

    if(t!=0.)fprintf(fp,"%g %g %g\n", 
      t+dt,100.*(s.sum-mass_ls_init)/mass_ls_init,
      100.*(s2.sum-mass_vof_init)/mass_vof_init); 
  }
}

/**
If we are using adaptivity, we also output the levels of refinement at
maximum stretching. */

#if TREE
event levels (t = T/2) {
  if (N == (1<<N_display)) {
    scalar l[];
    foreach()
      l[] = level;
    output_ppm (l, file = "levels.png", n = 400, min = 0, max = MAXLEVEL);
  }
}
#endif


/**
Level set reinitialization event. The number of iteration of the LS_reinit 
function is set to  $10$.
*/


event LS_reinitialization(i++,last){
  LS_reinit(dist,it_max = 10);
}

#if 0
event movies ( i++,last){  
  if((MAXLEVEL == N_display) && (i%10 == 1)){
    view(tx=  -0.5, ty = -0.5);
    draw_vof("f");
    squares("dist", min =-0.02, max = 0.02);
    save ("visu.mp4");
  }
}

#endif

/**
## Results

We use gnuplot to compute the convergence rate
of the error norms with adaptation. 
Issues with the VOF function. (unexplained)

~~~gnuplot Convergence rates for constant- and adaptive grids

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2

f3(x)=a3+b3*x
fit f3(x) 'log' u (log($1)):(log($7)) via a3,b3
f4(x)=a4+b4*x
fit f4(x) 'log' u (log($1)):(log($5)) via a4,b4

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [32:512]
set xtics 32,2,512
set grid ytics
set cbrange [1:1]
set format y "%.1e"
plot 'log' u 1:4 t 'VOF max (cart)'  , exp(f(log(x)))       t ftitle(a ,b ), \
     'log' u 1:2 t 'VOF norm1 (cart)', exp(f2(log(x))) lw 3 t ftitle(a2,b2), \
     'log' u 1:7 t 'LS max (cart)'   , exp(f3(log(x)))      t ftitle(a3,b3), \
     'log' u 1:5 t 'LS norm1 (cart)' , exp(f4(log(x))) lw 3 t ftitle(a4,b4)

~~~

The shapes of the interface at $t=0$, $t=T/4$, $t=T/2$, $t=3T/4$ and
$t=T$ are displayed below for both sets of simulations (constant and
adaptive), for $N=N_{display}$. The shapes for $t=T/4$ should be identical to
those for $t=3T/4$ and similarly for $t=0$ and $t=T$ (for which we
measure the error). Note that the errors for $t=3T/4$ seem to be much
larger than those for $t=T$.

~~~gnuplot LS_interface
reset
set size ratio -1
unset format y
set xrange  [0:1]
set yrange  [0:1]
plot  'out2' w l t "LS" 
~~~
LS-Interface is slightly deformed at the end of the simulation.

~~~gnuplot VOF_interface
reset
set size ratio -1
set xrange  [0:1]
set yrange  [0:1]
plot  'out' w l t "VOF" 
~~~
VOF-Interface is streched to the limit of resolution of the VOF Model.

~~~gnuplot Mass loss during calculation
reset
set ylabel 'Mass error (%)'
set xlabel 'Time'
plot 'mass' u 1:2 t 'LS' w l lw 3, \
     'mass' u 1:3 t 'VOF' w l lw 3
~~~
Mass loss between 1-2%
*/
