
/**
# Time-reversed advection of a level-set function

This classical test advects and stretches an initially circular
interface in a non-divergent vortical flow. The flow reverts in time
and the interface should come back to its original position. The
difference between the initial and final shapes is a measure of the
errors accumulated during advection. We are in particular interested in the
ability of the thinning methd to extract a skeleton the interface

We will need the advection solver combined with the VOF advection
scheme and the reinitialization function of the LS function.


![Animation of the interface + skeleton](reversed_skeleton/skeleton.mp4)(loop)  
*/

#define QUADRATIC 1
#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
#include "embed.h"
#include "utils.h"
#include "vof.h"
#include "advection.h"
#include "../LS_reinit.h"
#include "../basic_geom.h"
#include "view.h"
#include "../LS_curvature.h"
#include "../thinning.h"
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

/**
Here are the parameters for the simulation. We use a narrow band (NB) approach 
meaning that the level set function has meaning only in the direct vicinity of 
the 0 value of the level set function. For this test case, the NB is made of 
only 4 cells. 
 */

int     MAXLEVEL;
int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
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

  for (MAXLEVEL = 7; MAXLEVEL <= 7; MAXLEVEL++) {
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
  boundary({dist});

  scalar curveLS[];
}

event velocity (i++) {

  /**
  This event defines the velocity field.
  
  On trees we first adapt the grid so that the estimated error on
  the volume fraction is smaller than $5\times 10^{-3}$. We limit the
  resolution at `MAXLEVEL` and we refine the volume fraction field
  `f` and on the level set function `dist`. */

#if 0
  adapt_wavelet ({f,dist}, (double[]){1.e-3, 1.e-3}, MAXLEVEL, minlevel = 4 );
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
  // dump();

  scalar l[];
    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_final.png", n = 400,min = -NB_width,
      max = NB_width);

    scalar curveLS[];
  vertex scalar distn[];
  scalar cs1[];
  face vector fs1[];

  curvature_LS(dist, curveLS);
  foreach(){
    if(fabs(dist[]) > 0.9*NB_width){
      curveLS[] = 0.;
    }
    else{
      curveLS[]=  curveLS[];
    }
  }
  boundary({curveLS});
  cell2node(dist,distn);
  fractions(distn, cs1, fs1);
  dump();
  // exit(1);


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



#if 1
event movies ( i++,last){  
  if((MAXLEVEL == N_display) && (i%10 == 1)){
    view(tx=  -0.5, ty = -0.5);
    draw_vof("f");
    squares("dist", min =-0.02, max = 0.02);
    save ("visu.mp4");
  }
}

event curvature_movie(i+=4, last){
  vertex scalar distn[];
  scalar cs1[];
  face vector fs1[];
  cell2node(dist,distn);
  fractions(distn, cs1, fs1);

  view(tx = -0.5, ty = -0.5);

/**
Skeleton part
*/
  scalar skeleton[];
  foreach(){
    if(dist[]> 0) skeleton[] = 1;
    else skeleton[] = 0;
  }
  boundary({skeleton});

  thinning2D(skeleton);
  draw_vof("cs1");
  squares("skeleton", min = 0, max = 1);
  save("skeleton.mp4");
}
#endif

/**
##See also

[Accuracy of this test case with a level-set function](/src/test/reversed.c)

[Accuracy of this test case with a level-set function](reversed.c)

*/