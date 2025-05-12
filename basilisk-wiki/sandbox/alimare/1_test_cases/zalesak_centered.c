/**
# Zalesak's notched disk

This classical test advects a notched disk from [Zalesak (1979)](https://www.sciencedirect.com/science/article/pii/0021999179900512). The interface should come back to its original position.
The difference between the initial and final shapes is 
a measure of the errors accumulated during advection.

We will need the advection solver combined with the VOF advection
scheme. */

#define QUADRATIC 1
#define quadratic(x,a1,a2,a3)           \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
#define BGHOSTS 2
#include "utils.h"
#include "vof.h"
#include "advection.h"
#include "../alex_functions.h"
#include "../LS_reinit.h"
#include "../simple_discretization.h"
#include "../basic_geom.h"
#include "view.h"
#define Pi 3.14159265358979323846

#define Nb_tour 2
#define T 2*Nb_tour*3.14
#define N_display 1 << 8

double NB_width;
int     nb_cell_NB =  1 << 4 ;  // number of cells for the NB

/**
The volume fraction is stored in scalar field `f` which is listed as
an *interface* for the VOF solver. We do not advect any tracer with
the default (diffusive) advection scheme of the advection solver. */

scalar f[], dist[],dist2[], dist3[],dist4[];
scalar * interfaces = {f}, * tracers = NULL;
scalar * level_set = {dist};

vector uc[];

uc.n[left]   = dirichlet(Pi/3.14*(0.5-y));
uc.n[right]  = dirichlet(Pi/3.14*(0.5-y));
uc.n[top]    = dirichlet(Pi/3.14*(x-0.5));
uc.n[bottom] = dirichlet(Pi/3.14*(x-0.5));

double geometry(double x, double y) {

  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = center_rectangle.x = 0.5;
  center_circle.y = 0.75;

  center_rectangle.y = 0.725;

  size_rectangle.x = 0.05;
  size_rectangle.y = 0.25; 

  double s = -circle (x, y, center_circle, 0.15);
  double r = -rectangle (x, y, center_rectangle, size_rectangle);

  double zalesak = difference(s,r) ;

  return zalesak;
}



int MAXLEVEL;

/**
We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
  
  /**
  We then run the simulation for different levels of refinement. */

  for (MAXLEVEL = 8; MAXLEVEL <= 8; MAXLEVEL++) {
    init_grid (1 << MAXLEVEL);
    run();
  }
}

/**
We define the levelset function $\phi$ on each vertex of the grid and
compute the corresponding volume fraction field. */

event init (i = 0) {
  DT = T/(1280*Nb_tour);
  NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);

  foreach(){
    dist[] = clamp(geometry(x,y),-1.2*NB_width, 1.2*NB_width);
    dist2[] = clamp(geometry(x,y),-1.2*NB_width, 1.2*NB_width);
    dist3[] = clamp(geometry(x,y),-1.2*NB_width, 1.2*NB_width);
    dist4[] = clamp(geometry(x,y),-1.2*NB_width, 1.2*NB_width);
  }
  fraction(f, geometry(x,y));
  boundary({f,dist,dist2,dist3,dist4});
  restriction({f,dist,dist2,dist3,dist4});


}

/**
This event defines the velocity field. The disk rotates around the center of
the grid (0.5,0.5).

On trees we first adapt the grid so that the estimated error on
the volume fraction is smaller than $5\times 10^{-3}$. We limit the
resolution at `MAXLEVEL` and we only refine the volume fraction field
`f`.
*/

event velocity (i++) {
  foreach(){
    dist[] = clamp(dist[],-1.2*NB_width, 1.2*NB_width);
    dist2[] = clamp(dist2[],-1.2*NB_width, 1.2*NB_width);
    dist3[] = clamp(dist3[],-1.2*NB_width, 1.2*NB_width);
    dist4[] = clamp(dist4[],-1.2*NB_width, 1.2*NB_width);
  }
  boundary({dist,dist2,dist3,dist4});
  restriction({dist,dist2,dist3,dist4});
#if TREE
  adapt_wavelet ({f,dist,dist2,dist3,dist4}, (double[]){1.e-3, 1.e-4,1.e-4,1.e-4,1.e-4}, MAXLEVEL, 
    MAXLEVEL-2,list= {f,dist,dist2,dist3,dist4});
#endif
  trash ({u,uc});


  foreach_face(x)
    u.x[] = Pi/3.14*(0.5-y);

  foreach_face(y)
    u.y[] = Pi/3.14*(x-0.5);


  foreach(){
    uc.x[] = Pi/3.14*(0.5-y);
    uc.y[] = Pi/3.14*(x-0.5);
  }
  boundary ((scalar *){u,uc});
  restriction((scalar *){uc});
}

/**
At the start and end of the simulation we check the sum, min and max
values of the volume fraction field. The sum must be constant to
within machine precision and the volume fraction should be bounded by
zero and one. */

event logfile (t = {0,T}) {
  stats s = statsf (f);
  fprintf (stderr, "#REZ %f %.12f %.9f %g \n", t, s.sum, s.min, s.max);
  stats s2 = statsf (dist);
  fprintf (stderr, "#REZ %f %.12f %.9f %g\n", t, s2.sum, s2.min, s2.max);
}


/**
To compute the error, we reinitialise field `e` at the end of the
simulation with the initial shape and compute the difference with the
final shape. We output the norms as functions of the maximum
resolution `N`. */

event field (t = T) {
  scalar e[], e2[];

  foreach()
    e2[] = geometry(x,y);
  boundary ({e2});
  face vector s[];
  fractions (e2, e, s);

  boundary({dist});

  foreach(){
    e[]  -= f[];
    e2[] -= fabs(dist[])< NB_width/2. ? dist[] : e2[];
  }
  norm n  = normf (e);
  norm n2 = normf (e2);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, 
            n2.avg, n2.rms, n2.max);
}

event LS_advection(i++,last){
  advection (level_set, u, dt);
  RK2 (dist2, uc, dt, NB_width);

  scalar dist_copy[];
  foreach()
    dist_copy[] = dist3[];
  boundary({dist_copy});
  restriction({dist_copy});
  FE (dist3, dist_copy, uc, dt, NB_width);

  RK3_WENO5 (dist4, uc, dt, NB_width);
}

event LS_reinitialization(i++,last){
  LS_reinit(dist, it_max = 10);
  LS_reinit(dist2, it_max = 10);
  LS_reinit(dist3, it_max = 10);
  LS_reinit(dist4, it_max = 10);
}

/**
We output the interfaces of the vof variable.
*/
event shape (t = {0,T}) {
    output_facets (f);

    scalar cs[];
    face vector fs[];
    vertex scalar distn[];
    cell2node(dist,distn);
    fractions (distn, cs, fs);
    static FILE * fp2 = fopen ("out2","w");
    if(t>0)output_facets(cs, fp2);
    
    cell2node(dist2,distn);
    fractions (distn, cs, fs);
    static FILE * fp3 = fopen ("out3","w");
    if(t>0)output_facets(cs, fp3);

    cell2node(dist3,distn);
    fractions (distn, cs, fs);
    static FILE * fp4 = fopen ("out4","w");
    if(t>0)output_facets(cs, fp4);

    cell2node(dist4,distn);
    fractions (distn, cs, fs);
    static FILE * fp5 = fopen ("out5","w");
    if(t>0)output_facets(cs, fp5);
}

event movie(i+=20,last){
  view(tx = -0.5, ty = -0.5);
  squares("dist");
  save("dist.mp4");
  squares("dist2");
  save("dist2.mp4");
  squares("dist3");
  save("dist3.mp4");
  
  scalar cs[];
  face vector fs[];
  vertex scalar distn[];
  cell2node(dist4,distn);
  fractions (distn, cs, fs);
  squares("dist4");
  draw_vof("cs", "fs");
  save("dist4.mp4");
}

/**
## Results

~~~gnuplot Shapes of the interface - VOF variable
reset
set size ratio -1
plot [0.25:0.75][0.5:1]'out' w l t "VOF interface"
~~~

~~~gnuplot Shapes of the interface - LS variable
reset
set size ratio -1
plot [0.25:0.75][0.5:1]'out2' w l t "LS BCG",\
  'out3' w l t "LS RK2 - WENO2",\
  'out4' w l t "LS FE",\
  'out5' w l t "LS RK3 - WENO5"
~~~

Interface is slightly more deformed for level-set function.

 */
