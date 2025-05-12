/**# Complex 3D vascular model simulation
Our goal is to simulate blood flow from a real 3D model, in this case we will use a real aortic arch model.

ATTENTION: This is not a complete code, you may need to manually adjust many values
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "lambda2.h"

#define MIN_LEVEL 4
#define LEVEL 6
#define MAX_LEVEL 8
#define eq_r 1.

#define tmax   15*2.*M_PI
#define alpha_w  10.


/**
##Importing the geometry 
see details in [distance.c](http://basilisk.fr/src/examples/distance.c)*/
void fraction_from_stl (scalar f, FILE * fp, int maxlevel)
{
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
    maxl = max.x - min.x;
    scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){1e-3*L0}, MAX_LEVEL).nf);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f); 
  view (fov = 19.9677, quat = {0.614306,0.390179,0.155159,0.66807}, 
        tx = -0.120804, ty = -0.340923, width = 640, height = 480);
  isosurface ("d", 0, color = "level", min = 5, max = 10);
  save ("stl.png");
}

/**
## volume fraction field declaire */
scalar f0[];
int main()
{
  init_grid (64);
  size (8.);
  origin (-1*L0/4 ,-2*L0/5, 0.);
  run();
}

/**
We then set an inflow velocity , this is a little tricky, we need to use the fact f0 >0 while inside the fraction field. 

Z - front z =L0 back z=0 

Y - top bottom 

X - left right
*/

u.n[back] = ( y >= 0)? dirichlet(10.*f0[]):neumann(0);
u.t[back] = ( y >= 0)? dirichlet(0):neumann(0);
u.r[back] = ( y >= 0)? dirichlet(0):neumann(0);

p[back] = ( y >= 0)? dirichlet(0):neumann(0) ;
pf[back] = ( y >= 0)? dirichlet(0):neumann(0) ;

p[front] = neumann(0) ;
pf[front] = neumann(0) ;
/**## Inital event
We read the STL file, and set the domain center and size using the value of model. */
event init (t = 0) {
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("aortic_arch.stl", "r");
    fraction_from_stl (f0, fp, LEVEL);
    f0.refine = f0.prolongation = fraction_refine;
    fclose (fp);
    //  Finally we display the surface reconstructed from volume fractions
    clear();
    draw_vof ("f0", edges = true, lw = 0.5);
    save ("vof.png");
    double viscosity =  sq(eq_r) / sq(alpha_w);
    const face vector muc[] = {viscosity, viscosity, viscosity};
    mu = muc;
  }
}
/** ## boundary condition
We use a simple (but crude) imposition of u=0 outside the fraction field. */

event bc (t <= tmax; i += 1) {
  //BC : lateral velocity = 0
  foreach()
    foreach_dimension()
    u.x[] = f0[]*u.x[];
    boundary ((scalar *){u , p});
  //position du test point
  double px = 3.5;
  double py = 0.8;
  double pz = 2.;
  fprintf(stderr, "%d %g %g %g \n", i, t, dt,interpolate(u.x , px, py, pz))
}
/**
We save a dumpfile in every 1000 iteration*/
event snapshot (i += 1000)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  scalar l2[];
  lambda2 (u,l2);
  dump (file = name);
}

/**##
We output the values of velocity and pressure in a choosen position at the last period with time step \delta t=0.2T*/
int j = 0;
event profill (t += 0.4*M_PI)
{
  if (t > tmax - 2.*M_PI){
    j += 1;
    char name2[80];
    sprintf (name2, "test-%d", j);	
    FILE *fp2;
    fp2 = fopen (name2, "w");
    double uyy = -0.5 ;
    double uxx = 2.5;
    double velou,velov,velo;
    for (double uzz = 0.5; uzz <= 4.; uzz += 4./100. ){
      velou = interpolate(u.x , uxx, uyy, uzz);
      velov = interpolate(u.y , uxx, uyy, uzz);
      velo = velou * cos(M_PI/4.) + velov * sin(M_PI/4.);
      fprintf(fp2,"%g %g %g %g %g\n", t, uzz-2.5, velo, velou, velov);
    }
  }
}

/**
Mesh adaption*/
event adapt (i++) {
  double uemax = 0.01;
  adapt_wavelet ({f0,u}, (double[]){0.01,uemax,uemax,uemax}, MAX_LEVEL, 4);
}