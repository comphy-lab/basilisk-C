/**
# 3D aorta-tube-endosize womersley flow
Our goal is to simulate womersley flow from a real 3D internal-aortic tube, this model should be available for all almost straight models
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#define MIN_LEVEL 5//the control of refine level
#define LEVEL 6
#define MAX_LEVEL 8
#define tmax 10*2.*M_PI
#define alpha_w   10.

/**
##Importing the geometry 
see details in [distance.c](http://basilisk.fr/src/examples/distance.c) */
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
  while (adapt_wavelet ({d}, (double[]){1e-3*L0}, LEVEL).nf);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f);
  view (fov = 27.1723, quat = {-0.707795,-0.108707,-0.69171,0.0934929}, 
        tx = -0.545947, ty = -0.00320332, bg = {0.3,0.4,0.6}, width = 640, height = 320);
  isosurface ("d", 0, color = "level", min = 5, max = 10);
  //save ("endotubestl.png");
}


/**
##main fonction */
//fonction & field declaire
scalar f0[];
double cradius;

int main()
{
  //We set the domain center and it's size (self-defining or shrinking the STL model is needed)
  init_grid (32);
  size (10.);
  origin (-L0/2 ,-L0/2, 0.);
  run();
}

/**
the definition use in Basilisk:
Z - front z =L0 back z=0 &
Y - top bottom &
X - left right
*/

/**
## boundary conditions
the main axi direction here is z*/
//not sure it's useful
//u.n[front] = neumann(0);
//u.t[front]= dirichlet(0);
//u.n[back] = neumann(0);
//u.t[back]= dirichlet(0);
p[front] = dirichlet(0);


/**
##initial event */
event init (t = 0) {
  //read the STL file 
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("endotube.stl", "r");
    fraction_from_stl (f0, fp, LEVEL);
    f0.refine = f0.prolongation = fraction_refine;
    fclose (fp);

    /**
We display the surface reconstructed from volume fractions. */
    clear();
    draw_vof ("f0", edges = true, lw = 0.5);
    //save ("endotubevof.png");

    // initial velocity in order to avoid possible error(not sure really useful)
    foreach() 
      u.z[] = 0.01*f0[];
    boundary ({u.z});
  }

  /**
We calculate the total volume and the equivalent radius */
  //calculate volume
  double volume = statsf(f0).sum;
  fprintf(stderr, "total volume : %g \n", volume);
  //calculate caracteristic radius
  double eq_r = sqrt( volume/ (L0 *M_PI) ) ;
  cradius = 1.0 / eq_r;
  fprintf(stderr, "caracteristic radius :%g \n \n \n",cradius);
  /**
In order to compare with the theoretical results given by the adimensional equation, we need to inverse calculate the value of \ mu in our case. */
  double viscosity = sq(cradius) / sq(alpha_w);
  const face vector muc[] = {viscosity, viscosity, viscosity};
  mu = muc;
}

/**
##accelerate event
We ossilate the system periodicly. Here we take $T = 2 \Pi$ to match the adimensionalization.
*/
const face vector av[];
event acceleration (i++) {
  foreach()
    av.z[] = sin(t)*f0[] ;
  boundary ((scalar *){av});
  a = av;
}

/**
##Boundary Conditions
We force the velocity outside (fraction filed) equals zero, which is a simple way to match the no slip boundary conditions*/
event bc (t <= tmax; i += 1) {
  //no slip boundary conditions
  foreach()
    foreach_dimension()
    u.x[] = f0[]*u.x[];
    boundary ((scalar *){u , p});
  /**
We put a test point to verify the resulte is stational*/
  double px = 0.8;
  double py = 1.1;
  double pz = 5.0;
  fprintf(stderr, "%d %g %g %g \n", i, t, dt,interpolate(u.z , px, py, pz))
}

/**
## Dump files 
It's use to restore the data from each 100 interation*/
event snapshot (i += 100)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}

/**
## velocity profile output
We output the velocity profil in the middle of our tube, in order to compare with the theory value*/
int j =0;
event profill (t += 0.2*M_PI)
{
  if (t > tmax - 2.*M_PI){
    j += 1;
    char name3[80];
    sprintf (name3, "test-%d", j);	
    FILE *fp2;
    fp2 = fopen (name3, "w");
    double uyy = 1.1;
    double uzz = 5.;
    for (double uxx = -1.; uxx <= 2.; uxx += 1./100. ){
      fprintf(fp2,"%g %g %g\n", t, (uxx - 0.8), 1.5*interpolate(u.z , uxx, uyy, uzz));
    }
  }
}

/**
## Mesh adaption
mesh adaptation based both on volume fraction and velocity accuracy*/
event adapt (i++) {
  double uemax = 0.01;
  adapt_wavelet ({f0,u}, (double[]){0.01,uemax,uemax,uemax}, MAX_LEVEL, MIN_LEVEL);
}

/**
##fraction field
![Isosurface of the distance function coloured with level of refinement.](endotube/endotubestl.png)
![Reconstructed VOF surface.](endotube/endotubevof.png)

##gnuplot 
~~~gnuplot compare with theory value 
set xlabel "r"
plot [0:1.05][]'../2dwo/thwo10' us 1:2 t'theo' w l ,\
'test-1' us 2:3 t't1' w lp,\
'test-2' us 2:3 t't2' w lp,\
'test-3' us 2:3 t't3' w lp,\
'test-4' us 2:3 t't4' w lp,\
'test-5' us 2:3 t't5' w lp,\
'test-6' us 2:3 t't6' w lp,\
'test-7' us 2:3 t't7' w lp,\
'test-8' us 2:3 t't8' w lp,\
'test-9' us 2:3 t't9' w lp,\
'test-10' us 2:3 t't10' w lp
~~~
*/