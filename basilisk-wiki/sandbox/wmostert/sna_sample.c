/**
   This is a standard file to reproduce numerical simulation of bubble breakup in the publication:

   [Daniel J. Ruth, Wouter Mostert, Stéphane Perrard, and Luc Deike. Bubble pinch-off in turbulence. Proceedings of the National Academy of Sciences, 116(51):25412–25417, 2019.](https://www.pnas.org/content/116/51/25412)

   Source Basilisk code for numerical simulation of bubble
   breakup. This problem file is for the L=12 nonaxisymmetric bubbles
   which resolve on curvature - see Supplementary Information for
   description.

   Bubble pinchoff studied through a dumbbell-shaped STL file.  To run
   successfully, both a restart file and the bubble-STL should be
   present in the directory. The STL file can be generated by
   e.g. MATLAB. To generate the restart file, run the code on a single
   process until the very first dump file is produced (titled
   "dump0"). Rename that dump file to "restart". Then, the code can be
   re-run on any number of processes. */

#include "grid/octree.h"
#include "distance.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

/**
   Set the maximum and minimum refinement levels. */
#define LEVEL 12
#define MINLEVEL 5
double maxl;

/**
   Set various constants for the simulation. */
#define RATIO 1.0/850.0
#define MURATIO  17.4e-6/8.9e-5

/**
   Define the distance field for the initial condition. */
scalar d[];
vertex scalar psi[];
face vector s[];

/**
   Define curvature refinement terms. */
scalar curv[];
double DoC = 1./16.;
void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}

int main ()
{
  fprintf(ferr, "Beginning problem.\n");
  init_grid(1 << MINLEVEL);
  /** Periodic boundary conditions. */
  periodic (top);
  periodic (right);
  periodic (front);
  /**
     Set physical constants. */
  fprintf(ferr, "Setting physical constants...\n");
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = rho1 * 0.001;
  mu2 = mu1*MURATIO;
  f.sigma = 1.0;
  fprintf(ferr, "Getting coordinates...\n");
  coord * p =input_stl ( fopen ("dumbbell_nonaxi.stl", "r"));
  coord min, max;
  fprintf(ferr, "Setting bounding box...\n");
  bounding_box (p, &min, &max);
  maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  fprintf(ferr, "maxl=%g\n", maxl);
  size (1.2*maxl);
  fprintf(ferr, "parameters=%g %g %g %g %g %g\n", max.x, min.x, max.y, min.y, max.z, min.y);
  fprintf(ferr, "Corner, x %g\n", (max.x + min.x)/2.-L0/2.);
  fprintf(ferr, "Corner, y%g\n", (max.y + min.y)/2.);
  origin ((max.x + min.x)/2. - L0/2.,
	  (max.y + min.y)/2. - L0/2.,
	  (max.z + min.z)/2. - L0/2.);
  run();
  return 0;
}

//---------------------------INITIAL CONDITIONS---------------------------------------
event init (i=0)
{
  if (!restore("restart")){
    /**
       Input coordinate file for the dumbbell. */
    char filename[100];

    /**
       Initialize the distance field. */
    fprintf(ferr, "Initializing distance field...\n");
    coord * p =input_stl ( fopen ("dumbbell_nonaxi.stl", "r"));
    distance (d, p);
    while (adapt_wavelet ({d}, (double[]){1e-6}, LEVEL).nf);
    fprintf(ferr, "Initializing vertex distance field...\n");
    /**
       Construct a vertex field and interpolate from the centered field. */
    foreach_vertex()
      psi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.0;
    fprintf(ferr, "Initializing fractions field...\n");
    fractions (psi, f, s);
    while (adapt_wavelet ({f}, (double[]){1e-6}, LEVEL).nf);
    fprintf(ferr, "Complete.\n");
  }
}

//----------------------------ADAPTIVITY--------------------------------------------
//

event adapt (i++) {
  curvature(f, curv);
  foreach()
   if (abs(z) > 0.05)
    curv[] = nodata;
  curv.prolongation = prolongate_ratio;
  double femax = 3e-1;
  adapt_wavelet({f, curv}, (double[]){femax, DoC}, LEVEL );
}

event dumpfile (t += 0.01)
{
  char dname[100];
  sprintf(dname, "dump%g", t);
  fprintf(ferr, "t = %g\n", t);
  char fname[100];
  sprintf(fname, "f.ppm");
  static FILE * fp = fopen ("f.ppm", "w");
  output_ppm (f, fp, n=1024);
  static FILE * fu = fopen ("ux.ppm", "w");
  output_ppm (u.x, fu, n=1024);
  static FILE * fv = fopen ("uy.ppm", "w");
  output_ppm (u.y, fv, n=1024);
}

event dumpnext (i += 10)
{
  char dname[100];
  sprintf(dname, "dump%d", i);
  dump(dname);
}

event dumpfirst(i=0)
{
  dump ("dump0");
}

event end (i=30000)
{
  fprintf(ferr, "t=%g", t);
}
