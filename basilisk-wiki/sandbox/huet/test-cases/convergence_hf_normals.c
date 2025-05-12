/**
# Question about the convergence of the interface normals computed by the height functions

In this file we create a circular (respectively spherical) droplet for three different mesh sizes (levels 5 to 7) in order to compare the convergence of interface normals computed by the height functions in 2D and 3D.

We assume that the height function normals are provided at the center of the interfacial cells, while the are actually computed on the interface: it means this test should show a first order convergence at best. This is what is found in 2D, however in 3D the interfacial normals do not converge.
Since this result is surprising, normal vectors from a levelset function are also computed, just to make sure no careless mistake as been made when computing the error. The normals from the level-set show the expected order of convergence in both 2 and 3 dimensions

To reproduce these results, run this file twice, commenting the first and second line alternatively. Then, look at the files "hf_2d.txt", "hf_3d.txt", "levelset_2d.txt", "levelset_3d.txt".  In each file, the three columns are the level, the L2 norm and the L-infinity norm.

I used the previous version of qcc (pre-December 2021, so without automatic boundary conditions), although I wouldn't be surprised if this file compiles with the new version as well.
*/

#include "grid/quadtree.h"
// #include "grid/octree.h"
#include "run.h"
#include "fractions.h"
#include "curvature.h"
scalar f[], * interfaces = {f};

#define L0 1.
#define RADIUS (.25 + 1.e-8)
#define LEVEL_MIN 5
#define LEVEL_MAX 7

/** Define the spherical coordinates for the analytical solution */
#define RAD (sqrt(sq(x) + sq(y) + sq(z)) > 1.e-30 ? \
  sqrt(sq(x) + sq(y) + sq(z)) : 0.)
#if dimension > 2
  #define THETA (fabs(z) > 1.e-30 ? atan2(sqrt(sq(x) + sq(y)),z) : pi/2.)
#else
 #define THETA (pi/2.)
#endif
#define PHI (fabs(x) > 1.e-30 ? atan2(y,x) : ((y > 0 ? (pi/2.) : (3*pi/2.))))

/** Compute normal vectors with either the height-functions or the level-set*/
scalar levelset[];
vector normal[];
void impose_n(bool use_height_function) {
  if (use_height_function) {
    // Compute normals in interfacial cells using the height functions
    coord n_coord;
    vector fh = f.height, h = automatic (fh);
    if (!fh.x.i)
    heights (f, h);
    foreach() {
      if (interfacial(point, f)) {
        n_coord = height_normal(point, f, h);
        foreach_dimension() {
          normal.x[] = n_coord.x;
        }
      }
      else foreach_dimension() normal.x[] = 0.;
    }
  }
  else {
    // Compute normals in interfacial cells using a levelset function
    foreach() {
      if (interfacial(point, f)) {
        foreach_dimension() {
          normal.x[] = -(levelset[1] - levelset[-1])/(2.*Delta);
        }
        double normn = sqrt(sq(normal.x[]) + sq(normal.y[]) + sq(normal.z[]));
        if (normn > 1.e-30) {
          foreach_dimension() {
            normal.x[] /= normn;
          }
        }
      }
    }
  }
  boundary((scalar*){normal});
}

int my_level = 0.;
FILE* foutput = NULL;
bool use_height_function;
int main() {
  size(L0);
  origin(-.5*L0, -.5*L0, -.5*L0);
  DT = 1.;

  use_height_function = 1;
  #if dimension > 2
    foutput = fopen("hf_3d.txt","w");
  #else
    foutput = fopen("hf_2d.txt","w");
  #endif
  for (int level=LEVEL_MIN; level<=LEVEL_MAX; level++) {
    my_level = level;
    N = 1 << level;
    run();
  }
  fclose(foutput);

  use_height_function = 0;
  #if dimension > 2
    foutput = fopen("levelset_3d.txt","w");
  #else
    foutput = fopen("levelset_2d.txt","w");
  #endif
  for (int level=LEVEL_MIN; level<=LEVEL_MAX; level++) {
    my_level = level;
    N = 1 << level;
    run();
  }
  fclose(foutput);
}

event init (i = 0) {
  fraction(f, sq(RADIUS) - sq(x) - sq(y) - sq(z));
}

event do_the_work (i++) {
  fraction(f, sq(RADIUS) - sq(x) - sq(y) - sq(z));
  foreach() {
    levelset[] = RADIUS - RAD;
  }
  boundary({levelset});
  impose_n(use_height_function);
}

/** Compute and output the error associated with our normal vectors, assuming they are computed at the center of the cells (which means we expect a 1st order accuracy in the case of height functions).*/
vector n_exact[];
event compute_error (t = 2.) {
  vector n_error[];
  double l2_err = 0.;
  double linf_err = 0.;
  int nb_cells = 0.;
  foreach() {
    if (interfacial(point,f)) {
      nb_cells++;
      n_exact.x[] = cos(PHI)*sin(THETA);
      n_exact.y[] = sin(PHI)*sin(THETA);
      n_exact.z[] = cos(THETA);
      foreach_dimension() n_error.x[] = sqrt(sq(n_exact.x[] - normal.x[]));
      l2_err += (n_error.x[] + n_error.y[] + n_error.z[])/3.;
      foreach_dimension() if (n_error.x[] > linf_err) linf_err = n_error.x[];
    }
    else foreach_dimension() n_error.x[] = 0.;
  }
  if (nb_cells > 0) {
    l2_err /= nb_cells;  }
  fprintf(foutput, "%d %g %g", my_level, l2_err, linf_err);
  fprintf(foutput, "\n");
  fflush(fout);
}

event end (t = 2.) {
  return (1.);
}
