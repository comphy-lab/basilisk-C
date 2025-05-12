/**
# Field scale wave breaking (multilayer solver)
*/

#include "grid/multigrid.h"
#include "view.h"
#include "layered/hydro.h"
//#include "layered/hydro_test.h"
#include "layered/nh.h"
#include "layered/remap.h"
//#include "remap_test.h"
#include "layered/perfs.h"
// #include "input.h" // Used for older input method (matrices of quantities)
#include "spectrum.h" // Used for new input method (the spectrum info)
#include "output_mpi.h" // Antoon's function for MPI compatible matrix output

/**
Definition of some controlling parameters. */

#define g_ 9.8
double h_ = 10; // depth of the water
double gpe_base = 0; // gauge of potential energy
/* double ETAE = 0.1; // refinement criteria for eta */
/* int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function */
/* int MINLEVEL = 6; // min level */
double TEND = 50.; // t end
int NLAYER = 10; // number of layers
int LEVEL_data = 7;

/** 
## Parameters
*/

int main(int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]); // # of layers
  if (argc > 2)
    LEVEL_data = atoi(argv[2]); // Horizontal resolution
  if (argc > 3)
    TEND = atof(argv[3]); // End time
  if (argc > 4)
    nu = atof(argv[4]); // Viscosity for vertical diffusion
  else
    nu = 0.;
  if (argc > 5) 
    RANDOM = atoi(argv[5]); // An integer to seed the random number generator (don't use 0 or 1)
  if (argc > 6)
    L0 = atof(argv[6]); // Box size (not necessarily related to peak wave number)
  else
    L0 = 50.;
  if (argc > 7)
    theta_H = atof(argv[7]); // Numerical parameter to dump fast barotropic modes
  else
    theta_H = 0.5;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; 
  nl = NLAYER;
  G = g_;
  h_ = L0/5.; // set the water depth to be 1/5 the field size (ad hoc)
#if dimension == 2
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
#else
  gpe_base = -0.5*sq(h_)*L0*g_;
#endif
  CFL_H = 1; // Smaller time step
  // max_slope = 0.8; // Change the slope limiter threshold
  run();
}

/**
## Initialization
New method: read in the power spectrum and initialize the wave field. 
In addition we use the already written remapping function that defines a sigma vertical coordinate. See [breaking.c](http://basilisk.fr/src/examples/breaking.c) for details. 

TO-DO: test different remapping method and $\theta_H$ (default 0.51).

~~~~~~~
coeff = 0.05;
geometric_beta (1./3., true);
theta_H = 0.51;
~~~~~~~

*/

event init (i = 0)
{
  if (!restore ("restart")) {
    power_input();
    dkx_ = kx_[1] - kx_[0];
    dky_ = ky_[1] - ky_[0];
    fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_);
    geometric_beta (1./3., true); // Varying layer thickness
    foreach() {
      zb[] = -h_;
      eta[] = wave(x, y);
      double H = wave(x, y) - zb[];
      foreach_layer() {
	      h[] = H/nl;
      }
    }
    // remap?
    vertical_remapping (h, tracers);
    foreach() {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        u.x[] = u_x(x, y, z);
        u.y[] = u_y(x, y, z);
        w[] = u_z(x, y, z);
        z += h[]/2.;
      }
    }
    fprintf (stderr,"Done initialization!\n");
    dump("initial");
  }
}

/**
## Outputs
This is not necessary. It seems that the remapping does not change the energy.
*/

event energy_before_remap (i++, last)
{
  if (i==10) {
    fprintf(stderr, "energy output before remap!\n");
    fflush(stderr);
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer () {
      double norm2 = sq(w[]);
      foreach_dimension()
        norm2 += sq(u.x[]);
        ke += norm2*h[]*dv();
        gpe += (zc + h[]/2.)*h[]*dv();
        zc += h[];
    }
  }
  static FILE * fp = fopen("energy_before_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe - gpe_base);
  fflush (fp);
}

event energy_after_remap (i++, last)
{
  if (i==10) {
    fprintf(stderr, "energy output after remap!\n");
    fflush(stderr);
  }
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer () {
      double norm2 = sq(w[]);
      foreach_dimension()
        norm2 += sq(u.x[]);
        ke += norm2*h[]*dv();
        gpe += (zc + h[]/2.)*h[]*dv();
        zc += h[];
    }
  }
  static FILE * fp = fopen("energy_after_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe - gpe_base);
  fflush (fp);
}

/**
   Note that the movie generation below is very expensive. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)
#if 1
event movie (t += 0.1; t <= TEND)
{
  char s[80];
  view (fov = 20, quat = {0.475152,0.161235,0.235565,0.832313}, width = 800, height = 600);
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 30);
  sprintf (s, "u%d.x", nl-1);
  squares (s, linear = true, z = "eta", min = -2./7.*sqrt(L0), max = 2./7.*sqrt(L0));
  {
  static FILE * fp = POPEN ("ux", "a");
  save (fp = fp);
  }
  scalar slope[];
  foreach () {
    slope[] = (eta[1]-eta[-1])/(2.*Delta);
  }
  clear();
  squares ("slope", linear = true, z = "eta", min = -1./50.*L0, max = 1./50.*L0);
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 30);
  {
  static FILE * fp = POPEN ("slope", "a");
  save (fp = fp);
  }
  char filename1[50], filename2[50], filename3[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  sprintf (filename3, "surface/uy_matrix_%g", t);  
  FILE * feta = fopen (filename1, "w");
  // Might need to change to mpi function later
  output_matrix_mpi (eta, feta, N, linear = true);
  fclose (feta);
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi (u_temp.x, fux, N, linear = true);
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi (u_temp.y, fuy, N, linear = true);
  fclose (fuy);  
}
#endif

/**
   Output 3-D field (not just the surface laye) if needed for Paraview visualization or 
   other analyses. */
#if PARAVIEW
event paraview (t = 100; t += 0.2; t <= TEND) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50];
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_matrix_t%g_l%d", t, j);
    sprintf (filename2, "field/uy_matrix_t%g_l%d", t, j);  
    sprintf (filename3, "field/uz_matrix_t%g_l%d", t, j);  
    sprintf (filename4, "field/h_matrix_t%g_l%d", t, j);  
    sprintf (s, "u%d", j);
    vector u_temp = lookup_vector (s);
    FILE * fux = fopen (filename1, "w");
    output_matrix_mpi (u_temp.x, fux, N, linear = true);
    fclose (fux);
    FILE * fuy = fopen (filename2, "w");
    output_matrix_mpi (u_temp.y, fuy, N, linear = true);
    fclose (fuy);
    FILE * fuz = fopen (filename3, "w");
    sprintf (s, "w%d", j);
    scalar w_temp = lookup_field (s);
    output_matrix_mpi (w_temp, fuz, N, linear = true);
    fclose (fuz);
    FILE * fh = fopen (filename4, "w");
    sprintf (s, "h%d", j);
    scalar h_temp = lookup_field (s);
    output_matrix_mpi (h_temp, fh, N, linear = true);
    fclose (fh);    
  }
}
#endif

/** 
The mesh is not adaptive yet. */

#if QUADTREE
event adapt (i++) {
  /* fprintf(stderr, "Adapting start!\n"); */
  /* fflush(stderr); */
  my_adapt();
}
#endif

event endrun (t = TEND) {
  dump();
}