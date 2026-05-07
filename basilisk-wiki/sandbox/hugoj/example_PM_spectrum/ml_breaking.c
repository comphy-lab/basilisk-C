/**
# Wave breaking field decay

This example shows how to use [spectrum.h](../lib/spectrum.h).
There is no stratification, and no forcing (wave decay).

The initial sea state is chosen to have a azimuthal integrated spectrum close to
the Pierson-Moscowitz spectrum representative of a developped sea.

The shape of this spectrum shows a peak at $k_p$ that we set as $k_p=10\pi/L$
with $L$ the domain size. The peak Reynolds number is 
$Re_p = \sqrt{g\lambda_p}/\nu = 3.17 \times 10^7$. **SP: this is not a Reynolds number**

![Evolution of the surface elevation](ml_breaking/eta.mp4)(loop width="50%")

![Evolution of the horizontal velocity component of the top
 layer](ml_breaking/u14.mp4)(loop width="50%")
*/

#include "grid/gpu/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "profiling.h"

/** Let us define some functions */
#define g_ 9.81
#include "hugoj/lib/spectrum.h"

/**
## PARAMETERS

**SP: Some of your parameters are redundant/not used??**
*/

// -> Initial conditions
double P = 0.02;         // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10;     // kpL0 = coeff_kpL0 * pi
int N_mode = 32;         // Number of modes in wavenumber space
int N_power = 5;         // directional spreading coeff
// -> Domain definition
int N_grid = 9;          // 2^N_grid : number of x and y gridpoints
double L = 200.0;        // domain size
int N_layer = 15;        // number of layers
double kp;               // peak wave number
double h0 = 100.0;       // depth of water
// -> Runtime parameters
double tend = 200.0;        // end time of simulation
// -> physical properties

/**
 **SP: this is way too small. Did you compute the corresponding wave Reynolds number?**
 */
double nu0 = 0.000025;      // Viscosity for vertical diffusion

int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


// diag
double *u_profile;
double dt_mean = 5;  // s

/* Main program */
int main(int argc, char *argv[])  
{
  
  kp = pi * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength
  
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1;
  TOLERANCE = 1e-4;

  /** Boundary condition */
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  // allocate diag
  u_profile = (double*)calloc(nl, sizeof(double));
   
  fprintf (stderr, "Read in parameters!\n");
  run();
}

/**
## Initial conditions */

event init (i = 0) {
  geometric_beta (1./3., true); // Varying layer thickness

  /** Generate linearly spaced kx, ky according to specified # of modes, and
      then interpolated on a cartesian grid to get F(kx,ky). */
  T_Spectrum spectrum = spectrum_gen_linear (N_mode, N_power, L, P, kp);

  /** set eta and h */
  foreach (cpu) {
    zb[] = - h0;
    eta[] = wave (x, y, spectrum);
    double H = eta[] - zb[];
    foreach_layer()
      h[] = H*beta[point.l];
  }
  /** set currents */
  foreach (cpu) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      coord wu = wave_u (x, y, z, spectrum);
      u.x[] = wu.x;
      u.y[] = wu.y;
      w[]   = wu.z;
      z += h[]/2.;
    }
  }
  fprintf (stderr,"Done initialization!\n");
  free_spectrum (spectrum);
}

/**
## Log */

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  double gpe0 = - sq(L0)*G*sq(h0)/2.;
  fprintf (stderr, "%g %g %g\n", t, ke/2., G*gpe - gpe0);
}

/**
## Layer average

Diagnostic of the layer averaged zonal velocity.

**SP: I don't understand what you are doing here...**. */

event compute_horizontal_avg (t += dt_mean) {
  foreach_layer() {
    double sum = 0.;
    foreach (reduction(+:sum))
      sum += u.x[];
    u_profile[_layer] += sum / (N*N); // * dt / dt_mean;
  }
}

/** Writing the layer average diagnostic to a file  */
event write_diag (t = 0, t += dt_mean)
{
  for (int i = 0; i < nl; i++)
    printf ("%g %d %g\n", t, i, u_profile[i]);
  printf ("\n");
  fflush (stdout);
  
  // Reset the profile for all workers
  for (int i = 0; i < nl; i++)
    u_profile[i] = 0.0;
}

/**
## Movies */

event movies (t += 4./24.) { //  2.*24.
  output_ppm (eta,
#if SHOW  // do not use in "batch" mode: very slow somehow
              fps = 24,
#endif
              file = "eta.mp4", min = -2., max = 2., 
              n = max(N, 512), linear = true, map=gray);

  vector u = lookup_vector ("u14");
  output_ppm (u.x,
#if SHOW  // do not use in "batch" mode: very slow somehow
              fps = 24,
#endif
              file = "u14.mp4", min = -1., max = 4., 
              n = max(N, 512), linear = true, map=gray);
}

// Clean for my diag of layer avg
event cleanup (t = end) {
  free (u_profile);  
}

event ending (t = tend);

/**
## Results

~~~gnuplot Evolution of the kinetic, potential and total energies
set xlabel 'Time (s)'
set ylabel 'Energies'
plot [][0:]'log' u 1:2 w l t 'Kinetic', '' u 1:3 w l t 'Potential', '' u 1:($2+$3)/2. w l t 'Total/2'
~~~

The layer-averaged velocity shows the developement of a sub-surface jet with a maximum velocity of about 6 cm/s

~~~gnuplot Layer-averaged u.x profiles
# Plot the heatmap
set pm3d map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "u.x (m/s)"
set terminal pngcairo size 800,600 enhanced
set output 'u_profile.png'
unset key
splot "out" using 1:2:3 with pm3d
~~~

## Profiling

This simulation runs in 1h44 min (using a RTX PRO 1000 Blackwell Laptop).

~~~bash
calls    total     self   % total   function
    1296     4.47     4.44     81.2%   relax_nh():/home/jacqhugo/basilisk/src/layered/nh.h:223
       1   6252.83     0.39      7.2%   run():/home/jacqhugo/basilisk/src/run.h:37
      11     0.14     0.14      2.6%   vertical_remapping():/home/jacqhugo/basilisk/src/layered/remap.h:165
      22     0.12     0.11      2.0%   advect():/home/jacqhugo/basilisk/src/layered/hydro.h:404
      47     0.10     0.10      1.9%   residual_nh():/home/jacqhugo/basilisk/src/layered/nh.h:299
    7140     0.04     0.04      0.8%   setup_shader():/home/jacqhugo/basilisk/src/grid/gpu/grid.h:1956
      69     0.04     0.04      0.7%   multigrid_restriction():/home/jacqhugo/basilisk/src/grid/gpu/../multigrid-common.h:459
      11     4.67     0.03      0.6%   pressure_1():/home/jacqhugo/basilisk/src/layered/nh.h:460
     396     0.03     0.02      0.5%   relax_hydro():/home/jacqhugo/basilisk/src/layered/implicit.h:104
      47     4.54     0.02      0.4%   mg_cycle():/home/jacqhugo/basilisk/src/poisson.h:92
      11     0.03     0.02      0.4%   acceleration_0():/home/jacqhugo/basilisk/src/layered/implicit.h:206
      11     0.07     0.02      0.3%   pressure():/home/jacqhugo/basilisk/src/layered/hydro.h:461
~~~

This simulation runs in ~10 minutes on an RTX 4090 (for a resolution
of 512^2^ i.e. `N_grid = 9`).

Profiling output (`CFLAGS=-DTRACE=2`) gives:

~~~bash
   calls    total     self   % total   function
    1840     2.63     2.58     75.4%   relax_nh():/src/layered/nh.h:223
      65     0.15     0.15      4.4%   residual_nh():/src/layered/nh.h:300
      38     0.16     0.15      4.4%   advect():/src/layered/hydro.h:404
      19     0.11     0.11      3.2%   vertical_remapping():/src/layered/remap.h:165
      19     2.93     0.07      2.1%   pressure_1():/src/layered/nh.h:463
   11485     0.06     0.06      1.8%   setup_shader():/src/grid/gpu/grid.h:1956
     103     0.05     0.05      1.4%   multigrid_restriction():/src/grid/gpu/../multigrid-common.h:459
      19     0.05     0.04      1.3%   acceleration_0():/src/layered/implicit.h:206
      65     2.73     0.04      1.1%   mg_cycle():/src/poisson.h:92
      19     0.12     0.03      1.0%   pressure():/src/layered/hydro.h:461
~~~

## TODO

* add bview snapshot when using CPU
* crashes at after some time as passed ? (for `N_grid = 10`, this can be fixed with the new version of the RHS of the multilayer).
* Initialisation is very slow: use FFTs!
*/
