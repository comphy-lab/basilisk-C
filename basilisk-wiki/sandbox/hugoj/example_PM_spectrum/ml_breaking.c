/**
  # Wave breaking field decay

No stratification, wave decay (no forcing)
*/
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/perfs.h"
#include "layered/remap.h"


/** Let us define some functions */
#define g_ 9.81
#include "hugoj/lib/spectrum.h"

/*
 ## PARAMETERS

Dimensions : [Length, Time]
*/
// -> Initial conditions
double P = 0.2 [1, -1];     // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];     // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];         // Number of modes in wavenumber space
int N_power = 5 [];         // directional spreading coeff
int F_shape = 0 [];         // shape of the initial spectrum
// -> Domain definition
int N_grid = 8 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 15 [];         // number of layers
double kp = PI*10/200.0 [-1];// peak wave number
double h0 = 100.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 200.0;          // end time of simulation
// -> saving outputs
double dtout = 2.0;         // dt for output in netcdf
// -> physical properties
double nu0 = 0.000025;      // Viscosity for vertical diffusion
int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


// diag
double *u_profile;
double dt_mean = 5;  // s
static FILE * fp;

/* Main program */
int main(int argc, char *argv[])  
{
  
  kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength
  
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; 

  double omegap=sqrt(g_*kp);
  double cphase=omegap/kp;

  /** Boundary condition */
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  // allocate diag
  u_profile = (double*)calloc(nl, sizeof(double));
  fp  = fopen("u_profile.dat","w"); // reset file
  fclose(fp);
   
  fprintf (stderr, "Read in parameters!\n");
  run();
}

/* Initial conditions */
event init(i =  0) {
  geometric_beta (1./3., true); // Varying layer thickness

  /** Generate linealy spaced kx, ky according to specified # of modes, and
      then interpolated on a cartesian grid to get F(kx,ky). */
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);

  /** set eta and h */
  foreach(cpu) {
    zb[] = -h0;
    eta[] = wave(x, y, N_grid, spectrum, 0.);
    double H = wave(x, y, N_grid, spectrum, 0.) - zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
    } 
  }
  /** set currents */
  foreach(cpu) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      u.x[] = u_x(x, y, z, N_grid, spectrum, 0.);
      u.y[] = u_y(x, y, z, N_grid, spectrum, 0.);
      w[] = u_z(x, y, z, N_grid, spectrum, 0.);
      z += h[]/2.;
    }
  }
  fprintf (stderr,"Done initialization!\n");
}


// event film surface
event movies (t += 1./24.){ //  2.*24.
  output_ppm (eta, fps = 24, file = "eta.mp4", min = -2., max = 2., 
              n = 512, linear = true, map=gray);

  vector u = lookup_vector ("u14");
  output_ppm (u.x, fps = 24, file = "u14.mp4", min = -1., max = 5., 
              n = 512, linear = true, map=gray);

}


/** Diagnostic of the layer averaged zonal velocity */
This event compute layer average of u.x
event compute_horizontal_avg (i++; t<=tend+1e-10){
event compute_horizontal_avg (t+=dt_mean; t<=tend+1e-10){
   foreach(reduction(+:u_profile[:nl]))
     foreach_layer(){
       u_profile[point.l] += u.x[] / (N*N) * dt / dt_mean;

     }
 }

// This event writes to a file the layer average
event write_diag(t=0., t+=dt_mean){
     // main worker is writing the file
     if (pid()==0) {
       fp  = fopen("u_profile.dat","a");
       if (fp == NULL){
         fprintf(stderr, "Error opening file u_profile.dat");
         return 2;
       }
       for (int i=0; i<nl; ++i) {
         fprintf (fp, "%f %d %g\n", t, i, u_profile[i]);
       }
       fprintf(fp,"\n");
       fclose(fp);
     }
     // Reset the profile for all workers
     for (int i=0; i<nl; ++i) {
       u_profile[i] = 0.0;
     }
  }

// event snapshot (t = end)
// {
//   clear();
	
  // 3/4 view
  // view (quat = {0.567, 0.137, 0.196, 0.789},
  //     fov = 30, near = 0.01, far = 1000,
  //     tx = 0.048, ty = -0.001, tz = -2.096,
  //     width = 1398, height = 803);
  //
  //
  // char s[80];
  // sprintf (s, "t = %.2f", t);
  // draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x) {
  //     squares ("eta", linear = true, z = "eta*10", min = -1.0, max = 1.0 , map=gray);
  //   }
  // box ();
  //
  //
  //
  //
  //
  // colorbar(map=gray, label="eta (m)", min=-1.0,max=1.0, pos={-0.95,-0.5},
  //          levels=10);
  // save ("snap_side.png");
  //
  //
  // // top view
  // clear();
  // view (quat = {0.000, 0.000, 0.000, 1.000},
  //     fov = 30, near = 0.01, far = 1000,
  //     tx = 0.000, ty = 0.000, tz = -2.505,
  //     width = 1398, height = 803);
  //
  //
  // sprintf (s, "t = %.2f", t);
  // draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x) {
  //     squares ("eta", min = -1.0, max = 1.0 , map=gray);
  //   }
  // box ();
  //
  // colorbar(map=gray, label="eta (m)", min=-1.0,max=1.0, pos={-0.95,-0.5},
  //          levels=10);
  // save ("snap_top.png");


//}

// Clean for my diag of layer avg
event cleanup(t=tend){
  free(u_profile);
}

event ending (t = tend);




/**
Results: layer average u.x profile
~~~gnuplot 
# Plot the heatmap
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "u.x (m/s)"
set yrange [0:14]
set xrange [0:200]
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'u_profile.png'
set size 0.9, 0.9
splot "u_profile.dat" using 1:2:3 with pm3d
unset output
~~~
**/

