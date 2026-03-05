/**
  Field scale wave breaking (multilayer solver)

No stratification


USAGE

make

*/


#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "bderembl/libs/netcdf_bas.h"
#define g_ 9.81
#include "spectrum.h" // Initial conditions generation

/*
DEFAULT PARAMETERS
*/
char file_out[20] = "out.nc";         // file name of output
// -> Initial conditions
double P = 0.02 [1, -1];     // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];     // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];         // Number of modes in wavenumber space
int N_power = 5 [];         // directional spreading coeff
int F_shape = 1 [];         // shape of the initial spectrum
// -> Domain definition
int N_grid = 6 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 10 [];         // number of layers
double kp = PI*10/200.0 [-1];// peak wave number
double h0 = 40.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 10.0;          // end time of simulation
// -> saving outputs
double dtout = 1.0;        // dt for output in netcdf
// -> physical properties
double nu0 = 0.;            // Viscosity for vertical diffusion
int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


// diag
double *u_profile;
double dt_mean = 1.0;
double U=0.0;
int l;
static FILE * fp;




int main(int argc, char *argv[])  
{
 
  kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength

  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; // Smaller time step

  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  // allocate diag
  u_profile = (double*)calloc(nl, sizeof(double));
  l=nl-1;
  fprintf(stderr, "u%d ini = %f\n",l, U);
  fp  = fopen("u_profile.dat","w"); // reset file
  fp  = fopen("u_sfx.dat","w"); // reset file
  fclose(fp);

  run();
}

event init(i =  0) {

  geometric_beta (1./3., true); // Varying layer thickness
  T_Spectrum spectrum;
  if (F_shape==1){
    spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);
  }
  else{
    spectrum = read_spectrum("init_spectrum", N_mode); // TO DO: 
  }

  // step 1: set eta and h
  foreach() {
    zb[] = -h0;
    eta[] = wave(x, y, N_grid, spectrum);
    double H = wave(x, y, N_grid, spectrum) - zb[];
    foreach_layer() {
      h[] = H/nl;
    } 
  }
  
  // step 2: remap
  vertical_remapping (h, tracers);
  // step 3: set currents
  foreach() {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      u.x[] = u_x(x, y, z, N_grid, spectrum);
      u.y[] = u_y(x, y, z, N_grid, spectrum);
      w[] = u_z(x, y, z, N_grid, spectrum);
      z += h[]/2.;
    }
  }
  
  // step 4: save initial conditions
  create_nc({zb, h, u, w, eta}, file_out);
  
  // Free memory
  free(spectrum.kmod);
  free(spectrum.kx);
  free(spectrum.ky);  
  free(spectrum.F_kmod);
  free(spectrum.F_kxky);
  free(spectrum.phase);
  free(spectrum.omega);
  fprintf (stderr,"Done initialization!\n");
}


// COMPUTING U_PROFILE
event compute_horizontal_avg (i++; t<=tend+1e-10){
  // fprintf(stderr, "dt %f\n",dt);
  
  //#if 0
  foreach(reduction(+:u_profile[:nl])){
    foreach_layer(){
      // if (point.l==29) {
      //   if (u.x[0,0,_layer])
      //   fprintf(stderr, "t %f, i %d, value %f\n", t, i,  u.x[0,0,_layer] / (N*N) * dt / 1.0);
      // }
      u_profile[point.l] += u.x[] / (N*N) * dt / dt_mean;

    }
  }
  //#else
  fprintf(stderr, "   u_sfx pre avg = %f\n", U);
    foreach(reduction(+:U)){
      U += u.x[0,0,l] / (N*N) * dt / dt_mean;
    }
  fprintf(stderr, "   t=%f u_sfx= %f\n",t, U);
  //#endif
  //fprintf(stderr, "t %f, i %d, value %f\n", t, i,  u_profile[nl-2]);
  fprintf (stderr, "t %f u.max %f\n", t, statsf(u.x).max);
}

// WRITING U_PROFILE
event write_diag(t=0., t+=dt_mean, t<=tend+1e-10){
    //#if 0
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
    for (int i=0; i<nl; ++i) {
      u_profile[i] = 0.0;
    }
    //#else
    if (pid()==0) {
      fp  = fopen("u_sfx.dat","a");
      if (fp == NULL){
        fprintf(stderr, "Error opening file u_sfx.dat");
        return 2;
      }
      fprintf (fp, "%f %d %g\n", t, l, U);
     
      //fprintf(fp,"\n");
      fclose(fp);
    }
    fprintf(stderr, "writing   t=%f u_sfx= %f\n",t, U);
    U=0.0;
    //#endif
}


// OUTPUT NETCDF
event output(t = 0.; t<= tend+1e-10; t+=dtout){
  //fprintf(stdout, "output at t=%f, i=%d\n", t, i);
  write_nc();
}

event cleanup(t=end){
  free(u_profile);
}

/**
## Results

~~~gnuplot The profile of layer averaged u.x in time (serial)
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "Value"
splot "u_profile.dat" using 1:2:3 with pm3d
~~~

~~~gnuplot The profile of layer averaged u.x in time (MPI)
splot "../ml_mpi/u_profile.dat" using 1:2:3 with pm3d
~~~

~~~gnuplot The first layer average of u.x in time, MPI is clearly now correct with patch a983e576ec608cfaef43c343d30140e1a3e20c0c, thanks Stéphane.
reset
set xlabel "Time (s)"
set ylabel "U_{sfx} (m/s)"
set logscale y
plot "u_sfx.dat" using 1:3 w l t 'serial', \
     "../ml_openmp/u_sfx.dat" u 1:3 t 'openmp', \
     "../ml_mpi/u_sfx.dat" u 1:3 t 'MPI'
~~~

~~~gnuplot u.x.max
set ylabel "u.x.max"
unset logscale
plot "< grep u.max log" u 2:4 w lp t 'serial', \
     "< grep u.max ../ml_openmp/log" u 2:4 w lp t 'openmp', \
     "< grep u.max ../ml_mpi/log" u 2:4 w lp t 'mpi'
~~~
**/


