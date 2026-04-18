/**
 
 # Investigating the initial velocity profiles under a directional spectrum

**/
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files



#define g_ 9.81

#ifndef WAVE
    #define WAVE 0 // default is to produce monochromatic wave
#endif

// General parameters

double P = 0.02;
int N_mode = 32;
int N_power = 5;
double L = 200.0;
double kp;
double h0 = 40.0;        // depth of water

double * list_thetam;
double thetam=0.;
int Nthetam = 4;

int N_grid = 6;
int N_layer = 15;
char ncname[20] = "out";         // file name of output

int main(){
  
  // Common part, solver values
  L0 = L;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = 0.503;
  CFL_H = 1; 
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  

    
  
  
  // // loop through the thetam values
  // for (int kt=0; kt<Nthetam; ++kt){
  //   thetam = list_thetam[kt];
  //   fprintf(stderr, "run loop, thetam=%f\n", thetam);
  //   snprintf(ncname, sizeof(ncname), "out_%d.nc", kt);
  run();
  
}

event init(i =  0) {
  
  /*
  ## Monochromatic wave
  */
  #if WAVE == 0
  
  snprintf(ncname, sizeof(ncname), "monoc.nc");
  
  create_nc({zb, h, u, w, eta}, ncname);
  write_nc();

  /*
  ## Stokes wave
  */
  #elif WAVE == 1
  #include "examples/test/stokes.h"
  snprintf(ncname, sizeof(ncname), "stokes.nc");

  double a = 0.25;
  foreach() {
    zb[] = - 0.5 + a*sin(k_/2.*y);
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
    }
  }
  create_nc({zb, h, u, w, eta}, ncname);
  write_nc();

  /*
  ## Synthetic wave field (4 directions)
  */
  #elif WAVE == 2
  #include "hugoj/lib/spectrum.h" // Initial conditions generation

  kp=10*PI/L;
  // initialise thetam list 
  list_thetam = (double*)calloc(Nthetam, sizeof(double));
  for (int i=0; i<Nthetam; ++i){
    list_thetam[i] = i*PI;
  }

  // Varying layer thickness
  geometric_beta (0., true); 
  
  // loop through the thetam values
  for (int kt=0; kt<Nthetam; ++kt){
    thetam = list_thetam[kt];
    snprintf(ncname, sizeof(ncname), "synth_eta_%d.nc", kt);

    // step 0: generate a synthetic spectrum
    T_Spectrum spectrum;
    spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp, thetam);

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
    // initialise nc file
    create_nc({zb, h, u, w, eta}, ncname);
    write_nc();
  }
  
  #else
    fprintf(stderr, "the chosen wave flag should be 0 or 1 or 2, it is WAVE");
    return 2;
  #endif
  return 1;
}

/**
TODO: 
- thetam is correctly input in spectrum_gen_linear but the eta and velocity
fields are not correctly initialised. I think this comes to the fact the kx, ky
are not changed accordingly.

**/
