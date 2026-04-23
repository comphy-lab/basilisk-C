/**
 
# Investigating horizontal average in the multilayer framework 

The multilayer solver uses layer. Here we investigate how to do the horizontal
average for different cases, from a simple monochromatic linear wave to a
realistic synthetic wave field.

 */
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files

/**
 ## Parameters

General */

double g_=9.81;

#ifndef WAVE
    #define WAVE 0 // default is to produce monochromatic wave
#endif

double L = 200.0;         // Size of the domain
double h0 = 100.0;        // depth of water
double kp;                // peak wavelength

/** Solver parameters */
int N_grid = 6;           // Number of cells
int N_layer = 30;         // Number of layers
char ncname[20] = "out";  // file name of output

/** Stokes wave parameters */
double ak = 0.033;

/** Stynthetic field parameters */
double P = 0.02;          // Energy level
int N_mode = 32;          // Number of modes
int N_power = 5;          // power exponent of cos^n
double thetam=0.;         // default direction 
int Nthetam = 4;          // number of direction (dir=i/Nthetam for i in
                          //                                [0,Nthetam])

/** Initial conditions generation functions */
#if WAVE<2
#include "hugoj/lib/common_waves.h"
#elif WAVE==2
#include "hugoj/lib/spectrum.h" 
#endif


int main(){
  
  /** Common part, solver values */
  L0 = L;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = 0.503;
  CFL_H = 1; 
  origin (-L0/2., -L0/2.);
  periodic (top);
  periodic (left);

  run(); 
}



event init(i =  0) {
  
  /**
  ## Monochromatic wave
  */
  #if WAVE == 0
  
  snprintf(ncname, sizeof(ncname), "monoc.nc");
  double a = 1.;
  double k = 2*pi/L;
  geometric_beta (1./3., true);
  foreach() {
    zb[] = - h0;
    double H = wave_monolin(0, x, a, k) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x_monolin(0, x, z, a, k);
      w[] = u_y_monolin(0, x, z, a, k);
      eta[] += h[];
      z += h[]/2.;
    }
  }
  create_nc({zb, h, u, w, eta}, ncname);
  write_nc();

  /**
  ## Stokes wave
  */
  #elif WAVE == 1
  
  snprintf(ncname, sizeof(ncname), "stokes.nc");
  kp=2*pi/L;
  geometric_beta (1./3., true);
  foreach() {
     zb[] = - h0;
     double H = wave_stokes(x, 0, ak, kp, h0) - zb[];
     double z = zb[];
     foreach_layer() {
       h[] = H*beta[point.l];
       z += h[]/2.;
       u.x[] = u_x_stokes(x, z, ak, kp, h0);
       w[] = u_y_stokes(x, z, ak, kp, h0);
       z += h[]/2.;
       eta[] += h[];
     }
  }
  create_nc({zb, h, u, w, eta}, ncname);
  write_nc();

  /**
  ## Synthetic wave field (4 directions)
  */
  #elif WAVE == 2
  

  kp=10*PI/L;
  // initialise thetam list 
  double * list_thetam = (double*)calloc(Nthetam, sizeof(double));
  for (int i=0; i<Nthetam; ++i){
    list_thetam[i] = i*PI;
  }

  // Varying layer thickness
  geometric_beta (1./3., true); 
  
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
