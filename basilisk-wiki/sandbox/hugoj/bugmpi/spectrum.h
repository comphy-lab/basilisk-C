/*
 Code translation from python to C.
 Original code : Jiarong Wu,
 https://github.com/jiarong-wu/multilayer_breaking/blob/main/specgen/specgen.py

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "interpolate.h"

#define PI 3.14159265358979323846

typedef struct {
  int N_mode;
  double *kmod;
  double *kx;
  double *ky;
  double *F_kmod;
  double *F_kxky;
  double *phase;
  double *omega;
} T_Spectrum;

void cart2pol(double x, double y, double *rho, double *phi) {
  *rho = sqrt(x * x + y * y);
  *phi = atan2(y, x);
}

void pol2cart(double rho, double phi, double *x, double *y) {
  *x = rho * cos(phi);
  *y = rho * sin(phi);
}

double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}


// Some common spectra
double spectrum_PM(double P, double kp, double kmod) {
  // Note: P here is in fact P/sqrt(g) in the PM spectrum equation
  return P * pow(kmod, -2.5) * exp(-1.25 * pow(kp / kmod, 2.0));
}
double spectrum_JONSWAP(double alpha, double kp, double kmod) {
  return alpha * pow(kmod, -3.0) * exp(-1.25 * pow(kp / kmod, 2.0));
}
double spectrum_Gaussian(double G, double span, double kp, double kmod) {
  return (G / span) * exp(-0.5 * pow((kmod - kp) / span, 2.0));
}

T_Spectrum spectrum_gen_linear(int N_mode, int N_power, double L, double P,
                               double kp) 
{

  /* The function to generate a kx-ky spectrum based on uni-directional
      spectrum and a cos^N (theta) directional spreading.
        Arguments:
            shape: the spectrum shape (a function)
            N_mode: # of modes (Right now it's N_mode for kx and N_mode+1 for
                    ky; has to much what's hard-coded in the spectrum.h
                    headerfile.
            N_power: directional spectrum spreading coefficient
            L: physical domain size
  */

  int N_kmod = 64; // Uniform grid in kmod and ktheta, can be finer than N_mode
  int N_theta = 64;
  double thetam = 0.0; // midline direction
  double theta[N_theta];
  double dtheta;
  double Dtheta[N_theta];               // for directional spectrum
  double sum = 0.0;                     // for Dtheta normalization
  double F_kmodtheta[N_kmod * N_theta] ; // directional spectrum

  T_Spectrum spectrum;
  spectrum.N_mode = N_mode;
  spectrum.kmod = (double *)calloc(N_kmod, sizeof(double));
  spectrum.kx = (double *)calloc(N_mode, sizeof(double));
  spectrum.ky = (double *)calloc(N_mode + 1, sizeof(double));  
  spectrum.F_kmod = (double *)calloc(N_kmod, sizeof(double));
  spectrum.F_kxky = (double *)calloc(N_mode * (N_mode + 1), sizeof(double));
  spectrum.phase = (double *)calloc(N_mode * (N_mode + 1), sizeof(double));
  spectrum.omega = (double *)calloc(N_mode * (N_mode + 1), sizeof(double));


  for (int i = 0; i < N_kmod; ++i) {
    spectrum.kmod[i] =
        2 * PI / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * PI / L;
  }
  for (int i = 0; i < N_theta; ++i) {
    theta[i] = -0.5 * PI + 1.0 * i / (N_theta - 1) * PI;
    Dtheta[i] = fabs(pow(cos(theta[i] - thetam), N_power));
  }

  // Normalizing Dtheta
  dtheta = theta[1] - theta[0];
  for (int i = 0; i < N_theta - 1; ++i) {
    sum = sum + dtheta * 0.5 * (Dtheta[i] + Dtheta[i+1]);
  }
  for (int i = 0; i < N_theta; ++i) {
    Dtheta[i] = Dtheta[i] / sum;
  }

  // Pick the spectrum shape : for now PM only
  // TO DO: add more shapes ?
  for (int i = 0; i < N_kmod; ++i) {
    spectrum.F_kmod[i] = spectrum_PM(P, kp, spectrum.kmod[i]);
  }
  for (int ik = 0; ik < N_kmod; ++ik) {
    for (int itt = 0; itt < N_theta; ++itt) {
      F_kmodtheta[ik * N_kmod + itt] =
          spectrum.F_kmod[ik] * Dtheta[itt] / spectrum.kmod[ik];
      //printf("F_kmodtheta (ik=%d,itt=%d) = %f\n", ik,itt,F_kmodtheta[ik * N_kmod + itt]);
      // Notice the normalize by kmod !
    }
  }

  // Uniform grid in kx ky
  for (int i = 0; i < N_mode; ++i) {
    spectrum.kx[i] = 2 * PI / L * (i + 1);
  }
  for (int i = 0; i < N_mode + 1; ++i) {
    spectrum.ky[i] = 2 * PI / L * (i - N_mode / 2);
  }

  // interp F_kmodtheta on kx,ky grid
  double rho, phi;
  
  for (int ix = 0; ix < N_mode; ++ix) {
    for (int iy = 0; iy < (N_mode + 1); ++iy) {
      // first we get polar coords
      cart2pol(spectrum.kx[ix], spectrum.ky[iy], &rho, &phi);
      // then interp at these coords
      spectrum.F_kxky[ix * N_mode + iy] = interp_lin(
        spectrum.kmod, theta, N_kmod, N_theta, rho, phi, F_kmodtheta);
    }
  }

  // Random phase
  int RANDOM=2;
  srand(RANDOM); // We can seed it differently for different runs
  int index = 0;
  double kmod = 0;
  for (int i=0; i<N_mode; i++) {
    for (int j=0; j<N_mode+1; j++) {
      index = j*N_mode + i;
      kmod = sqrt(sq(spectrum.kx[i]) + sq(spectrum.ky[j]));
      spectrum.omega[index] = sqrt(g_*kmod);
      spectrum.phase[index] = randInRange (0, 2.*PI);
    }
  }

  return spectrum;
}

T_Spectrum read_spectrum(char path[], int N_mode) {
  /* This function reads a file and extract the spectrum from it.
      The spectra can be generated by spec_gen_linear or another function,
      given that it can be read by 'read_spectrum'
   */

  int N_kmod = 64;
  


  T_Spectrum spectrum;
  spectrum.N_mode = N_mode;
  spectrum.kmod = (double *)calloc(N_kmod, sizeof(double));
  spectrum.kx = (double *)calloc(N_mode, sizeof(double));
  spectrum.ky = (double *)calloc((N_mode + 1), sizeof(double));  
  spectrum.F_kmod = (double *)calloc(N_kmod, sizeof(double));
  spectrum.F_kxky = (double *)calloc(N_mode, (N_mode + 1) * sizeof(double));
  spectrum.phase = (double *)calloc(N_mode, (N_mode + 1) * sizeof(double));
  spectrum.omega = (double *)calloc(N_mode, (N_mode + 1) * sizeof(double));
  //fprintf(stderr, "read_spectrum WIP\n");
  return spectrum;
}




// Surface elevation following the linear wave theory
// (strictly speaking we look for a solution that is the sum of linear modes, no
// need for the  hypothesis of linear theory)
double wave(double x, double y, int N_grid, T_Spectrum spec) {
  double eta = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  int index = 0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  for (int i = 0; i < N_mode; i++) {
    for (int j = 0; j < N_mode + 1; j++) {
      index = j * N_mode + i;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      eta += ampl * cos(a);
    }
  }
  // printf("wave(): ampl %f, a %f, eta %f\n", ampl, a, eta);
  return eta;
}

// Velocities following the linear wave theory
double u_x(double x, double y, double z, int N_grid, T_Spectrum spec) 
{
  int index = 0;
  double u_x = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  double z_actual = 0.0;
  double kmod = 0.0;
  double theta = 0.0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  for (int i = 0; i < N_mode; i++) {
    for (int j = 0; j < N_mode + 1; j++) {
      index = j * N_mode + i;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      z_actual = (z < ampl ? (z) : ampl);
      // fprintf(stderr, "z = %g, ampl = %g, z_actual = %g\n", z, ampl,
      // z_actual);
      kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      theta = atan(spec.ky[j] / spec.kx[i]);
      a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      u_x +=
          sqrt(g_ * kmod) * ampl * exp(kmod * z_actual) * cos(a) * cos(theta);
    }
  }
  return u_x;
}
double u_y(double x, double y, double z, int N_grid, T_Spectrum spec) 
{
  int index = 0;
  double u_y = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  double z_actual = 0.0;
  double kmod = 0.0;
  double theta = 0.0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  for (int i = 0; i < N_mode; i++) {
    for (int j = 0; j < N_mode + 1; j++) {
      index = j * N_mode + i;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      theta = atan(spec.ky[j] / spec.kx[i]);
      a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      u_y +=
          sqrt(g_ * kmod) * ampl * exp(kmod * z_actual) * cos(a) * sin(theta);
    }
  }
  return u_y;
}
double u_z(double x, double y, double z, int N_grid, T_Spectrum spec) 
{
  int index = 0;
  double u_z = 0.0;
  double ampl = 0.0;
  double a = 0.0;
  double z_actual = 0.0;
  double kmod = 0.0;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  for (int i = 0; i < N_mode; i++) {
    for (int j = 0; j < N_mode + 1; j++) {
      index = j * N_mode + i;
      ampl = sqrt(2. * spec.F_kxky[index] * dkx * dky);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      a = (spec.kx[i] * x + spec.ky[j] * y + spec.phase[index]);
      u_z += sqrt(g_ * kmod) * ampl * exp(kmod * z_actual) * sin(a);
    }
  }
  return u_z;
}
// TO DO: u_y and w !

