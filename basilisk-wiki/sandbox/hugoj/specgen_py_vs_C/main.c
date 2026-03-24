#define g_ 9.81
#include "spectrum.h" // Initial conditions generation

double P = 0.02;
int N_mode = 32;
int N_power = 5;
double L = 200.0;
double *kmod;
double *F_kmod;
double *x, *y;
int N_cells = 1024;
int N_kmod = 64;
double dx;
double *eta;
int index_a;

int main(){
  double kp=10*PI/L;
  
  
  // F(k)
  // Note: this is a copy of whats inside spectrum_gen_linear  
  F_kmod = (double *)malloc(N_kmod * sizeof(double));
  kmod = (double *)malloc(N_kmod * sizeof(double));
  for (int i = 0; i < N_kmod; ++i) {
    kmod[i] = 2 * PI / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * PI / L;
  }
  for (int i = 0; i < N_kmod; ++i) {
    F_kmod[i] = spectrum_PM(P, kp, kmod[i]);
  }


  
  // F(kx,ky)
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);
  
  // eta
  dx = L/(1.0*N_cells);
  x = (double *)malloc(N_cells * sizeof(double));
  y = (double *)malloc(N_cells * sizeof(double));
  for (int i=0; i<N_cells; i++) {
    x[i] = L/2 + i*dx;
    y[i] = x[i];
  }
  eta = (double *)malloc(N_cells*N_cells * sizeof(double));
  for (int i=0; i<N_cells; i++) {
    for (int j=0; j<N_cells; j++){
      index_a = i*N_cells + j;
      eta[index_a] = wave(x[i], y[j], N_cells, spectrum);
    }
  }
  //eta[] = wave(x, y, N_grid, spectrum);

  // Writting to file
  FILE *fptr0 = fopen("F_k_C", "wb");
  fwrite(F_kmod, sizeof(double), N_kmod, fptr0);
  fclose(fptr0);

  FILE *fptr1 = fopen("F_kxky_C", "wb");
  fwrite(spectrum.F_kxky, sizeof(double), N_mode*(N_mode+1), fptr1);
  fclose(fptr1);

  FILE *fptr2 = fopen("kx_C", "wb");
  fwrite(spectrum.kx, sizeof(double), N_mode, fptr2);
  fclose(fptr2);

  FILE *fptr3 = fopen("ky_C", "wb");
  fwrite(spectrum.ky, sizeof(double), N_mode+1, fptr3);
  fclose(fptr3);

  FILE *fptr4 = fopen("kmod_C", "wb");
  fwrite(kmod, sizeof(double), N_kmod, fptr4);
  fclose(fptr4);
  printf("ok");
  FILE *fptr5 = fopen("eta_C", "wb");
  fwrite(eta, sizeof(double), N_cells*N_cells, fptr5);
  fclose(fptr5);


}
