/** 
 
The goal of this script is to test the ability of the spectrum.h file to
generate the same spectrum but with different directions (in the physical
space).

We compute the variance from eta in the python script. The values arent stricly
the same but this is expected because of the truncature done on the modes in the
wave() function in spectrum.h

*/
#define g_ 9.81
#include "hugoj/lib/spectrum.h"
double P = 0.02;
int N_mode = 32;
int N_power = 5;
double L = 200.0;
int N_cells = 1024;
double *eta;
double *u;
double *v;
double dx;
double x;
double y;
int index_a;

int main(){
  double dir;
  int Ndir=2;
  double kp=10*PI/L;
  dx = L/(1.0*N_cells);

  /** We generate a 2D spectrum on the half plane $k_x$ > 0 */
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);

  /** Giving a direction to the wave field at the computation of the surface
     step using the 'dir' argument. We test dir=[0,pi/4,pi/2,3pi/4] */
  fprintf(stderr, "N_cells=%d, Ndir=%d, eta=%p\n", N_cells, Ndir, (void*)eta);
  eta = (double *)malloc(N_cells*N_cells*Ndir * sizeof(double));
  u = (double *)malloc(N_cells*N_cells*Ndir * sizeof(double));
  v = (double *)malloc(N_cells*N_cells*Ndir * sizeof(double));

  for (int i=0; i<N_cells; i++) {
    x = L/2 + i*dx;
    for (int j=0; j<N_cells; j++){
      y =  L/2 + j*dx;
      for (int d=0; d<Ndir; ++d){
        dir = 1./4. * pi * d;
          //
          // 1.0*d/(1.0*Ndir);
        index_a = i*N_cells*Ndir + j*Ndir + d;
        eta[index_a] = wave(x, y, N_cells, spectrum, dir*pi); // 
        u[index_a] = u_x(x, y, eta[index_a], N_cells, spectrum, dir*pi);
        v[index_a] = u_y(x, y, eta[index_a], N_cells, spectrum, dir*pi);

      }
    }
  }

  /** Saving the file */
  FILE *fptr = fopen("eta_C", "wb");
  fwrite(eta, sizeof(double), N_cells*N_cells*Ndir, fptr);
  fclose(fptr);

  FILE *fptr2 = fopen("u_C", "wb");
  fwrite(u, sizeof(double), N_cells*N_cells*Ndir, fptr2);
  fclose(fptr2);

  FILE *fptr3 = fopen("v_C", "wb");
  fwrite(v, sizeof(double), N_cells*N_cells*Ndir, fptr3);
  fclose(fptr3);


  /** ![direction = 0.](eta_dir0.png) */

  /** ![direction = pi/4](eta_dir1.png) */

}
