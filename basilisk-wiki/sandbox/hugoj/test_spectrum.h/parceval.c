/**

# Test of spectrum.h: parceval equality

Reminder: Parceval equalities

$$
\begin{aligned}
	<\eta^{2}> = \int_{0}^{\infty} \int_{-\pi}^{\pi} E_1(k,\theta) k dk d\theta =
	\int_{0}^{\infty} \int_{0}^{\infty} E_2 (k_x,k_y) d k_x d k_y
\end{aligned}
$$

The omnidirectionnal wavenumber spectrum is:

$$
\begin{aligned}
	\phi(k) = \int_{- \pi}^{\pi} E(k,\theta) k d \theta
\end{aligned}
$$


 */


#define g_ 9.81
#include "hugoj/lib/spectrum.h" // Initial conditions generation

double P = 0.02;
int N_mode = 32;
int N_power = 5;
double L = 200.0;
double *kmod;
double *F_kmod;
int N_cells = 512;
int N_kmod = 128;
double dx;
double x;
double y;
double *eta;
int index_a;

int main(){
  double kp=10*PI/L;
  
  
  /** F(k)
  Note: this is a copy of whats inside spectrum_gen_linear  */
  F_kmod = (double *)malloc(N_kmod * sizeof(double));
  kmod = (double *)malloc(N_kmod * sizeof(double));
  for (int i = 0; i < N_kmod; ++i) {
    kmod[i] = 2 * PI / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * PI / L;
    F_kmod[i] = spectrum_PM(P, kp, kmod[i]);
  }

  /** F(kx,ky) */
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);
  int Ntmode = N_mode*2+1;

  // Writting to file
  FILE *fptr0 = fopen("F_k", "wb");
  fwrite(F_kmod, sizeof(double), N_kmod, fptr0);
  fclose(fptr0);

  FILE *fptr1 = fopen("kmod", "wb");
  fwrite(kmod, sizeof(double), N_kmod, fptr1);
  fclose(fptr1);

  write_spectrum(spectrum);
  free_spectrum(&spectrum);


  /** ## Verification
   We open the files to check that everything went well
   */
  T_Spectrum spectrum2;
  spectrum2 = read_spectrum(N_mode);

  /** 
   Computing eta*/
  dx = L/(1.0*N_cells);
  eta = (double *)malloc(N_cells*N_cells * sizeof(double));
  for (int i=0; i<N_cells; i++) {
    x = L/2 + i*dx;
    for (int j=0; j<N_cells; j++){
      index_a = i*N_cells + j;
      y =  L/2 + j*dx;
      eta[index_a] = wave(x, y, spectrum2);
    }
  }

  /** Printing the variances */
  double sum=0.;
  for (int i=0; i<N_kmod; ++i){
    sum += F_kmod[i]*(kmod[1]-kmod[0]);
  }
  fprintf(stderr,"%f\n", sum);

  sum = 0.;
  double dkx = spectrum2.kx[1]-spectrum2.kx[0];
  double dky = spectrum2.ky[1]-spectrum2.ky[0];

  for (int i=0; i<2*spectrum2.N_mode+1; ++i){
    for (int k=0; k<2*spectrum2.N_mode+1; ++k){
      sum += spectrum2.F_kxky[i*spectrum2.N_mode + k]*dkx*dky;
    }
  }
  fprintf(stderr,"%f\n", sum);

  /** Writing eta to plot it later */
  FILE *fptr7 = fopen("eta_C", "wb");
  fwrite(eta, sizeof(double), N_cells*N_cells, fptr7);
  fclose(fptr7);
  
  // free memory
  free(eta);
  free(F_kmod);
  free(kmod);
  free_spectrum(&spectrum2);


  /** Executing the python code to work on the binary files */
  // TODO:

  //system ("python3 parceval.py");

  // Problem for now: the custom import for fftlib isnt recognized (module not
  // found)
  // Solution: make parceval.tst then py parceval.py
}

/**
 ![Synthetic wave spectrum](F_kxky_C.pdf)

![Comparison between the original PM spectrum (dashed line) and the azimuthal
integrated 1D spectra derived from the synthetic wave spectrum](comparison_Fk.pdf)

![Synthetic wave field](eta.png)

*/

