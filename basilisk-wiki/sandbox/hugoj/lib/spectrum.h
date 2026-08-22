/**
 
This script is used to generate a synthetic sea state from a 1D spectrum
following [Wu et al., 2023](#Wu2023) and their python implementation 
 ([here](https://github.com/jiarong-wu/multilayer_breaking/blob/main/specgen/specgen.py)).
We give a direction to a azimuthal integrated spectrum

$$
\begin{aligned}
  E(k,\theta) = \frac{\phi(k)}{k} \psi(\theta)
\end{aligned}
$$

we choose $\psi$ to be

$$
\begin{aligned}
  \psi(\theta) = cos^N(\theta-\theta_m)/\int_{-\pi/2}^{\pi/2} cos^N(\theta-\theta_m)d\theta
\end{aligned}
$$

with $\theta_m$ the main direction (positive anti-clockwise, zero when aligned with $x$). 
As of now, $\phi$ is a Pierson-Moscowitz spectrum

$$
\begin{aligned}
  \phi(k) = P g^{-1/2}k^{-2.5} exp(-1.25 (k_p/k)^2)
\end{aligned}
$$

From $E(k,\theta)$ we use the Airy theory to build the surface elevation $\eta$
and currents ($u$,$v$,$w$).

Note: g_ has to be defined before #include spectrum.h !

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_fft_complex.h>
//#include <string.h>
#include "hugoj/lib/interpolate.h"
#define PI 3.14159265358979323846

#ifndef REAL
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#endif

#pragma autolink -lgsl

#define CHECK_ALLOC(ptr, name) \
    if ((ptr) == NULL) { \
        fprintf(stderr, "Error: allocation failed for %s\n", name); \
        exit(EXIT_FAILURE); \
    }
#define CHECK_FILE(fp, filename) \
    if (fp == NULL) { \
        fprintf(stderr, "Error: could not open file %s\n", filename); \
        exit(EXIT_FAILURE); \
      }
#define CHECK_NDOUBLE(nread, fp, N, filename) \
    if (nread != (size_t)N) { \
        fprintf(stderr, "Error: expected %d doubles, got %zu in file %s\n", N, nread, filename); \
        fclose(fp); \
        exit(EXIT_FAILURE); \
    }



typedef struct {
  int N_mode;
  double thetam;
  double *kx;
  double *ky;
  double *F_kxky;
  double *phase;
  double *omega;
} T_Spectrum;

void free_spectrum(T_Spectrum *s) {
    free(s->kx);
    free(s->ky);
    free(s->F_kxky);
    free(s->phase);
    free(s->omega);
    
    // set to NULL so that I can check later if it's allocated
    s->kx = NULL;
    s->ky = NULL;
    s->F_kxky = NULL;
    s->phase = NULL;
    s->omega = NULL;
}

void alloc_spectrum(T_Spectrum *s, int N_mode) {
    int Ntmode = 2*N_mode+1;
    s->N_mode = N_mode;
    s->kx = (double *)calloc(Ntmode, sizeof(double));
    CHECK_ALLOC(s->kx, "kx");
    s->ky = (double *)calloc(Ntmode, sizeof(double));  
    CHECK_ALLOC(s->kx, "ky");
    s->F_kxky = (double *)calloc(Ntmode*Ntmode, sizeof(double));
    CHECK_ALLOC(s->kx, "F_kxky");
    s->phase = (double *)calloc(Ntmode*Ntmode, sizeof(double));
    CHECK_ALLOC(s->kx, "phase");
    s->omega = (double *)calloc(Ntmode*Ntmode, sizeof(double));
    CHECK_ALLOC(s->kx, "omega");
}

void cart2pol(double x, double y, double *rho, double *phi) {
  *rho = sqrt(x * x + y * y);
  *phi = atan2(y, x);
}

void pol2cart(double rho, double phi, double *x, double *y) {
  *x = rho * cos(phi);
  *y = rho * sin(phi);
}


double randInRange(double min, double max)
{
  return min + (rand() / (RAND_MAX+1.0) * (max - min));
}




/** ## Some common spectra */

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

/** ## Generate a spectrum */
trace
T_Spectrum spectrum_gen_linear(int N_mode, int N_power, double L, double P,
                               double kp, double thetam=0.) 
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
            P: energy of the wave field (dimension=velocity)
            kp: peak wavenumber
            thetam: midline direction (rad, positive anticlockwise, along x = 0.)
  */

  int N_kmod = 128; // Uniform grid in kmod and ktheta, can be finer than N_mode
  int N_theta = 128;
  //double thetam = 0.; // midline direction
  double theta[N_theta];
  double dtheta;
  double Dtheta[N_theta];               // for directional spectrum
  double sum = 0.0;                     // for Dtheta normalization
  double F_kmodtheta[N_kmod * N_theta] ; // directional spectrum
  double kmod[N_kmod];
  double F_kmod[N_kmod];
  int Ntmode = 2*N_mode + 1;

  T_Spectrum spectrum;
  alloc_spectrum(&spectrum, N_mode);

  // building kmod
  for (int i = 0; i < N_kmod; ++i) {
    kmod[i] =
        2 * pi / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * pi / L;
    // fixme: what are these weird constants?? 1.41 * 100 * 2 - 2
  }
  spectrum.thetam = thetam;
  // build Dtheta
  for (int i = 0; i < N_theta; ++i) {
    theta[i] = -pi + 2.0*pi*i / (N_theta - 1);
    Dtheta[i] = fabs(pow(cos(theta[i] - thetam), N_power));  
  }
  // Normalizing Dtheta
  dtheta = theta[1] - theta[0];
  for (int i = 0; i < N_theta - 1; ++i) {
    sum = sum + dtheta * 0.5 * (Dtheta[i] + Dtheta[i+1]); // trapezoid integ
  }
  for (int i = 0; i < N_theta; ++i) {
    Dtheta[i] = Dtheta[i] / sum;
  }

  // Pick the spectrum shape : for now PM only
  // TODO: add more shapes ?
  for (int i = 0; i < N_kmod; ++i) {
    F_kmod[i] = spectrum_PM(P, kp, kmod[i]);
  }
  for (int ik = 0; ik < N_kmod; ++ik) {
    for (int itt = 0; itt < N_theta; ++itt) {
      F_kmodtheta[ik * N_theta + itt] =
          F_kmod[ik] * Dtheta[itt] / kmod[ik];
      // Notice the normalize by kmod !
    }
  }

  // Uniform grid in kx ky
  for (int i = 0; i < Ntmode; ++i) {
    spectrum.kx[i] = -2*pi*N_mode/L + 2*pi/L*i;
  }
  for (int i = 0; i < Ntmode; ++i) {
    //spectrum.ky[i] = 2 * pi / L * (i - N_mode / 2);
    spectrum.ky[i] = -2*pi*N_mode/L + 2*pi/L*i;
  }
  
  

  // interp F_kmodtheta on kx,ky grid
  double rho, phi;
  double localkx;
  for (int ix = 0; ix < Ntmode; ++ix) {
    for (int iy = 0; iy < Ntmode; ++iy) {
      // first we get polar coords
      cart2pol(spectrum.kx[ix], spectrum.ky[iy], &rho, &phi);

      // then interp at these coords
      spectrum.F_kxky[ix*Ntmode + iy] = 2*interp_lin(
        kmod, theta, N_kmod, N_theta, rho, phi, F_kmodtheta);
      
      // we remove the negative kx values
      // this removes the tilted half plane 
      // Fixme: just not generate these ! faster :)
      localkx = spectrum.kx[ix]*cos(thetam) + spectrum.ky[iy]*sin(thetam);
      if (localkx < 0.){
        spectrum.F_kxky[ix*Ntmode + iy] = 0.; // and x2 the half plane
      }
      // we remove the center point too, not defined behavior
      if (ix==N_mode && iy==N_mode){
        spectrum.F_kxky[ix*Ntmode + iy] = 0.;
      }
      // Uncomment this if F_kxky is < 0.
      // if (spectrum.F_kxky[ix*Ntmode + iy] < 0.){
      //   fprintf(stderr,"i=%d j=%d %f \n", ix,iy, spectrum.F_kxky[ix*Ntmode + iy]);
      // }
      
    }
  }

  // Random phase
  int RANDOM=2;
  srand(RANDOM); // We can seed it differently for different runs
  int index = 0;
  double k = 0;
  for (int i=0; i<Ntmode; i++) {
    for (int j=0; j<Ntmode; j++) {
      index = i*Ntmode + j;
      k = sqrt(sq(spectrum.kx[i]) + sq(spectrum.ky[j]));
      spectrum.omega[index] = sqrt(g_*k); // we use linear dispersion relation 
      spectrum.phase[index] = randInRange (0, 2.*pi); // random phase in [0,2pi)
    }
  }
  return spectrum;
}

/**
## Write a spectrum to file 

Name of the files are:
- 2D spectrum : F_kxky
- kx
- ky
- omega
- phase
*/

trace
void write_spectrum(T_Spectrum spectrum) {

  int Ntmode = spectrum.N_mode*2+1;

  FILE *fptr1 = fopen("F_kxky", "wb");
  fwrite(spectrum.F_kxky, sizeof(double), Ntmode*Ntmode, fptr1);
  fclose(fptr1);

  FILE *fptr2 = fopen("kx", "wb");
  fwrite(spectrum.kx, sizeof(double), Ntmode, fptr2);
  fclose(fptr2);

  FILE *fptr3 = fopen("ky", "wb");
  fwrite(spectrum.ky, sizeof(double), Ntmode, fptr3);
  fclose(fptr3);

  FILE *fptr4 = fopen("omega", "wb");
  fwrite(spectrum.omega, sizeof(double), Ntmode*Ntmode, fptr4);
  fclose(fptr4);

  FILE *fptr5 = fopen("phase", "wb");
  fwrite(spectrum.phase, sizeof(double), Ntmode*Ntmode, fptr5);
  fclose(fptr5);
}


/** 
## Reading a spectrum from file
*/
trace
T_Spectrum read_spectrum(int N_mode) {
  /* This function reads a file and extract the spectrum from it.
      The spectra can be generated by spec_gen_linear or another function,
      given that it can be read by 'read_spectrum'

      We look for files like
      F_kxky, kx, ky, omega, phase

      Note: 

      - the length of the array in the files has to match with N_mode ! 

      - if the main program is run using mpi, this read_spectrum function will be
      used by all procs. As its read only, no problem. 

      - we read/write omega and phase. Doing this ensure that the same random number has been used.
      We could give this task to the main proc. 
      
   */

  T_Spectrum spectrum;
  int Ntmode=2*N_mode+1;
  alloc_spectrum(&spectrum, N_mode);
  int length1D = Ntmode;
  int length2D = Ntmode*Ntmode;
  char filename[100];
  size_t nread;

  // Read F_kxky
  sprintf (filename, "F_kxky");
  FILE * fp = fopen (filename, "rb");
  CHECK_FILE(fp, filename);
  nread = fread(spectrum.F_kxky, sizeof(double), length2D, fp);
  CHECK_NDOUBLE(nread, fp, length2D, filename)
  fclose (fp);

  // Read kx
  sprintf (filename, "kx");
  FILE * fp2 = fopen (filename, "rb");
  CHECK_FILE(fp2, filename);
  nread = fread(spectrum.kx, sizeof(double), length1D, fp2);
  CHECK_NDOUBLE(nread, fp2, length1D, filename)
  fclose (fp2);

  // Read ky
  sprintf (filename, "ky");
  FILE * fp3 = fopen (filename, "rb");
  CHECK_FILE(fp3, filename);
  nread = fread(spectrum.ky, sizeof(double), length1D, fp3);
  CHECK_NDOUBLE(nread, fp3, length1D, filename)
  fclose (fp3);
 
  // Read omega
  sprintf (filename, "omega");
  FILE * fp4 = fopen (filename, "rb");
  CHECK_FILE(fp4, filename);
  nread = fread(spectrum.omega, sizeof(double), length2D, fp4);
  CHECK_NDOUBLE(nread, fp4, length2D, filename)
  fclose (fp4);

  // Read phase
  sprintf (filename, "phase");
  FILE * fp5 = fopen (filename, "rb");
  CHECK_FILE(fp5, filename);
  nread = fread(spectrum.phase, sizeof(double), length2D, fp5);
  CHECK_NDOUBLE(nread, fp5, length2D, filename)
  fclose (fp5);

  return spectrum;
}
/** ## Inverse FFT 
 see the original version from Andrés [here](https://basilisk.fr/sandbox/acastillo/input_fields/initial_conditions_dimonte_fft2.h)
*/
void ifft2D(double *data, int NI, int NJ){  
  
  // Inverse FFT along rows 
  for (int i = 0; i < NI; ++i){
    gsl_fft_complex_radix2_backward(data + 2 * i * NJ, 1, NJ);
  }

  // Inverse FFT along columns
  double *column = malloc(2 * NI * sizeof(double));
  for (int j = 0; j < NJ; ++j){
    for (int i = 0; i < NI; ++i){
      REAL(column,i) = REAL(data, i*NJ + j);
      IMAG(column,i) = IMAG(data, i*NJ + j);
    }
    gsl_fft_complex_radix2_backward(column, 1, NI);
    for (int i = 0; i < NI; ++i)
    {
      REAL(data, i*NJ + j) = REAL(column,i);
      IMAG(data, i*NJ + j) = IMAG(column,i);
    }
  }
  free(column);
}

/**

## Surface elevation

$$
\begin{aligned}
  \eta(x,y) = \sum_{ij} a_{ij} cos(k_{x,i}x+k_{y,j}y + \phi_{rand})
\end{aligned}
$$

with $a_{ij} = \sqrt{2F(k_{x,i},k_{y,j})dk_x dk_y}$

Strictly speaking we look for a solution that is the sum of linear modes, no
need for the  hypothesis of linear theory)
*/

void transpose_NxN_inplace(double *A, int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {

            double tmp = A[j*N + i];
            A[j*N + i] = A[i*N + j];
            A[i*N + j] = tmp;
        }
    }
}

/**
 This function scatter T_Spectrum onto an FFT-ordered, zero-padded N*N grid
 ready to be used with ifft2D
 */
trace
static void eta_spectrum_scatter (double *data, T_Spectrum spec, int N)
{
  int N_mode = spec.N_mode;
  int Ntmode = 2*N_mode + 1;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];

  

  /**
  GSL FFT ordering is not the usual spectral space.
  We need to go from a space centered on (i,j)=(Nmode/2,Nmode/2)
  to a space centered on (bi,bj)=(0,0)

  example for N_mode = 2
  
  index F_kxky space from 'spec'
  |----|----|----|----|----|----|
            0    1    2

  index GSL space
  |----|----|----|----|----|----|
  0    1    2    3    4    5    6
  
  shift is done like this:
  if i>N_mode, bi = i
  else bi = N-i

  (and expand to 2D)
   */
  for (int i = 0; i < Ntmode; i++) {
    int m  = i - N_mode;
    int bi = (m >= 0) ? m : N + m;
    for (int j = 0; j < Ntmode; j++) {
      int n  = j - N_mode;
      int bj = (n >= 0) ? n : N + n;
      int index = i*Ntmode + j;
      double ampl = sqrt(2*spec.F_kxky[index]*dkx*dky);
      int o = bi*N + bj;
    
      REAL(data, o) = ampl*cos(spec.phase[index]);
      IMAG(data, o) = ampl*sin(spec.phase[index]);
      // Note: Ermitian symmetry is not enforced, so var = var(real) + var(imag)
    }
  }

}

/** spectrum -> FFT grid -> Basilisk scalar field */
trace
void initial_condition_wave_fft (scalar eta, T_Spectrum spec, int N)
{
  
  double *zdata = malloc(N*N*sizeof(double));
    if (pid() == 0) {
    double *data = malloc(2*N*N*sizeof(double));
    memset(data, 0, 2*N*N*sizeof(double));

    eta_spectrum_scatter(data, spec, N);
    
    #if CHECK_PARSEVAL
    // Verify Parseval's theorem: 
    // 0) Compute target variance
    int Ntmode = 2*spec.N_mode +1;
    double dkx = spec.kx[1]-spec.kx[0];
    double dky = spec.ky[1]-spec.ky[0];
    double variance_target = 0.;
    for (int i=0;i<Ntmode*Ntmode; ++i){
      variance_target += spec.F_kxky[i]*dkx*dky;
    }
    // 1) compute variance in spectral space
    double variance_spectral = 0.;  
    for (int i = 0; i < N*N; ++i){
      variance_spectral += sq(REAL(data, i)) + sq(IMAG(data, i));
    }
    
    // Note: I do not enforce Hermitian symmetry in 'eta_spectrum_scatter'
    // because my spectrum is tilted (I would need to impose F(-k) = F(k)*)
    // As my modes are independant, I get the energy split between
    // real and imaginary part EXACTLY. So here I divid by 2 the variance computed from the 
    // complex spectrum.
    variance_spectral /= 2.0;
    #endif // CHECK_PARSEVAL

    // spectral space -> physical space, change 'data' in place 
    ifft2D(data, N, N);
    for (int n = 0; n < N*N; n++)
      zdata[n] = REAL(data, n);

    #if CHECK_PARSEVAL
    // 2) compute variance in real space
    double mean_real = 0.;
    for (int i = 0; i < N * N; i++){
      mean_real += zdata[i];
    }
    mean_real /= 1.0*N*N;
    double variance_real = 0.;
    for (int i = 0; i < N * N; i++){
      variance_real += sq(zdata[i] - mean_real);
    }
    variance_real /= 1.0*N*N;
    fprintf(stdout, "wave: Spectral target variance (Fkxky): %f \n", variance_target);
    fprintf(stdout, "wave: Spectral variance (Parseval): %f \n", variance_spectral);
    fprintf(stdout, "wave: Real space variance (after IFFT): %f\n", variance_real);
    fprintf(stdout, "wave: Parseval check: spectral vs real space energy ratio = %f\n", 
            variance_real / variance_spectral);
    #endif // CHECK_PARSEVAL
    
    free(data);
  }
  
  #if _MPI
  // Broadcast to other mpi processes
  MPI_Bcast(zdata, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif // _MPI

  double dx = L0/N;
  double dy = dx;
  // Assign to eta
  foreach(cpu) {
    int i = (int)floor(x/dx);
    int j = (int)floor(y/dy);
    i = (i % N + N) % N ;
    j = (j % N + N) % N ;
    eta[] = zdata[i*N + j];
    // Fixme: there is still a X,Y shift of 1/2 cells when compared to wave_v1
    // I think this is because my spectrum is at every dx starting 0. but the Basilisk
    // grid is every dx starting at dx...
  }

  free(zdata);
}

/** Legacy function, use 'initial_condition_wave_fft' */
trace
double wave_v1 (double x, double y, T_Spectrum spec)
{
    double eta = 0.;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  int Ntmode = N_mode*2 + 1;

  for (int i = 0; i < Ntmode; i++)
    for (int j = 0; j < Ntmode; j++) {
      int index = i*Ntmode + j;
      double ampl = sqrt(2. * spec.F_kxky[index]*dkx*dky);
      double a = spec.kx[i]*x + spec.ky[j]*y + spec.phase[index];
      eta += ampl*cos(a);
    }

  return eta;
}

/** 
 
## Currents

$$
\begin{aligned}
  u(x,y,z) = \sum_{ij} a_{ij} \sqrt{gk} e^{kz} cos(\beta_{ij}) cos(k_{x,i}x+k_{y,j}y + \phi_{rand})
\end{aligned}
$$

with $\beta_{ij} = arctan2(k_{y,i}/k_{x,i})$ defined in $[-\pi,\pi)$, $g$ the gravity,
$k=\sqrt{k_{x,i}^2+k_{y,j}^2}$ and $\phi_{rand}$ and random phase in $[-\pi,\pi)$.

*/


/** Tools functions and struct */
typedef struct {
  int N, Nz;
  double *z;      // [Nz], physical z of each plane
  double *Ux, *Uy, *Uz;  // [Nz*N*N] each
} T_UStack;

void ustack_alloc(T_UStack *S, int N, int Nz)
{
  S->N = N;
  S->Nz = Nz;
  S->z  = malloc(Nz*sizeof(double));
  S->Ux = malloc(Nz*N*N*sizeof(double));
  S->Uy = malloc(Nz*N*N*sizeof(double));
  S->Uz = malloc(Nz*N*N*sizeof(double));
}

void ustack_free (T_UStack *S)
{
  free(S->z); 
  free(S->Ux); 
  free(S->Uy); 
  free(S->Uz); 

  S->z = NULL;
  S->Ux = NULL;
  S->Uy = NULL;
  S->Uz = NULL;
}

/**
Vertical interpolation at 2D index 'p' and altitude 'zq'.
Input data is 'field' of size Nz*N*N, 'zlev' of size Nz
 */
double column_interp (const double *field, int Nz, int N, int p,
                                const double *zlev, double zq)
{
  if (zq <= zlev[0])       
    return field[0*N*N + p];
  if (zq >= zlev[Nz - 1])  
    return field[(Nz - 1)*N*N + p];

  int k_below, k_above;
  find_ibounds(zlev, Nz, zq, &k_below, &k_above);

  double z0 = zlev[k_below],       z1 = zlev[k_above];
  double v0 = field[k_below*N*N + p], v1 = field[k_above*N*N + p];

  double t = (zq - z0)/(z1 - z0);
  return v0 + t*(v1 - v0);
}

/** This function is similar to 'eta_spectrum_scatter': 
  it makes the original F_kxky array (Ntmode**2 with Ntmode=2*N_mode+1)
  in the right shape and indexing for use with ifft2D 
  */
trace
static void u_spectrum_scatter (double *datau, 
                                double *datav, 
                                double *dataw, 
                                double z, 
                                T_Spectrum spec, 
                                int N)
{
  int N_mode = spec.N_mode;
  int Ntmode = 2*N_mode + 1;
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  
  for (int i = 0; i < Ntmode; i++) {
    int m  = i - N_mode;
    int bi = (m >= 0) ? m : N + m;
    for (int j = 0; j < Ntmode; j++) {
      int n  = j - N_mode;
      int bj = (n >= 0) ? n : N + n;
      int index = i*Ntmode + j;
      int o = bi*N + bj;

      double kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      double kmod_safe = (kmod > 0.) ? kmod : 1.;   // never zero, safe in all SIMD lanes (HPC optimisation)
      if (kmod > 0.) {
        double ampl = sqrt(2*spec.F_kxky[index]*dkx*dky); 
        // We check that z is in the water. If outside, we set the current to it's surface value
        double z_actual = (z < ampl ? (z) : ampl);
        double uampl = sqrt(g_*kmod_safe)*ampl*exp(kmod_safe*z_actual);
        double cosa = cos(spec.phase[index]);
        double sina = sin(spec.phase[index]);


        REAL(datau, o) = uampl*spec.kx[i]/kmod_safe*cosa;
        IMAG(datau, o) = uampl*spec.kx[i]/kmod_safe*sina;

        REAL(datav, o) = uampl*spec.ky[j]/kmod_safe*cosa;
        IMAG(datav, o) = uampl*spec.ky[j]/kmod_safe*sina;

        REAL(dataw, o) = uampl*cosa;
        IMAG(dataw, o) = uampl*sina;

      }
    }
  }
}

/** This function builds a 3D cartesian grid with currents. 
   This is usefull because we can then use ifft2D on z planes, then interp on layers*/
trace
T_UStack ustack_build (T_Spectrum spec, int N, double zmin, double zmax, int Nz)
{
  T_UStack S;
  ustack_alloc(&S, N, Nz);

  if (pid() == 0) {

    double dz = (zmax-zmin)/((double)Nz-1);
    for (int l = 0; l < Nz; l++) {
      S.z[l] = zmin + l*dz;
    }

    double *datau = malloc(2*N*N*sizeof(double));
    double *datav = malloc(2*N*N*sizeof(double));
    double *dataw = malloc(2*N*N*sizeof(double));

    for (int l = 0; l < Nz; l++) {
      memset(datau, 0, 2*N*N*sizeof(double));
      memset(datav, 0, 2*N*N*sizeof(double));
      memset(dataw, 0, 2*N*N*sizeof(double));

      u_spectrum_scatter(datau, datav, dataw, S.z[l], spec, N);
      ifft2D(datau, N, N);
      ifft2D(datav, N, N);
      ifft2D(dataw, N, N);

      for (int n = 0; n < N*N; n++) {
        S.Ux[l*N*N + n] = REAL(datau, n);
        S.Uy[l*N*N + n] = REAL(datav, n);
        S.Uz[l*N*N + n] = IMAG(dataw, n);
      }
    }
    free(datau); free(datav); free(dataw);
  }

  #if _MPI
    MPI_Bcast(S.Ux, Nz*N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(S.Uy, Nz*N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(S.Uz, Nz*N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif

  return S;
}


/**  main driver : scatter spectrum on iFFT grid -> interp columns -> Basilisk scalar field */
trace
void initial_condition_u_fft (vector u, 
                              T_Spectrum spec, 
                              double zmin, 
                              double zmax, 
                              int N, 
                              int Nz=10)
{
  // Generate currents on 3D cartesian grid
  T_UStack S;
  S = ustack_build (spec, N, zmin, zmax, Nz);

  double dx = L0/N;
  double dy = dx;
  
  foreach(cpu) {
    int i = (int)floor(x/dx)+1;
    int j = (int)floor(y/dy)+1;
    i = (i % N + N) % N; // X cyclic
    j = (j % N + N) % N; // Y cyclic

    // p is the position in the flattened array (N*N) of (i,j) in the 2D array (N,N)
    int p = i*N + j; 
    double z = zb[];
    foreach_layer() {
      z += h[]/2;
      u.x[] = column_interp(S.Ux, Nz, N, p, S.z, z);  
      u.y[] = column_interp(S.Uy, Nz, N, p, S.z, z); 
      w[] =  column_interp(S.Uz, Nz, N, p, S.z, z);
      z += h[]/2;
    }
  }
  ustack_free(&S);
}





/** Legacy function, use 'wave_u' */
trace
coord wave_u_v1 (double x, double y, double z, T_Spectrum spec) 
{
  coord u = {0.,0.,0.};
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  int Ntmode = 2*N_mode + 1;

  for (int i = 0; i < Ntmode; i++)
    for (int j = 0; j < Ntmode; j++) {
      double kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      // kmod == 0 corresponds to the (kx,ky) = (0,0) grid point (the zero/mean
      // mode of the discretized spectrum). This mode carries no wave direction
      // and is physically excluded from the orbital velocity sum.
      double kmod_safe = (kmod > 0.) ? kmod : 1.;   // never zero, safe in all SIMD lanes (HPC optimisation)
      if (kmod > 0.) {
        int index = i*Ntmode + j;
        double ampl = sqrt(2.*spec.F_kxky[index]*dkx*dky);
        double z_actual = (z < ampl ? (z) : ampl); // fixme: why?
        double a = spec.kx[i]*x + spec.ky[j]*y + spec.phase[index];
        double uampl = sqrt(g_*kmod_safe)*ampl*exp(kmod_safe*z_actual);
        double cosa = cos(a);
        u.x += uampl*cosa*spec.kx[i]/kmod_safe;
        u.y += uampl*cosa*spec.ky[j]/kmod_safe;
        u.z += uampl*sin(a);
      }
    }
  
  return u;
}


/**
 
## Todos

- use noise() from Basilisk to generate $\phi_{rand}$
- add the possiblity to switch to other form of spectrum
- use Basilisk's interpolation instead of mine ?
- make T_spectrum GPU compatible

## References

~~~bib
@article{Wu2023,
	title = {Breaking wave field statistics with a multi-layer model},
	volume = {968},
	issn = {0022-1120, 1469-7645},
	url = {https://www.cambridge.org/core/product/identifier/S0022112023005220/type/journal_article},
	doi = {10.1017/jfm.2023.522},
	abstract = {The statistics of breaking wave ﬁelds are characterised within a novel multi-layer framework, which generalises the single-layer Saint-Venant system into a multi-layer and non-hydrostatic formulation of the Navier–Stokes equations. We simulate an ensemble of phase-resolved surface wave ﬁelds in physical space, where strong nonlinearities, including directional wave breaking and the subsequent highly rotational ﬂow motion, are modelled, without surface overturning. We extract the kinematics of wave breaking by identifying breaking fronts and their speed, for freely evolving wave ﬁelds initialised with typical wind wave spectra. The Λ(c) distribution, deﬁned as the length of breaking fronts (per unit area) moving with speed c to c + dc following Phillips (J. Fluid Mech., vol. 156, 1985, pp. 505–531), is reported for a broad range of conditions. We recover the Λ(c) ∝ c−6 scaling without wind forcing for sufﬁciently steep wave ﬁelds. A scaling of Λ(c) based solely on the root-mean-square slope and peak wave phase speed is shown to describe the modelled breaking distributions well. The modelled breaking distributions are in good agreement with ﬁeld measurements and the proposed scaling can be applied successfully to the observational data sets. The present work paves the way for simulations of the turbulent upper ocean directly coupled to a realistic breaking wave dynamics, including Langmuir turbulence, and other sub-mesoscale processes.},
	pages = {A12},
	journaltitle = {Journal of Fluid Mechanics},
	shortjournal = {J. Fluid Mech.},
	author = {Wu, Jiarong and Popinet, Stéphane and Deike, Luc},
	urldate = {2024-08-23},
	date = {2023-08-10},
	langid = {english},
}
~~~
*/


