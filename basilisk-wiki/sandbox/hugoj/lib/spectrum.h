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
#include "hugoj/lib/interpolate.h"

typedef struct {
  int N_mode;
  double *kx;
  double *ky;
  double *F_kxky;
  double *phase;
  double *omega;
} T_Spectrum;

void free_spectrum(T_Spectrum s) {
    free(s.kx);
    free(s.ky);
    free(s.F_kxky);
    free(s.phase);
    free(s.omega);
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
  spectrum.N_mode = N_mode;
  spectrum.kx = (double *)calloc(Ntmode, sizeof(double));
  spectrum.ky = (double *)calloc(Ntmode, sizeof(double));  
  spectrum.F_kxky = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.phase = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.omega = (double *)calloc(Ntmode*Ntmode, sizeof(double));

  // building kmod
  for (int i = 0; i < N_kmod; ++i) {
    kmod[i] =
        2 * pi / L + 1.0 * i / (N_kmod - 1) * (1.41 * 100 * 2 - 2) * pi / L;
    // fixme: what are these weird constants?? 1.41 * 100 * 2 - 2
  }
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

      // Log out-of-range coordinates to diagnose extrapolation
      // if (rho < kmod[0] || rho > kmod[N_kmod-1])
      //   fprintf(stderr, "rho out of range: %f (kmod: %f to %f)\n",
      //           rho, kmod[0], kmod[N_kmod-1]);
      // if (phi < theta[0] || phi > theta[N_theta-1])
      //   fprintf(stderr, "phi out of range: %f (theta: %f to %f)\n",
      //           phi, theta[0], theta[N_theta-1]);

      // then interp at these coords
      spectrum.F_kxky[ix*Ntmode + iy] = 2*interp_lin(
        kmod, theta, N_kmod, N_theta, rho, phi, F_kmodtheta);

      // we remove the negative kx values
      localkx = spectrum.kx[ix]*cos(thetam) + spectrum.ky[iy]*sin(thetam);
      if (localkx < 0.){
        spectrum.F_kxky[ix*Ntmode + iy] = 0.; // and x2 the half plane
      }
      //
      // we remove the center point too, not defined behavior
      if (ix==N_mode && iy==N_mode){
        //fprintf(stderr, "i'm in ! ix=%d iy=%d\n", ix, iy);
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
      //index = j*N_mode + i;
      index = i*Ntmode + j;
      k = sqrt(sq(spectrum.kx[i]) + sq(spectrum.ky[j]));
      spectrum.omega[index] = sqrt(g_*k); // we use linear dispersion relation 
      spectrum.phase[index] = randInRange (0, 2.*pi); // random phase in [0,2pi]
    }
  }
  return spectrum;
}

/** 
## Reading a spectrum from file
*/

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
  spectrum.N_mode = N_mode;
  spectrum.kx = (double *)calloc(Ntmode, sizeof(double));
  spectrum.ky = (double *)calloc(Ntmode, sizeof(double));  
  spectrum.F_kxky = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.phase = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  spectrum.omega = (double *)calloc(Ntmode*Ntmode, sizeof(double));
  int length1D;
  int length2D;
  length1D = Ntmode;
  length2D = Ntmode*Ntmode;
  char filename[100];

  // Read F_kxky
  double * a = (double*) malloc (sizeof(double)*length2D);
  sprintf (filename, "F_kxky");
  FILE * fp = fopen (filename, "rb");
  fread (a, sizeof(double), length2D, fp);
  for (int i=0;i<length2D;i++) {
    spectrum.F_kxky[i] = a[i];
  }
  fclose (fp);

  // Read kx
  double * b1 = (double*) malloc (sizeof(double)*length1D);
  sprintf (filename, "kx");
  FILE * fp2 = fopen (filename, "rb");
  fread (b1, sizeof(double), length1D, fp2);
  for (int i=0;i<length1D;i++) {
    spectrum.kx[i] = b1[i];
  }
  fclose (fp2);

  // Read ky
  double * b2 = (double*) malloc (sizeof(double)*length1D);
  sprintf (filename, "ky");
  FILE * fp3 = fopen (filename, "rb");
  fread (b2, sizeof(double), length1D, fp3);
  for (int i=0;i<length1D;i++) {
    spectrum.ky[i] = b2[i];
  }
  fclose (fp3);

  // Read omega
  double * b3 = (double*) malloc (sizeof(double)*length2D);
  sprintf (filename, "omega");
  FILE * fp4 = fopen (filename, "rb");
  fread (b3, sizeof(double), length2D, fp4);
  for (int i=0;i<length2D;i++) {
    spectrum.omega[i] = b3[i];
  }
  fclose (fp4);

  // Read phase
  double * b4 = (double*) malloc (sizeof(double)*length2D);
  sprintf (filename, "phase");
  FILE * fp5 = fopen (filename, "rb");
  fread (b4, sizeof(double), length2D, fp5);
  for (int i=0;i<length2D;i++) {
    spectrum.phase[i] = b4[i];
  }
  fclose (fp5);

  return spectrum;
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

double wave (double x, double y, T_Spectrum spec)
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
 
## U current

$$
\begin{aligned}
  u(x,y,z) = \sum_{ij} a_{ij} \sqrt{gk} e^{kz} cos(\beta_{ij}) cos(k_{x,i}x+k_{y,j}y + \phi_{rand})
\end{aligned}
$$

with $\beta_{ij} = arctan2(k_{y,i}/k_{x,i})$ defined in $[-\pi,\pi)$, $g$ the gravity,
$k=\sqrt{k_{x,i}^2+k_{y,j}^2}$ and $\phi_{rand}$ and random phase in $[-\pi,\pi)$.

*/

// Velocities following the linear wave theory
coord wave_u (double x, double y, double z, T_Spectrum spec) 
{
  coord u = {0.,0.,0.};
  double dkx = spec.kx[1] - spec.kx[0];
  double dky = spec.ky[1] - spec.ky[0];
  int N_mode = spec.N_mode;
  int Ntmode = 2*N_mode + 1;

  for (int i = 0; i < Ntmode; i++)
    for (int j = 0; j < Ntmode; j++) {
      double kmod = sqrt(sq(spec.kx[i]) + sq(spec.ky[j]));
      if (kmod > 0.) { // fixme: why is kmod sometimes zero!!
        int index = i*Ntmode + j;
        double ampl = sqrt(2.*spec.F_kxky[index]*dkx*dky);
        double z_actual = (z < ampl ? (z) : ampl); // fixme: why?
        double a = spec.kx[i]*x + spec.ky[j]*y + spec.phase[index];
        double uampl = sqrt(g_*kmod)*ampl*exp(kmod*z_actual);
        double cosa = cos(a);
        u.x += uampl*cosa*spec.kx[i]/kmod;
        u.y += uampl*cosa*spec.ky[j]/kmod;
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


