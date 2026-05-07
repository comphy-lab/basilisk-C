/** 

# Direction.c

The goal of this script is to test the ability of the spectrum.h file to
generate the same spectrum but with different directions (in the physical
space).

We give a direction to a azimuthal integrated spectrum following
$$
\begin{aligned}
  E(k,\theta) = \frac{\phi(k)}{k} \psi(\theta)
\end{aligned}
$$

we choose $\psi$ following [Wu et al., 2023](#Wu2023)

$$
\begin{aligned}
  \psi(\theta) = cos^N(\theta)/\int_{-\pi/2}^{\pi/2} cos^N(\theta)d\theta
\end{aligned}
$$

and $\phi$ to be a Pierson-Moscowitz spectrum
$$
\begin{aligned}
  \phi(k) = P g^{-1/2}k^{-2.5} exp(-1.25 (k_p/k)^2)
\end{aligned}
$$

We compute the variance from eta in the python script. The values arent stricly
the same for case with direction at pi/2 ? 

How to use this script:
make test_direction

*/
#define g_ 9.81
#include "hugoj/lib/spectrum.h"

/** First we define some parameters */

// Spectrum related
double P = 0.02;
int N_mode = 32;
int N_power = 5;
// Size of domain
double L = 200.0;
int N_cells = 256; // 1024
// Direction related
int Ndir=4;
double base_angle=1./4.; // *pi
// Other declaration
double dir;
double *eta;
double *u;
double *v;
double *w;
double dx;
double x;
double y;
int index_a;

/** ## Generate 2D spectrum and surface fields */

int main(){
  
  double kp=10*PI/L;
  dx = L/(1.0*N_cells);

  /** We generate a 2D spectrum on the half plane $k_x$ > 0 */
  T_Spectrum spectrum[Ndir];
  //T_Spectrum spectrum;
  for (int d=0; d<Ndir; ++d){
    dir = base_angle*pi*d;
    spectrum[d] = spectrum_gen_linear(N_mode, N_power, L, P, kp, dir);
  }
  //spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp, dir);

  /** Giving a direction to the wave field at the computation of the surface
     step using the 'dir' argument. We test dir=[0,pi/4,pi/2,3pi/4] */
  eta = (double *)calloc(N_cells*N_cells*Ndir , sizeof(double));
  u = (double *)calloc(N_cells*N_cells*Ndir , sizeof(double));
  v = (double *)calloc(N_cells*N_cells*Ndir , sizeof(double));
  w = (double *)calloc(N_cells*N_cells*Ndir , sizeof(double));

  /** Then, given a 2D spectrum, we build eta and currents */
  for (int i=0; i<N_cells; i++) {
    x = L/2 + i*dx;
    for (int j=0; j<N_cells; j++){
      y =  L/2 + j*dx;
      for (int d=0; d<Ndir; ++d){
        dir = base_angle * pi * d;
        index_a = i*N_cells*Ndir + j*Ndir + d;
        eta[index_a] = wave(x, y, N_cells, spectrum[d]); // 
        u[index_a] = u_x(x, y, eta[index_a], N_cells, spectrum[d]);
        v[index_a] = u_y(x, y, eta[index_a], N_cells, spectrum[d]);
        w[index_a] = u_z(x, y, eta[index_a], N_cells, spectrum[d]);

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

  free(eta);
  free(u);
  free(v);
  free(w);
  for (int d=0; d<Ndir; ++d){
    free_spectrum(spectrum[d]);
  }
}

/** ![eta for direction = pi/4](eta_dir1.png) */

/** ![u.x for direction = pi/4](u_dir1.png) */

/** ![F_kxky for direction = pi/4](F_kxky_dir1.png) */


/** TODO: make this more Basilisk like */

/**
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


