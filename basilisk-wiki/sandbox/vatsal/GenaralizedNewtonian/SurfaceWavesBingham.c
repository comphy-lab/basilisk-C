/**
#  Surface wave in a Viscoplastic fluid

Using the Generalized Newtonian Fluid implementation presented in [Couette_NonNewtonian.c](/sandbox/vatsal/GenaralizedNewtonian/Couette_NonNewtonian.c),
Surface waves in a Viscoplastic medium are studied. The problem has been studied by
[Prosperetti (1976)](https://doi.org/10.1063/1.861446) for Newtonian fluid. A more general version of the problem for
The Newtonian fluid has been simulated using Basilisk (as a test case) [(here)](http://basilisk.fr/src/test/capwave.c).<br/>

# Theory
<figure>
<p align="center">
  <img src="https://www.dropbox.com/s/95jxa147spkjgj0/Schematic.png?dl=1" width="50%">
  <caption><p align="center">Schematic of the problem: Free Surface Wave in a Viscoplastic Fluid.
  We track the temporal variation of the amplitude of these oscillations $\left(\zeta\right)$.</caption>
</figure>

## Newtonian Case
Following [Prosperetti (1976)](https://doi.org/10.1063/1.861446), the equation of motion for the interface can be written as:
$$\frac{\partial^2 \zeta}{\partial t^2} + 4\epsilon\frac{\partial \zeta}{\partial t} + \zeta - 4\epsilon^2\int_0^t\left\{\frac{e^{-\epsilon(t-\theta)}}{\sqrt{\pi\epsilon(t-\theta)}} - erfc\left(\sqrt{\epsilon(t-\theta)}\right)\right\}\frac{d}{dt}\left(\zeta(\theta)\right)d\theta = 0$$
The above equation is a simplified form of the equation presented in [Prosperetti (1981)](https://doi.org/10.1063/1.863522),
under the assumption: $\rho_{\text{upper}} \to 0\:\&\:\rho_{\text{lower}} = \rho, \nu_{\text{lower}} = \nu$ and can be solved
using Laplace Transformation. If $\hat{\zeta}(s)$ is the Laplace Transformation of $\zeta(t)$,
then the solution of the above equation in the Laplace Space with $\zeta(0) = 1.0$ and $\dot{\zeta}(0) = 0.0$ is
$$ \hat{\zeta}(s) = \frac{1}{s}\left(1 - \frac{1}{1 + s^2 + 4\epsilon s + 4\epsilon^2 - 4\epsilon\sqrt{\left(1+s/\epsilon\right)}}\right) $$
Here, $\omega^2 = gk + \frac{\sigma}{\rho}k^3$ and $\epsilon = \frac{\nu k^2}{\omega}$. The above equation can be inverted using the MATLAB code by
[McClure (2016)](https://in.mathworks.com/matlabcentral/fileexchange/39035-numerical-inverse-laplace-transform) which is based on the numerical method described by
[Abate & Whitt (2006)](https://doi.org/10.1287/ijoc.1050.0137).

## Bingham Fluid

We are still working on this. It would be interesting to study the surface waves in Viscoplastic medium theoretically. Notably, we are looking for a method to
get the stopping time of these waves as a function of the Yield Stress. For now, we only have the numerical simulations.

# Code

We solve the full Navier-Stokes equation, using VOF method for interface tracking and reconstruction.
*/
#include "navier-stokes/centered.h"
/**
Smear density and viscosity field because of high ratios
*/
#define FILTERED
/**
[two-phaseVP.h](/sandbox/vatsal/GenaralizedNewtonian/two-phaseVP.h) contains the Generalized Newtonian Fluid implementation for two-phase flows.
*/
#include "two-phaseVP.h"
#include "tension.h"
/**
The density and viscosity of the upper fluid is so low that it does not
effect the flow inside the first one, resulting in a surface wave.
*/
#define Rhor (1e-3)
#define MUr (1e-3)

#define Oh 0.02

#define tmax (1.30)
#define tsnap (tmax/200.)
/**
$$ \omega_0^2 = gk + \frac{\sigma}{\rho}k^3 $$
*/
#define w0 15.75
#define fErr (5e-2)
#define VelErr (2e-1)
#define KAPPAErr (1e-3)
#define OmegaErr (1e-3)

char nameOut[80], amplitudeFile[80], name[80];
int counter;
static FILE * fp1 = NULL;
static FILE * fp2 = NULL;
#define MAXlevel 7
#define MINlevel 0
int LEVEL, counter;
/**
We make sure that the boundary conditions for the face-centered
velocity field are consistent with the centered velocity field (this
affects the advection term). */

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

/**
The initial condition is a small amplitude plane wave of wavelength ($\lambda$)
unity. This wavelength is the relevant length scale for this problem. Note that $k = 2\pi$ */
event init (t = 0) {
  fraction (f, - y + 0.01*cos (2.*pi*x));
}

int main() {
  /**
  The domain is 2x2 to minimise finite-size effects. The surface
  tension is one. */
  L0 = 2.0;
  Y0 = -L0/2.;
  f.sigma = 1.0;
  TOLERANCE = 1e-6;

  // // Uncomment the following lines to run Grid Independent Study.
  // for (LEVEL = MAXlevel-3; LEVEL < MAXlevel+1; LEVEL++){
  //   counter = 0;
  //   init_grid(1 << LEVEL);
  //   rho1 = 1.0; rho2 = Rhor;
  //   mu1 = Oh; mu2 = MUr*Oh;
  //   mumax = (1e0)*Oh;
  //   sprintf (name, "kLEVEL-%d", LEVEL);
  //   sprintf (amplitudeFile, "AmpLEVEL-%d", LEVEL);
  //   fp1 = fopen (name, "w");
  //   fp2 = fopen (amplitudeFile, "w");
  //   char comm[80];
  //   sprintf (comm, "mkdir -p intermediate-%d", counter);
  //   system(comm);
  //   run();
  //   fclose (fp1);
  //   fclose (fp2);
  // }
  /**
  Cases for variation of Yield Stress
  */
  for (counter = 1; counter < 5; counter++){
    LEVEL = MAXlevel;
    init_grid(1 << LEVEL);
    rho1 = 1.0; rho2 = Rhor;
    mu1 = Oh; mu2 = MUr*Oh;

    switch (counter) {
    case 1: // Newtonian Prosperetti
      tauy = 0.0;
      mumax = (1e0)*Oh;
      sprintf (name, "k-%d", counter);
      sprintf (amplitudeFile, "Amp-%d", counter);
      break;
    case 2: // Bingham 1
      tauy = (1e-3);
      mumax = (1e5)*Oh;
      sprintf (name, "k-%d", counter);
      sprintf (amplitudeFile, "Amp-%d", counter);
    break;
    case 3: // Bingham 2
      tauy = (1e-2);
      mumax = (1e5)*Oh;
      sprintf (name, "k-%d", counter);
      sprintf (amplitudeFile, "Amp-%d", counter);
    break;
    case 4: // Bingham 3
      tauy = (1e-1);
      mumax = (1e5)*Oh;
      sprintf (name, "k-%d", counter);
      sprintf (amplitudeFile, "Amp-%d", counter);
    break;
    }
    fprintf(ferr, "mu0 = %g, tauy = %g, mumax = %g\n",mu1, tauy, mumax);
    char comm[80];
    sprintf (comm, "mkdir -p intermediate-%d", counter);
    system(comm);
    fp1 = fopen (name, "w");
    fp2 = fopen (amplitudeFile, "w");
    run();
    fclose (fp1);
    fclose (fp2);
  }

}


event adapt (i++) {
  scalar KAPPA[], omega[];
  curvature(f, KAPPA);
  vorticity (u, omega);
  boundary ((scalar *){omega});
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, omega},
     (double[]){fErr, VelErr, VelErr, KAPPAErr, OmegaErr},
      maxlevel = MAXlevel, minlevel = MINlevel);
}

/**
## Writing Output files
*/

event writingFiles (t += tsnap; t <= tmax) {
  dump (file = "dump");
  sprintf (nameOut, "intermediate-%d/snapshot-%6.5f", counter, t);
  dump(file=nameOut);
}

/**
The calculation of amplitude of the surface wave is same as done [(here)](http://basilisk.fr/src/test/capwave.c).<br/>
By default tracers are defined at $t-\Delta t/2$. We use the *first*
keyword to move VOF advection before the *amplitude* output i.e. at
$t+\Delta/2$. This improves the results. */

event vof (i++, first);

/**
We output the amplitude of the standing surface wave.
*/

event amplitude (i++; t <= tmax) {

  /**
  To get an accurate amplitude, we reconstruct interface position
  (using height functions) and take the corresponding maximum. */

  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;

  /**
  We output the corresponding evolution in a file indexed with the
  case number. */

  fprintf (fp2, "%g %g\n", t*w0, max);
  fflush (fp2);

}

event logfile (i++; t <= tmax) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*RHO(sf[]);
  }
  fprintf (fp1, "%g %g %d\n", t*w0, ke, mgp.i);
  fprintf (ferr, "%d %g %g %d\n", i, t*w0, ke, mgp.i);
  fflush (fp1);
}

/**
## Running the code

Use the following `run.sh` script

~~~bash
#!/bin/bash
qcc -O2 -Wall SurfaceWavesBingham.c -o SurfaceWavesBingham -lm
./SurfaceWavesBingham
~~~

# Output and Results
The post-processing codes and simulation data are available at: [PostProcess](https://www.dropbox.com/sh/x27mvdqc6725pwf/AACOfh5tVijv-Q0Givh_nD0Ga?dl=0)
<figure>
<p align="center">
  <img src="https://www.dropbox.com/s/x2dd39ag8fcwgyw/Amplitude.png?dl=1" width="50%">
  <figcaption><p align="center">Variation of the amplitude of oscillations with time. Initially, all fluids yield because of the high-stress initial condition.
  In time, oscillation is arrested as the fluid reaches the plastic limit. This stopping time decreases with increasing yield stress $\tau_y$.</figcaption>
</figure>
## Newtonian Surface Waves
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/1fhothyjhwgy22k/Case1.mp4?dl=1" type="video/mp4">
  <caption><p align="center">Surface wave oscillation for Newtonian Case. The oscillations are damped because of viscous forces. We can also notice the
  diffusion of vortices from the free surface to the bulk as discussed in [Prosperetti (1976)](https://doi.org/10.1063/1.861446) and
  [Prosperetti (1981)](https://doi.org/10.1063/1.863522).</caption>
</video>
## Surface Waves: Bingham $\tau_y = 0.01$
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/ppb584yyb3jt28f/Case3.mp4?dl=1" type="video/mp4">
  <caption><p align="center">Interfacial oscillations for $\tau_y = 0.01$. The oscillation ceases after three cycles. We show the yield surface on the left
   and vorticity contour on the right. Notice the oscillation of unyielded surface with the surface waves. For visualization of the unyielded surface,
   $log_{10}\left(\|D_{ij}\|\right)$ is plotted, with threshold at -4.</caption>
</video>
## Surface Waves: Bingham $\tau_y = 0.10$
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/wj8h6gyh2mtwq45/Case4.mp4?dl=1" type="video/mp4">
  <caption><p align="center">Interfacial oscillations for $\tau_y = 0.1$. The oscillation ceases before the completion of the first cycle. For visualization of the unyielded surface,
   $log_{10}\left(\|D_{ij}\|\right)$ is plotted, with threshold at -4. Vorticies freeze near the interface.</caption>
</video>
# Cautionary Note
We should be careful in selecting the regularisation parameter, $\mu_{\text{max}}$. It should not affect the results of the simulation. Therefore, we carried out a
sensitivity study to ensure that the temporal variation of the amplitude is unaffected on changing $\mu_{\text{max}}$.
<figure>
<p align="center">
  <img src="https://www.dropbox.com/s/xqrpgzogkpdzkzk/AmplitudeSensitivity.png?dl=1" width="50%">
  <figcaption><p align="center">Sensitivity Test: $\zeta(t)$ saturates after $\mu_{\text{max}} = (1e5)Oh$. For these simulations, $\tau_y$ = 0.10</figcaption>
</figure>
*/
