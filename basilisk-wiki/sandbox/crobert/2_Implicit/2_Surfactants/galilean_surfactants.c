/**
# Advection of a rippled interface with surfactants

We test the ability of the [multilayer solver](/src/layered/hydro.h)
to transport a corrugated interface and the adsorbed specis on the surface
with a constant velocity profile, i.e. we check that the solver preserves 
the property of
[Galilean invariance](https://en.wikipedia.org/wiki/Galilean_invariance).

![Advection of a rippled interface](galilean_surfactants/wave.gif)

In the absence of body forces, the interface is initialised with a
 sinusoidal shape and a sinusoidal population of surfactants.
 The wave is then swept away with a uniform velocity profile $U$ (this
 moving interface is therefore not to be confused with a 
 [inertial-gravity propagating wave](/src/examples/tsunami.c)).

The resulting solid-body like motion is an exact solution of the
[multilayer set of equations with
the non-hydrostatic corrections](/src/layered/nh.h), with $\phi = 0$. */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "crobert/2_Implicit/surface.h"
#include "layered/remap.h"

/**
We non-dimensionalise the problem with the advection velocity $U$
and the wavelength $\lambda$, and we set the amplitude of the wave
to one tenth of the channel depth. */

#define wavelength 1.
#define U 1.
#define wavenumber (2.*pi/wavelength)
#define depth 0.5
#define amplitude (depth/10.)
#define endoftime (1.*wavelength/U)

double maxwaveerror;
double maxsurferror;
double mass_init;
/**
## Main function

The boundary conditions are set to periodic. The test is done in a
zero-gravity environment. The default minmod slope limiter is turned
off to avoid blunting off the wave crests. */

int main()
{
  periodic (right);
  G = 0.;
  gradient = NULL;

  breaking = 0.; // fixme: without this the computation diverges for N = 256

  for (N = 64; N <= 256; N *= 2)
    for (nl = 1; nl <= 8; nl *= 2) {
      mass_init = 0;
      run();
    }
}

/**
## Initialisation

The initial shape of the interface is a simple sine, and the velocity field
is set to a constant velocity $U = 1$ and $w = 0$ everywhere. */

event init (i = 0)
{
  foreach() {
    double H = depth + amplitude*cos(wavenumber*x);
    foreach_layer() {
      h[] = H/nl;
      u.x[] = U;
    }
    M[] = H;
  }
  maxwaveerror = 0.;
  maxsurferror = 0.;
}

/**
## Outputs

We keep track of the amplitude of the wave vs time, and for a typical
case ($N = 128$ and $\text{nl} = 4$) we export deeper checks on the
non-hydrostatic pressure $\phi$ and velocity $u$. */

event monitoring (t += endoftime/100; t = 0.; t <= endoftime)
{
  stats s = statsf(eta);
  double etaamp = fabs((s.max - s.min)/(2*amplitude) - 1);
  scalar Gamma[];
  foreach ()
    Gamma[] = M[]/(area(point));      
  s = statsf(Gamma);
  double Gamma_amp = fabs((s.max - s.min)/(2*amplitude) - 1);
  double mass_tot = 0;
  foreach ()
    mass_tot += Delta*M[];
  if (mass_init == 0)
    mass_init = mass_tot;
  
  if (N == 128) {
    char name[80];
    sprintf (name, "output-N-%d-nl-%d", N, nl);
    static FILE * fpout = fopen (name, "w");

    double absdevU = 0., absPhi = 0.;
    foreach()
      foreach_layer() {
        if (fabs(u.x[] - U) > absdevU)
	  absdevU = fabs(u.x[] - U);
	if (fabs(phi[]) > absPhi)
	  absPhi = fabs(phi[]);
      }
    fprintf (fpout, "%g %g %g %g %g %g\n", t, absPhi, absdevU/U, etaamp,
     Gamma_amp, mass_tot - mass_init);
    assert (absPhi < 2e-13 && absdevU/U < 1e-14);
  }
  maxwaveerror = max(maxwaveerror, etaamp);
  maxsurferror = max(maxsurferror, Gamma_amp);
}

/**
The reference file of this test is just based on the overall maximum wave
amplitude error. */

event errorlog (t = end)
{
  fprintf (stderr, "%d %d %.3e %.3e\n", N, nl, maxwaveerror, maxsurferror);
}
