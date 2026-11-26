/**
# Impinging Planar Jets — Bertsch *et al.* (2020) & Bongarzone *et al.* (2021)

This simulation reproduces the **planar impinging jets configuration** studied in  
- *Bertsch et al., "Feedback-free microfluidic oscillator with impinging jets", Phys. Rev. Fluids, 2020*  
- *Bongarzone et al., "Impinging planar jets: hysteretic behaviour and origin of the self-sustained oscillations", J. Fluid Mech., 2021*  

In this configuration, two **opposing laminar jets** impinge at the mid-plane of a rectangular domain, forming a stagnation region and recirculation cells.  
At moderate Reynolds numbers, the steady symmetric solution can lose stability and lead to periodic oscillations even **without any external feedback** — a purely hydrodynamic oscillator.

This case studies that geometry and basic setup in a 2D viscous, incompressible flow using the *Basilisk* framework. **Many thanks to Timothée Salamon for the collaboration**

## Simulation Overview

We solve the incompressible Navier–Stokes equations with a passive tracer field `s` to visualize the jet interaction and recirculating structures.

Parameters:
- $Re = 20$ (based on jet radius $R_d$ and inlet velocity $U_0$)
- $R_d = 0.03$
- $U_0 = 1$
- Domain size $L_0 = 1$

Slip/outflow conditions are imposed on the left and right sides; two opposing inflows are prescribed at the bottom and top boundaries.

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"
#include "harmonic.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;
double Re;

face vector muv[];
FILE * fpmax;

// Harmonic analysis of the oscillation
double t_harm = 70.;           // start of the harmonic analysis
double omegas[2] = {0., 0.};   // measured pulsation, 0-terminated for harmonic.h
double t_cross[1000];          // upward zero crossings of u.x at the centre
int n_cross = 0;
double ux_prev = 0., t_prev = 0.;

int main() {

  Re = 20;         
  R_d = 0.03;      // half jet width
  L0  = 1.0;       // domain size
  U0  = 1.0;       // reference velocity

  TOLERANCE = 1e-3 [*];
  
  //Each jet is defined by a rectangular inlet of half-width `R_d`.

  u.n[bottom] = dirichlet ( U0 * (x > -R_d && x < R_d) );
  u.t[bottom] = dirichlet (0.);
  s[bottom]   = dirichlet ( 2 * (x > -R_d && x < R_d) );

  u.n[top] = dirichlet ( -U0 * (x > -R_d && x < R_d) );
  u.t[top] = dirichlet (0.);
  s[top]   = dirichlet ( -2 * (x > -R_d && x < R_d) );

  u.n[left] = neumann(0);
  p[left]   = dirichlet(0);

  u.n[right] = neumann(0);
  p[right]   = dirichlet(0);

  origin(-L0/2, 0);
  N = 128;
  init_grid(N);
  
  char param_dim[80];
  sprintf(param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf(fparam, "%g %g %g %d\n", L0, R_d, U0, N);
  fclose(fparam);

  mu = muv;
  fpmax = fopen("log.dat", "w"); 

  run();
}

/**
## Initialization

We initialize the tracer `s` to distinguish the top and bottom flow regions visually.
*/

event init (t = 0) {
  foreach()
    s[] = (y < L0/2); 
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[] * U0 * R_d / Re;
}

event logfile (i++) { 
  fprintf(stderr, "%d %g\n", i, t); 
  fprintf(fpmax, "%d %g\n", i, t);
}

/**
## Visualization Outputs

We output the vertical velocity component and the tracer field for visualization.
*/

event ppm_output (t = 0; t += 0.1; t <= 200) {

  char name1[80];
  sprintf(name1, "uY.mp4");
  output_ppm(u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf(name2, "s.mp4");
  output_ppm(s, file = name2, n = 512, min = -2, max = 2, linear = true);
}

/**
![Passive tracer field](imp_jets/s.mp4)
![Vertical velocity field](imp_jets/uY.mp4)

We observe the interaction region where the opposing jets collide at mid-height,
forming recirculating cells and a stagnation point.
For Reynolds numbers near the onset threshold ($Re \approx 20$),
small asymmetries can grow, leading to a **self-sustained oscillation** of the impingement point,
as reported by Bongarzone *et al.* (2021) and Bertsch *et al.* (2020).
*/

/**
## Mode structure (harmonic.h)

The pulsation $\omega$ of the oscillation is measured from the upward zero
crossings of $u_x$ at the centre $(0, L_0/2)$ during $[t_\text{harm}-20, t_\text{harm}]$.
From $t_\text{harm}$ to the end, the velocity is decomposed with
[harmonic.h](/src/harmonic.h):
$\mathbf{u} \simeq \mathbf{Z} + \mathbf{A}\cos(\omega t) + \mathbf{B}\sin(\omega t)$,
with $\mathbf{Z}$ the mean flow and $\sqrt{|\mathbf{A}|^2 + |\mathbf{B}|^2}$ the
local amplitude of the oscillation. $t_\text{harm}$ must be in the saturated regime.
*/

event harmonic_analysis (t += 0.01; t <= 200) {
  double ux = interpolate (u.x, 0., L0/2);

  if (t > t_harm - 20. && t <= t_harm && ux_prev < 0. && ux >= 0. && n_cross < 1000)
    t_cross[n_cross++] = t_prev + (t - t_prev)*(-ux_prev)/(ux - ux_prev);
  ux_prev = ux;
  t_prev = t;

  if (t > t_harm && omegas[0] == 0. && n_cross >= 2) {
    omegas[0] = 2.*pi*(n_cross - 1)/(t_cross[n_cross - 1] - t_cross[0]);
    fprintf (stderr, "Oscillation: omega = %g\n", omegas[0]);
  }

  if (omegas[0] > 0.) {
    harmonic_decomposition (u.x, t, omegas);
    harmonic_decomposition (u.y, t, omegas);
  }
}

event mode_output (t = end) {
  if (omegas[0] > 0. && u.x.harmonic.invertible) {
    scalar Zx = u.x.harmonic.Z, Ax = u.x.harmonic.A[0], Bx = u.x.harmonic.B[0];
    scalar Zy = u.y.harmonic.Z, Ay = u.y.harmonic.A[0], By = u.y.harmonic.B[0];
    FILE * fp = fopen ("mode.txt", "w");
    fprintf (fp, "# omega = %g\n", omegas[0]);
    foreach()
      fprintf (fp, "%g %g %g %g %g %g %g %g\n",
               x, y, Zx[], Zy[], Ax[], Bx[], Ay[], By[]);
    fclose (fp);
  }
}

event profile (t = end) {
  printf("-----END-----\n");
}


/**
~~~pythonplot Mean flow and mode structure
import numpy as np
import matplotlib.pyplot as plt

with open('mode.txt') as fh:
    omega = float(fh.readline().split('=')[1])
m = np.loadtxt('mode.txt')
xu, ix = np.unique(m[:, 0], return_inverse=True)
yu, iy = np.unique(m[:, 1], return_inverse=True)

def grid(col):
    G = np.zeros((len(yu), len(xu)))
    G[iy, ix] = m[:, col]
    return G

Zx, Zy, Ax, Bx, Ay, By = [grid(c) for c in range(2, 8)]
ext = [xu[0], xu[-1], yu[0], yu[-1]]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

im = ax[0].imshow(np.hypot(Zx, Zy), origin='lower', extent=ext, cmap='viridis')
xs = np.linspace(xu[0], xu[-1], len(xu))       # exact uniform grid for streamplot
ys = np.linspace(yu[0], yu[-1], len(yu))
ax[0].streamplot(xs, ys, Zx, Zy, color='w', linewidth=0.5, density=1.2)
ax[0].set_title(r'Mean flow $|\mathbf{Z}|$')
fig.colorbar(im, ax=ax[0])

im = ax[1].imshow(np.sqrt(Ax**2 + Bx**2 + Ay**2 + By**2), origin='lower', extent=ext, cmap='magma')
ax[1].set_title(r'Mode amplitude $\sqrt{|\mathbf{A}|^2 + |\mathbf{B}|^2}$')
fig.colorbar(im, ax=ax[1])

vmax = np.abs(Ax).max()
im = ax[2].imshow(Ax, origin='lower', extent=ext, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
ax[2].set_title(r'Mode $A_x$ ($\omega = %.3f$, $f = %.3f$)' % (omega, omega / (2 * np.pi)))
fig.colorbar(im, ax=ax[2])

for a_ in ax:
    a_.set_xlabel(r'$x$')
    a_.set_ylabel(r'$y$')
    a_.set_xlim([ext[0], ext[1]])
    a_.set_ylim([ext[2], ext[3]])

plt.tight_layout()
plt.savefig('mode.png')
~~~
*/

