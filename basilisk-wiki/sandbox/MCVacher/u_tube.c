/**
# Oscillations of a liquid column in a U-tube (NS, two-phase, embed)

U-tube: two vertical arms of width $w$, centred at $x=\pm d$,
joined by a semicircular bend of centreline radius $d$ centred at
$(0, y_c)$. The tube is open at the top (atmospheric pressure).

The liquid is initially at rest with a level $h_0 + A$ in the left arm
and $h_0 - A$ in the right arm, and is then released.

## Simulation
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "navier-stokes/perfs.h"
#include "tag.h"

/**
Set `NOSLIP` to 0 for free-slip walls: this removes the wall friction
and is a good check of the inviscid frequency $\omega_0$. */

#ifndef NOSLIP
# define NOSLIP 1
#endif

/**
### Geometry and physical parameters */

double d  = 0.25;  // half distance between the arm axes (= bend radius)
double w  = 0.10;  // tube width
double yc = 0.35;  // height of the bend centre (start of the straight arms)
double h0 = 0.60;  // mean liquid level
double A  = 0.10;  // initial level perturbation
double tend = 10.;

double omega0, gamma_th, omegad;
FILE * fplog;

double utube (double x, double y) //U-tube geometry
{
  if (y > yc)
    return w/2. - fabs (fabs(x) - d);
  double r = sqrt (sq(x) + sq(y - yc));
  return min ((d + w/2.) - r, r - (d - w/2.));
}

int main()
{
  L0 = 1.;
  origin (-L0/2., 0.);
  N = 256;               // ~26 cells across the tube

  rho1 = 1000., rho2 = 1.3;
  mu1  = 0.05,  mu2  = 0.;  // nu = 5e-5: Stokes layer ~ 1-2 cells at N = 256
  G.y  = -9.81;

  TOLERANCE = 1e-4 [*];

  /**
  Theoretical frequency and damping. */

  double Lc = 2.*(h0 - yc) + M_PI*d;
  double nu = mu1/rho1;
  omega0   = sqrt (2.*9.81/Lc);
  gamma_th = NOSLIP ? sqrt (nu*omega0/2.)/w : 0.;
  omegad   = sqrt (sq(omega0) - sq(gamma_th));
  fprintf (stderr, "# Lc = %g  omega0 = %g  T0 = %g  gamma = %g  delta = %g\n",
	   Lc, omega0, 2.*M_PI/omega0, gamma_th, sqrt (2.*nu/omega0));

  fplog = fopen ("log.dat", "w");
  fprintf (fplog, "# t zL zR eta eta_theory\n");

  init_grid (N);
  run();
  fclose (fplog);
}

/**
### Boundary conditions

Walls of the tube (embedded boundary): no-slip (or free-slip).
Top of the domain: open, atmospheric pressure. */

u.n[embed] = dirichlet (0.);
#if NOSLIP
u.t[embed] = dirichlet (0.);
#else
u.t[embed] = neumann (0.);
#endif

u.n[top] = neumann (0.);
p[top]   = dirichlet (0.);

/**
### Initial condition

The liquid fills the tube up to $h_0 + A$ on the left and $h_0 - A$ on
the right. The step at $x = 0$ lies inside the solid (the bend is far
below these levels), so it is harmless. */

event init (t = 0)
{
  solid (cs, fs, utube (x, y));
  fraction (f, (x < 0. ? h0 + A : h0 - A) - y);}

/**
To get rid of the small bubbles at the wall: the solid cells are never advected, so they would keep the initial levels forever, and the interface would see a fixed step there as the column moves. After each step, the solid cells next to the fluid take
the (cut-cell weighted) average of their fluid neighbours. This acts roughly as a neutral (about 90°) contact condition. */

event solid_fraction (i++)
{
  foreach()
    if (cs[] == 0.) {
      double fsum = 0., wsum = 0.;
      foreach_neighbor (1)
	if (cs[] > 0.)
	  fsum += cs[]*f[], wsum += cs[];
      if (wsum > 0.)
	f[] = fsum/wsum;
    }
}

/**
### Diagnostics

The liquid area in each straight arm (above $y_c$) gives the arm level
$z = y_c + \text{area}/w$. The antisymmetric part
$\eta = (z_L - z_R)/2$ is insensitive to the small bias coming from
cells cut by the line $y = y_c$. */

event logfile (i++)
{
  double aL = 0., aR = 0.;
  foreach (reduction(+:aL) reduction(+:aR))
    if (y > yc && cs[] > 0.) {
      double a = f[]*cs[]*dv();
      if (x < 0.) aL += a;
      else        aR += a;
    }
  double zL = yc + aL/w, zR = yc + aR/w;
  double eta = (zL - zR)/2.;
  double eta_th = A*exp(-gamma_th*t)*(cos(omegad*t) +
				       gamma_th/omegad*sin(omegad*t));
  fprintf (fplog, "%g %g %g %g %g\n", t, zL, zR, eta, eta_th);
  fflush (fplog);
  if (i % 100 == 0)
    fprintf (stderr, "%d %g %g\n", i, t, eta);
}

event remove_droplets (i++)
{
  remove_droplets (f, threshold = 0.05, bubbles = true);
}

/**
### Outputs */

int isave = 1;
event interfaces (t += 0.02; t <= tend)
{
  char name[80];
  sprintf (name, "interface-%d.txt", isave++);
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

event movies (t = 0; t += 0.02; t <= tend)
{
  scalar m[];
  foreach()
    m[] = cs[] - 0.5;   // hide the solid in the images
  output_ppm (f,   file = "f.mp4",   n = 512, min = 0, max = 1,
	      linear = true, mask = m);
  output_ppm (u.x, file = "u_x.mp4", n = 512, min = -1, max = 1,
	      linear = true, mask = m);
  output_ppm (u.y, file = "u_y.mp4", n = 512, min = -1, max = 1,
	      linear = true, mask = m);
}

event end (t = tend)
{
  printf ("-----END-----\n");
}

/**
## Theory

Let $\eta$ be the level displacement in the left arm
($z_L = h_0 + \eta$, $z_R = h_0 - \eta$). For an inviscid liquid column
of centreline length
$$
L_c = 2\,(h_0 - y_c) + \pi d,
$$
the restoring force $2\rho g \eta$ acts on a mass $\rho L_c$
(per unit cross-section), so
$$
\ddot{\eta} + 2\gamma\dot{\eta} + \omega_0^2\,\eta = 0,
\qquad \omega_0 = \sqrt{\frac{2g}{L_c}}.
$$

For no-slip walls with a thin Stokes layer
($\delta = \sqrt{2\nu/\omega_0} \ll w$), the time-averaged dissipation
in the two wall layers of a 2D channel gives the amplitude decay rate
$$
\gamma \simeq \frac{1}{w}\sqrt{\frac{\nu\,\omega_0}{2}}.
$$
This neglects bend losses, the free-surface/contact-line dissipation
and the non-uniform velocity profile, so it is only an estimate. The
solution is
$\eta(t) = A\, e^{-\gamma t} \left[\cos \omega_d t + \frac{\gamma}{\omega_d}
\sin\omega_d t\right]$, with $\omega_d = \sqrt{\omega_0^2 - \gamma^2}$.

## Results

![Interface](u_tube/f.mp4)
![Vertical velocity](u_tube/u_y.mp4)
![Horizontal velocity](u_tube/u_x.mp4)

~~~gnuplot Level displacement in the left arm
set xlabel 't (s)'
set ylabel '{/Symbol h} (m)'
set grid
plot 'log.dat' u 1:4 w l lw 2 t 'Basilisk', \
     'log.dat' u 1:5 w l dt 2 lw 2 t 'damped oscillator (theory)'
~~~
*/