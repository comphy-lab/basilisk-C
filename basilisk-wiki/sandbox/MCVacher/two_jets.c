/**
# Jet/jet interaction

Jet 1 (left,  x = -d/2): velocity U1, radius R1, tracer s = +1

Jet 2 (right, x = +d/2): velocity U2, radius R2, tracer s = -1

Both jets are issued from two submerged needles, each made of two thin solid
walls (`embed.h`), so that the bath can be entrained from all around the
nozzle exits, as in the experiments. The flow is
initialised with two Schlichting/Bickley profiles starting at the needle exits.
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

face vector muv[];

scalar s[];                // +1 jet 1, -1 jet 2 (visualisation)
scalar c2[];               // pre-fiber fluid only: 1 in jet 2, 0 elsewhere
scalar * tracers = {s, c2};

double U1, U2;
double R1, R2;
double Re, Re2;             // Reynolds numbers of the two jets (U R / nu)
double d;
double x1j, x2j;           // jet centres

int maxlevel = 8, minlevel = 4;   // adaptive grid: finest cell L0/256 (proof of concept)
double uemax = 0.01;               // adaptation threshold on u
double cemax = 0.01;               // threshold on the tracers

double nu_val;
double H1, H2;             // needle lengths (depth of the exits below the top)
double tw;                 // needle wall thickness
double A1, B1, r01;        // Bickley constants (initial condition)
double A2, B2, r02;

FILE * fpmax;

/** ## Initial profiles definition

Bickley jet along $r$ (distance from the nozzle), lateral coordinate $\xi = x - x_\text{jet}$:
$$\psi = \frac{A}{B}\,r^{1/3}\tanh(B\eta), \qquad \eta = \frac{\xi}{r^{2/3}},$$
$$V = \frac{\partial \psi}{\partial \xi} = A\,r^{-1/3}\,\mathrm{sech}^2(B\eta), \qquad
u = -\frac{\partial \psi}{\partial r} = -\frac{A}{3B}\,r^{-2/3}\tanh(B\eta) + \frac{2A}{3}\,\xi\,r^{-4/3}\,\mathrm{sech}^2(B\eta),$$
so that the field is exactly divergence-free, with $u$ pointing towards the jet (entrainment).
The jets point downwards from the needle exits: $r = (L_0 - H) - y + r_0$ and $u_y = -V$. */

static inline double schlichting_v (double xpos, double r, double x_jet,
                                    double A, double B)
{
  double eta = (xpos - x_jet) / pow(r, 2./3.);
  double ch  = cosh(B * eta);
  return A * pow(r, -1./3.) / (ch * ch);
}

static inline double schlichting_u (double xpos, double r, double x_jet,
                                    double A, double B)
{
  double xi    = xpos - x_jet;
  double eta   = xi / pow(r, 2./3.);
  double ch    = cosh(B * eta);
  double sech2 = 1./(ch*ch);
  return - (A / (3.*B)) * pow(r, -2./3.) * tanh(B * eta)
         + (2./3.) * A * xi * pow(r, -4./3.) * sech2;
}

/**
Constants of the Bickley jet (Schlichting's book), from the kinematic momentum
flux $K = 2RU^2$:
$$A = \left(\frac{3K^2}{32\nu}\right)^{1/3}, \qquad B = \left(\frac{K}{48\nu^2}\right)^{1/3}.$$
The virtual origin $r_0 = (A/U)^3$ gives a centreline velocity $U$ at the nozzle. */

static inline void schlichting_constants (double U, double R, double nu,
                                          double *A, double *B, double *r0)
{
  double K = 2. * R * U * U;
  *A  = pow(3. * K * K / (32. * nu), 1./3.);
  *B  = pow(K / (48. * nu * nu), 1./3.);
  *r0 = pow(*A / U, 3.);
}

/** ## Needles

Needle centred at $x_j$, from the top down to the depth $H_j$: two walls of
thickness $t_w$ on each side of the channel of half-width $R_j$ (positive
inside the walls). */

static inline double needle (double x, double y, double xj, double Rj, double Hj)
{
  return min (min (y - (L0 - Hj), fabs(x - xj) - Rj), Rj + tw - fabs(x - xj));
}

/** ## Properties */

int main ()
{
  Re  = 20.;      // assisting jet (paper: 60; lower here to stay steady)
  Re2 = 0.16;     // pre-fiber jet (paper: 0.16), sets U2 below
  L0 = 0.3;      // box size
  U1 = 1.0;
  R1 = 0.005;
  R2 = 0.015;
  d  = 0.1;       // about 17 R1
  H2 = 0.1/3.;    // depth of the pre-fiber needle exit (about 6 R1)
  H1 = 2.*H2;     // assisting needle twice longer (height difference)

  /**
  The nozzles are put on whole cells of the finest grid: centres on cell faces,
  radii a whole number of cells. The injected flux and momentum are then exact. */

  double dmin = L0/(1 << maxlevel);
  R1  = max(1., round(R1/dmin))*dmin;
  R2  = max(1., round(R2/dmin))*dmin;
  x1j = -L0/2. + round((-d/2. + L0/2.)/dmin)*dmin;
  x2j = -L0/2. + round(( d/2. + L0/2.)/dmin)*dmin;
  d   = x2j - x1j;
  H1  = round(H1/dmin)*dmin;
  H2  = round(H2/dmin)*dmin;
  tw  = 3.*dmin; // needle walls: 3 cells at the finest level

  nu_val = U1 * R1 / Re;
  U2 = Re2 * nu_val / R2;

  schlichting_constants(U1, R1, nu_val, &A1, &B1, &r01);
  schlichting_constants(U2, R2, nu_val, &A2, &B2, &r02);

  TOLERANCE = 1e-3 [*];

  /**
  Boundary conditions:

  * top: inflow inside the two needles, slip wall elsewhere (far above the exits);
  * needle walls: no slip;
  * left, right and bottom: open ($p = 0$), so that the jets can entrain fluid
  from the sides and leave at the bottom, as in a large bath. */

  u.n[top] = dirichlet( - U1 * (fabs(x - x1j) < R1) - U2 * (fabs(x - x2j) < R2) );
  u.t[top] = (fabs(x - x1j) < R1 || fabs(x - x2j) < R2) ? dirichlet(0.) : neumann(0.);

  s[top] = dirichlet( 1. * (fabs(x - x1j) < R1) - 1. * (fabs(x - x2j) < R2) );
  c2[top] = dirichlet( 1. * (fabs(x - x2j) < R2) );

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  u.n[bottom] = neumann(0.);
  p[bottom]   = dirichlet(0.);

  u.n[left]  = neumann(0.);
  p[left]    = dirichlet(0.);
  u.n[right] = neumann(0.);
  p[right]   = dirichlet(0.);
  
  origin(-L0/2., 0.);
  init_grid(1 << 7);

  FILE * fparam = fopen("param_dim.txt", "w");
  fprintf(fparam, "%g %g %g %g %g %g %g %d %g %g %g\n",
          L0, R1, R2, U1, U2, d, nu_val, maxlevel, H1, H2, tw);
  fclose(fparam);

  mu = muv;
  fpmax = fopen("log.dat", "w");
  run();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[] * nu_val;
}

/** ## Initial conditions

The needles are refined at the finest level, then defined as embedded solids.
Below each exit, the jet starts as a Bickley profile (virtual origin $r_0$);
inside each needle, as a uniform inflow. */

event init (t = 0)
{
  s.gradient  = minmod2;     // limited advection: tracers stay within [0, 1]
  c2.gradient = minmod2;

  refine (y > L0 - max(H1, H2) - 0.03
          && x > x1j - 3.*R1 - tw && x < x2j + 3.*R2 + tw
          && level < maxlevel);
  solid (cs, fs, - max (needle (x, y, x1j, R1, H1), needle (x, y, x2j, R2, H2)));

  foreach() {
    double ux = 0., uy = 0.;
    if (y < L0 - H1) {
      double r = (L0 - H1) - y + r01;
      ux += schlichting_u (x, r, x1j, A1, B1);
      uy -= schlichting_v (x, r, x1j, A1, B1);
    }
    else if (fabs(x - x1j) < R1)
      uy -= U1;
    if (y < L0 - H2) {
      double r = (L0 - H2) - y + r02;
      ux += schlichting_u (x, r, x2j, A2, B2);
      uy -= schlichting_v (x, r, x2j, A2, B2);
    }
    else if (fabs(x - x2j) < R2)
      uy -= U2;
    u.x[] = cs[] > 0. ? ux : 0.;
    u.y[] = cs[] > 0. ? uy : 0.;
  }
}

/** ## Convergence monitoring

RMS change of the velocity between two iterations, and bounds of the tracers
($c_2$ must stay in $[0,1]$, $s$ in $[-1,1]$). */

scalar ux_prev[], uy_prev[];

event logfile (i++) {
  double res_u = 0.;
  int n = 0;
  foreach (reduction(+:res_u) reduction(+:n)) {
    res_u += sq(u.x[] - ux_prev[]) + sq(u.y[] - uy_prev[]);
    ux_prev[] = u.x[];
    uy_prev[] = u.y[];
    n++;
  }
  res_u = i > 0 ? sqrt(res_u/n) : 0.;

  stats sc = statsf (c2), ss = statsf (s);
  fprintf(stderr,"%d %g %g %g %g %g %g\n",i,t,res_u,sc.min,sc.max,ss.min,ss.max);
  fprintf(fpmax,"%d %g %g %g %g %g %g\n",i,t,res_u,sc.min,sc.max,ss.min,ss.max);
  fflush(fpmax);
}
 
/**
~~~pythonplot Residuals and tracer bounds
import matplotlib.pyplot as plt
import numpy as np

it, t, res, cmin, cmax, smin, smax = np.loadtxt('log.dat', unpack=True)
keep = it >= 10                    # skip the initial projection (first iterations)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

ax[0].semilogy(t[keep], res[keep], 'r-', label=r'$\mathcal{R}^{(n)}(u)$')
ax[0].set_xlabel('t')
ax[0].legend()

ax[1].plot(t, cmin, 'b-', label=r'min $c_2$')
ax[1].plot(t, cmax, 'b--', label=r'max $c_2$')
ax[1].plot(t, smin, 'r-', label=r'min $s$')
ax[1].plot(t, smax, 'r--', label=r'max $s$')
for b in (-1, 0, 1):
    ax[1].axhline(b, color='k', linewidth=0.5)
ax[1].set_xlabel('t')
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.savefig('residuals.png')
~~~
*/

/** ## Adaptive mesh

The grid is adapted on the solid, the velocity and, with a tighter criterion,
on the tracers (so that the thin pre-fiber jet stays at the finest level); the
needles and the region around their exits are always kept at the finest level
(otherwise the slow pre-fiber fluid leaving its needle is smeared by the
numerical diffusion of coarse cells). */

event adapt (i++) {
  adapt_wavelet ({cs, u, s, c2}, (double[]){1e-3, uemax, uemax, cemax, cemax},
                 maxlevel, minlevel);
  refine (y > L0 - max(H1, H2) - 0.03       // needles and region around the exits
          && x > x1j - 3.*R1 - tw && x < x2j + 3.*R2 + tw
          && level < maxlevel);
}

/** ## Movies */

event ppm_output (t = 0; t += 0.01; t <= 25) {   // end of the run
  scalar m[];                // needles in black
  foreach()
    m[] = cs[] - 0.5;
  output_ppm(u.y, file = "uY.mp4", n = 512,
             min = -1.5*max(U1,U2), max = 0., linear = true, mask = m);
  output_ppm(s, file = "s.mp4", n = 512,
             min = -1., max = 1., linear = true, mask = m);
  output_ppm(c2, file = "c2.mp4", n = 512,
             min = 0., max = 1., linear = true, mask = m);

  scalar l[];
  foreach()
    l[] = level;
  output_ppm(l, file = "level.mp4", n = 512, min = minlevel, max = maxlevel);
}

event profile (t = end) {
  printf("-----END-----\n");
  fclose(fpmax);
}

/** ## Thinning of the pre-fiber jet

At the end, on each row of the finest grid, we integrate the pre-fiber
concentration $c_2$ across the domain:

* equivalent thickness $w(y) = \int c_2\,\mathrm{d}x$ (robust to numerical smearing);
* flux $q(y) = \int c_2\,(-u_y)\,\mathrm{d}x$, which must stay equal to the injected
flux $Q_2 = 2R_2U_2$ in the steady state (conservation check);
* trajectory $x_c(y)$: centre of the cells where $c_2$ exceeds half its maximum
on the row (insensitive to the numerical noise of $c_2$ elsewhere);
* lateral transport $q_x(y) = \int c_2\,u_x\,\mathrm{d}x$.

Close to its nozzle, the pre-fiber jet runs almost horizontally towards the
assisting jet: a row then cuts it lengthwise and $w$, $x_c$ are meaningless.
They are only kept where the pre-fiber jet goes more downwards than sideways,
$|q_x| < q$ (slope below 45°). */

event thinning (t = end) {
  FILE * fp = fopen ("thinning.txt", "w");
  int n = 1 << maxlevel;
  double dmin = L0/n;
  double cr[n], xr[n];
  for (double yy = L0 - dmin/2.; yy > 0.; yy -= dmin) {
    double w = 0., q = 0., qx = 0., cmax = 0.;
    for (int k = 0; k < n; k++) {
      xr[k] = -L0/2. + (k + 0.5)*dmin;
      cr[k] = interpolate (c2, xr[k], yy);   // all cells: small over/undershoots cancel
      w  += cr[k]*dmin;
      q  += - cr[k]*interpolate (u.y, xr[k], yy)*dmin;
      qx += cr[k]*interpolate (u.x, xr[k], yy)*dmin;
      cmax = max (cmax, cr[k]);
    }
    double xs = 0., csum = 0.;                 // centre of the core (c2 > cmax/2)
    for (int k = 0; k < n; k++)
      if (cr[k] > 0.5*cmax) {
        xs += cr[k]*xr[k];
        csum += cr[k];
      }
    fprintf (fp, "%g %g %g %g %g\n", yy, w, q, csum > 0. ? xs/csum : 0., qx);
  }
  fclose (fp);
}

/**
~~~pythonplot Trajectory of the pre-fiber jet
import numpy as np
import matplotlib.pyplot as plt

L0, R1, R2, U1, U2, d, nu, maxlevel, H1, H2, tw = np.loadtxt('param_dim.txt')
y, w, q, xc, qx = np.loadtxt('thinning.txt', unpack=True)
z = L0 - y                                # depth below the top
ok = (w > 0) & (np.abs(qx) < q)           # pre-fiber jet going more down than sideways

plt.figure(figsize=(5.5, 5.5))
plt.plot(xc[ok], z[ok], 'b-', label='Pre-fiber jet')
for xj, Rj, Hj, col, lab in [(-d / 2, R1, H1, 'r', 'Assisting nozzle'),
                             (d / 2, R2, H2, 'b', 'Pre-fiber nozzle')]:
    plt.axvline(xj - Rj, color=col, linestyle='--', linewidth=1, label=lab)
    plt.axvline(xj + Rj, color=col, linestyle='--', linewidth=1)
    for xw in (xj - Rj - tw, xj + Rj):    # needle walls
        plt.gca().add_patch(plt.Rectangle((xw, 0), tw, Hj, color='k'))
plt.xlabel(r'Trajectory $x_c$')
plt.ylabel(r'Depth $z$')
plt.xlim([-L0 / 2, L0 / 2])
plt.ylim([L0, 0])
plt.gca().set_aspect('equal')             # same scale in x and z
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('trajectory.png')
~~~

~~~pythonplot Conservation of the pre-fiber flux
import numpy as np
import matplotlib.pyplot as plt

L0, R1, R2, U1, U2, d, nu, maxlevel, H1, H2, tw = np.loadtxt('param_dim.txt')
Q2 = 2 * R2 * U2
y, w, q, xc, qx = np.loadtxt('thinning.txt', unpack=True)
z = L0 - y

plt.figure(figsize=(5.5, 4.5))
plt.plot(q[w > 0] / Q2, z[w > 0], 'b-')
plt.axvline(1, color='k', linestyle='--')
plt.xlabel(r'Flux $q/Q_2$')
plt.ylabel(r'Depth $z$')
plt.xlim([0, 1.5])
plt.gca().invert_yaxis()
plt.grid(True)
plt.tight_layout()
plt.savefig('flux.png')
~~~
*/

/**
# Movies

![Passive tracer (blue=jet1, red=jet2)](two_jets/s.mp4)
![Pre-fiber fluid c2](two_jets/c2.mp4)
![Vertical velocity](two_jets/uY.mp4)
![Refinement level](two_jets/level.mp4)
*/