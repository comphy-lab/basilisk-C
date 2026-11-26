/**
# Low Reynolds jet confined in a box (bottom aspiration)

This simulation models a 2D planar jet confined inside a rectangular box. A low Reynolds instability of the jet takes place, reaching a limit cycle.
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"
#include "harmonic.h"

// Passive tracer
scalar s[];
scalar * tracers = {s};

// Physical and geometrical parameters

double U0;        // Jet inlet velocity
double R_d;       // Jet radius
double ceiling;   // Height of the solid ceiling

// Domain extents used for visualization
double xmin, ymin;
double xmax, ymax;

double Re;        // Reynolds number

FILE * fpmax; 

face vector muv[];

/**
Probes on the jet axis ($x = 0$), and harmonic analysis of the limit cycle.
The pulsation of the limit cycle is not known in advance: it is measured
from the upward zero crossings of $u_x$ at the middle probe during
$[t_\text{harm}-20, t_\text{harm}]$, then used for the harmonic decomposition
from $t_\text{harm}$ to the end. $t_\text{harm}$ must be in the limit cycle:
check it on the probe signals. */

#define NP 3
double yprobe[NP] = {0.2, 0.4, 0.6};
FILE * fprobe;

double t_harm = 70.;              // start of the harmonic analysis (or set omegas[0] by hand)
double omegas[2] = {0., 0.};      // measured pulsation, 0-terminated for harmonic.h
double t_cross[1000];             // upward zero crossings of u.x (middle probe)
int n_cross = 0;
double ux_prev = 0., t_prev = 0.;

int main() {

  // Flow parameters
  Re      = 26;      // Reynolds number
  
  N       = 256;     // Grid resolution
  int nR  = 3;       // Jet radius in cells (use nR = 6 with N = 512: same R_d)
  L0      = 1.;      // Box size
  R_d     = nR*L0/N; // Jet radius: a whole number of cells (= 0.0117)
  U0      = 1.;      // Inlet velocity
  ceiling = 0.8;     // Ceiling height

  TOLERANCE = 1e-3 [*];

  
  // Boundary conditions
  
  // Bottom inlet: top-hat velocity profile
  u.n[bottom] = dirichlet(U0*(x > -R_d && x < R_d));
  u.t[bottom] = dirichlet(0.);

  // No-slip condition on the embedded solid ceiling
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  // Tracer injected with the jet
  s[bottom] = dirichlet(U0*(x > -R_d && x < R_d));

  // Lateral walls: no-slip walls except for the lower part where there is imposed outflow
  u.n[left]  = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.t[left]  = y < R_d ? neumann(0.)  : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.)  : dirichlet(0.);
  
  // Grid initialization (N is set above)
  origin (-L0/2, 0);
  init_grid(N);
  
  char param_dim[80];
  sprintf (param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf (fparam, "%g %g %g %d\n",L0,R_d,U0,N);
  fclose (fparam);

  mu=muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
}

/**
Viscosity definition: $\mu = U_0 \times R_d / Re$ inside the fluid.
*/
event properties (i++)
{
  foreach_face() {
    muv.x[] = fm.x[]*U0*R_d/Re;
  }
}

/**
Initialization: domain extents and embedded ceiling.
*/
event init (t = 0) {

  xmin = -L0/2;
  xmax = L0/2;
  ymin = 0.0;
  ymax = ceiling;

  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
    fprobe = fopen ("probes.dat", "w");
  } else {
    fprintf(stderr, "Restarting from previous dump\n");
    fprobe = fopen ("probes.dat", "a");
  }

  solid (cs, fs, y<ceiling); // Embedded ceiling using embed.h
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

/**
Movie outputs for velocity, tracer, and pressure fields.
*/
event ppm_output (t = 0; t += 0.05; t <= 100) {
  char name[80];
  sprintf (name, "uX.mp4");
  output_ppm (u.x, file = name, n = 512, min = -U0, max = U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true , box = {{xmin, ymin}, {xmax, ymax}});

  // the pressure is defined up to a constant: we remove its mean
  stats sp = statsf (p);
  scalar pr[];
  foreach()
    pr[] = p[] - sp.sum/sp.volume;
  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (pr, file = name3, n = 512, min = -0.05*U0*U0, max = 0.05*U0*U0, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  // vorticity (saturated close to the nozzle, where it reaches ~ U0/R_d)
  scalar vort[];
  vorticity (u, vort);
  char name4[80];
  sprintf (name4, "vorticity.mp4");
  output_ppm (vort, file = name4, n = 512, min = -10.*U0/L0, max = 10.*U0/L0, linear = true, map = cool_warm, box = {{xmin, ymin}, {xmax, ymax}});
}

/**
## Probes and harmonic analysis
*/

event probes (t += 0.01; t <= 100) {
  double ux[NP], uy[NP];
  for (int k = 0; k < NP; k++) {
    ux[k] = interpolate (u.x, 0., yprobe[k]);
    uy[k] = interpolate (u.y, 0., yprobe[k]);
  }
  fprintf (fprobe, "%g", t);
  for (int k = 0; k < NP; k++)
    fprintf (fprobe, " %g %g", ux[k], uy[k]);
  fprintf (fprobe, "\n");
  fflush (fprobe);

  // pulsation of the limit cycle, from the upward zero crossings of u.x
  if (t > t_harm - 20. && t <= t_harm && ux_prev < 0. && ux[1] >= 0. && n_cross < 1000)
    t_cross[n_cross++] = t_prev + (t - t_prev)*(-ux_prev)/(ux[1] - ux_prev);
  ux_prev = ux[1];
  t_prev = t;

  if (t > t_harm && omegas[0] == 0. && n_cross >= 2) {
    omegas[0] = 2.*pi*(n_cross - 1)/(t_cross[n_cross - 1] - t_cross[0]);
    fprintf (stderr, "Limit cycle: omega = %g\n", omegas[0]);
  }

  // harmonic decomposition at the measured pulsation
  if (omegas[0] > 0.) {
    harmonic_decomposition (u.x, t, omegas);
    harmonic_decomposition (u.y, t, omegas);
  }
}

/**
At the end, we save the mean flow $Z$ and the unstable mode $(A, B)$:
$u \simeq Z + A\cos(\omega t) + B\sin(\omega t)$. */

event mode_output (t = end) {
  if (omegas[0] > 0. && u.x.harmonic.invertible) {
    scalar Zx = u.x.harmonic.Z, Ax = u.x.harmonic.A[0], Bx = u.x.harmonic.B[0];
    scalar Zy = u.y.harmonic.Z, Ay = u.y.harmonic.A[0], By = u.y.harmonic.B[0];
    FILE * fp = fopen ("mode.txt", "w");
    fprintf (fp, "# omega = %g\n", omegas[0]);
    foreach()
      fprintf (fp, "%g %g %g %g %g %g %g %g %g\n",
               x, y, cs[], Zx[], Zy[], Ax[], Bx[], Ay[], By[]);
    fclose (fp);
  }
}


/**
Periodic dump for restart and post-processing.
*/
double t_prev_dump = 0.;

event dump_state (t = 0; t += 10; t <= 100) { //t_max
  dump("restart");
  t_prev_dump = t;
  fprintf(stderr, "Dumped state at t = %g\n", t);

  char name[80];
  sprintf(name, "res_t%.1f.txt", t);
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf(fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[], s[]);
  fclose(fpres);
}

event profile (t = end) {
  printf ("----END----\n");
}

/**
## Movies

![Passive tracer](box_jet/s.mp4)
![Vertical velocity](box_jet/uY.mp4)
![Horizontal velocity](box_jet/uX.mp4)
![Pressure (mean removed)](box_jet/p.mp4)
![Vorticity](box_jet/vorticity.mp4)
*/
/**
## Limit cycle: probes and spectrum

The signals of $u_x$ on the jet axis at $y = 0.2, 0.4, 0.6$ show the growth
and saturation of the instability. The spectrum (computed on the last 40 % of
the run) gives the frequency $f$ of the limit cycle, made dimensionless with
the jet width, $St_D = f\,2R_d/U_0$, or with the jet length, $St_H = f\,H/U_0$
($H$ = ceiling height). A closed orbit in the phase portrait
$(u_x(y_1), u_x(y_3))$ is the signature of a limit cycle.

~~~pythonplot Probes, spectrum and phase portrait
import numpy as np
import matplotlib.pyplot as plt

R_d = 3 / 256          # jet radius (nR*L0/N)
U0 = 1.
H = 0.8                # ceiling
yprobe = [0.2, 0.4, 0.6]

d = np.loadtxt('probes.dat')
# after a restart, keep only the latest data for each time
keep, tmin = [], np.inf
for i in range(len(d) - 1, -1, -1):
    if d[i, 0] < tmin - 1e-9:
        keep.append(i)
        tmin = d[i, 0]
d = d[keep[::-1]]
t, ux = d[:, 0], d[:, 1::2]

sat = t > t[0] + 0.6 * (t[-1] - t[0])          # last 40 % of the run
s = ux[sat, 1] - np.mean(ux[sat, 1])
F = np.abs(np.fft.rfft(s * np.hanning(len(s))))
f = np.fft.rfftfreq(len(s), t[1] - t[0])
f_peak = f[1 + np.argmax(F[1:])]
print("f = %.4f, St_D = %.4f, St_H = %.4f" % (f_peak, f_peak * 2 * R_d / U0, f_peak * H / U0))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

for k in range(3):
    ax[0].plot(t, ux[:, k], linewidth=0.8, label=r'$y = %.1f$' % yprobe[k])
ax[0].set_xlabel(r'$t$')
ax[0].set_ylabel(r'$u_x$ on the jet axis')
ax[0].legend()
ax[0].grid(True)

ax[1].semilogy(f, F / F.max(), 'b-')
ax[1].axvline(f_peak, color='r', linestyle='--',
              label=r'$f = %.3f$, $St_D = %.3f$, $St_H = %.2f$'
              % (f_peak, f_peak * 2 * R_d / U0, f_peak * H / U0))
ax[1].set_xlim([0, 5 * f_peak])
ax[1].set_ylim([1e-4, 2])
ax[1].set_xlabel(r'$f$')
ax[1].set_ylabel('Spectrum (normalised)')
ax[1].legend()
ax[1].grid(True, which='both', alpha=0.4)

ax[2].plot(ux[sat, 0], ux[sat, 2], 'b-', linewidth=0.8)
ax[2].set_xlabel(r'$u_x(y = 0.2)$')
ax[2].set_ylabel(r'$u_x(y = 0.6)$')
ax[2].set_title('Phase portrait (limit cycle)')
ax[2].grid(True)

plt.tight_layout()
plt.savefig('probes.png')
~~~

## Unstable mode (harmonic.h)

Over the limit cycle ($t \ge t_\text{harm}$), the velocity is decomposed with
[harmonic.h](/src/harmonic.h) at the measured pulsation $\omega$:
$\mathbf{u} \simeq \mathbf{Z} + \mathbf{A}\cos(\omega t) + \mathbf{B}\sin(\omega t)$.
$\mathbf{Z}$ is the mean flow and $|\mathbf{A}|^2 + |\mathbf{B}|^2$ the local
amplitude of the oscillation. On the jet axis, writing
$u_x = Z_x + R\cos(\omega t - \varphi)$ with $\varphi = \arctan(B_x/A_x)$, the phase
$\varphi(y)$ tells how the disturbances propagate: a linear increase
$\varphi = k\,y$ is a wave travelling along the jet at the phase speed
$c = \omega/k$, whereas jumps of $\pi$ indicate a standing wave.

~~~pythonplot Mean flow and unstable mode
import numpy as np
import matplotlib.pyplot as plt

with open('mode.txt') as fh:
    omega = float(fh.readline().split('=')[1])
m = np.loadtxt('mode.txt')
x, y, cs = m[:, 0], m[:, 1], m[:, 2]
xu, ix = np.unique(x, return_inverse=True)
yu, iy = np.unique(y, return_inverse=True)

def grid(col):
    G = np.full((len(yu), len(xu)), np.nan)
    G[iy, ix] = np.where(cs > 0, m[:, col], np.nan)   # solid ceiling masked
    return G

Zx, Zy, Ax, Bx, Ay, By = [grid(c) for c in range(3, 9)]
ext = [xu[0], xu[-1], yu[0], yu[-1]]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax = list(ax)

im = ax[0].imshow(np.hypot(Zx, Zy), origin='lower', extent=ext, cmap='viridis')
ax[0].set_title(r'Mean flow $|\mathbf{Z}|$')
fig.colorbar(im, ax=ax[0])

im = ax[1].imshow(np.sqrt(Ax**2 + Bx**2 + Ay**2 + By**2), origin='lower', extent=ext, cmap='magma')
ax[1].set_title(r'Mode amplitude $\sqrt{|\mathbf{A}|^2 + |\mathbf{B}|^2}$')
fig.colorbar(im, ax=ax[1])

for a_ in ax[:2]:
    a_.set_xlabel(r'$x$')
    a_.set_ylabel(r'$y$')
    a_.set_ylim([0, 0.8])

# phase of u_x along the jet axis (column closest to x = 0)
j0 = np.argmin(np.abs(xu))
ok = np.isfinite(Ax[:, j0])
y_ax, A_ax, B_ax = yu[ok], Ax[ok, j0], Bx[ok, j0]
R_ax = np.hypot(A_ax, B_ax)
sig = R_ax > 0.05 * R_ax.max()                    # phase only where the mode exists
phi = np.unwrap(np.arctan2(B_ax[sig], A_ax[sig]))
k, phi0 = np.polyfit(y_ax[sig], phi, 1)
print("k = %.3f, phase speed c = omega/k = %.4f" % (k, omega / k))

# consistency: mode amplitude at (0, 0.4) vs probe amplitude (sqrt(2) x rms)
t_harm = 70.
p = np.loadtxt('probes.dat')
sel = p[:, 0] >= t_harm
A_probe = np.sqrt(2) * np.std(p[sel, 3])
R_mode = R_ax[np.argmin(np.abs(y_ax - 0.4))]
print("amplitude at (0, 0.4): mode = %.4f, probe = %.4f" % (R_mode, A_probe))

ax[2].plot(phi, y_ax[sig], 'b.', markersize=4, label=r'$\varphi(y)$')
ax[2].plot(phi0 + k * y_ax[sig], y_ax[sig], 'r--', linewidth=2,
           label=r'$k = %.2f$, $c = \omega/k = %.3f$' % (k, omega / k))
ax[2].set_xlabel(r'Phase $\varphi$ of $u_x$ on the axis (rad)')
ax[2].set_ylabel(r'$y$')
ax[2].set_ylim([0, 0.8])
ax[2].set_title(r'$\omega = %.3f$; at $(0, 0.4)$'
                % (omega))
ax[2].legend(loc='lower right')
ax[2].grid(True)

plt.tight_layout()
plt.savefig('mode.png')
~~~
*/