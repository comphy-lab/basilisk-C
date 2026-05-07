/**
# Flow around a 2D cylinder with Navier slip condition on embedded boundaries.
Navier boundary condition is imposed _iteratively_ on an embedded boundary
as if it was a dirichlet BC.

Test case somehow derived from [karman.c](https://basilisk.fr/src/examples/karman.c)

### The idea:
I can not solve directly for
$$ \mathbf{u}_{t} +\lambda \dfrac{\partial \mathbf{u}_t}{\partial n} = 0.$$ 
Compute the normal component of the gradient of the tangential velocity
and impose it as a BC in the embedded boundary, i.e.
$$ \mathbf{u}_{t} = dirichlet \left( -\lambda \dfrac{\partial \mathbf{u}_t}{\partial n} \right) $$

This is true only if I can get the right $\mathbf{u}_t$ at first iteration.

But for a steady case, since I've got to use some iterative solver, the BC at iteration $i$ will read as: 
$$ \mathbf{u}_{t}^{i} = dirichlet \left( -\lambda \dfrac{\partial \mathbf{u}_t^{i-1}}{\partial n} \right) $$

Here I prove that for the flow around a circular cylinder @ $Re=20$, it fits quite well with available literature.

I compare with the work by [Legendre, Lauga & Magnaudet, 2009](https://arxiv.org/abs/0905.0648).

### Note:
It is not perfect, but it kinda does the job.
*/


/**
~~~pythonplot Force coefficients on a 2D cylinder with slippery surfaces.
import numpy as np
import matplotlib.pyplot as plt

# Dataset from Legendre Lauga Magnaudet, JFM 2009
# doi:10.1017/S0022112009008015
# DATA FOR Re = 20
# data from figure 2(a)
data = np.array([
    [1.01092e-2, 9.85294e-1],
    [2.02467e-2, 9.66667e-1],
    [5.04290e-2, 9.19608e-1],
    [1.00933e-1, 8.42157e-1],
    [2.01905e-1, 7.14706e-1],
    [5.01859e-1, 4.81373e-1],
    [9.96740e-1, 3.03922e-1],
    [1.99386e+0, 1.76471e-1],
    [4.96338e+0, 7.84314e-2],
    [9.93846e+0, 4.01961e-2],
    [1.99042e+1, 1.96078e-2],
])
Cd_inf = 1.33
Cd_0   = 2.045

x = data[:,0]
y = data[:,1]

# get the computed values from LEGENDRE 2009

Cd_ref = (Cd_0-Cd_inf)*y + Cd_inf;
Kn_ref = x;

### Data from my simulation
# # MAXLEVEL MINLEVEL Delta_min Reynolds Re_cell Kn lambda Cx_stop converged
data = np.loadtxt("Cx_vs_Kn_MAXLEVEL.dat", comments="#")

MAXLEVEL = data[:, 0].astype(int)
Kn       = data[:, 5]
Cd       = data[:, 7]
conv     = data[:, 8].astype(int)

plt.figure()
plt.semilogx(Kn_ref, Cd_ref, 'o-', label='Legendre et al. 2009')

for lev in np.unique(MAXLEVEL):
    m = (MAXLEVEL == lev)
    order = np.argsort(Kn[m])
    plt.semilogx(Kn[m][order], Cd[m][order], 's-', label=f'present work, LEVEL={lev}')

plt.xlabel(r'$Kn$')
plt.ylabel(r'$C_d(Kn)$')
plt.grid(True, which='both')
plt.legend()
plt.tight_layout()
plt.savefig('Cd_Vs_Kn.svg')

~~~
*/

/**
~~~pythonplot streamwise velocity profile, measured at equator

import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

files = sorted(glob.glob("ux_profile_LEVEL*.dat"))

for fname in files:
    data = np.loadtxt(fname, comments="#")

    # columns:
    # MAXLEVEL Reynolds Kn lambda x y ux
    level = data[:, 0].astype(int)
    Kn    = data[:, 2]
    y     = data[:, 5]
    ux    = data[:, 6]

    lev = level[0]

    fig, ax = plt.subplots()

    kn_values = np.unique(Kn)
    kn_positive = kn_values[kn_values > 0.]

    norm = LogNorm(vmin=kn_positive.min(), vmax=kn_positive.max())
    cmap = plt.cm.viridis

    for kn in kn_values:
        m = np.isclose(Kn, kn)
        order = np.argsort(y[m])

        color = "k" if kn == 0. else cmap(norm(kn))

        ax.plot(ux[m][order], y[m][order], color=color)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r"$Kn$")

    ax.set_xlabel(r"$u_x(x=0,y)$")
    ax.set_ylabel(r"$y$")
    ax.set_title(fr"Velocity profile, MAXLEVEL={lev}")
    ax.grid(True)
    fig.tight_layout()

    fig.savefig(f"ux_profile_LEVEL{lev}_colorbar.svg")
    #plt.close(fig)

~~~
*/



#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double R0 = 0.50;
#define xc 0.
#define yc 0.
#define endTime 50.

face vector muv[];

int MAXLEVEL;
int MINLEVEL;

double Reynolds = 20.;

double Kn = 0.;
double lambda = 0.;

scalar dnut[];

double relax_dnut = 0.005;
int ndnut = 1;

// stopping criterion
double Cx_tol = 1e-7;
double t_min_stop = 20.;
int i_min_stop = 100;
double Cx_old, Cx_stop;
int converged;

static double dnut_embed (Point point, vector u)
{
  coord b, n;
  double area = embed_geometry (point, &b, &n);

  if (area <= 0.)
    return 0.;

#if dimension == 2
  coord t = {-n.y, n.x};

  double ubx = embed_interpolate (point, u.x, b);
  double uby = embed_interpolate (point, u.y, b);

  double coefx = 0., coefy = 0.;

  double dudn = dirichlet_gradient (point, u.x, cs, n, b, ubx, &coefx);
  double dvdn = dirichlet_gradient (point, u.y, cs, n, b, uby, &coefy);

  if (coefx)
    dudn += coefx*u.x[];
  if (coefy)
    dvdn += coefy*u.y[];

  return dudn*t.x + dvdn*t.y;
#else
  return 0.;
#endif
}

static coord ub_navier_embed (Point point, double ut)
{
  coord b, n;
  double area = embed_geometry (point, &b, &n);

  coord ub = {0., 0.};

  if (area <= 0.)
    return ub;

  double nn = sqrt(sq(n.x) + sq(n.y)) + 1e-30;
  n.x /= nn;
  n.y /= nn;

  coord t = {-n.y, n.x};

  ub.x = ut*t.x;
  ub.y = ut*t.y;

  return ub;
}

// inlet/outlet
u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// embedded Navier-like BC
u.x[embed] = dirichlet(ub_navier_embed(point, -lambda*dnut[]).x);
u.y[embed] = dirichlet(ub_navier_embed(point, -lambda*dnut[]).y);

int main()
{
  size (64.);
  origin (-L0/2., -L0/2.);

  mu = muv;

  FILE * fsweep = fopen ("Cx_vs_Kn_MAXLEVEL.dat", "w");
  fprintf (fsweep,
           "# MAXLEVEL MINLEVEL Delta_min Reynolds Re_cell Kn lambda Cx_stop converged\n");
  fclose (fsweep);

  int level_min = 9;
  int level_max = 9;

  double Kn_values[10];
  Kn_values[0] = 0.;
  for (int k = 1; k < 10; k++)
    //Kn_values[k] = pow(10., -3. + (k - 1)*(5./8.)); // 1e-3 ... 1e2
  	Kn_values[k] = pow(10., -2. + (k - 1)*(3./8.)); // 1e-2, 1e1

  for (MAXLEVEL = level_min; MAXLEVEL <= level_max; MAXLEVEL++) {
    MINLEVEL = max(4, MAXLEVEL - 4);

    for (int k = 0; k < 10; k++) {
      Kn = Kn_values[k];
      lambda = Kn*R0; // since Kn = lambda/R0

      Cx_old = HUGE;
      Cx_stop = 0.;
      converged = 0;

      init_grid (1 << MAXLEVEL);
      run();

      double Delta_min = L0/(1 << MAXLEVEL);
      double Re_cell = Reynolds*Delta_min;

      FILE * fsweep = fopen ("Cx_vs_Kn_MAXLEVEL.dat", "a");
      fprintf (fsweep,
               "%d %d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %d\n",
               MAXLEVEL, MINLEVEL, Delta_min,
               Reynolds, Re_cell, Kn, lambda, Cx_stop, converged);
      fclose (fsweep);
    }
  }
}

event init (t = 0)
{
  solid (cs, fs, sq(x - xc) + sq(y - yc) - sq(R0));

  foreach()
    dnut[] = 0.;
  boundary ({dnut});

  Cx_old = HUGE;
  Cx_stop = 0.;
  converged = 0;
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/Reynolds;
}

event update_dnut (i++)
{
  if (i % ndnut == 0) {
    foreach() {
      double newdnut = (cs[] > 0. && cs[] < 1.) ?
        dnut_embed(point, u) : 0.;
      dnut[] = (1. - relax_dnut)*dnut[] + relax_dnut*newdnut;
    }
    boundary ({dnut});
  }
}

event adapt (i++)
{
  adapt_wavelet ({cs,u},
                 (double[]){1e-2, 1e-2, 1e-2},
                 MAXLEVEL, MINLEVEL);
}

event compute_forces (i += 1; t <= endTime)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double Cx = 2.*(Fp.x + Fmu.x);
  double Cy = 2.*(Fp.y + Fmu.y);
  Cx_stop = Cx;

  double residual = fabs(Cx - Cx_old)/max(fabs(Cx), 1e-30);

  static FILE * fevol = NULL;

  if (pid() == 0) {
    if (!fevol) {
      fevol = fopen ("all_temporal_evolution.dat", "w");
      fprintf (fevol, "# i t dt MAXLEVEL Reynolds Kn lambda Cx Cy residual\n");
    }

    fprintf (fevol,
             "%06d %+6.5e %+6.5e %03d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
             i, t, dt, MAXLEVEL, Reynolds, Kn, lambda, Cx, Cy, residual);

    fflush (fevol);
  }

  if (i > i_min_stop && t > t_min_stop && residual < Cx_tol) {
    converged = 1;
    return 1;
  }

  Cx_old = Cx;
}

event output_profile (t = end)
{
  char name[256];
  sprintf (name, "ux_profile_LEVEL%d.dat", MAXLEVEL);

  FILE * fp = fopen (name, Kn == 0. ? "w" : "a");

  if (Kn == 0.)
    fprintf (fp, "# MAXLEVEL Reynolds Kn lambda x y ux\n");

  double xprof = 0.;
  int Np = 401;

  for (int j = 0; j < Np; j++) {
    double yp = 0.0 + j*(3 - 0.5)/(Np - 1);
    double ux = interpolate (u.x, xprof, yp);

    fprintf (fp, "%d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
             MAXLEVEL, Reynolds, Kn, lambda, xprof, yp, ux);
  }

  fprintf (fp, "\n");
  fclose (fp);
}
