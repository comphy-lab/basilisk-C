/**
# Schlichting 2D Jet

We try to recover the famous problem of the laminar 2D jet in an immersed domain, also called the Schlichting jet or the Bickley jet. This analytical solution can be found in H. Schlichting's book "Boundary layer theory" (McGraw-Hill) in section IX.f. It was derived by H. Schlichting (1933) himself and W. Bickley (1939).

It is a monophasic problem depending only on the Reynolds number $Re$. Here we take $Re = 10$.

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#define BVIEW 1


scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;

FILE * fpmax; 

face vector muv[];

/**
We choose to impose left-right symmetry by simulating only one half of the domain (using slip boundary conditions at the right). Left and top boundaries are outflow conditions, and imposed flow at the bottom. 
*/

int main() {

  Re=10;

  R_d=0.01; //Initial radius for the jet
  L0=1;  
  U0=1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (U0*(x > L0/2-R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);

  s[bottom] = dirichlet (U0*(x > L0/2-R_d));

  u.n[left] = neumann(0);
  p[left] = dirichlet(0);

  u.n[right] = dirichlet(0.);
  p[right]=neumann(0.);

  N=256;
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

face vector u_old[];
scalar p_old[];

event init (t = 0)
{
  foreach_face(){
    u_old.x[] = u.x[];
    u_old.y[] = u.y[];
  }

  foreach()
    p_old[] = p[];
}

/**
We choose the viscosity so that it matches the Re.
*/

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*U0*R_d/Re;
}

/**
At each iteration, we calculate the RMS difference between
two successive velocity and pressure fields to monitor
numerical convergence. The residuals are defined as

$$
R^{(n)}(u) =
\sqrt{
\frac{1}{N^2}
\sum_{i=1}^{N}
\sum_{j=1}^{N}
\left[
\left(u_{i,j}^{(n+1)}-u_{i,j}^{(n)}\right)^2
+
\left(v_{i,j}^{(n+1)}-v_{i,j}^{(n)}\right)^2
\right]
}
$$

and

$$
R^{(n)}(p) =
\sqrt{
\frac{1}{N^2}
\sum_{i=1}^{N}
\sum_{j=1}^{N}
\left(p_{i,j}^{(n+1)}-p_{i,j}^{(n)}\right)^2
}
$$

The residuals are saved in "log.dat".
*/

double res_u = 0.;
double res_p = 0.;

event logfile (i++)
{
  double sum_res_u = 0.;
  double sum_res_p = 0.;

  foreach_face()
    sum_res_u += sq(u.x[] - u_old.x[]) + sq(u.y[] - u_old.y[]);

  foreach()
    sum_res_p += sq(p[] - p_old[]);

  res_u = sqrt(sum_res_u/(N*N));
  res_p = sqrt(sum_res_p/(N*N));

  fprintf (stderr, "%d %g %g %g\n",
           i, t, res_u, res_p);

  fprintf (fpmax, "%d %g %g %g\n",
           i, t, res_u, res_p);

  foreach_face(){
    u_old.x[] = u.x[];
    u_old.y[] = u.y[];
  }

  foreach()
    p_old[] = p[];
}

/**
~~~pythonplot Residuals
import matplotlib.pyplot as plt

list_t = []
list_res_u = []
list_res_p = []

with open('log.dat', encoding='utf8') as file:
    for line in file:
        if line.strip():
            temp = line.split()

            list_t.append(float(temp[1]))
            list_res_u.append(float(temp[2]))
            list_res_p.append(float(temp[3]))

plt.figure()

plt.semilogy(
    list_t,
    list_res_u,
    'r-',
    label=r'$R^{(n)}(u)$'
)

plt.semilogy(
    list_t,
    list_res_p,
    'g-',
    label=r'$R^{(n)}(p)$'
)

plt.xlabel(r'$t$')
plt.ylabel('Residual')
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig('residuals.png')
~~~

Residuals are converged around t=800.
*/

/**
We plot the vertical velocity as well as a passive tracer to observe the establishment of the flow.
*/

event ppm_output (t = 0; t += 0.05; t <= 200) {
  
    char name1[80];
    sprintf (name1, "uY.mp4");
    output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

    // optionally tracer
    char name2[80];
    sprintf (name2, "s.mp4");
    output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

event movie_streamlines (t += 0.05; t <= 200)
{
  scalar omega[];
  vertex scalar stream[];
  scalar psi[];

  vorticity(u, omega);

  psi[bottom] = dirichlet (0.);
  psi[top]    = neumann (0.);
  psi[left]   = neumann (0.);
  psi[right]  = dirichlet (0.);

  poisson (psi, omega);
  boundary ({psi});

  foreach_vertex()
    stream[] = (psi[0,-1] + psi[-1,-1] + psi[] + psi[-1])/4.;

  clear();
  view(
    width = L0/2,
    height = L0/2,
    tx     = 0,       // center horizontally
    ty     = -1.5*L0/4     // shift up so bottom aligns with y=0
);
  isoline("stream", n = 30);
  box();
  save ("streamlines.mp4");
}

/**
![Passive tracer](schlichting/s.mp4)
![Vertical velocity](schlichting/uY.mp4)
![Streamlines](schlichting/streamlines.mp4)

We can clearly see that the Neumann condition at the top for the velocity creates a slow down of the flow, because the condition does not respect the flow solution : *equations to be written soon*.

*/

/**
We save the different fields at the last iteration.
*/

event profile (t = end) {
  char name[80];
  sprintf (name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[],s[]);
  fclose(fpres);
  
  printf ("-----END-----\n");
}

/**
## Comparison with theory

We consider the boundary layer equations (also called Prandtl model) :

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0\\
        u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = \nu\frac{\partial^{2}v}{\partial x^{2}} $$

The flow has to recover zero velocity at infinity. We write this set of equation with the stream function as varaible:

$$u(x,y) = \frac{\partial \psi}{\partial y} \\ v(x, y) = -\frac{\partial \psi}{\partial x}$$

Mass conservation:

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = \frac{\partial^{2} \psi}{\partial x \partial y} - \frac{\partial^{2} \psi}{\partial y \partial x} = 0$$

Momentum:

$$-\frac{\partial\psi}{\partial y}\frac{\partial^{2} \psi}{\partial x^{2}} + \frac{\partial\psi}{\partial x}\frac{\partial^{2} \psi}{\partial y \partial x} = -\nu \frac{\partial^{3}\psi}{\partial x^{3}}$$

We look for a self-similar solution for the stream function $\psi$ that would be satisified far from the nozzle: 

$$\psi(x,y) = y^{p}F(\frac{x}{y^{q}})$$

To find $p$ and $q$, we refer to the two following conditions. First, the total y-momentum fluw has to be constant (see Schlichting's book for details):

$$ J = \rho \int_{-\infty}^{\infty} v^{2}\mathrm{d}x = cst~~~~(=\mathrm{2R_dU_0^{2}~in~this~case})$$
$$\Rightarrow y^{2(p-q)+q } \sim 1\\ \Rightarrow 2p - q = 0$$

Second, inertial is balanced by viscous forces:

$$v\frac{\partial v}{\partial y} \sim \nu\frac{\partial^{2}v}{\partial x^{2}} $$

$$\Rightarrow p+q = 1$$

Hence:

$$\psi(x,y) = y^{1/3}F(\frac{x}{y^{2/3}})$$

This is equivalent to: 

$$\psi(x,y) = \nu^{1/2}y^{1/3}F(\frac{1}{3\nu^{1/2}}\eta)$$

with $\eta = \frac{x}{y^{2/3}}$. Using this into the previous momentum equation, it gives that F is solution of the following ODE:

$$ (F')^{2} + FF'' + F''' = 0$$

It can be shown that the solution of this ODE is: 

$$ F'(\frac{1}{3\nu^{1/2}}\eta) = 2\alpha \times \mathrm{sech^{2}}(\alpha\frac{1}{3\nu^{1/2}} \eta)$$

where $\alpha$ is deduced from $J$ (see Schlichting's book). With $K = J/\rho$ the kinematic momentum flux, this gives the closed form:

$$ v(x,y) = \left(\frac{3K^2}{32\,\nu\,y}\right)^{1/3}\mathrm{sech}^2\xi, \qquad \xi = \left(\frac{K}{48\,\nu^2}\right)^{1/3}\frac{x}{y^{2/3}}$$

i.e. $y^{1/3}\,v = A\,\mathrm{sech}^2(B\eta)$ with $A = (3K^2/32\nu)^{1/3}$ and $B = (K/48\nu^2)^{1/3}$: no free constant once $K$ is known.

### Effective momentum flux and virtual origin

As in the [axisymmetric case](/sandbox/MCVacher/schlichting_axi.c), the simulated jet does not come from a point source and does not carry exactly the nozzle momentum $K = 2R_dU_0^2$ (inlet discretised on a few cells, suction along the bottom wall, outflow at the top). We therefore measure, in the window $10 \le y/D \le 40$:

* $K_\text{eff} = \langle \int v^2\,\mathrm{d}x \rangle$, the momentum flux actually carried by the jet;
* the virtual origin $y_0$ and a second estimate $K_\text{slope}$, from a linear fit of the centreline velocity $v_c^{-3} = 32\nu\,(y-y_0)/(3K^2)$.

The theory is then plotted with $K_\text{eff}$ and $y - y_0$, without any fitted constant. With the half-width $b_{1/2}$ ($v = v_c/2$), the normalised profile $v/v_c = \mathrm{sech}^2(0.8814\,x/b_{1/2})$ has no parameter at all.
*/



/**
~~~pythonplot Self-similarity check
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Simulation parameters
# ============================================================

Re = 10
V_0 = 1
R_d = 0.01
D = 2 * R_d
nu = V_0 * R_d / Re
K0 = 2 * R_d * V_0**2        # nozzle momentum flux (full jet, per unit depth)

def A(K): return (3 * K**2 / (32 * nu))**(1/3)
def B(K): return (K / (48 * nu**2))**(1/3)

# ============================================================
# Load numerical data, structured array V[y, x], jet axis first
# ============================================================

data = np.loadtxt("res_end.txt")
x_unique, ix = np.unique(data[:, 0], return_inverse=True)
y_unique, iy = np.unique(data[:, 1], return_inverse=True)
V = np.full((len(y_unique), len(x_unique)), np.nan)
V[iy, ix] = data[:, 3]

Delta = x_unique[1] - x_unique[0]
x_axis = x_unique[-1] + Delta / 2           # symmetry plane (right boundary)
xs = (x_axis - x_unique)[::-1]              # distance from the jet axis
V = V[:, ::-1]

yD_min, yD_max = 10, 40                     # avoid inlet & top boundary
win = (y_unique / D >= yD_min) & (y_unique / D <= yD_max)

# ============================================================
# Effective momentum flux and virtual origin
# ============================================================

K_y = 2 * np.sum(V**2, axis=1) * Delta      # int v^2 dx over both halves
K_eff = np.mean(K_y[win])

vc = V[:, 0]                                # centreline velocity
slope, intercept = np.polyfit(y_unique[win], vc[win]**-3, 1)
y0 = -intercept / slope                     # vc^-3 = slope (y - y0)
K_slope = np.sqrt(32 * nu / (3 * slope))

print("K_eff/K0 = %.3f, K_slope/K0 = %.3f, y0/D = %.2f"
      % (K_eff / K0, K_slope / K0, y0 / D))

# ============================================================
# Self-similar collapse
# ============================================================

eta_all, f_all, s_all, g_all, yD_all = [], [], [], [], []

for i in np.where(win)[0]:
    v = V[i, :]
    Y = y_unique[i] - y0                    # distance from virtual origin

    # half-width: first x where v < vc/2, linear interpolation
    j = np.argmax(v < vc[i] / 2)
    b12 = xs[j-1] + (vc[i]/2 - v[j-1]) * (xs[j] - xs[j-1]) / (v[j] - v[j-1])

    eta_all.extend(xs / Y**(2/3))
    f_all.extend(Y**(1/3) * v)
    s_all.extend(xs / b12)
    g_all.extend(v / vc[i])
    yD_all.extend([y_unique[i] / D] * len(xs))

eta_all, f_all = np.array(eta_all), np.array(f_all)
s_all, g_all = np.array(s_all), np.array(g_all)
yD_all = np.array(yD_all)

eta_th = np.linspace(-1, 1, 500)
s_th = np.linspace(-3, 3, 500)

# ============================================================
# Plot (mirrored for symmetry)
# ============================================================

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

for sgn in (1, -1):
    sc = ax[0].scatter(sgn * eta_all, f_all, c=yD_all, cmap='viridis',
                       vmin=yD_min, vmax=yD_max, s=5, alpha=0.6)
    ax[1].scatter(sgn * s_all, g_all, c=yD_all, cmap='viridis',
                  vmin=yD_min, vmax=yD_max, s=5, alpha=0.6)

ax[0].plot(eta_th, A(K_eff) / np.cosh(B(K_eff) * eta_th)**2, 'r', linewidth=2,
           label=r'Theory with $K_\mathrm{eff}$')
ax[0].plot(eta_th, A(K0) / np.cosh(B(K0) * eta_th)**2, 'k--', linewidth=1.5,
           label=r'Theory with nozzle $K$')
ax[0].set_xlabel(r'$\eta = x/(y-y_0)^{2/3}$')
ax[0].set_ylabel(r'$(y-y_0)^{1/3}\,v(x,y)$')
ax[0].set_xlim([-0.5, 0.5])
ax[0].legend(loc='upper left')
ax[0].grid(True)
ax[0].text(0.97, 0.6,
           r'$K_\mathrm{eff}/K = %.2f$' '\n' r'$K_\mathrm{slope}/K = %.2f$' '\n' r'$y_0/D = %.1f$'
           % (K_eff / K0, K_slope / K0, y0 / D),
           transform=ax[0].transAxes, ha='right', va='top',
           bbox=dict(facecolor='white', alpha=0.8))

ax[1].plot(s_th, 1 / np.cosh(0.8814 * s_th)**2, 'r', linewidth=2,
           label=r'$\mathrm{sech}^2(0.8814\,x/b_{1/2})$')
ax[1].set_xlabel(r'$x/b_{1/2}$')
ax[1].set_ylabel(r'$v/v_c$')
ax[1].set_xlim([-3, 3])
ax[1].legend(loc='upper left')
ax[1].grid(True)

cbar = fig.colorbar(sc, ax=ax.ravel().tolist())
cbar.set_label(r'Vertical position $y/D$')

plt.savefig('self_sim.png')
~~~

~~~pythonplot Centreline decay and momentum flux
import numpy as np
import matplotlib.pyplot as plt

Re = 10
V_0 = 1
R_d = 0.01
D = 2 * R_d
nu = V_0 * R_d / Re
K0 = 2 * R_d * V_0**2

data = np.loadtxt("res_end.txt")
x_unique, ix = np.unique(data[:, 0], return_inverse=True)
y_unique, iy = np.unique(data[:, 1], return_inverse=True)
V = np.full((len(y_unique), len(x_unique)), np.nan)
V[iy, ix] = data[:, 3]
V = V[:, ::-1]                              # jet axis first
Delta = x_unique[1] - x_unique[0]

yD_min, yD_max = 10, 40
win = (y_unique / D >= yD_min) & (y_unique / D <= yD_max)

vc = V[:, 0]
K_y = 2 * np.sum(V**2, axis=1) * Delta
K_eff = np.mean(K_y[win])
slope, intercept = np.polyfit(y_unique[win], vc[win]**-3, 1)
y0 = -intercept / slope

yD = y_unique / D

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

ax[0].plot(yD, vc**-3, 'b-', label='Simulation')
ax[0].plot(yD, 32 * nu * y_unique / (3 * K0**2), 'k--', linewidth=1.5,
           label=r'Point source at nozzle, nozzle $K$')
ax[0].plot(yD, 32 * nu * (y_unique - y0) / (3 * K_eff**2), 'r', linewidth=2,
           label=r'Point source at $y_0$, $K_\mathrm{eff}$')
ax[0].axvspan(yD_min, yD_max, color='gray', alpha=0.15, label='Comparison window')
ax[0].set_xlabel(r'$y/D$')
ax[0].set_ylabel(r'$v_c^{-3}$')
ax[0].legend()
ax[0].grid(True)

ax[1].plot(yD, K_y / K0, 'b-', label=r'Simulation: $\int v^2\,\mathrm{d}x$')
ax[1].axhline(1, color='k', linestyle='--', linewidth=1.5,
              label=r'Nozzle: $K = 2R_dU_0^2$')
ax[1].axhline(K_eff / K0, color='r', linewidth=2,
              label=r'$K_\mathrm{eff}/K = %.2f$' % (K_eff / K0))
ax[1].axvspan(yD_min, yD_max, color='gray', alpha=0.15)
ax[1].set_xlabel(r'$y/D$')
ax[1].set_ylabel(r'$K(y)/K$')
ax[1].set_ylim([0, 1.5])
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.savefig('centreline.png')
~~~
*/

/**
## See also

See also the axisymmetric case [here](/sandbox/MCVacher/schlichting_axi.c).

*/