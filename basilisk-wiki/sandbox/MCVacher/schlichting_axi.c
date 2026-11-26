/**
# Schlichting Jet (axisymmetric case)

We try to recover the famous problem of the laminar 2D jet in an immersed domain, also called the Schlichting jet. We look at the axisymmetric case. This analytical solution can be found in H. Schlichting's book "Boundary layer theory" (McGraw-Hill) in section XI.a.2. It was derived by H. Schlichting himself.

It is a monophasic problem depending only on the Reynolds number $Re$. Here we take $Re = 30$.

**Reminder** : *"For problems with a symmetry of revolution around the z-axis of a cylindrical coordinate system : The longitudinal coordinate (z-axis) is x and the radial coordinate (r-axis) is y"*.

Here the videos are rotated for convenience (using the most infamous manner I found ^^).
*/

/**
## Simulation
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "axistream.h"
#include "view.h"
#define BVIEW 1

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;

FILE * fpmax; 

face vector muv[];

int main() {

  Re=30;

  R_d=0.025; //Rayon du jet initial
  L0=1.28; //Taille de la boite (R_d = 5 cells with N = 256)
  U0=1;

  TOLERANCE = 1e-3 [*];

  // Nozzle (y < R_d): imposed top-hat jet. Elsewhere: slip wall
  // (an open boundary there lets a spurious through-flow develop,
  //  since a uniform stream would satisfy all the boundary conditions).
  u.n[left] = dirichlet (U0*(y < R_d));
  u.t[left] = y < R_d ? dirichlet (0.) : neumann (0.);
  s[left] = dirichlet (U0*(y < R_d));

  u.n[right] = neumann(0);
  p[right] = dirichlet(0.);

  u.n[top] = neumann(0);
  p[top] = dirichlet(0);

  N=256;
  origin (L0, 0);
  init_grid(N);

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
*/

/**
We plot the vertical velocity as well as a passive tracer to observe the establishment of the flow.
*/

//Ask Antoon what is going on, but it works !
#define vertex_value(b) ((b[] + b[-1] + b[0,-1] + b[-1,-1])/4.)

event movie_streamlines (t += 0.05; t <= 50) {
  scalar omega[];
  vertex scalar stream[];
  scalar psi[];

  psi[left] = dirichlet (0.);
  psi[right]    = neumann (0.);
  psi[top]   = neumann (0.);
  
  vorticity(u, omega);  
  axistream(u, psi);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = vertex_value(psi);
  boundary ({omega, phi});

  clear();
  view(fov=24, 
       width = L0/2,
       height = L0/2,
       tx     = -1.45*L0,      
       ty     = -0.4*L0);
  
  isoline("phi", n = 30);

  // Draw velocity vectors optionally
  // vect("u", scale=0.1);

  box();
  save("streamlines.mp4");
}

event ppm_output (t = 0; t += 0.05; t <= 50) {

    char name1[80];
    sprintf (name1, "uX.mp4");
    output_ppm (u.x, file = name1, n = 512, min = -U0, max = +U0, linear = true);

    char name2[80];
    sprintf (name2, "s.mp4");
    output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

event profile (t = end) {
  char name[80];
  sprintf (name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g \n", x, y, u.x[], u.y[]);
  fclose(fpres);
  
  printf ("-----END-----\n");
}


/**
<div style="display:flex; flex-direction:column; gap:60px; align-items:center;">

  <video src="schlichting_axi/s.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(90deg);">
  </video>

  <video src="schlichting_axi/uX.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(90deg);">
  </video>
  
  <video src="schlichting_axi/streamlines.mp4"
         controls
         loop
         muted
         playsinline
         style="height:400px; transform: rotate(90deg);">
  </video>

</div>

*/

/**
## Comparison with theory

Far enough from the nozzle, the jet is slender and the boundary layer
equations apply. In cylindrical coordinates, with $x$ the axial coordinate
(`x` in the code, measured from the nozzle), $r$ the radial coordinate
(`y` in the code), $u$ the axial and $v$ the radial velocity:

$$\frac{\partial (r u)}{\partial x} + \frac{\partial (r v)}{\partial r} = 0, \qquad
u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial r} = \frac{\nu}{r}\frac{\partial}{\partial r}\left(r\frac{\partial u}{\partial r}\right)$$

The pressure is uniform in the surrounding fluid, so the (kinematic) momentum
flux of the jet is conserved:

$$K = \frac{J}{\rho} = 2\pi\int_0^{\infty} u^2\, r\,\mathrm{d}r = \mathrm{cst}
\qquad (= \pi R_d^2 U_0^2~\mathrm{in~this~case})$$

### Scaling

Inertia is balanced by viscous forces, $u^2/x \sim \nu u/\delta^2$, and the
momentum flux is constant, $u^2\delta^2 \sim K$. Hence

$$u \sim \frac{K}{\nu\,x}, \qquad \delta \sim \frac{\nu\,x}{\sqrt{K}}.$$

Unlike the plane jet ($v \sim y^{-1/3}$, $\delta \sim y^{2/3}$), the round jet
spreads linearly and its local Reynolds number $u\delta/\nu \sim \sqrt{K}/\nu$
is constant along the jet.

### Self-similar solution

With the stream function ($u = \frac{1}{r}\frac{\partial\psi}{\partial r}$,
$v = -\frac{1}{r}\frac{\partial\psi}{\partial x}$), the scaling above suggests

$$\psi = \nu\,x\,F(\xi), \qquad \xi = \gamma\,\frac{r}{x},$$

with $\gamma$ a constant to be fixed by $K$. The momentum equation reduces to

$$F F' = F' - \xi F''.$$

Since $F' - \xi F'' = \frac{\mathrm{d}}{\mathrm{d}\xi}\left(2F - \xi F'\right)$,
it integrates once, with $F(0) = 0$, into $\xi F' - 2F + \frac{1}{2}F^2 = 0$, whose solution is

$$F(\xi) = \frac{\xi^2}{1 + \xi^2/4}.$$

The velocity field is then

$$u = \frac{2\gamma^2\nu}{x}\,\frac{1}{\left(1+\xi^2/4\right)^2}, \qquad
v = \frac{\gamma\,\nu}{x}\,\frac{\xi\left(1-\xi^2/4\right)}{\left(1+\xi^2/4\right)^2}.$$

The radial velocity changes sign at $\xi = 2$: the jet entrains the
surrounding fluid. The constant $\gamma$ follows from the momentum flux,
$K = \frac{16\pi}{3}\gamma^2\nu^2$, so that finally

$$\boxed{u(x,r) = \frac{3K}{8\pi\nu\,x}\,\frac{1}{\left(1+\xi^2/4\right)^2}, \qquad
\xi = \sqrt{\frac{3K}{16\pi}}\,\frac{1}{\nu}\,\frac{r}{x}}$$

### Effective momentum flux $K_\text{eff}$ and virtual origin $x_0$

Schlichting's solution knows the jet through **one number only**, its
momentum flux $K$, and assumes the jet comes from a **point source** at
$x = 0$ in an **unbounded fluid at rest**. A simulated jet satisfies none of
these exactly, so the comparison must be done with the values the jet
*actually* has far from the nozzle:

* **$K_\text{eff}$**, the momentum flux really carried in the self-similar region,
$$K_\text{eff} = \left\langle 2\pi\int_0^{\infty} u^2\, r\,\mathrm{d}r \right\rangle_{6 \le x/D \le 20};$$
* **$x_0$**, the virtual origin: the position of the point source that would
produce the same far field. Since $1/u_c = 8\pi\nu\,(x-x_0)/(3K)$, a linear fit
of $1/u_c(x)$ gives both $x_0$ (where the line crosses zero) and a second,
independent estimate of the momentum flux from its slope,
$K_\text{slope} = 8\pi\nu/(3\times\text{slope})$.

If the jet is a Schlichting jet, $K_\text{eff}$ and $K_\text{slope}$ must agree.
Then the theory is used with $x \to x - x_0$ and $K \to K_\text{eff}$:

$$u(x,r) = \frac{3K_\text{eff}}{8\pi\nu\,(x-x_0)}\,\frac{1}{\left(1+\xi^2/4\right)^2}, \qquad
\xi = \sqrt{\frac{3K_\text{eff}}{16\pi}}\,\frac{1}{\nu}\,\frac{r}{x-x_0}.$$

### Why can the jet lose momentum?

Take a control volume between the inlet plane and a section $x$, extending
radially to infinity. Its axial momentum balance reads

$$\underbrace{2\pi\int_0^\infty u^2\,r\,\mathrm{d}r}_{K(x)}
= K_\text{nozzle}
+ \underbrace{\frac{2\pi}{\rho}\int_\text{inlet plane} p\,r\,\mathrm{d}r}_{\text{pressure on the inlet plane}}
- \underbrace{\frac{2\pi}{\rho}\int_{\text{section }x} p\,r\,\mathrm{d}r}_{\text{pressure in the jet}}
+ \underbrace{(\ldots)}_{\text{viscous stresses}}
- \underbrace{(\ldots)}_{\text{momentum leaving through the outer boundaries}}$$

Schlichting's $K = \mathrm{cst}$ assumes that all the extra terms vanish. 

### What we check

1. **Measured $K_\text{eff}$ and $x_0$**, and the agreement between $K_\text{eff}$
(from $u^2$) and $K_\text{slope}$ (from the centreline decay).
2. **Self-similarity**: $(x-x_0)\,u$ depends only on $r/(x-x_0)$, and must
follow the theory computed with $K_\text{eff}$ (the theory with the nozzle
$K$ is shown for reference).
3. **Normalised profile**: with the centreline velocity $u_c$ and the
half-width $r_{1/2}$ (where $u = u_c/2$, i.e. $\xi_{1/2} = 2\sqrt{\sqrt{2}-1}$),
the profile has no free parameter at all, neither $K$ nor $x_0$:
$$\frac{u}{u_c} = \frac{1}{\left(1 + (\sqrt{2}-1)\,(r/r_{1/2})^2\right)^2}.$$

With $Re = 30$, $R_d = 0.025$, $U_0 = 1$: $\nu = U_0 R_d/Re = 8.3\cdot10^{-4}$,
$K = \pi R_d^2 U_0^2 = 1.96\cdot10^{-3}$, $\gamma \approx 13.0$, so
$u_c \approx 0.28/x$ and $r_{1/2} \approx 0.099\,x$: the jet is slender.
We keep the positions $6 \le x/D \le 20$ ($D = 2R_d$), away from the nozzle
and from the outflow boundary ($x/D = 25.6$).

~~~pythonplot Self-similarity check
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Simulation parameters
# ============================================================

Re = 30
U0 = 1
R_d = 0.025
D = 2 * R_d
nu = U0 * R_d / Re
K0 = np.pi * R_d**2 * U0**2                 # momentum flux of the nozzle

# ============================================================
# Load numerical data and build structured 2D array U[x, r]
# ============================================================

data = np.loadtxt("res_end.txt")
x_unique, ix = np.unique(data[:, 0], return_inverse=True)
r_unique, ir = np.unique(data[:, 1], return_inverse=True)

U = np.full((len(x_unique), len(r_unique)), np.nan)
U[ix, ir] = data[:, 2]                      # axial velocity u.x

Delta = x_unique[1] - x_unique[0]
xs = x_unique - (x_unique[0] - Delta / 2)   # distance from the nozzle

xD_min, xD_max = 6, 20
win = (xs / D >= xD_min) & (xs / D <= xD_max)

# ============================================================
# Effective momentum flux and virtual origin
# ============================================================

K_x = 2 * np.pi * np.sum(U**2 * r_unique, axis=1) * Delta
K_eff = np.mean(K_x[win])

uc = U[:, 0]                                # centreline (first cell)
slope, intercept = np.polyfit(xs[win], 1 / uc[win], 1)
x0 = -intercept / slope                     # 1/u_c = slope (x - x0)
K_slope = 8 * np.pi * nu / (3 * slope)

print("K_eff/K0 = %.3f, K_slope/K0 = %.3f, x0/D = %.2f"
      % (K_eff / K0, K_slope / K0, x0 / D))

# ============================================================
# Self-similar collapse
# ============================================================

eta_all, xu_all, r_r12_all, u_uc_all, xD_all = [], [], [], [], []

for i in np.where(win)[0]:
    u = U[i, :]
    X = xs[i] - x0                          # distance from virtual origin

    # half-width: first r where u < uc/2, linear interpolation
    j = np.argmax(u < uc[i] / 2)
    r12 = r_unique[j-1] + (uc[i]/2 - u[j-1]) * (r_unique[j] - r_unique[j-1]) / (u[j] - u[j-1])

    eta_all.extend(r_unique / X)
    xu_all.extend(X * u)
    r_r12_all.extend(r_unique / r12)
    u_uc_all.extend(u / uc[i])
    xD_all.extend([xs[i] / D] * len(r_unique))

# ============================================================
# Theory
# ============================================================

def xu_theory(eta, K):
    gamma = np.sqrt(3 * K / (16 * np.pi)) / nu
    return 3 * K / (8 * np.pi * nu) / (1 + gamma**2 * eta**2 / 4)**2

eta_th = np.linspace(0, 0.5, 500)
s_th = np.linspace(0, 4, 500)
u_uc_th = 1 / (1 + (np.sqrt(2) - 1) * s_th**2)**2

# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

sc = ax[0].scatter(eta_all, xu_all, c=xD_all, cmap='viridis',
                   vmin=xD_min, vmax=xD_max, s=5, alpha=0.6)
ax[0].plot(eta_th, xu_theory(eta_th, K_eff), 'r', linewidth=2,
           label=r'Theory with $K_\mathrm{eff}$')
ax[0].plot(eta_th, xu_theory(eta_th, K0), 'k--', linewidth=1.5,
           label=r'Theory with nozzle $K$')
ax[0].set_xlabel(r'$r/(x-x_0)$')
ax[0].set_ylabel(r'$(x-x_0)\,u(x,r)$')
ax[0].set_xlim([0, 0.4])
ax[0].legend()
ax[0].grid(True)
ax[0].text(0.97, 0.55,
           r'$K_\mathrm{eff}/K = %.2f$' '\n' r'$K_\mathrm{slope}/K = %.2f$' '\n' r'$x_0/D = %.1f$'
           % (K_eff / K0, K_slope / K0, x0 / D),
           transform=ax[0].transAxes, ha='right', va='top',
           bbox=dict(facecolor='white', alpha=0.8))

ax[1].scatter(r_r12_all, u_uc_all, c=xD_all, cmap='viridis',
              vmin=xD_min, vmax=xD_max, s=5, alpha=0.6)
ax[1].plot(s_th, u_uc_th, 'r', linewidth=2,
           label=r'$\left(1+(\sqrt{2}-1)(r/r_{1/2})^2\right)^{-2}$')
ax[1].set_xlabel(r'$r/r_{1/2}$')
ax[1].set_ylabel(r'$u/u_c$')
ax[1].set_xlim([0, 4])
ax[1].legend()
ax[1].grid(True)

cbar = fig.colorbar(sc, ax=ax.ravel().tolist())
cbar.set_label(r'Axial position $x/D$')

plt.savefig('self_sim.png')
~~~

~~~pythonplot Centreline decay and momentum flux
import numpy as np
import matplotlib.pyplot as plt

Re = 30
U0 = 1
R_d = 0.025
D = 2 * R_d
nu = U0 * R_d / Re
K0 = np.pi * R_d**2 * U0**2

data = np.loadtxt("res_end.txt")
x_unique, ix = np.unique(data[:, 0], return_inverse=True)
r_unique, ir = np.unique(data[:, 1], return_inverse=True)
U = np.full((len(x_unique), len(r_unique)), np.nan)
U[ix, ir] = data[:, 2]

Delta = x_unique[1] - x_unique[0]
xs = x_unique - (x_unique[0] - Delta / 2)

xD_min, xD_max = 6, 20
win = (xs / D >= xD_min) & (xs / D <= xD_max)

uc = U[:, 0]                                                # centreline velocity
K_x = 2 * np.pi * np.sum(U**2 * r_unique, axis=1) * Delta   # 2 pi int u^2 r dr
K_eff = np.mean(K_x[win])
slope, intercept = np.polyfit(xs[win], 1 / uc[win], 1)
x0 = -intercept / slope

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

ax[0].plot(xs / D, 1 / uc, 'b-', label='Simulation')
ax[0].plot(xs / D, 8 * np.pi * nu * xs / (3 * K0), 'k--', linewidth=1.5,
           label=r'Point source at nozzle, nozzle $K$')
ax[0].plot(xs / D, 8 * np.pi * nu * (xs - x0) / (3 * K_eff), 'r', linewidth=2,
           label=r'Point source at $x_0$, $K_\mathrm{eff}$')
ax[0].axvspan(xD_min, xD_max, color='gray', alpha=0.15, label='Comparison window')
ax[0].set_xlabel(r'$x/D$')
ax[0].set_ylabel(r'$1/u_c$')
ax[0].legend()
ax[0].grid(True)

ax[1].plot(xs / D, K_x / K0, 'b-', label=r'Simulation: $2\pi\int u^2 r\,\mathrm{d}r$')
ax[1].axhline(1, color='k', linestyle='--', linewidth=1.5,
              label=r'Nozzle: $K = \pi R_d^2 U_0^2$')
ax[1].axhline(K_eff / K0, color='r', linewidth=2,
              label=r'$K_\mathrm{eff}/K = %.2f$' % (K_eff / K0))
ax[1].axvspan(xD_min, xD_max, color='gray', alpha=0.15)
ax[1].set_xlabel(r'$x/D$')
ax[1].set_ylabel(r'$K(x)/K$')
ax[1].set_ylim([0, 1.5])
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.savefig('centreline.png')
~~~

Outside the nozzle, the inlet plane is a slip wall: the entrained fluid is
sucked in radially along it, the pressure on the wall drops below ambient
and the pressure term on the inlet plane pulls the jet backwards. Hence
$K_\text{eff} < K$ is expected; $K(x)$ should however be flat in the window. If $K_\text{eff}$ and
$K_\text{slope}$ disagree, the jet is not yet self-similar in the window (or
the run is not converged).
*/