/**
# Natural convection along a heated vertical plate

Numerical simulation of the natural convection boundary layer
along a vertical isothermal plate ($x = 0$), maintained at $T_0 > T_\infty$.

## Non-dimensionalisation
The Navier-Stokes equations are non-dimensionalised by the thermal diffusivity
$\alpha$ and the domain length $L_0$:

$$x^* = \frac{x}{L_0}, \quad t^* = \frac{\alpha \, t}{L_0^2}, \quad
\mathbf{u}^* = \frac{L_0}{\alpha}\,\mathbf{u}, \quad
\theta = \frac{T - T_\infty}{T_0 - T_\infty}$$

The non-dimensional momentum equation reads:

$$\frac{\partial \mathbf{u}^*}{\partial t^*} + (\mathbf{u}^* \cdot \nabla)\mathbf{u}^* = -\nabla p^* + Pr \, \nabla^2 \mathbf{u}^* + Ra \cdot Pr \; \theta \, \mathbf{e}_y$$

and the energy equation:

$$\frac{\partial \theta}{\partial t^*} + \mathbf{u}^* \cdot \nabla \theta = \nabla^2 \theta$$

with $Ra = \dfrac{g \beta (T_0 - T_\infty) L_0^3}{\nu \alpha}$ and $Pr = \dfrac{\nu}{\alpha}$.
The buoyancy term `Ra * Pr * θ` appears directly in the `acceleration` event.

## Simulation parameters
 *   Ra = 1e5    (Rayleigh number based on $L_0$)
 *   Pr = 0.71   (air)

With $Ra = 10^5$, the thermal boundary layer thickness at the top of the
domain is $\delta \approx 6\,Ra^{-1/4} \approx 0.34\,L_0$: the boundary
layer is thin enough to be compared with the self-similar solution
(see below), and it is resolved by about 90 cells with $N = 256$.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#define Ra 1e5
#define Pr 0.71

double tmax     = 0.5;
double t_output = 0.005;

scalar T[];
scalar * tracers = {T};

mgstats mgT;

face vector av[], muc[], D[];

int main()
{
  L0 = 1.;
  origin (0., 0.);

  N = 256;

  DT = 1e-4;
  TOLERANCE = 1e-4;
  NITERMIN = 4;

  a  = av;
  mu = muc;

  /**
  ## Boundary conditions

  The left boundary is the hot no-slip wall. The other boundaries are open:
  the pressure is set to zero there, so that the fluid can be entrained
  from the right and bottom boundaries and leave at the top. At the top, the
  outflow condition is switched to a no-flow condition wherever the fluid
  would re-enter the domain. */

  // Left: hot wall
  u.n[left] = dirichlet(0.);
  u.t[left] = dirichlet(0.);

  // Right: open boundary (entrainment)
  u.n[right]  = neumann(0.);
  u.t[right]  = neumann(0.);
  p[right]    = dirichlet(0.);
  pf[right]   = dirichlet(0.);

  // Bottom: open boundary (entrainment)
  u.n[bottom] = neumann(0.);
  u.t[bottom] = neumann(0.);
  p[bottom]   = dirichlet(0.);
  pf[bottom]  = dirichlet(0.);

  // Top: outflow, no backflow
  u.n[top] = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
  u.t[top] = neumann(0.);
  p[top]   = dirichlet(0.);
  pf[top]  = dirichlet(0.);

  // Temperature
  T[left]   = dirichlet(1.);
  T[right]  = dirichlet(0.);
  T[bottom] = dirichlet(0.);
  T[top]    = neumann(0.);

  run();
}

event init (i = 0)
{
  foreach() {
    T[] = 0.;
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

/**
Non-dimensional viscosity ($Pr$) and thermal diffusivity (1). */

event properties (i++)
{
  foreach_face() {
    muc.x[] = fm.x[] * Pr;
    D.x[]   = fm.x[];
  }
}


/**
Thermal diffusion
*/

event tracer_diffusion (i++)
{
  mgT = diffusion (T, dt, D);
}


/**
Buoyancy force
*/

event acceleration (i++)
{
  face vector av = a;

  foreach_face(y)
    av.y[] += Ra * Pr * (T[] + T[0,-1]) / 2.;
}

event logfile (i++)
{
  double umax = 0.;

  foreach(reduction(max:umax))
    umax = max(umax, norm(u));

  fprintf(stderr,
          "i=%d t=%g umax=%g\n",
          i, t, umax);
}


/**
Save final solution
*/

event res_save (t = end)
{
  FILE * fp = fopen("res_final.txt", "w");

  foreach()
    fprintf(fp, "%g %g %g %g %g %g\n",
            x, y, u.x[], u.y[], p[], T[]);

  fclose(fp);
}


/**
Movies
*/

event movie (t = 0; t += t_output; t <= tmax)
{
  output_ppm(T,
             file = "T.mp4",
             n = 512,
             min = 0.,
             max = 1.,
             linear = true,
             map = jet);

  output_ppm(u.x,
             file = "ux.mp4",
             n = 512,
             spread = -1,
             linear = true,
             map = cool_warm);

  output_ppm(u.y,
             file = "uy.mp4",
             n = 512,
             spread = -1,
             linear = true,
             map = cool_warm);
}


event end (t = tmax)
{
  printf("-----END-----\n");
}

/**
# Outputs
![Temperature](heated_plate/T.mp4)
![Vertical velocity](heated_plate/uy.mp4)
![Horizontal velocity](heated_plate/ux.mp4)
*/

/**
# Theory: self-similar solution

The non-dimensional temperature is $\theta = \dfrac{T - T_\infty}{T_0 - T_\infty} \in [0,1]$.

The similarity variable is built from the thermal boundary layer scale
$\delta_T \sim y \, Ra_y^{-1/4}$:

$$\eta = \frac{x}{\delta_T} = \frac{x}{y} Ra_y^{1/4}, \qquad Ra_y = \frac{g \beta \Delta T \, y^3}{\nu \alpha}$$

The stream function is introduced as:

$$\psi = \alpha \, Ra_y^{1/4} \, F(\eta), \qquad v = \frac{\partial \psi}{\partial x} = \frac{\alpha}{y}\,Ra_y^{1/2}\,F'(\eta)$$

Non-dimensional equations (valid for any $Pr$):

$$\boxed{\begin{cases} F''' + \dfrac{1}{Pr}\!\left(\dfrac{3}{4} F F'' - \dfrac{1}{2} F'^2\right) + \theta = 0 \\[8pt] \theta'' + \dfrac{3}{4} F \theta' = 0 \end{cases}}$$

Boundary conditions:

$$\eta = 0 : \quad F = 0, \quad F' = 0, \quad \theta = 1$$
$$\eta \to \infty : \quad F' \to 0, \quad \theta \to 0$$

This boundary value problem is solved by a shooting method on the two
unknowns $F''(0)$ and $\theta'(0)$. At $Pr = 0.71$:
$F''(0) = 0.879$ and $\theta'(0) = -0.387$.

## Self-similar variables in simulation units

With the non-dimensionalisation above, $Ra_y = Ra\,y^{*3}$, so that

$$\eta = x^*\,Ra^{1/4}\,y^{*-1/4}, \qquad
\theta = \theta(\eta), \qquad
\frac{v^*}{Ra^{1/2}\,y^{*1/2}} = F'(\eta).$$

If the boundary layer is self-similar, all the profiles taken at different
heights $y^*$ must collapse on a single curve when plotted against $\eta$.
We keep only the heights $0.3 \le y^* \le 0.9$, to stay away from the
leading edge (bottom) and from the outflow boundary (top).

The local Nusselt number gives a second check:

$$Nu_y = \frac{h_y\,y}{\lambda} = -\theta'(0)\,Ra_y^{1/4}
\quad\Rightarrow\quad
\frac{Nu_y}{Ra_y^{1/4}} = -\left.\frac{\partial \theta}{\partial x^*}\right|_{x^*=0} \frac{y^{*1/4}}{Ra^{1/4}} = -\theta'(0) = 0.387.$$

The wall gradient is computed with a second-order one-sided formula using the
wall value $\theta_w = 1$ and the first two cells:
$\partial_x \theta|_0 \approx (-8\,\theta_w + 9\,\theta_1 - \theta_2)/(3\Delta)$.

~~~pythonplot Self-similarity check
import numpy as np
import matplotlib.pyplot as plt

Ra = 1e5
Pr = 0.71

# ============================================================
# Self-similar solution (shooting method, RK4)
# ============================================================

eta_max, h = 10., 0.01
eta_th = np.arange(0., eta_max + h/2, h)

def rhs(Y):
    F, Fp, Fpp, th, thp = Y
    return np.array([Fp, Fpp,
                     -th - (0.75*F*Fpp - 0.5*Fp**2)/Pr,
                     thp,
                     -0.75*F*thp])

def integrate(s):
    Y = np.array([0., 0., s[0], 1., s[1]])
    sol = [Y]
    for _ in eta_th[1:]:
        k1 = rhs(Y); k2 = rhs(Y + h/2*k1)
        k3 = rhs(Y + h/2*k2); k4 = rhs(Y + h*k3)
        Y = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
        sol.append(Y)
    return np.array(sol)

def miss(s):
    # conditions at infinity: F'(eta_max) = 0, theta(eta_max) = 0
    Y = integrate(s)[-1]
    return np.array([Y[1], Y[3]])

s = np.array([0.9, -0.4])          # initial guess for F''(0), theta'(0)
for _ in range(20):                # Newton iterations
    r = miss(s)
    if np.max(np.abs(r)) < 1e-10:
        break
    J = np.zeros((2, 2))
    for k in range(2):
        ds = np.zeros(2); ds[k] = 1e-6
        J[:, k] = (miss(s + ds) - r)/1e-6
    s = s - np.linalg.solve(J, r)

sol = integrate(s)
Fp_th, theta_th = sol[:, 1], sol[:, 3]
dtheta0 = s[1]
print("F''(0) = %.4f, theta'(0) = %.4f" % (s[0], s[1]))

# ============================================================
# Load numerical data and build structured 2D arrays [y, x]
# ============================================================

data = np.loadtxt("res_final.txt")
x_unique, ix = np.unique(data[:, 0], return_inverse=True)
y_unique, iy = np.unique(data[:, 1], return_inverse=True)

V = np.full((len(y_unique), len(x_unique)), np.nan)
T = np.full((len(y_unique), len(x_unique)), np.nan)
V[iy, ix] = data[:, 3]
T[iy, ix] = data[:, 5]

# ============================================================
# Self-similar collapse, away from leading edge and outflow
# ============================================================

y_min, y_max = 0.3, 0.9
rows = np.where((y_unique >= y_min) & (y_unique <= y_max))[0]

eta_all, theta_all, fp_all, y_all = [], [], [], []
for i in rows:
    yi = y_unique[i]
    eta_all.extend(x_unique * Ra**0.25 * yi**(-0.25))
    theta_all.extend(T[i, :])
    fp_all.extend(V[i, :] / (Ra**0.5 * yi**0.5))
    y_all.extend([yi]*len(x_unique))

eta_all = np.array(eta_all)
y_all = np.array(y_all)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

sc = ax[0].scatter(eta_all, theta_all, c=y_all, cmap='viridis', s=3, alpha=0.6)
ax[0].plot(eta_th, theta_th, 'r', linewidth=2, label='Theory')
ax[0].set_xlabel(r'$\eta = x^*\,Ra^{1/4}\,y^{*-1/4}$')
ax[0].set_ylabel(r'$\theta$')
ax[0].set_xlim([0, 8])
ax[0].legend()
ax[0].grid(True)

ax[1].scatter(eta_all, fp_all, c=y_all, cmap='viridis', s=3, alpha=0.6)
ax[1].plot(eta_th, Fp_th, 'r', linewidth=2, label='Theory')
ax[1].set_xlabel(r'$\eta = x^*\,Ra^{1/4}\,y^{*-1/4}$')
ax[1].set_ylabel(r"$F'(\eta) = v^* / (Ra^{1/2}\,y^{*1/2})$")
ax[1].set_xlim([0, 8])
ax[1].legend()
ax[1].grid(True)

cbar = fig.colorbar(sc, ax=ax.ravel().tolist())
cbar.set_label(r'Vertical position $y^*$')

plt.savefig('self_sim.png')
~~~

~~~pythonplot Local Nusselt number
import numpy as np
import matplotlib.pyplot as plt

Ra = 1e5
theta_prime_0 = -0.387   # self-similar solution at Pr = 0.71

data = np.loadtxt("res_final.txt")
x_unique, ix = np.unique(data[:, 0], return_inverse=True)
y_unique, iy = np.unique(data[:, 1], return_inverse=True)
T = np.full((len(y_unique), len(x_unique)), np.nan)
T[iy, ix] = data[:, 5]

Delta = x_unique[1] - x_unique[0]

# second-order one-sided wall gradient (theta_w = 1)
dTdx = (-8.*1. + 9.*T[:, 0] - T[:, 1]) / (3.*Delta)
nu = -dTdx * y_unique**0.25 / Ra**0.25      # Nu_y / Ra_y^(1/4)

plt.figure(figsize=(6, 4.5))
plt.plot(y_unique, nu, 'b-', label='Simulation')
plt.axhline(-theta_prime_0, color='r', linestyle='--', linewidth=2,
            label=r"Theory: $-\theta'(0) = 0.387$")
plt.axvspan(0.3, 0.9, color='gray', alpha=0.15,
            label='Range used for the collapse')
plt.xlabel(r'$y^*$')
plt.ylabel(r'$Nu_y / Ra_y^{1/4}$')
plt.ylim([0, 1])
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('nusselt.png')
~~~
*/