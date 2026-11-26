/**
# Blasius boundary layer

Semi-infinite flat plate along $y = 0$ with no-slip at $x > X_0$.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#define U_INF  1.0
#define RE     500.0        
#define NU     (U_INF/RE)
#define X0     0.2          // leading edge of the plate

#define LEVEL  9            
#define T_END  3.0          // a few L/U crossing times to reach steady state

double t_output = 0.01;

// Left inlet: uniform profile
u.n[left]   = dirichlet(U_INF);
u.t[left]   = dirichlet(0.);

// Right outlet: free outflow
u.n[right]  = neumann(0.);
u.t[right]  = neumann(0.);
p[right]    = dirichlet(0.);

// Bottom wall: slip before X0, no-slip after
u.t[bottom] = x > X0 ? dirichlet(0.) : neumann(0.);
u.n[bottom] = dirichlet(0.);

// Top boundary: open, so the fluid displaced by the boundary layer can leave
u.n[top]    = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
u.t[top]    = neumann(0.);
p[top]      = dirichlet(0.);

int main() {
  size(1.0);
  origin(0., 0.);
  init_grid(1 << LEVEL); //2^LEVEL

  const face vector muv[] = {NU, NU};
  mu = muv;

  run();
}

event init (t = 0) {
  foreach() {
    u.x[] = U_INF;
    u.y[] = 0.;
  }
}

event logfile (i += 10) {
  fprintf(stderr, "%d %g %g\n", i, t, dt);
}

event movie (t = 0; t += t_output; t <= T_END) {
  output_ppm(u.x, file = "ux.mp4", n = 512,
             min = 0., max = U_INF, linear = true, map = cool_warm);
  output_ppm(u.y, file = "uy.mp4", n = 512,
             spread = -1, linear = true, map = cool_warm);
}


event profile (t = end) {
  FILE * fp = fopen("res_end.txt", "w");
  foreach()
    fprintf(fp, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  fclose(fp);
}

event end (t = T_END) {}

/**
# Video results

![Streamwise velocity u](blasius/ux.mp4)
![Vertical velocity v](blasius/uy.mp4)

## Self-similarity check

We extract cell columns at several distances $x - X_0$ from the leading
edge. We avoid the immediate vicinity of the leading edge, where the solution
is singular, and the outlet.

**Normalisation by the local outer velocity.** The outer flow is not exactly
$U_\infty$: the fluid displaced by the boundary layer slightly accelerates it.
Blasius theory must be compared with the velocity actually seen at the edge of
the layer. The flow above the layer is not uniform in $y$
either: $u$ reaches a maximum just above the layer, then decreases slowly
towards the top of the domain. The displacement thickness acts like a growing
bump on the wall, and the flow over a bump is fastest close to its surface.
At leading order (Prandtl model) it vanishes and $\partial_y p = 0$. For each column we therefore define the edge velocity as the maximum of the profile, $U_e(x) = \max_{y < y_e} u(x, y)$ with $y_e = 0.5$. We then use:

$$
\delta(x) = \sqrt{\frac{\nu (x-X_0)}{U_e}}, \qquad
\eta = \frac{y}{\delta(x)}, \qquad
Re_x = \frac{U_e (x-X_0)}{\nu}
$$

and plot $u/U_e$ and $v\sqrt{Re_x}/U_e$ as functions of $\eta$, to compare with
$f'(\eta)$ and $\frac{1}{2}(\eta f' - f)$.

**Pressure check (Blasius or Falkner–Skan?).** Outside the layer, the Prandtl
model gives $-\frac{1}{\rho}\frac{\mathrm{d}p}{\mathrm{d}x} = U_e\frac{\mathrm{d}U_e}{\mathrm{d}x}$.
Blasius requires $\mathrm{d}p/\mathrm{d}x = 0$, i.e. a constant $U_e$. If instead
$U_e \propto (x-X_0)^m$, the self-similar solutions are the Falkner–Skan
profiles, with $\beta = 2m/(m+1)$:

* $m > 0$: accelerating flow, favourable pressure gradient, which gives a
thinner layer and a steeper wall slope than Blasius;
* $m < 0$: decelerating flow, adverse pressure gradient, which gives a thicker
layer.

The last plot shows $p(x)$ and $u(x, y_e)$ along $y = y_e$, together with
$p + \rho u^2/2$, which should be constant (Bernoulli). The exponent $m$ is
fitted by least squares on $\ln U_e$ against $\ln(x - X_0)$ over the
measurement stations (grey band), with $U_e$ the edge velocity defined above.
Because $U_e(x)$ is not exactly a power law, $m$ is only an average over the
stations. Blasius is recovered if $|m| \ll 1$. Note
that a virtual origin for the plate slightly biases this fit.

~~~pythonplot Self-similarity
import numpy as np
import matplotlib.pyplot as plt

# Parameters (must match the #define in the C code)
U_INF, RE, X0 = 1.0, 500.0, 0.2
NU = U_INF/RE
Y_E = 0.5                              # height where the outer velocity is taken

# ---------- Blasius solution: 2f''' + f f'' = 0 (RK4) ----------
def blasius(fpp0=0.332057, eta_max=10., n=4001):
    eta = np.linspace(0., eta_max, n)
    h = eta[1] - eta[0]
    rhs = lambda Y: np.array([Y[1], Y[2], -0.5*Y[0]*Y[2]])
    F = np.zeros((n, 3))
    F[0] = [0., 0., fpp0]
    for k in range(n - 1):
        Y = F[k]
        k1 = rhs(Y); k2 = rhs(Y + h/2*k1)
        k3 = rhs(Y + h/2*k2); k4 = rhs(Y + h*k3)
        F[k+1] = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return eta, F[:, 0], F[:, 1]

eta_th, f_th, fp_th = blasius()
v_th = 0.5*(eta_th*fp_th - f_th)

# ---------- Numerical data ----------
data = np.loadtxt("res_end.txt")
x, y, ux, uy, p = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
x_cells = np.unique(x)
y_cells = np.unique(y)

stations = X0 + np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])   # x - X0
colors = plt.cm.viridis(np.linspace(0, 0.9, len(stations)))

fig, ax = plt.subplots(2, 2, figsize=(11, 9))
ax = ax.ravel()

Ue_st, xl_st = [], []
for xs, c in zip(stations, colors):
    xc = x_cells[np.argmin(np.abs(x_cells - xs))]   # closest cell column
    m = (x == xc)
    order = np.argsort(y[m])
    yy, uu, vv = y[m][order], ux[m][order], uy[m][order]

    xl = xc - X0                                    # distance from the leading edge
    Ue = uu[yy < Y_E].max()                         # edge velocity = max of the profile
    delta = np.sqrt(NU*xl/Ue)
    Re_x = Ue*xl/NU
    eta = yy/delta
    lab = r'$x-X_0 = %.2f$' % xl
    Ue_st.append(Ue); xl_st.append(xl)

    ax[0].plot(uu/U_INF, yy, '.', ms=3, color=c, label=lab)
    ax[1].plot(eta, uu/Ue, '.', ms=3, color=c, label=lab)
    ax[2].plot(eta, vv/Ue*np.sqrt(Re_x), '.', ms=3, color=c, label=lab)

ax[0].set_xlabel(r'$u/U_\infty$'); ax[0].set_ylabel(r'$y$')
ax[0].set_title('Raw profiles')

ax[1].plot(eta_th, fp_th, 'r-', lw=2, label=r"Blasius $f'(\eta)$")
ax[1].set_xlabel(r'$\eta = y/\delta(x)$'); ax[1].set_ylabel(r'$u/U_e$')
ax[1].set_xlim(0, 8)
ax[1].set_title('Streamwise velocity')

ax[2].plot(eta_th, v_th, 'r-', lw=2, label=r"Blasius $\frac{1}{2}(\eta f'-f)$")
ax[2].set_xlabel(r'$\eta = y/\delta(x)$'); ax[2].set_ylabel(r'$v\sqrt{Re_x}/U_e$')
ax[2].set_xlim(0, 8)
ax[2].set_ylim(-0.05, 1.2)
ax[2].set_title('Vertical velocity')

# ---------- Outer flow: pressure and U_e along y = Y_E ----------
yr = y_cells[np.argmin(np.abs(y_cells - Y_E))]
m = (y == yr)
order = np.argsort(x[m])
xr, ur, pr = x[m][order], ux[m][order], p[m][order]
bern = pr + 0.5*ur**2

# Falkner-Skan exponent: U_e ~ (x - X0)^m fitted over the stations
m_fs = np.polyfit(np.log(xl_st), np.log(Ue_st), 1)[0]
beta = 2*m_fs/(m_fs + 1)

ax[3].plot(xr, pr, 'k-', lw=2, label=r'$p$')
ax[3].plot(xr, bern - bern[0], 'k:', lw=1.5, label=r'$p + u^2/2$ (shifted)')
ax[3].axvline(X0, color='grey', ls='--', lw=1)
ax[3].axvspan(stations[0], stations[-1], color='grey', alpha=0.15)
ax[3].set_xlabel(r'$x$'); ax[3].set_ylabel(r'$p$')
ax3b = ax[3].twinx()
ax3b.plot(xr, ur, 'b-', lw=2, label=r'$u(x, y_e)$')
ax3b.set_ylabel(r'$u(x, y_e)$', color='b')
ax3b.tick_params(axis='y', colors='b')
ax[3].set_title(r'Outer flow at $y = %.2f$:  $m = %.3f$,  $\beta = %.3f$'
                % (Y_E, m_fs, beta))
h1, l1 = ax[3].get_legend_handles_labels()
h2, l2 = ax3b.get_legend_handles_labels()
ax[3].legend(h1 + h2, l1 + l2, fontsize=8, loc='best')

for a in ax:
    a.grid(True, alpha=0.4)
ax[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig('self_sim.png')
~~~
*/

/**
## Theory

An order-of-magnitude analysis of the Navier–Stokes equations (Prandtl
approximation, $\delta \ll x$) gives $\partial_y p = 0$. Since the external
flow is uniform, we also have $\partial_x p = 0$. This leads to the Prandtl
model:

$$
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0, \qquad
u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}
= \nu\frac{\partial^{2}u}{\partial y^{2}}
$$

with $u = v = 0$ at $y = 0$ and $u \to U_\infty$ as $y \to \infty$.

The inertia/viscosity balance $U_\infty^2/x \sim \nu U_\infty/\delta^2$ sets
the boundary layer thickness and the similarity variable:

$$
\delta(x) = \sqrt{\frac{\nu x}{U_\infty}}, \qquad \eta = \frac{y}{\delta(x)}
$$

We look for a self-similar stream function
$\psi = U_\infty\,\delta(x)\,f(\eta)$, with $u = \partial_y\psi$ and
$v = -\partial_x\psi$. This gives the **Blasius equation**:

$$
2f''' + f f'' = 0, \qquad f(0) = f'(0) = 0, \qquad f'(\infty) = 1
$$

It has no analytical solution. It is integrated numerically (RK4) from
$\eta = 0$, using the known wall value $f''(0) \simeq 0.332$ that satisfies
$f'(\infty) = 1$. The velocity profiles read:

$$
\frac{u}{U_\infty} = f'(\eta), \qquad
\frac{v}{U_\infty}\sqrt{Re_x} = \frac{1}{2}\left(\eta f' - f\right),
\qquad Re_x = \frac{U_\infty x}{\nu}
$$
*/