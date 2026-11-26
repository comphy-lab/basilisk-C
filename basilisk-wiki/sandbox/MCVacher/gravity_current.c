/**
# Gravity current

An analog of the [splouch](/sandbox/M1EMN/Exemples/column_viscous.c), here with twophase.h.

## Simulation
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"

// I know understood that embed.h + twophase.h + "gravity" => need for reduced.h to impose gravity.


double h_1;
double l_1;
double mu_b;
double rho_b;
double R;

// Height of the solid ceiling
double ceiling;

// Domain extents used for visualization
double xmin, ymin;
double xmax, ymax;

FILE * fpmax; //

int main() {

  L0=4; //size of the box
  rho2 = 1000;
  rho1=rho2*1.011;
  mu1 = 0.001;
  mu2 = mu1;
  h_1=0.3;
  l_1=0.5;
  ceiling=1.5*h_1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.n[embed] = dirichlet(0.);
  u.t[embed] = neumann(0.);

  u.n[left] = dirichlet(0.);
  u.t[left] = neumann(0.);

  u.n[right] = dirichlet(0.);
  u.t[right] = neumann(0.);

  G.y = -9.81;

  N=1024;
  origin (0, 0);
  init_grid(N);

  fpmax =  fopen("log.dat", "w");

  run();
}

event init (t = 0) {
  mu_b=mu2/mu1;
  rho_b=rho2/rho1;
  R=h_1/l_1;

  xmin = 0;
  xmax = L0;
  ymin = 0.0;
  ymax = ceiling;

  char dim_adim[80];
  sprintf (dim_adim, "dim_adim.txt"); //to compute the non-dim. parameters
  FILE * fparam = fopen(dim_adim, "w");
  fprintf (fparam, "%g %g %g\n",mu_b,rho_b,R);
  fclose (fparam);

  fraction (f, intersection((l_1-x),(h_1-y)));
  solid (cs, fs, y<ceiling);
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 0.01; t <= 20){

  char name[80];

  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);

  isave1++;
}

/**
We generate a video of the fraction ...
*/

event ppm_output (t = 0; t += 0.01; t <= 20) {

  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 2048, min = 0, max = 1, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  //u.y
  char name_u[80];
  sprintf (name_u, "u_y.mp4");
  output_ppm (u.y, file = name_u, n = 2048,min = -0.3, max = 0.3, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

  //u.y
  char name_u_2[80];
  sprintf (name_u_2, "u_x.mp4");
  output_ppm (u.x, file = name_u_2, n = 2048,min = -0.3, max = 0.3, linear = true, box = {{xmin, ymin}, {xmax, ymax}});

}

/**
![Gravity current](gravity_current/f.mp4)
![Vertical velocity](gravity_current/u_y.mp4)
![Horizontal velocity](gravity_current/u_x.mp4)
*/

/**
## Comparison with theory : Huppert's box model

We follow Huppert & Simpson (*J. Fluid Mech.* 99, 1980) and Huppert (*J. Fluid Mech.* 121, 1982), and only use the position $x_f(t)$ of the front.

**Slumping phase.** Before the box-model regime, Huppert & Simpson (1980) observe an initial *slumping* phase, lasting a few lock lengths, during which the front moves at a constant speed. The head then has roughly half the lock depth, $h_s \simeq h_1/2$, so that the front condition gives

$$U_s \simeq Fr(\phi_s)\,\sqrt{g' h_s}, \qquad \phi_s = \frac{h_1}{2H} = \frac{1}{3} \quad\Rightarrow\quad U_s \approx 0.092~\mathrm{m/s}.$$

With $x_f/l_1 \le 5$ at $t = 20$ s, the simulation is expected to stay mostly in this phase.

**Box model.** The current is a rectangle of length $x_f$ and height $h$. Volume conservation (2D, per unit width) gives

$$h\,x_f = A = h_1 l_1,$$

and the front velocity is set by the front condition

$$U = \frac{\mathrm{d}x_f}{\mathrm{d}t} = Fr\,\sqrt{g' h}, \qquad g' = g\,\frac{\rho_1-\rho_2}{\rho_2},$$

with the empirical Froude number of Huppert & Simpson, which depends on the depth ratio $\phi = h/H$ ($H$ = depth of the channel, here the ceiling):

$$Fr = \frac{1}{2}\,\phi^{-1/3} \quad (0.075 \le \phi \le 1), \qquad Fr = 1.19 \quad (\phi < 0.075).$$

This gives an ODE for $x_f(t)$, integrated numerically from $x_f(0) = l_1$.

**Inertial self-similar regime.** For a thin current ($\phi < 0.075$), $Fr$ is constant and $\mathrm{d}x_f/\mathrm{d}t = Fr\sqrt{g'A/x_f}$ integrates into

$$x_f = \left(\frac{3Fr}{2}\right)^{2/3}(g'A)^{1/3}\,t^{2/3} \simeq 1.47\,(g'A)^{1/3}\,t^{2/3}, \qquad U \propto t^{-1/3}.$$

**Viscous self-similar regime** (Huppert 1982), when viscosity balances buoyancy:

$$x_f = 1.411\left(\frac{g'A^3}{\nu}\right)^{1/5} t^{1/5}, \qquad U \propto t^{-4/5}.$$

Equating both laws, the viscous regime starts at $t_v \sim (A^4/g'^2\nu^3)^{1/7}$.

**Measurements.** The front $x_f$ is the last position where the current is thicker than $h_\text{min} = 1$ cm (a thin film of heavy fluid can run ahead of the head along the free-slip bottom, and is ignored). The head height $h_\text{head}$ is the maximum height of the interface within 0.3 m (about two head heights) behind the front, smoothed by a running median over 3 s (about one billow period) to filter the Kelvin–Helmholtz billows passing behind the head; the first instants, when this window still overlaps the collapsing lock ($x_f < l_1 + 0.3$ m), are discarded. The front velocity is a linear fit of $x_f(t)$ over a sliding window of 1 s. The front condition is then checked with the measured head: $Fr = U/\sqrt{g' h_\text{head}}$ against $h_\text{head}/H$. Note that the box-model height $A/x_f$ is not the head height during slumping.

~~~pythonplot Front position and velocity
import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import re

# ============================================================
# Parameters (same as in main())
# ============================================================

L0 = 4.
N = 1024
h_1 = 0.3
l_1 = 0.5
H = 1.5 * h_1                 # ceiling
drho = 0.011                  # (rho1 - rho2)/rho2, as in main()
gp = 9.81 * drho              # g' = g (rho1 - rho2)/rho2
A = h_1 * l_1
dt_save = 0.01
t_end = 20.                   # end of the res_save event

h_min = 0.01                  # minimal thickness of the current at the front
head_len = 0.3                # length behind the front used for the head height

def Fr(phi):                  # Huppert & Simpson (1980)
    return np.where(phi < 0.075, 1.19, 0.5 * phi**(-1/3))

def front(data):
    x, y = data[:, 0], data[:, 1]
    dx = L0 / N
    bins = np.floor(x / dx).astype(int)
    hmax = np.zeros(bins.max() + 1)
    np.maximum.at(hmax, bins, y)                  # interface height per column
    xf = (np.where(hmax > h_min)[0].max() + 1) * dx
    hh = y[(x > xf - head_len) & (x <= xf)].max()  # head height
    return xf, hh

# ============================================================
# Front position and head height
# ============================================================

files = gb.glob('interface-*.txt')
files.sort(key=lambda f: int(re.search(r'\d+', f).group()))

t_list, xf_list = [], []
for name in files:
    tn = (int(re.search(r'\d+', name).group()) - 1) * dt_save
    if tn > t_end + 1e-9:
        continue
    data = np.loadtxt(name, ndmin=2)
    if data.size == 0:
        continue
    xf, hh = front(data)
    t_list.append(tn)
    xf_list.append(xf)
t = np.array(t_list[:-1])       # the last file (t = t_end) is discarded
xf = np.array(xf_list[:-1])

keep = xf < L0 - 0.3          # stop before the end wall
t, xf = t[keep], xf[keep]

# front velocity: sliding linear fit over 1 s
w = int(round(0.5 / dt_save))
t_U, U = [], []
for i in range(w, len(t) - w):
    U.append(np.polyfit(t[i-w:i+w+1], xf[i-w:i+w+1], 1)[0])
    t_U.append(t[i])
t_U, U = np.array(t_U), np.array(U)

# ============================================================
# Theory
# ============================================================

# slumping phase: constant speed, head of height h_1/2
U_s = Fr(h_1 / (2 * H)) * np.sqrt(gp * h_1 / 2)

# box model: dx/dt = Fr(h/H) sqrt(g' h), h = A/x
dt = 1e-3
t_box = np.arange(0, t[-1] + dt, dt)
x_box = np.zeros_like(t_box)
x_box[0] = l_1
for n in range(len(t_box) - 1):
    hb = A / x_box[n]
    x_box[n+1] = x_box[n] + dt * Fr(hb / H) * np.sqrt(gp * hb)
U_box = Fr(A / (x_box * H)) * np.sqrt(gp * A / x_box)

# inertial self-similar asymptote (Fr = 1.19)
x_inert = 1.47 * (gp * A)**(1/3) * t_box**(2/3)
U_inert = 2/3 * 1.47 * (gp * A)**(1/3) * t_box[1:]**(-1/3)

print("U_s = %.4f m/s, mean measured U for t > 2 s = %.4f m/s"
      % (U_s, np.mean(U[t_U > 2])))

# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

ax[0].plot(t, xf, 'b-', label='Simulation')
ax[0].plot(t_box, l_1 + U_s * t_box, 'g-.', linewidth=2, label=r'Slumping: $U_s$')
ax[0].plot(t_box, x_box, 'r--', linewidth=2, label=r'Box model, $Fr(h/H)$')
ax[0].plot(t_box, x_inert, 'k:', linewidth=1.5,
           label=r'$1.47\,(g^\prime A)^{1/3}\,t^{2/3}$')
ax[0].set_xlabel(r'$t$ (s)')
ax[0].set_ylabel(r'$x_f$ (m)')
ax[0].legend()
ax[0].grid(True)

ax[1].loglog(t_U, U, 'b-', label='Simulation')
ax[1].axhline(U_s, color='g', linestyle='-.', linewidth=2, label=r'Slumping: $U_s$')
ax[1].loglog(t_box, U_box, 'r--', linewidth=2, label=r'Box model, $Fr(h/H)$')
ax[1].loglog(t_box[1:], U_inert, 'k:', linewidth=1.5, label=r'$\propto t^{-1/3}$ (thin current)')
ax[1].set_xlabel(r'$t$ (s)')
ax[1].set_ylabel(r'$U = \mathrm{d}x_f/\mathrm{d}t$ (m/s)')
ax[1].set_xlim([0.3, t[-1]])
ax[1].set_ylim([0.02, 0.3])
ax[1].legend()
ax[1].grid(True, which='both', alpha=0.4)

plt.tight_layout()
plt.savefig('front.png')
~~~

~~~pythonplot Front condition
import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import re

L0 = 4.
N = 1024
h_1 = 0.3
l_1 = 0.5
H = 1.5 * h_1
drho = 0.011
gp = 9.81 * drho
dt_save = 0.01
t_end = 20.
h_min = 0.01
head_len = 0.3

def front(data):
    x, y = data[:, 0], data[:, 1]
    dx = L0 / N
    bins = np.floor(x / dx).astype(int)
    hmax = np.zeros(bins.max() + 1)
    np.maximum.at(hmax, bins, y)
    xf = (np.where(hmax > h_min)[0].max() + 1) * dx
    hh = y[(x > xf - head_len) & (x <= xf)].max()
    return xf, hh

files = gb.glob('interface-*.txt')
files.sort(key=lambda f: int(re.search(r'\d+', f).group()))
t_list, xf_list, hh_list = [], [], []
for name in files:
    tn = (int(re.search(r'\d+', name).group()) - 1) * dt_save
    if tn > t_end + 1e-9:
        continue
    data = np.loadtxt(name, ndmin=2)
    if data.size == 0:
        continue
    xf, hh = front(data)
    t_list.append(tn)
    xf_list.append(xf)
    hh_list.append(hh)
t, xf, hh = np.array(t_list[:-1]), np.array(xf_list[:-1]), np.array(hh_list[:-1])   # last file discarded
keep = xf < L0 - 0.3
t, xf, hh = t[keep], xf[keep], hh[keep]

w = int(round(0.5 / dt_save))
wm = int(round(1.5 / dt_save))   # half-window of the running median (3 s)
t_U, U, h_U, x_U = [], [], [], []
for i in range(w, len(t) - w):
    U.append(np.polyfit(t[i-w:i+w+1], xf[i-w:i+w+1], 1)[0])
    t_U.append(t[i])
    x_U.append(xf[i])
    h_U.append(np.median(hh[max(0, i-wm):i+wm+1]))   # running median: filters the billows
t_U, U, h_U, x_U = np.array(t_U), np.array(U), np.array(h_U), np.array(x_U)

ok = x_U > l_1 + head_len      # the head window no longer overlaps the lock

Fr_num = U / np.sqrt(gp * h_U)
hh_med = np.array([np.median(hh[max(0, i-wm):i+wm+1]) for i in range(len(hh))])

phi = np.linspace(0.03, 1, 500)
Fr_HS = np.where(phi < 0.075, 1.19, 0.5 * phi**(-1/3))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

ax[0].plot(t, hh / H, color='lightgray', label='Max height (raw)')
ax[0].plot(t, hh_med / H, 'b-', label='Running median (3 s)')
ax[0].axvspan(0, t_U[ok][0] if ok.any() else 0, color='gray', alpha=0.15, label='Window overlaps the lock')
ax[0].axhline(h_1 / (2 * H), color='g', linestyle='-.', linewidth=2, label=r'$h_1/2$ (slumping)')
ax[0].set_xlabel(r'$t$ (s)')
ax[0].set_ylabel(r'$h_\mathrm{head}/H$')
ax[0].set_ylim([0, 1])
ax[0].legend()
ax[0].grid(True)

sc = ax[1].scatter(h_U[ok] / H, Fr_num[ok], c=t_U[ok], cmap='viridis', s=8)
ax[1].plot(phi, Fr_HS, 'r--', linewidth=2, label='Huppert & Simpson (1980)')
ax[1].set_xlabel(r'$h_\mathrm{head}/H$')
ax[1].set_ylabel(r'$Fr = U/\sqrt{g^\prime h_\mathrm{head}}$')
ax[1].set_xlim([0, 0.7])
ax[1].set_ylim([0, 1.5])
ax[1].legend()
ax[1].grid(True)
cbar = fig.colorbar(sc, ax=ax[1])
cbar.set_label(r'$t$ (s)')

plt.tight_layout()
plt.savefig('front_condition.png')
~~~
*/