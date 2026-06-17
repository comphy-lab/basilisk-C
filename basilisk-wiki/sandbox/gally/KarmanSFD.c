/**
# Selective Frequency Damping method

The low-pass filtered version of the initial flow is calculated at each iteration using an Euler explicite scheme ($\Delta$ being to control parameter of the method) :
$$
    \mathbf{\bar{q}}^{n+1}_i \simeq
    \left( \mathbf{q}^n_i-\mathbf{\bar{q}}^n_i \right)
    \frac{\delta t}{\Delta}
    +\mathbf{\bar{q}}^n_i
$$

The flow is then forced by adding to the acceleration term ($\chi$ the second control parameter) : 
$$
-\chi(\mathbf{q}^n_i-\mathbf{\bar{q}}^n_i)
$$


# Bénard–von Kármán Vortex Street's stabilisation for flow around a cylinder at Re = 160

See the [example code](https://basilisk.fr/src/examples/karman.c) for the visualisation of the initial flow.

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
int maxlevel = 9;
face vector muv[];

/**
We set the constants of the flow.*/
double Reynolds = 160.;
double D = 0.125, U0 = 1.;

/**
We create all variables necessary for the SFD.
The Strouhal number has to be changed according to the Reynolds chosen. */

double Strouhal = 0.1871;
float SFD_startingtime = 10.;
float SFD_stoppingtime = 25.;

double SFD_delta, SFD_chi;
vector ubar[];
scalar pbar[];
face vector av[];
double SFD_res, res;

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

  /** We set the SFD control terms. 
  $\Delta$ and $\chi$ not yet optimised but functional.*/
  a = av;
  SFD_delta = D/(M_PI*Strouhal*U0);
  SFD_chi = 1/SFD_delta;

  /**
  When using bview we can interactively control the Reynolds number
  and maximum level of refinement. */
 
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);

  run();
}

/**
We set a constant viscosity based on the Reynolds number, the cylinder
diameter $D$ and the inflow velocity $U0$. */

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

/**
The fluid is injected on the left boundary with velocity $U0$. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The cylinder is no-slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

/** We use the same boundary conditions for the filtered flow.*/

ubar.n[left]  = dirichlet(U0);
pbar[left]    = neumann(0.);

ubar.n[right] = neumann(0.);
pbar[right]   = dirichlet(0.);

ubar.n[embed] = dirichlet(0.);
ubar.t[embed] = dirichlet(0.);

event init (t = 0) {
 
  solid (cs, fs, sqrt(sq(x) + sq(y)) - D/2.);
 
  foreach() {
    u.x[] = cs[] ? U0 : 0.;
    ubar.x[] = cs[] ? U0 : 0.;
    ubar.y[] = 0.;
    pbar[] = 0.;
  }
}

/**
We add the SFD correction term to NS. */

event acceleration (i++) {
    foreach_face()
      av.x[] = ((t >= SFD_startingtime) & (t <= SFD_stoppingtime)) ? - SFD_chi * ((u.x[] - ubar.x[]) + (u.x[-1] - ubar.x[-1]))/2. : 0.;
}

/**
We calculate the next step's filtered field.*/

event bar (i++) {

  int n = 0;
  double res = 0;
 
  foreach() {

    res += sqrt( sq((u.x[] - ubar.x[])/U0) + sq((u.y[] - ubar.y[])/U0) + sq((p[] - pbar[])/sq(U0)) );
    n++;
   
    ubar.x[] = (u.x[] - ubar.x[]) * (dt/SFD_delta) + ubar.x[];
    ubar.y[] = (u.y[] - ubar.y[]) * (dt/SFD_delta) + ubar.y[];
    pbar[] = (p[] - pbar[]) * (dt/SFD_delta) + pbar[];
  }

  SFD_res = res / n;
 
  static FILE * fbar = NULL;
  if (!fbar) {
    char filenamebar[25];
    snprintf(filenamebar, sizeof(filenamebar), "res_Re%.0f.dat", Reynolds);
    fbar = fopen(filenamebar, "w");
  }
  fprintf(fbar, "%g %g\n", t, SFD_res);
}

/**
We check the number of iterations of the Poisson and viscous problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We follow the evolution of a point situated after the cylinder. */

double xp, yp;

event probe (i++) {
 
  xp = 0.5*D + 3*D;
  yp = 0.;
  
  double vp = interpolate(u.y, xp, yp);

  static FILE * fp = NULL;
  if (!fp) {
    char filename[20];
    snprintf(filename, sizeof(filename), "probe_Re%.0f.dat", Reynolds);
    fp = fopen(filename, "w");
  }
  fprintf(fp, "%g %g\n", t, vp);
}

/**
We produce animations of the vorticity field... */

event movies (i += 4; t <= 35.) {
  
  scalar omega[], m[];
  
  vorticity (u, omega);
 
  foreach() {
    m[] = cs[] - 0.5;
    if (((x < -0.4) && (y > 0.4)) && ((t >= SFD_startingtime) && (t <= SFD_stoppingtime)))
	m[] = -0.5; // SFD on/off indicator
  }

  char name_video[30];
  snprintf(name_video, sizeof(name_video), "vort_Re%.0f.mp4", Reynolds);
  output_ppm (omega, file = name_video, box = {{-0.5,-0.5},{7.5,0.5}},
              min = -10, max = 10, linear = true, mask = m);
  
  char name_videof[30];
  snprintf(name_videof, sizeof(name_videof), "f_Re%.0f.mp4", Reynolds);
  output_ppm (f, file = name_videof, box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
  
}

/**
We adapt according to the error on the embedded geometry, velocity and tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, {1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}

/**
# Visualisations

The SFD is activated at t = 10 and deactivated at t = 25. A dot appears on the upper-left corner when the SFD is active.

![Animation of the vorticity field](KarmanSFD/vort_Re160.mp4)(loop)

![Animation of the tracer field](KarmanSFD/f_Re160.mp4)(loop)

# Residual's evolution through time

~~~pythonplot
import numpy as np
import matplotlib.pyplot as plt

# Parameters
Re = 160

# Load data: t  res
data = np.loadtxt(f"res_Re{Re}.dat")
t   = data[:, 0]
res = data[:, 1]

# (De-)activation SFD
t_activation = 10
t_deactivation = 25

plt.plot(t, res, 'k')
plt.axvline(t_activation, linestyle = '--', color = 'r', label = 'SFD activation')
plt.axvline(t_deactivation, linestyle = '-.', color = 'r', label = 'SFD deactivation')
plt.xlabel("Time")
plt.ylabel("Residual")
plt.xlim([0, 35])
plt.ylim([0, 0.2])
plt.legend()
plt.tight_layout()
plt.savefig("residual.png")

~~~
*/

/**
# Vertical velocity behind the cylinder through time

~~~pythonplot
import numpy as np
import matplotlib.pyplot as plt

# Parameters
Re = 160

# Load data: t  vp
data = np.loadtxt(f"probe_Re{Re}.dat")
t  = data[:, 0]
vp = data[:, 1]

# (De-)activation SFD
t_activation = 10
t_deactivation = 25

# Plot
plt.figure()
plt.plot(t, vp, 'k')
plt.axvline(t_activation, linestyle = '--', color = 'r', label = 'SFD activation')
plt.axvline(t_deactivation, linestyle = '-.', color = 'r', label = 'SFD deactivation')
plt.xlabel("Time")
plt.ylabel("Velocity")
#plt.xlim([0, 35])
plt.legend()
plt.tight_layout()
plt.savefig("velocity.png")

~~~
*/