/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re = 160

Re-use of the [example code](https://basilisk.fr/src/examples/karman.c)

Modifications to extract velocity and pressure after the cylinder.*/

#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

/**
The domain is eight units long, centered vertically. */

int main()
{
  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

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

double D = 0.125, U0 = 1.;

event properties (i++)
{
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

event init (t = 0)
{

  /**
  The domain contains a cylinder of diameter 0.125. */

  solid (cs, fs, sqrt(sq(x) + sq(y)) - D/2.);

  /**
  We set the initial velocity field. */
  
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}

/**
We check the number of iterations of the Poisson and viscous problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
   We follow the evolution of velocity and pressure at a point situated after the cylinder. */

double xp, yp;

event probe (i++){

  xp = 0.5*D + 3*D;
  yp = 0.;
  
  double up = interpolate(u.x, xp, yp);
  double vp = interpolate(u.y, xp, yp);
  double pp = interpolate(p,   xp, yp);

  static FILE * fp = NULL;
  if (!fp) {
    char filename[20];
    snprintf(filename, sizeof(filename), "probe_Re%.0f.dat", Reynolds);
    fp = fopen(filename, "w");
  }
  fprintf(fp, "%g %g %g %g\n", t, up, vp, pp);
}

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 4; t <= 20.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach() {
    m[] = cs[] - 0.5;
    if ((fabs(x - xp) < 3.*Delta) && (fabs(y - yp) < 3.*Delta)){
      omega[] = 0;
    }
  }
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
}

/**
We adapt according to the error on the embedded geometry, velocity and tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, {1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}

/**
# Strouhal number at Re = 160

~~~pythonplot Fourier transform
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq

# Define Reynolds
Re = 160

# Load data: t  up  vp  pp
data = np.loadtxt(f"probe_Re{Re}.dat")
t  = data[:, 0]
vp = data[:, 2]

# Temporal plot
plt.figure()
plt.plot(t, vp)
plt.xlabel("Time")
plt.ylabel("Velocity")
plt.title("Vertical velocity behind the cylinder")
plt.tight_layout()
plt.savefig("temp160.png")
plt.close()

vp = vp - np.mean(vp)

# Mask (wait until instabilities are established)
mask = (t > 10.)
vp = vp[mask]
t = t[mask]

# Time-step
dt = t[1] - t[0]
N = len(vp)

# FFT
spectrum = rfft(vp)
freq = rfftfreq(N, d = dt)

# Strouhal number
f_peak = freq[np.argmax(spectrum)]
St = f_peak * 0.125    # St = f*D/U0, with D=0.125, U0=1
print(f"Vortex-shedding frequency: {f_peak:.4f}")
print(f"Strouhal number: {St:.4f}")

# Plot
plt.figure()
plt.plot(freq, np.real(spectrum))
plt.xlabel("Frequency")
plt.ylabel("Spectrum")
plt.xlim(0, 20)
plt.title("FFT vertical speed behind the cylinder")
plt.suptitle(f"Vortex-shedding frequency: {f_peak:.4f} Hz"
             f"\n Strouhal number:    {St:.4f}"
             f"\n\n (With Re = 160, U0 = 1, D = 0.125)",
             x = 0.6, y = 0.4)
plt.tight_layout()
plt.savefig("strouhal160.png")
plt.close()
~~~
*/

/**
# Strouhal number for different Reynolds

After running the simulation for different Reynolds, we can plot the following graph: */