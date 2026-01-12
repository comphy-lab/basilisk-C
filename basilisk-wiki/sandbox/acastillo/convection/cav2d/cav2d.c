/**
# Differentially heated cavity in 2D

For this example we solve the [Boussinesq
equations](https://basilisk.fr/sandbox/acastillo/convection/convection_boussinesq.h).

In this example, we fix $\text{Ra}=10^3$ and $\text{Pr}=0.71$ which leads to a steady-state
solution. 

~~~pythonplot Streamfunction and temperature fields at the steady state
import numpy as np
import matplotlib.pyplot as plt

import numpy as np

def input_matrix(file_path):
  with open(file_path, 'rb') as file:
    # Read n as a single-precision float
    n = np.fromfile(file, dtype=np.float32, count=1)[0]
    n = int(n)  # Convert to integer for indexing

    # Initialize arrays
    y = np.zeros(n, dtype=np.float32)
    x = np.zeros(n, dtype=np.float32)
    f = np.zeros((n, n), dtype=np.float32)

    # Read y vector
    for j in range(n):
      y[j] = np.fromfile(file, dtype=np.float32, count=1)[0]

    # Read x vector and f matrix
    for i in range(n):
      x[i] = np.fromfile(file, dtype=np.float32, count=1)[0]
      for j in range(n):
        f[j, i] = np.fromfile(file, dtype=np.float32, count=1)[0]

  return f, n, x, y


T, _, _, _ = input_matrix("temperature.bin")
psi, n, x, y = input_matrix("streamfunction.bin")

fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=False, figsize=(8.5,3.5))

im1 = ax1.pcolormesh(x, y, psi, cmap='gnuplot_r', rasterized=True)
#ax1.contour(x, y, psi, levels=10, colors='k')
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$')
ax1.set_title('Stream-function $\\psi(x,y)$')
fig.colorbar(im1, ax=ax1)

im2 = ax2.pcolormesh(x, y, T, cmap='PuOr', rasterized=True)
#ax2.contour(x, y, T, levels=10, colors='k')
ax2.set_aspect('equal', adjustable='box')
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$')
ax2.set_title('Temperature $\\theta(x,y)$')
fig.colorbar(im2, ax=ax2)

plt.tight_layout()
plt.savefig('steady_fields.svg')
~~~

~~~pythonplot Convergence
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('cav2d.asc', delimiter=' ', usecols=[2,3,4,10,11], skiprows=1)

psimax = np.maximum(np.abs(data[:,3]), np.abs(data[:,4]))

fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=False, figsize=(8.5,3.5))

ax1.scatter(data[:,0], data[:,1], marker='s', label='left wall');
ax1.scatter(data[:,0], data[:,2], marker='.', label='right wall');
ax1.set_xscale('log', base=2)
ax1.legend()
ax1.set_xlabel('$N$')
ax1.set_ylabel('Nusselt')

ax2.scatter(data[:,0], psimax, marker='s');
ax2.set_xscale('log', base=2)
ax2.set_xlabel('$N$')
ax2.set_ylabel(r'$|\psi|_{\text{max}}$')
plt.tight_layout()
plt.savefig('convergence.svg')
~~~

*/

#define MINLEVEL 5
#define MAXLEVEL 8
#define RAYLEIGH 1e3
#define PRANDTL 0.71
#define MAXTIME 100.
#include "grid/multigrid.h"
#include "acastillo/convection/convection_boussinesq.h"
#include "acastillo/convection/global_nusselt.h"

scalar omega[];
scalar psi[];

int main() {

  // Geometric Parameters
  L0 = 1.;
  X0 = Y0 = -0.5;

  // Physical Parameters
  Ra = RAYLEIGH; 
  Pr = PRANDTL; 	

  // Numerical Parameters  
  DT = 0.05;  
  TOLERANCE = 1e-4;
  NITERMIN = 4;

  FILE * fp = fopen("cav2d.asc", "w");
  fprintf (fp, "[1]Ra ");
  fprintf (fp, "[2]Pr ");
  fprintf (fp, "[3]N ");
  fprintf (fp, "[4]nuleft ");
  fprintf (fp, "[5]nuright ");
  fprintf (fp, "[6]nuvol ");
  fprintf (fp, "[7]umin ");
  fprintf (fp, "[8]umax ");
  fprintf (fp, "[9]vmin ");
  fprintf (fp, "[10]vmax ");
  fprintf (fp, "[11]psimin ");
  fprintf (fp, "[12]psimax \n");
  fclose (fp);

  for (int l = (MINLEVEL); l <= (MAXLEVEL); l++) {
    N = 1<<l;
    run();
  }
}

/**
## Boundary conditions

The left and right walls are iso-thermal $\theta = \mp 1/2$ at $x = \pm 1/2$,
while all other boundaries are adiabatic. In addition to no-slip walls on all
boundaries.

*/
T[left]  = dirichlet(-0.5);
T[right] = dirichlet(0.5);

u.n[top]    = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right]  = dirichlet(0.);
u.n[left]   = dirichlet(0.);

u.t[top]    = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right]  = dirichlet(0.);
u.t[left]   = dirichlet(0.);

/**
## Initial conditions
Initial conditions correspond to a linear temperature profile and no motion.
*/

event init (i = 0){
  foreach(){
    T[] = x;
    foreach_dimension()
      u.x[] = 0.;
  }
}

/**
We also initialize a velocity field at time $t_n$, $u_n$. This way, we may
follow the maximum change in the horizontal velocity over two consecutive time
units to identify if the steady-state is reached.
*/
scalar un[];
event init_un (i = 0)
  foreach()
    un[] = u.x[];

/**
## Outputs

We store the residuals from the Poisson solvers

~~~pythonplot 0D Quantities
import numpy as np
import matplotlib.pyplot as plt

filename1 = 'log'

data = np.genfromtxt(filename1, usecols=[0,2,6,12,18], invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]

series_id = data[:,0]
t = data[:,1]
residual_p = data[:,2]
residual_u = data[:,3]
residual_T = data[:,4]

fig, axes = plt.subplots(1,3, figsize=(4*3, 3))
ax = axes.ravel()

# Plot each series individually with different colors
markers = ['.', 'o', 's', '^', 'v', 'D', 'p', '*']
for i, sid in enumerate(np.unique(series_id)):
  mask = (series_id == sid) & (t > 1.0)
  marker = markers[i]
  ax[0].semilogy(t[mask], residual_p[mask], marker, label=f'N={sid:g}', alpha=0.7)
  ax[1].semilogy(t[mask], residual_u[mask], marker, label=f'N={sid:g}', alpha=0.7)
  ax[2].semilogy(t[mask], residual_T[mask], marker, label=f'N={sid:g}', alpha=0.7)

xlabels = ['Pressure', 'Velocity', 'Temperature']
for i in range(3):
  ax[i].set_xlabel('Time')
  ax[i].set_ylabel('Residual '+xlabels[i])
ax[0].legend()

plt.tight_layout()
plt.savefig('plot_residuals.svg', dpi=150)
~~~

It is also possible to plot the performance metrics from the out file.

~~~pythonplot Performance metrics from out file
import numpy as np
import matplotlib.pyplot as plt
import re

filename2 = 'out'

# Parse the out file
N_values = []
real_time = []
mpi_max_pct = []

with open(filename2, 'r') as f:
  lines = f.readlines()
  i = 0
  while i < len(lines):
      line = lines[i].strip()
      if line.endswith('var'):
        match = re.search(r'([\d.]+)\s+real', line)
        if match:
          real_time.append(float(match.group(1)))          

          # Check next line for MPI data
          if i + 1 < len(lines):
              mpi_line = lines[i + 1].strip()
              mpi_match = re.search(r'max\s+[\d.]+\s+\(([\d.]+)%\)', mpi_line)
              if mpi_match:
                mpi_max_pct.append(float(mpi_match.group(1)))
      i += 1

MINLEVEL = 5
N_values = [2**(MINLEVEL + i) for i in range(len(real_time))]

# Convert to numpy arrays
N_values = np.array(N_values)
real_time = np.array(real_time)
mpi_max_pct = np.array(mpi_max_pct)

# Create plots
fig, axes = plt.subplots(1, 2, figsize=(4*2, 3))
ax = axes.ravel()
ax[0].loglog(N_values, real_time, 'o-', markersize=8)
ax[0].set_xlabel('Grid size N')
ax[0].set_ylabel('Real time (s)')
ax[0].set_xscale('log', base=2)

ax[1].semilogx(N_values, mpi_max_pct, '^-', markersize=8)
ax[1].set_xlabel('Grid size N')
ax[1].set_ylabel('MPI max (%)')
ax[1].set_xscale('log', base=2)

plt.tight_layout()
plt.savefig('plot_performance.svg', dpi=150)
~~~

*/

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, " \t - \t %d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

event log (i++) {
  fprintf (stderr, " %d %d %g ", N, i, t);
  mg_print (mgp);
  mg_print (mgu);
  mg_print (mgT);
  fprintf (stderr, "\n");
  fflush (stderr);
}
/**
We are going to define two auxiliary functions: 
 
- Compute the stream-function $\psi(x,y)$ by solving a poisson problem for the
  vorticity.

*/

void streamfunction (vector u, scalar psi) {
  vorticity (u, omega);
  poisson (psi, omega);
}

/** 
- Compute (and save) the heat-fluxes and peak velocities
*/

void compute_and_store_fluxes (vector u, scalar T){
  
  streamfunction (u, psi);
  stats spsi = statsf (psi); 
  
  stats svelx = statsf (u.x); 
  stats svely = statsf (u.y);

  double nu1 = nusselt_right(T);
  double nu2 = nusselt_left(T);
  double nu3 = nusselt_vol(T,u);

  if (pid()==0){
    FILE * fp = fopen("cav2d.asc", "a");
    fprintf (fp, "%.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n", Ra, Pr, N,
      nu1, nu2, nu3, 
      svelx.min, svelx.max, 
      svely.min, svely.max, 
      spsi.min, spsi.max);
    fclose (fp);
  }
}



/**

- Save the temperature and streamfunction fields to create the figures

*/

void save_fields (scalar T, vector u, int nf){
  FILE * fp ;
  fp = fopen("temperature.bin", "w");
  output_matrix (T, fp, N, linear = true);
  fclose (fp);

  streamfunction (u, psi);

  fp = fopen("streamfunction.bin", "w");
  output_matrix (psi, fp, N, linear = true);
  fclose (fp);
}

/**

Quantities are stored at the end of the simulation (`t=MAXTIME`) or when a
steady-state is reached.

*/

event is_steady (t += 1.0; t <= MAXTIME) {
  double deltau = change (u.x, un);
  if (deltau < 1e-6 && t > 10.0){
    compute_and_store_fluxes (u, T);
    return 1;
  }
}

event finalize(t = end){
  save_fields(T,u,0);
}
	
/**

# References

~~~bib

@article{de1983natural,
  title={Natural convection of air in a square cavity: a bench mark numerical solution},
  author={de Vahl Davis, Graham},
  journal={International Journal for numerical methods in fluids},
  volume={3},
  number={3},
  pages={249--264},
  year={1983},
  publisher={Wiley Online Library}
}

~~~
*/