/**
# Differentially heated cavity in 2D

For this example we solve the Boussinesq equations. We fix $\text{Ra}=10^3$ and
$\text{Pr}=0.71$ which leads to a steady-state solution. We use a $256^2$ grid.

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

plt.savefig('convergence.svg')
~~~

*/

#define MINLEVEL 6
#define MAXLEVEL 8
#define RAYLEIGH 1e3
#define PRANDTL 0.71
#define MAXTIME 100.
#include "convection_boussinesq.h"

scalar omega[];
scalar psi[];

double tend= MAXTIME;
int main() {

  // Geometric Parameters
  L0 = 1.;
  X0 = Y0 = -0.5;

  // Physical Parameters
  Ra = RAYLEIGH; 
  Pr = PRANDTL; 	

  // Numerical Parameters  
  DT = 0.05;  
  TOLERANCE = 1e-3;
  NITERMIN = 4;

  FILE * fp = fopen("cav2d.asc", "w");
  fprintf (fp, "[1]Ra [2]Pr [3]N [4]nuleft [5]nuright [6]nuvol"
      " [7]umin [8]umax [9]vmin [10]vmax [11]psimin [12]psimax \n");
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

*/

event logfile(t+=5.0; t <= tend){
  fprintf(stderr, " residuals: %.5f \t %g %d %d \t %g %d %d \t %g %d %d \n", t, 
    mgp.resa, mgp.i, mgp.nrelax, 
    mgu.resa, mgu.i, mgu.nrelax, 
    mgT.resa, mgT.i, mgT.nrelax);
}

/**
We are going to  define two auxiliary functions: 
 
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

#include "global_nusselt.h"
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

Quantities are stored at the end of the simulation (`t=tend`) or when a
steady-state is reached.

*/

event is_steady (t += 1.0; t <= tend) {
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