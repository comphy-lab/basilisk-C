/**
# Simulation of a Jet Impact on a Slip Surface

Only half of the domain is simulated to capture the stationnary behaviour.

Parameters:
- $Re_j$: Reynolds number (default 40)
- $R_d$: Initial jet radius
- $U_0$: Jet velocity
- $L_0$: Domain size
- N: Grid resolution

Output:
- "log.dat": Residuals for velocity and pressure
- "res_end.txt": Final profiles
- "uY.mp4": Vertical velocity field over time
- "s.mp4": Passive tracer over time
- "res_t*.txt": State dumps for restart and analysis
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;
double Re_j;

FILE * fpmax; 

face vector muv[];

int main() {

  Re_j = 40;
  R_d = 0.01;   // Initial jet radius
  L0  = 1;      // Domain size
  U0  = 1;      // Jet velocity

  TOLERANCE = 1e-3 [*]; 

  // Boundary conditions
  u.n[bottom] = dirichlet (U0*(x > -R_d)); 
  u.t[bottom] = dirichlet(0.);
  s[bottom]   = dirichlet (U0*(x > -R_d));

  u.n[top] = dirichlet(0);
  p[top]   = neumann(0);

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);

  u.n[right] = dirichlet(0.);
  p[right]   = neumann(0.);

  //Initialize grid
  N = 256;  
  origin(-L0, 0);
  init_grid(N);
  
  //Save domain parameters to file
  char param_dim[80];
  sprintf(param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf(fparam, "%g %g %g %d \n", L0, R_d, U0, N);
  fclose(fparam);

  mu = muv;

  fpmax = fopen("log.dat", "a"); 

  run();
}

event init (t = 0) {
  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
  } else {
    fprintf(stderr, "Restarting from previous dump\n");
  }
}

/** Viscosity assignment for each face, ensuring Reynolds */
event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[] * U0 * R_d / Re_j;
}

/** Variables to track residuals */
double sum_u = 0.;
double sum_u_t = 0.;
double res_u = 0.;
double sum_p = 0.;
double sum_p_t = 0.;
double res_p = 0.;

/** Log residuals for velocity and pressure
- Computes L1 norms and appends to log file */
event logfile (i++) { 
  sum_u_t = 0.;
  sum_p_t = 0.;

  foreach_face() {
    sum_u_t += sqrt(sq(u.x[]) + sq(u.y[]));
  }
  foreach() {
    sum_p_t += p[];
  }

  res_u = fabs(sum_u_t - sum_u);
  res_p = fabs(sum_p_t - sum_p);

  sum_u = sum_u_t;
  sum_p = sum_p_t;

  fprintf(stderr, "%d %g %g %g \n", i, t, res_u, res_p); 
  fprintf(fpmax, "%d %g %g %g \n", i, t, res_u, res_p);
}

/** Save final profiles at the end of simulation */
event profile (t = end) {
  char name[80];
  sprintf(name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf(fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[], s[]);
  fclose(fpres);
  
  printf("-----END-----\n");
}

double t_prev_dump = 0.;

/** Output movie frames of velocity and tracer fields */
event movie_frames (t = 0; t += 0.5; t <= 500) {
  char ufile[80], sfile[80];
  sprintf(ufile, "uY.mp4");
  sprintf(sfile, "s.mp4");
  output_ppm(u.y, file = ufile, n = 512, min = -U0, max = +U0, linear = true);
  output_ppm(s,   file = sfile, n = 512, min = 0.,  max = U0,  linear = true);
}

/** Periodically dump simulation state for restart */
event dump_state (t = 0; t += 10; t <= 500) {
  dump("restart");
  t_prev_dump = t;
  fprintf(stderr, "Dumped state at t = %g\n", t);

  char name[80];
  sprintf(name, "res_t%.1f.txt", t);
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf(fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[], s[]);
  fclose(fpres);
}

/**
## Videos

![Passive tracer](box_jet_tim/s.mp4)
![Vertical velocity](box_jet_tim/uY.mp4)
*/

/**
## Residuals

We plot the residuals at each time :

~~~pythonplot Residuals
import matplotlib.pyplot as plt
import numpy as np

list_t=[]
list_res_u=[]
list_res_p=[]

with open('log.dat', encoding='utf8') as File:
    for line in File:
        if line.isspace()==0:
            temp=line.split()
            t=float(temp[1])
            list_t.append(t)
            norm_u=float(temp[2])
            list_res_u.append(norm_u)
            p=float(temp[3])
            list_res_p.append(p)
    File.close()

plt.figure()
plt.semilogy(list_t, list_res_p, 'g-', label=r'Res. $p$')
plt.semilogy(list_t, list_res_u, 'r-', label=r'Res. $|u|$')
plt.legend()
plt.savefig('residuals.png')
~~~

Residuals converge around t = 1800.
*/
