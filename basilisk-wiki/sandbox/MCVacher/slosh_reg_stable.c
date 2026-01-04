/**
# Stable jet under a free-surface

Description : *soon*
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h" 
#include "navier-stokes/conserving.h"
#include "tracer.h"
#include "tag.h"
#include "view.h"
#define BVIEW 1


double h;
double U0;
double R_d;

/**
A passive tracer is injected with the jet to follow lagrangian paths.
*/

scalar s[];
scalar * tracers = {s};

FILE * fpmax; //

int main() {

  R_d=0.003; 
  L0=0.4; 
  rho2 = 1.3;
  rho1=1000;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.45;
  h=0.15;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
  p[top] = dirichlet(0.);
  s[bottom] = dirichlet (f[]*U0*(x > -R_d && x <R_d));

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.n[right] = y < R_d ? dirichlet(U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);
  u.t[right] = y < R_d ? neumann(0.) : dirichlet(0.);
 
  N=256;
  origin (-L0/2, 0);
  init_grid(N);
  
  fpmax =  fopen("log.dat", "w");
  
  f.sigma = 0.072;

  run();
}

/**
We initiate gravity, which is opposing to the inertia of the jet.
*/

event init (t = 0) {
  fraction (f, y<h);
}

#if 1
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] = -9.81;
}
#endif

/**
At each iteration, we calculate residuals to see if the solution has converged. They are saved in "log.dat".
*/

double sum_u = 0.;
double sum_u_t = 0.;
double res_u = 0.;
double sum_p = 0.;
double sum_p_t = 0.;
double res_p = 0.;

event logfile (i++) { 
  sum_u_t=0.;
  sum_p_t=0.;

  foreach_face(){
    sum_u_t += sqrt(sq(u.x[]) + sq(u.y[]));
  }
  foreach(){
    sum_p_t += p[];
  }
  res_u = fabs(sum_u_t - sum_u);
  res_p = fabs(sum_p_t - sum_p);

  sum_u=sum_u_t;
  sum_p=sum_p_t;

  fprintf (stderr, "%d %g %g %g\n", i, t, res_u,res_p); 
  fprintf (fpmax, "%d %g %g %g\n", i, t, res_u,res_p);
}

/**
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
plt.semilogy(list_t,list_res_p,'g-',label= r'Res. $p$')
plt.semilogy(list_t,list_res_u,'r-',label= r'Res. $|u|$')
plt.legend()
plt.savefig('residuals.png')
~~~
*/

/**
We kill eventual numerical bubbles (not real bubbles because there is no surface tension here.
*/

event remove_droplets (i++) {
  remove_droplets (f, threshold=0.05, bubbles=true);
}


event profile (t = end) {
  printf ("-----END-----\n");
}

int isave1 = 1;
event res_save (t += 1; t <= 100) {
  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);
  
  isave1++;
}

/**
To visualize both the surface and the jet, we can follow lagrangian trajectories using a passive tracer. We generate videos:
*/
event ppm_output (t = 0; t += 0.5; t <= 100) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  // optionally tracer
  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

event movie_streamlines (t += 0.5; t <= 100)
{
  scalar omega[];
  vertex scalar stream[];
  scalar psi[];

  vorticity(u, omega);

  psi[bottom] = dirichlet (0);
  psi[top]    = neumann (0);
  psi[left]   = dirichlet (0);
  psi[right]  = dirichlet (0);

  poisson (psi, omega);
  boundary ({psi});

  foreach_vertex()
    stream[] = (psi[0,-1] + psi[-1,-1] + psi[] + psi[-1])/4.;

  clear();
  view(
    width = L0,
    height = L0,
    tx = 0,           // center horizontally
    ty = -L0        // center vertically
);
  isoline ("stream", n = 30);
  box();
  save ("streamlines.mp4");
}

/**
![Free-surface](slosh_reg_stable/f.mp4)
![Passive tracer](slosh_reg_stable/s.mp4)
![Vertical velocity](slosh_reg_stable/uY.mp4)
![Streamlines](slosh_reg_stable/streamlines.mp4)
*/


