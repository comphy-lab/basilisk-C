/**
# Schlichting 2D Jet

We try to recover the famous problem of the laminar 2D jet in an immersed domain, also called the Schlichting jet or the Bickley jet. This analytical solution can be found in H. Schlichting's book "Boundary layer theory" (McGraw-Hill) in section IX.f. It was derived by H. Schlichting (1933) himself and W. Bickley (1939).

It is a monophasic problem depending only on the Reynolds number $Re$. Here we take $Re = 10$.

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;
double grav;

FILE * fpmax; //

face vector muv[];

/**
We choose to impose left-right symmetry by simulating only one half of the domain (using slip boundary conditions at the right). Left and top boundaries are outflow conditions, and imposed flow at the bottom. 
*/

int main() {

  Re=10;

  R_d=0.01; //Rayon du jet initial
  L0=1; //Taille de la boite  
  U0=1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (U0*(x > -R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);

  s[bottom] = dirichlet (U0*(x > -R_d));

  u.n[left] = neumann(0);
  p[left] = dirichlet(0);

  u.n[right] = dirichlet(0.);
  p[right]=neumann(0.);

  N=256; //256 works even better
  origin (-L0, 0);
  init_grid(N);
  
  char param_dim[80];
  sprintf (param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf (fparam, "%g %g %g %d\n",L0,R_d,U0,N);
  fclose (fparam);

  mu=muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
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

Residuals are converged around t=800.
*/

/**
We plot the vertical velocity as well as a passive tracer to observe the establishment of the flow.
*/

event ppm_output (t = 0; t += 0.1; t <= 550) {
  
    char name1[80];
    sprintf (name1, "uY.mp4");
    output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

    // optionally tracer
    char name2[80];
    sprintf (name2, "s.mp4");
    output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

/**
![Passive tracer](schlichting/s.mp4)
![Vertical velocity](schlichting/uY.mp4)

We can clearly see that the Neumann condition at the top for the velocity creates a slow down of the flow, because the condition does not respect the flow solution : *equations to be written soon*.

*/

/**
We save the different fields at the last iteration.
*/

event profile (t = end) {
  char name[80];
  sprintf (name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g %g %g", x, y, u.x[], u.y[], p[],s[]);
  fclose(fpres);
  
  printf ("-----END-----\n");
}

/**
## Comparison with theory

*soon*
*/