/**
# Jet/jet interaction

Jet 1 (left,  x = -d/2): velocity U1, radius R1, tracer s = +1

Jet 2 (right, x = +d/2): velocity U2, radius R2, tracer s = -1

We initiate the simulation with the sum of two theoretical Schlichting/Bickley profiles, applied only where r > r_min to avoid the near-field singularity (r^{-1/3} -> inf).
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#define BVIEW 1

face vector muv[];

scalar s[];
scalar * tracers = {s};

double U1, U2;
double R1, R2;
double Re;
double d;

double nu_val;
double A1, B1;
double A2, B2;

FILE * fpmax;

/** ## Initial profiles definition */

static inline double schlichting_v (double xpos, double r, double x_jet,
                                    double A, double B)
{
  if (r <= 0.) return 0.;
  double eta = (xpos - x_jet) / pow(r, 2./3.);
  double ch  = cosh(B * eta);
  return A * pow(r, -1./3.) / (ch * ch);
}

static inline double schlichting_u (double xpos, double r, double x_jet,
                                    double A, double B)
{
  if (r <= 0.) return 0.;
  double eta      = (xpos - x_jet) / pow(r, 2./3.);
  double ch       = cosh(B * eta);
  double sh       = sinh(B * eta);
  double tanh_eta = sh / ch;
  double sech2    = 1./(ch*ch);
  return (A / B) * ((1./3.) * pow(r, -2./3.) * tanh_eta
                   - (2./3.) * (xpos - x_jet) * pow(r, -5./3.) * sech2);
}

static inline void schlichting_constants (double U, double R,
                                          double nu, double *A, double *B)
{
  double J     = U * U * 2. * R;
  double alpha = 0.8255 * pow(J / sqrt(nu), 1./3.);
  *A = 0.68*(2./3.) * alpha * alpha;
  *B = 0.9*alpha / (3. * sqrt(nu));
}

/** ## Properties */

int main ()
{
  Re = 10.;
  L0 = 1.;
  U1 = 1.0;
  U2 = 0.05;
  R1 = 0.005;
  R2 = 0.015;
  d  = 0.36;

  nu_val = U1 * R1 / Re;

  schlichting_constants(U1, R1, nu_val, &A1, &B1);
  schlichting_constants(U2, R2, nu_val, &A2, &B2);

  TOLERANCE = 1e-3 [*];

  u.n[top] = dirichlet( - U1 * (x > -d/2. - R1 && x < -d/2. + R1)
                         - U2 * (x >  d/2. - R2 && x <  d/2. + R2) );
  u.t[top] = dirichlet(0.);

  s[top] = dirichlet(   1. * (x > -d/2. - R1 && x < -d/2. + R1)
                      - 1. * (x >  d/2. - R2 && x <  d/2. + R2) );

  //u.n[bottom] = dirichlet(-U1*2.*R1/L0 - U2*2.*R2/L0);
  u.n[bottom] = u.n[] < 0. ? neumann(0.) : dirichlet(0.);
  p[bottom]   = dirichlet(0.);

  u.n[left]  = u.n[] < 0. ? neumann(0.) : dirichlet(0.);
  u.n[right] = u.n[] > 0. ? neumann(0.) : dirichlet(0.);
  
  N = 256;
  origin(-L0/2., 0.);
  init_grid(N);

  FILE * fparam = fopen("param_dim.txt", "w");
  fprintf(fparam, "%g %g %g %g %g %g %g %d\n", L0, R1, R2, U1, U2, d, nu_val, N);
  fclose(fparam);

  mu = muv;
  fpmax = fopen("log.dat", "w");
  run();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[] * nu_val;
}

/** ## Initial conditions

The Schlichting solution diverges as r^{-1/3} near the nozzle (r -> 0).
We therefore apply the theoretical profile only for r > r_min = 10*max(R1,R2),
and use a simple top-hat inlet profile in the near-field region. */

event init (t = 0)
{
  double x1 = -d / 2.;
  double x2 =  d / 2.;
  double r_min = 1.*max(R1, R2);

  foreach() {
    double r = L0 - y;
    int in1 = (x > x1 - R1 && x < x1 + R1);
    int in2 = (x > x2 - R2 && x < x2 + R2);

    if (r > r_min) {
      u.x[] =  schlichting_u(x, r, x1, A1, B1) + schlichting_u(x, r, x2, A2, B2);
      u.y[] = -(schlichting_v(x, r, x1, A1, B1) + schlichting_v(x, r, x2, A2, B2));
    } else {
      u.x[] = 0.;
      u.y[] = -(U1 * in1 + U2 * in2);
    }
    p[] = 0.;
  }
  boundary({u, p});
}

/** ## Convergence monitoring */
 
double ux_prev[256*256], uy_prev[256*256];
int res_init = 0;
 
event logfile (i++) {
  double res_u = 0.;
  int k = 0;
 
 if (!res_init) {
  foreach() {
    ux_prev[k] = u.x[];
    uy_prev[k] = u.y[];
    k++;
  }
  res_init = 1;
}
else {
  foreach() {
    res_u += sq(u.x[] - ux_prev[k]) + sq(u.y[] - uy_prev[k]);
    k++;
  }

  res_u = sqrt(res_u/(double)k);

  k = 0;
  foreach() {
    ux_prev[k] = u.x[];
    uy_prev[k] = u.y[];
    k++;
  }
}

fprintf(stderr,"%d %g %g\n",i,t,res_u);
fprintf(fpmax,"%d %g %g\n",i,t,res_u);
fflush(fpmax);
}
 
/**
~~~pythonplot Residuals
import matplotlib.pyplot as plt
import numpy as np
 
list_t=[]
list_res_u=[]
 
with open('log.dat', encoding='utf8') as File:
    for line in File:
        if line.isspace()==0:
            temp=line.split()
            list_t.append(float(temp[1]))
            list_res_u.append(float(temp[2]))
    File.close()
 
plt.figure()
plt.semilogy(list_t, list_res_u, 'r-', label=r'$\mathcal{R}^{(n)}(u)$')
plt.xlabel('t')
plt.legend()
plt.savefig('residuals.png')
~~~
*/

/** ## Movies */

event ppm_output (t = 0; t += 0.05; t <= 60) {
  output_ppm(u.y, file = "uY.mp4", n = 512,
             min = -max(U1,U2), max = 0., linear = true);
  output_ppm(s, file = "s.mp4", n = 512,
             min = -1., max = 1., linear = true);
}

event profile (t = end) {
  FILE * fpres = fopen("res_end.txt", "w");
  foreach()
    fprintf(fpres, "%g %g %g %g %g %g\n",
            x, y, u.x[], u.y[], p[], s[]);
  fclose(fpres);
  printf("-----END-----\n");
  fclose(fpmax);
}

/**
![Passive tracer (blue=jet1, red=jet2)](two_jets/s.mp4)
![Vertical velocity](two_jets/uY.mp4)
*/