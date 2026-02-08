/**
# 3D drop oscillation

This setup was adapted from [Shin](#Shin_2002_180).
According to Lamb's analytic solution, the oscillation frequency
of an inviscid droplet is
$$ \omega_n^2 = \frac{n (n + 1) (n - 1) (n + 2) \sigma}{[(n + 1) \rho_1 + n \rho_2] r^3} $$
and amplitude of a viscous droplet would decay as
$$ a_n(t) = a_0 \exp^{-t/\tau}, \tau = \frac{r_0^2}{(n - 1)(2n + 1)\nu}$$ */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase-ebit.h"
#include "tension.h"

uf.n[bottom] = 0.;
uf.n[top] = 0.;
uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[back] = 0.;
uf.n[front] = 0.;

int level = 5, maxlvl, minlvl;
const double R = 1.;
const double EPSR = 0.005;
double la = 1. [0];

FILE * fp;
double tp_max = 1. [0,1];
double t_step = 1. [0,1];
int main(int argc, char * argv[]) {
  if (argc > 1)
    level = atoi (argv[1]);
  maxlvl = level;
  minlvl = max(level - 5, 3);

  rho1 = 10. [-3,0,1];
  rho2 = 0.1;
  mu1 = 5.e-2;
  mu2 = 5.e-4;

  f.sigma = 10.;
  // f.sigma = 0.1;

  int no = 2;
  double omega2 = no*(no-1)*(no+1)*(no+2)*f.sigma/((no+1)*rho1 + no*rho2)/cube(R);
  double omega = sqrt(omega2);
  la = 2.*rho1*R*f.sigma/sq(mu1);

  tp_max = 4.*2.*pi/omega;

  TOLERANCE = 1e-6 [*];
  CFL = 0.1;

  size (4. [1]);
  origin (-L0/2., -L0/2., -L0/2.);
  init_grid (1 << maxlvl);

  DT = 2.*pi/omega/(1 << 10);
  t_step = DT;
  run();

  fclose(fp);
}

/**
The initial radial position of the droplet interface is
$$ r(\theta) = r_0 + \epsilon P_n(\cos \theta)$$
where $P_n$ is the nth-order Legendre polynomial,
and the second-order mode is investigated here. */

event init (i = 0) {
  vertex scalar phi[];
  foreach_vertex() {
    double cth, sth, dr2, rth, pn;
    dr2 = sq(z) + sq(y) + sq(x);
    cth = z/(sqrt(dr2) + 1.e-32);
    pn = 0.5*(3.*sq(cth) - 1.);
    rth = R + EPSR*pn;
    phi[] = rth - sqrt(dr2);
  }

  char name[80];
  int npr = (int) ((1 << maxlvl)*R/L0);

  char method_name[] = "EBIT";
  init_markers (phi);
  sprintf (name, "%sdroplet_3d_dr_ebit_%d_%d_La%d.dat",
    OUTPUTPATH, npr, (int) L0, (int) ceil(la));

  fp = fopen (name, "w");

  if (pid() == 0)
    printf ("%s R:%g, Sigma:%g, La:%g DT:%.5e\n", method_name, R, f.sigma, la, DT);
}

event logfile (t += t_step; t <= tp_max) {
  double zmax = 0., zmin = HUGE;
  double ke = 0., ke_total = 0.;
  foreach (reduction(+:ke) reduction(+:ke_total)) {
    ke += cube(Delta)*(sq(u.x[]) + sq(u.y[]) + sq(u.z[]))*rho(f[])*f[];
    ke_total += cube(Delta)*(sq(u.x[]) + sq(u.y[]) + sq(u.z[]))*rho(f[]);
  }

  foreach_vertex (reduction(max:zmax) reduction(min:zmin)) {
    if (with_marker.z[] > 0) {
      double zm = z + s.z[]*Delta;
      zmax = max(zmax, zm);
      zmin = min(zmin, zm);
    }
  }
  zmax -= R;
  zmin += R;

  fprintf (fp, "%.8e %.12e %.12e %.12e %.12e\n", t, zmax, zmin, volume, ke);

  fflush (fp);
}

event monitor (i++, last) {
  if (pid() == 0 && i % 1000 == 0) {
    printf ("i: %d Time: %g\n", i, t);
  }
}

/**
## References

~~~bib
@article{Shin_2002_180,
  title={Modeling Three-Dimensional Multiphase Flow Using a Level Contour Reconstruction
    Method for Front Tracking without Connectivity},
  author={Shin, Seungwon and Juric, Damir},
  journal={Journal of Computational Physics},
  volume={180},
  pages={427-470},
  year={2002},
  publisher={Elsevier}
}
~~~
*/