/**
# Three-dimensional bubble merging

This setup was adapted from [Unverdi](#Unverdi_1992_100).
See also [Sussman](#Sussman_2000_162) for the same setup.

The flow is characterized by the Morton number (M) and Eotvos number
(Eo, or Bond number Bo) 
$$ \mathrm{Eo} = \rho_o g D^2 / \sigma $$
$$ \mathrm{M} = g \mu_o^4 / \rho_o / \sigma^3 $$
The Galileo number (Ga) is an alternative of Morton number
$$ \mathrm{Ga} = \rho_o^2 g D^3 / \mu_o^2 = (\mathrm{Eo}^3/\mathrm{M})^{1/2}$$. */

face vector av[];

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "two-phase-ebit.h"
#include "tension.h"
#include "tag.h"

/**
Free-slip boundary conditions are imposed at all domain boundaries. */

uf.n[top] = 0.;
uf.n[bottom] = 0.;
uf.n[front] = 0.;
uf.n[back] = 0.;

int level = 6;
int case_id = 1;
double R = 0.25, xry = 0.5, width = 2. [1], R2;
double eo = 1.[0], mo = 1.[0], r_rho = 1.[0], r_mu = 1.[0];
double gra = 0.98;
FILE *fp_vel;

int main (int argc, char * argv[]) {
  if (argc > 1)
    level = atoi (argv[1]);

  DT = 1.e-3 [0,1];
  TOLERANCE = 1e-6 [*];
  rho1 = 1000. [-3,0,1];
  dimensions(nx = 2, ny = 1, nz = 1);

  eo = 50.;
  mo = 1.;
  r_rho = 20.;
  r_mu = 26.;

  f.sigma = rho1*gra*sq(2.*R)/eo;
  mu1 = sqrt(mo*rho1*cube(f.sigma)/gra);
  mu1 = sqrt(mu1);
  rho2 = rho1/r_rho;
  mu2 = mu1/r_mu;
  
  xry = 0.5 [1];

  size (width);
  init_grid (1 << level);
  a = av;

  run();
}

event init (i = 0) {
  R2 = R;
  coord xc1 = {0.5, 0.375, 0.5}, xc2 = {1.075, 0.625, 0.5};

  vertex scalar phi[];
  foreach_vertex() {
    double phi1, phi2;
    phi1 = sq(x - xc1.x) + sq(y - xc1.y) + sq(z - xc1.z) - sq(R);
    phi2 = sq(x - xc2.x) + sq(y - xc2.y) + sq(z - xc2.z) - sq(R2);
    phi[] = min(phi1, phi2);
  }

  init_markers (phi);

  eo = rho1*gra*sq(2.*R)/f.sigma;
  mo = gra*sq(mu1)*sq(mu1)/rho1/cube(f.sigma);
  if (pid() == 0) {
    printf ("Physical properties:\n");
    printf ("Density: %.4e %.4e, viscosity: %.4e %.4e\n", rho1, rho2, mu1, mu2);
    printf ("Surface tension: %.4e, gravity:%.4e\n", f.sigma, gra);
    printf ("Morton num.: %g  Eotvos num.: %g\n", mo, eo);
  }
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  foreach_face(x)
    av.x[] -= gra;
  boundary ((scalar *) {av});
}

/**
We log the position of the center of mass of the bubbles and their velocities. */

bool has_merged = false;

event logfile (i+=4) {
  scalar d[];
  double yb = 0., vb = 0., sb = 0.;
  double yb_s[2] = {0., 0.}, vb_s[2] = {0., 0.}, sb_s[2] = {0., 0.};
  if (!has_merged) {
    foreach() {
      if (f[] < 1.e-4)
        d[] = 1.;
      else
        d[] = 0.;
    }
    int n = tag (d);

    if (n != 2)
      has_merged = true;
  }

  if (has_merged) {
    foreach(reduction(+: yb) reduction(+: vb) reduction(+: sb)) {
      double dv = (1. - f[])*dv();
      yb += x*dv;
      vb += u.x[]*dv;
      sb += dv;
    }
    yb_s[0] = yb_s[1] = yb;
    vb_s[0] = vb_s[1] = vb;
    sb_s[0] = sb_s[1] = sb;
  }
  else {
    foreach(reduction(+: yb) reduction(+: vb) reduction(+: sb)
      reduction(+: yb_s[:2]) reduction(+: vb_s[:2]) reduction(+: sb_s[:2])) {
      double dv = (1. - f[])*dv();
      yb += x*dv;
      vb += u.x[]*dv;
      sb += dv;
      int ind = (int) d[] - 1;
      if (ind >= 0) {
        yb_s[ind] += x*dv;
        vb_s[ind] += u.x[]*dv;
        sb_s[ind] += dv;
      }
    }
  }

  if (pid() == 0) {
    fprintf (stderr, "%.8e %.8e %.8e %.8e %.8e\n",
      t, yb_s[0]/sb_s[0], vb_s[0]/sb_s[0], yb_s[1]/sb_s[1], vb_s[1]/sb_s[1]);
  }
}

/**
Output the shape of the bubble at diffrent time instants. */

event interface (t = {0., 0.5, 1., 1.5, 1.7, 2.}) {
  // char name[80], name_con[80];
  // if (pid() == 0)
  //   printf ("Save interface (mesh) at time step:%d\n", i);

  // sprintf (name, "%srising_merge_3d_ebit_case%d_%d_%.1f.dat", \
  //   OUTPUTPATH, case_id, N, t);
  // output_facets_semushin (name);

  // sprintf (name, "%srising_merge_3d_ebit_case%d_%d_%.1f_fem.dat", \
  //   OUTPUTPATH, case_id, N, t);
  // sprintf (name_con, "%srising_merge_3d_ebit_case%d_%d_%.1f_con.dat", \
  //   OUTPUTPATH, case_id, N, t);
  // output_facets_semushin (name, tri = true, file_con = name_con);
}

/**
## Results

~~~gnuplot Rise velocity as a function of time for test case 1.
set term pop
reset
set grid
set xlabel 't'
set key bottom right
plot [0:2.][0:] 'log' u 1:3 w l t 'Bubble 1', 'log' u 1:5 w l dt 2 t 'Bubble 2'
~~~

## References

~~~bib
@article{Sussman_2000_162,
title = {A Coupled Level Set and Volume-of-Fluid Method for Computing {3D} and Axisymmetric Incompressible Two-Phase Flows},
journal = {Journal of Computational Physics},
volume = {162},
number = {2},
pages = {301-337},
year = {2000},
doi = {https://doi.org/10.1006/jcph.2000.6537},
author = {Mark Sussman and Elbridge Gerry Puckett}
}

@article{Unverdi_1992_100,
  title={A front-tracking method for viscous, incompressible, multi-fluid flows},
  author={Unverdi, Salih Ozen and Tryggvason, Gr\'{e}tar},
  journal={Journal of Computational Physics},
  volume={100},
  pages={25-37},
  year={1992},
  publisher={Elsevier}
}
~~~
*/
