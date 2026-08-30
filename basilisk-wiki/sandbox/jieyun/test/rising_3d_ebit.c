/**
# Three-dimensional rising bubble

This setup was adapted from [Unverdi](#Unverdi_1992_100).
See also [Shin](#Shin_2002_180) for different setups.

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

/**
Free-slip boundary conditions are imposed at all domain boundaries. */

uf.n[top] = 0.;
uf.n[bottom] = 0.;
uf.n[front] = 0.;
uf.n[back] = 0.;

int level = 6;
int case_id = 1;
double R = 0.25, width = 2. [1];
double eo, mo;
double gra = 0.98;

int main(int argc, char * argv[]) {
  if (argc > 1)
    level = atoi (argv[1]);

  DT = 1.e-3 [0, 1];
  TOLERANCE = 1e-4 [*];
  rho1 = 1000. [-3,0,1];
  dimensions(nx = 2, ny = 1, nz = 1);

  if (case_id == 1) {
    // Eo: 10, M: 0.1
    // mu1 = 35. * sqrt(0.1);
    mu1 = 35.;
    rho2 = 25.;
    mu2 = mu1/150.;
    f.sigma = 24.5;
  }
  else if (case_id == 2) {
    // Eo: 100, M: 1
    mu1 = 11.07;
    rho2 = 25.;
    mu2 = mu1/50.;
    f.sigma = 2.45;
  }

  size (width);
  init_grid (1 << level);
  a = av;

  run();
}

event init (i = 0) {
  coord xc = {0.5, 0.5, 0.5};

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x - xc.x) + sq(y - xc.y) + sq(z - xc.z) - sq(R);

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
We log the position of the center of mass of the bubble, its velocity
and volume. */

event logfile (i++) {
  double yb = 0., vb = 0., sb = 0.;
  foreach(reduction(+: yb) reduction(+: vb) reduction(+: sb)) {
    double dv = (1. - f[])*dv();
    yb += x*dv;
    vb += u.x[]*dv;
    sb += dv;
  }

  if (pid() == 0) {
    fprintf (stderr, "%.8e %.8e %.8e %.8e\n", t, sb, yb/sb, vb/sb);
  }
}

/**
Output the shape of the bubble at diffrent time instants. */

event interface (t = {0., 1., 1.5}) {
  // char name[80], name_con[80];
  // if (pid() == 0)
  //   printf ("Save interface (mesh) at time step:%d\n", i);

  // sprintf (name, "%srising_3d_ebit_case%d_%d_%.1f.dat",
  //   OUTPUTPATH, case_id, N, t);
  // output_facets_semushin (name);

  // sprintf (name, "%srising_3d_ebit_case%d_%d_%.1f_fem.dat",
  //   OUTPUTPATH, case_id, N, t);
  // sprintf (name_con, "%srising_3d_ebit_case%d_%d_%.1f_con.dat",
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
plot [0:1.5][0:] 'log' u 1:4 w l t 'EBIT'
~~~

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
