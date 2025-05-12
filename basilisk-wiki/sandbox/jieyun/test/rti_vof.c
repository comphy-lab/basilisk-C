/**
# Rayleigh-Taylor instability with the VOF method

The [Rayleigh-Taylor instability, RTI](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Taylor_instability)
occurs when a heavy fluid is on top of a lighter one.
 */

#define LEVEL 8
#define ADAPT 1

face vector av[];

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

const double MEANPOS = 2.00;
const double DY = 0.1;
const double lref = 1. [1];
double gra, Reynolds;
FILE *fp_amp;

/** 
This setup with 0 surface tension was adapted from 
[Tryggvason](#tryggvason1988), which has been widely investigated in
several studies. The flow can be charaterized by two dimensionless numbers:
Atwood number $At$ and Reynolds number $Re$
$$ At = \frac{\rho_1 - \rho_2}{\rho_1 + \rho_2} = 3 $$
$$ Re = \frac{\rho_1 g^{1/2} d^{3/2}}{\mu_1}  = 3000, d = 1 $$
*/

int main() {
  Reynolds = 3000. [0];
  rho1 = 3. [0];
  rho2 = 1.;
  f.sigma = 0.;
  gra = 9.81 [1, -2];
  mu1 = rho1*sqrt(gra)/Reynolds*sqrt(cube(lref));
  mu2 = mu1;

  CFL = 0.05;
  DT = 2.e-4 [0, 1];
  TOLERANCE = 1e-4 [*];

  size (4. [1]);
  init_grid (1 << LEVEL);

  a = av;
  
  run();
  fclose (fp_amp);
}

/** 
The initial interface shape is
$$ y(x) = 2d + 0.1 d \cos(kx), \quad k = \frac{2 \pi}{d}$$
*/

event init (i = 0) {
  u.t[bottom] = dirichlet(0);
  u.t[top] = dirichlet(0);

  uf.n[right] = 0.;
  uf.n[left] = 0.;

  mask (x > 1. ? right : x < 0. ? left : none);

  vertex scalar phi[];
  foreach_vertex()
    phi[] = -(MEANPOS - y + DY*cos(2. [-1]*pi*x));
  
  fractions (phi, f);

  char name[80];
  sprintf (name, "rt_vof_dis.dat");
  fp_amp = fopen (name, "w");
}

/**
Output the time evolution of the hightest and lowest positions of the interface.
Time is normalized by reference time $\tau = t / t_{ref}$.

$$ t_{ref} = \sqrt{d / At g}$$
*/
const double TREF = 0.451523641;
event amplitude (i++) {
  scalar pos[];
  position (f, pos, {0, 1});
  stats pos_stats = statsf(pos);
  double ymax = pos_stats.max, ymin = pos_stats.min;

  fprintf (fp_amp, "%g %g %g\n", t/TREF, ymax - MEANPOS, ymin - MEANPOS);
  fflush (fp_amp);
}

/** The vertical acceleration is added here. */

event acceleration (i++) {
  foreach_face(y)
    av.y[] -= gra;
  boundary ((scalar *){av});
}

/** Ouput the interfaces at different time instants. */
event interface (t = {1.*TREF, 1.5*TREF, 1.75*TREF, \
  2.*TREF, 2.25*TREF, 2.5*TREF}) {
  char name[80];

  sprintf (name, "rt_intf_vof_%.2f.dat", t/TREF);

  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fflush (fp);
  fclose (fp);
}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){5e-4, 1e-3, 1e-3}, maxlevel = LEVEL, minlevel = LEVEL - 3);
}
#endif

/**
## Results

Evolution of the hightest and lowest positions of the interface as
a function of dimensionless time $\tau$

~~~gnuplot Time evolution of interface position.
reset
set grid
set xlabel 'tau'
set ylabel 'Amplitude'
set key bottom left
plot [0.:2.5][-2.:1.] 'rt_vof_dis.dat' u 1:2 w l lw 3 t "EBIT, upper", \
  'rt_vof_dis.dat' u 1:3 w l lw 3 t "EBIT, lower"
~~~

~~~gnuplot Shapes of the interface at $\tau = 1$.
reset
set size ratio -1
plot [0.:1.][0.5:3.] 'rt_intf_vof_1.00.dat' w l t "tau = 1.00"
~~~

~~~gnuplot Shapes of the interface at $\tau = 1.5$.
plot [0.:1.][0.5:3.] 'rt_intf_vof_1.50.dat' w l t "tau = 1.50"
~~~

~~~gnuplot Shapes of the interface at $\tau = 1.75$.
plot [0.:1.][0.5:3.] 'rt_intf_vof_1.75.dat' w l t "tau = 1.75"
~~~

~~~gnuplot Shapes of the interface at $\tau = 2.0$.
plot [0.:1.][0.5:3.] 'rt_intf_vof_2.00.dat' w l t "tau = 2.00"
~~~

~~~gnuplot Shapes of the interface at $\tau = 2.25$.
plot [0.:1.][0.5:3.] 'rt_intf_vof_2.25.dat' w l t "tau = 2.25"
~~~

~~~gnuplot Shapes of the interface at $\tau = 2.5$.
plot [0.:1.][0.5:3.] 'rt_intf_vof_2.50.dat' w l t "tau = 2.50"
~~~

## See also

* [Rayleigh-Taylor instability with the EBIT method](./rti_ebit.c)

## References

~~~bib
@article{tryggvason1988,
  title={Numerical simulations of the Rayleigh-Taylor instability},
  author={Tryggvason, Gr\'{e}tar},
  journal={Journal of Computational Physics},
  volume={75},
  pages={253-282},
  year={1988},
  publisher={Elsevier}
}
~~~
*/
