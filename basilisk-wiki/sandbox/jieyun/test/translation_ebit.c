/**
# Translation of a circular interface with the EBIT method

A uniform velocity field $(u, v) = (1, -1)$ is imposed and the
flow direction is reversed at $t = T/2$. */

#define ADAPT 1

/**
 The volume fraction (f) is only used for statistic purpose. */

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

#include "advection-ebit.h"
event stability (i++) {} // this ensures that time step is computed before marker advection
#include "ebit-2d.h"

const char *OUTNAME = "translation";
const double xcenter = 0.25, ycenter = 0.75, R = 0.15;
const double Cfl = 0.125 [0];
double EndT = 1. [0, 1];
double ddt;
int IT, level;

int main() {
  for (level = 5; level < 8; level++) {
    init_grid (1 << level);
    ddt = 1. [0, 1]*Cfl/N;
    IT = (int) (EndT*N/Cfl);
    run();
  }
}

event init (i = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(R) - (sq(x - xcenter) + sq(y - ycenter));

  init_markers (phi);
  init_circle (xcenter, ycenter, R, f, s, itmax = 20);

  semu2vof();
  area0 = area;
}

/** The timestep `dt` and the velocity field are set. */
event stability (i++, i < IT, first) {
  dt = dtnext (ddt);

  coord dir = {1., -1.};
  double reversed = (i >= N/Cfl/2.) ? -1. [0]: 1. [0];
  foreach()
    foreach_dimension()
      u.x[] = dir.x*reversed;

  tTime += dt;
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({mask_intf}, (double[]){0.02}, maxlevel = level, minlevel = level - 3);
}
#endif

event interface_out (i++, last) {
  if (2*(i + 1) % max(IT, 1) == 0 || i == 0) {
    int ii = 2*(i + 1)/max(IT, 1);
    char name[80];
    sprintf (name, "%s_ebit_%d_%d.dat", OUTNAME, N, ii);
    output_facets_ebit (name);
  }
}

/**
We can compute the shape error ($E_{shape}$) and area error ($E_{area}$).

$$ E_{shape}=\max_{i}| \mathrm{dist} (\boldsymbol{x}_i)| $$.
$$ \mathrm{dist}(\boldsymbol{x}_i)=\sqrt{(x_i - x_c)^2 + (y_i - y_c)^2} - R$$
where the reference solution is a circle centered in $(x_c,y_c)$ and
with radius $R$.

$$E_{area} = (A(T) - A(0)) / A(0)$$.
*/
event calc_infty_norm (t = end) {
  double l_inf = 0.;
  coord dir = {0., 1.};

  foreach_face(reduction(max:l_inf)) {
    if (with_marker.x[] > 1.e-6) {
      double ss = (s.x[] - 0.5)*Delta, xx, yy;
      xx = x + ss*dir.x;
      yy = y + ss*dir.y;
      double dist = fabs(sqrt(sq(xx - xcenter) + sq(yy - ycenter)) - R);
      if (dist > l_inf ) l_inf = dist;
    }
  }

  // shape error and area error
  printf ("%d %e %e %e %e\n", N, area0, area, fabs(area0 - area)/area0, l_inf);
}

/**
## Results

The shapes of the interface at $t = T/2$ and
$t = T$ are displayed below.

~~~gnuplot Shapes of the interface ($N = 128$).
reset
set size ratio -1
plot [0.:1.][0.:1.]'translation_ebit_128_1.dat' w l lw 3 t "EBIT, t = T/2", \
  'translation_ebit_128_2.dat' w l lw 3 t "EBIT, t = T"
~~~

## See also

* [Translation with the VOF method](./translation_vof.c)

*/
