/**
# Single Vortex with the EBIT method (unsplit scheme)

The single vortex test with a highly streched and deformed interface
was proposed by [Rider, 1998](#rider1998). A divergence-free velocity field
$(u, v) = (\partial \phi \big/ \partial y, -\partial \phi \big/ \partial x)$ described
by the stream function
$$\phi = \pi^{-1} \sin^2(\pi x) \sin^2(\pi y) \cos(\pi t / T)$$
is imposed. */

#define UNSPLIT 1
#define ADAPT 1

/**
 The volume fraction (f) is only used for statistic purpose. */

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

#include "advection-ebit.h"
event stability (i++) {} // this ensures that time step is computed before marker advection
#include "ebit-2d.h"

const char *OUTNAME = "vortex";
const double xcenter = 0.5, ycenter = 0.75, R = 0.15;
const double Cfl = 0.125 [0];
double EndT = 2. [0, 1];
double ddt;
int IT, level;

int main() {
  int orders[3] = {1, 2, 4};
  EndT = 8.;
  level = 7;
  init_grid (1 << level);
  ddt = 1. [0, 1]*Cfl/N;
  IT = (int) (EndT*N/Cfl);

  for (int i = 0; i < 3; i++) {
    ebit_order = orders[i];
    run();
  }
}

event init (i = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(R) - (sq(x - xcenter) + sq(y - ycenter));

  init_markers (phi);
  init_circle (xcenter, ycenter, R, f, s);

  semu2vof();
  area0 = area;
}

/** The timestep `dt` and the velocity field are set. */
event stability (i++, i < IT, first) {
  dt = dtnext (ddt);

  double u0 = 1. [1, -1], cdt, cdtp1, cdtp2;
  cdt = u0*cos(pi*tTime/EndT);
  cdtp1 = u0*cos(pi*(tTime + 0.5*dt)/EndT);
  cdtp2 = u0*cos(pi*(tTime + dt)/EndT);

  foreach() {
    double x0 = x/L0, y0 = y/L0;
    coord uc = {sq(sin(x0*pi))*sin(2.*y0*pi), -sin(2.*x0*pi)*sq(sin(y0*pi))};

    foreach_dimension() {
      u.x[] = cdt*uc.x;
      up.x[] = cdtp2*uc.x;
      urk1.x[] = cdtp1*uc.x;
      urk2.x[] = cdtp2*uc.x;
    }
  }

  tTime += dt;
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({mask_intf}, (double[]){0.02}, \
    maxlevel = level, minlevel = level - 4);
}
#endif

event interface_out (i++, last) {
  if (2*(i + 1) % max(IT, 1) == 0 || i == 0) {
    int ii = 2*(i + 1)/max(IT, 1);
    char name[80];
    sprintf (name, "%s_ebit_%d_%d_%d_order_%d.dat",\
      OUTNAME, N, (int) EndT, ii, ebit_order);
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

  // reference file
  output_facets_ebit ("", stderr);
}

/**
## Results

The shapes of the interface at $t = T/2$ and
$t = T$ are displayed below ($T = 8, N = 128$).

~~~gnuplot Shapes of the interface at $t = T/2$.
reset
set size ratio -1
plot [0.:1.][0.:1.]'../vortex_ebit/vortex_ebit_128_8_1.dat' w l lw 3 t "Split", \
  'vortex_ebit_128_8_1_order_1.dat' w l lw 3 t "Unsplit-Euler", \
  'vortex_ebit_128_8_1_order_2.dat' w l lw 3 t "Unsplit-PC", \
  'vortex_ebit_128_8_1_order_4.dat' w l lw 1.5 t "Unsplit-RK4", \
  '../vortex_ana_8_4.dat' w l dt 2 lt -1 t "Ref."
~~~

~~~gnuplot Shapes of the interface at $t = T$.
reset
set size ratio -1
plot [0.25:0.75][0.5:1.]'../vortex_ebit/vortex_ebit_128_8_2.dat' w l lw 3 t "Split", \
  'vortex_ebit_128_8_2_order_1.dat' w l lw 3 t "Unsplit-Euler", \
  'vortex_ebit_128_8_2_order_2.dat' w l lw 3 t "Unsplit-PC", \
  'vortex_ebit_128_8_2_order_4.dat' w l lw 1.5 t "Unsplit-RK4", \
  '../vortex_ana_8_8.dat' w l dt 2 lt -1 t "Ref."
~~~

## See also

* [Single vortex with the EBIT method (Split scheme)](./vortex_ebit.c)
* [Single vortex with the VOF method](./vortex_vof.c)

## References

~~~bib
@article{rider1998,
  title={Reconstructing volume tracking},
  author={Rider, William J. and Kothe, Douglas B.},
  journal={Journal of Computational Physics},
  volume={141},
  pages={112-152},
  year={1998},
  publisher={Elsevier}
}
~~~
*/
