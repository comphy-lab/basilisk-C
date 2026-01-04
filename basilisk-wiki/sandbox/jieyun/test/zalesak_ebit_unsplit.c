/**
# Zalesak's disk with the EBIT method (unsplit scheme)

In the [Zalesak's disk](#zalesak1979) test case, a notched circular
interface is rotated by a constant velocity field
$$(u, v) = (2 \pi (0.5 - y), 2 \pi (x - 0.5))$$ */

#define UNSPLIT 1
#define ADAPT 1

/**
 The volume fraction (f) is only used for statistic purpose. */

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

#include "advection-ebit.h"
event stability (i++) {} // this ensures that time step is computed before marker advection
#include "ebit-2d.h"

const char *OUTNAME = "zalesak";
const double xcenter = 0.5, ycenter = 0.75, R = 0.15;
double EndT = 1. [0, 1];
double ddt;
int IT, level;

int main() {
  int orders[3] = {1, 2, 4};
  level = 6;
  init_grid (1 << level);
  ddt = 1. [0, 1]/N/16.;
  IT = 16*N;

  for (int i = 0; i < 3; i++) {
    ebit_order = orders[i];
    run();
  }
}

/** Function used for a notched circle. */
double rectangle (double x, double y, coord center, coord size) {
  double P1_Plus = x - size.x/2.  - center.x;
  double P1_Minus = x + size.x/2. - center.x;
  double P1 = max (P1_Plus, -P1_Minus);
  
  double P2_Plus = y - size.y/2.  - center.y;
  double P2_Minus = y + size.y/2. - center.y;
  double P2 = max (P2_Plus, -P2_Minus);
  
  double c = max (P1, P2);
  return c;
}

double circle (double x, double y,  coord center, double radius) {
  double R2 = sq(x - center.x) + sq (y - center.y);
  return sqrt(R2) - radius;
}

double geometry (double x, double y) {
  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = center_rectangle.x = 0.5;
  center_circle.y = 0.75;
  
  center_rectangle.y = 0.725;
  
  size_rectangle.x = 0.05;
  size_rectangle.y = 0.25; 
  
  double s = -circle (x, y, center_circle, 0.15);
  double r = -rectangle (x, y, center_rectangle, size_rectangle);
  
  double zalesak = difference(s, r) ;
  
  return zalesak;
}

event init (i = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = geometry (x,y);

  init_markers (phi);

  semu2vof();
  area0 = area;
}

/** The timestep `dt` and the velocity field are set. */
event stability (i++, i < IT, first) {
  dt = dtnext (ddt);
  double u0 = 2. [1, -1];
  foreach() {
    double xx = x/L0, yy = y/L0;
    coord ur = {u0*pi*(0.5 - yy), u0*pi*(xx - 0.5)};

    foreach_dimension() {
      u.x[] = ur.x;
      up.x[] = ur.x;
      urk1.x[] = ur.x;
      urk2.x[] = ur.x;
    }
  }
  tTime += dt;
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({mask_intf}, (double[]){0.02},\
    maxlevel = level, minlevel = level - 3);
}
#endif

event interface_out (i++, last) {
  if (2*(i + 1) % max(IT, 1) == 0 || i == 0) {
    int ii = 2*(i + 1)/max(IT, 1);
    char name[80];
    sprintf (name, "%s_ebit_%d_%d_order_%d.dat", OUTNAME, N, ii, ebit_order);
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

The shapes of the interface at $t = T$ are displayed below ($N = 64$).

~~~gnuplot Shapes of the interface at $t = T$.
reset
set size ratio -1
plot [0.25:0.75][0.5:1.] '../zalesak_ebit/zalesak_ebit_64_2.dat' w l lw 3 t "Split", \
  'zalesak_ebit_64_2_order_1.dat' w l lw 3 t "Unsplit-Euler", \
  'zalesak_ebit_64_2_order_2.dat' w l lw 3 t "Unsplit-PC", \
  'zalesak_ebit_64_2_order_4.dat' w l lw 1.5 t "Unsplit-RK4", \
  '../zalesak_ana_4.dat' w l dt 2 lt -1 t "Ref."
~~~

## See also

* [Zalesak's disk with the EBIT method (Split scheme)](./zalesak_ebit.c)
* [Zalesak's disk with the VOF method](./zalesak_vof.c)

## References

~~~bib
@article{zalesak1979,
  title={Fully multidimensional flux-corrected transport algorithms for fluids},
  author={Zalesak, Steven T.},
  journal={Journal of Computational Physics},
  volume={31},
  pages={335-362},
  year={1979},
  publisher={Elsevier}
}
~~~
*/
