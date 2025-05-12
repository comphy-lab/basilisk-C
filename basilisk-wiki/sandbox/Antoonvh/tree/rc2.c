/**
# Enforcing refinement based on a function of the radius

   This page looks a bit like [this one](../rc.c).
*/
#include "fractions.h"
#include "curvature.h"
#include "view.h"

/**
Say, for the viscous boundary layer in a flow around a cylinder,
$\Delta_{min} = C R \times \mathrm{Re}^{-0.5} = C \sqrt{\frac{\nu R}{U}}$,
with dimensionless tuning constant $C \approx 1.25$ ([see
here:](../cyl2d.c)), altough $C = 0.2$ has also been proposed
[here](/src/test/starting.c). 
*/

double C = 1.25; 
double U = 1, nu = 1./500.;
/**
`prolongate_ratio()` is defined for this purpose only.
 */
void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}
/**
   We test the method with a circle of radius $R$.
 */

#define CIRC(R) (sq(x - R + 1.01) + sq(y) - sq(R))

scalar f[], res[];

int main() {
  f.prolongation = f.refine = fraction_refine; //f is the fraction field...
  L0 = 4;
  X0 = Y0 = -L0/2;
  init_grid (1 << 6);
  view (width = 700, height = 700);
  //Loop over radii
  for (double R = 2; R > 0.05; R /= 1.01) {
    fraction (f, CIRC(R));                 //Compute drop
    boundary ({f});
    double davg = 0;
    int i = 0;
    foreach() {
      res[] = nodata;
      if (f[] > 0 && f[] < 1) {
	res[] = U/sqrt(R*nu);
	davg += Delta;
	i++;
      }
    }
    printf ("%g %g\n", R, davg/(double)i);
    res.prolongation = prolongate_ratio;
    
    adapt_wavelet ({res}, (double[]){C}, maxlevel = 99);
    //Make a movie
    draw_vof ("f", lc = {1., 0., 1.}, lw = 2);
    cells();
    save ("plot.mp4");
  }
}

/**
## Results;

~~~gnuplot It worked
set key bottom right
set log x
set log y 
plot 'out' t 'size', sqrt (x/500) * 1.25 t 'request'
~~~

![movie](rc2/plot.mp4)

*/
