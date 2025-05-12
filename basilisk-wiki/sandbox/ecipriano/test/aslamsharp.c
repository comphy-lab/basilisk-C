/**
# Aslam Extrapolations

We want to apply [Aslam extrapolations](#aslam2004partial)
on a field which is non-null just at the interface,
described using a VOF field. */

#include "grid/multigrid.h"
#include "utils.h"
#include "aslam.h"
#include "redistance.h"
#include "view.h"

/**
We define the volume fraction field `f`, the function that we want to
extrapolate `m`, and the Heaviside functions `H` which defines the
extrapolation region (just for visualization). */

scalar f[], m[], H[];

int main (void) {

  /**
  We set the domain geometry and we initialize
  the grid and the volume fraction. */

  X0 = -0.5, Y0 = -0.5;
  init_grid (1 << 6);
  fraction (f, (sq(0.4) - sq(x) - sq(y)));

  /**
  We create two copies of the volume fraction fields. */

  scalar fL[], fG[];
  foreach() {
    fL[] = f[];
    fG[] = 1. - f[];
  }

  /**
  We convert the VOF field to a level set. */

  scalar ls1[], ls2[];
  vof_to_ls (f, ls1);
  foreach()
    ls2[] = -ls1[];

  /**
  The function to extrapolate is initialize to `x*y` just at the
  interface. */

  foreach()
    m[] = (f[] > 1.e-10 && f[] < 1.-1.e-10) ? x*y : 0.;

  /**
  We extrapolate `m` from the interface toward the gas phase. */

  clear();
  draw_vof ("f", lw = 1.5);
  squares ("m", spread = -1);
  save ("initial.png");

  /**
  We extrapolate `m` from the interface toward the liquid phase. */

  constant_extrapolation (m, ls1, 0.5, 10, c=fL, nl=0,
      nointerface=true);

  clear();
  draw_vof ("f", lw = 1.5);
  squares ("m", spread = -1);
  save ("extrapolated1.png");
  fprintf (stderr, "constant = %g\n", statsf(m).sum);

  /**
  We extrapolate `m` from the interface toward the liquid phase. */

  foreach()
    m[] = (f[] > 1.e-10 && f[] < 1.-1.e-10) ? x*y : 0.;

  constant_extrapolation (m, ls2, 0.5, 10, c=fG, nl=0,
      nointerface=true);

  clear();
  draw_vof ("f", lw = 1.5);
  squares ("m", spread = -1);
  save ("extrapolated2.png");
  fprintf (stderr, "linear = %g\n", statsf(m).sum);

  /**
  We print the Heaviside function used for these extrapolations. This
  function is non-null in the region where `m` must be extrapolated.  */

  mapregion (H, fL, nl=0, nointerface=true);

  clear();
  draw_vof ("f", lw = 1.5);
  squares ("H", spread = -1);
  save ("heaviside.png");
  fprintf (stderr, "heaviside = %g\n", statsf(H).sum);
}

/**
## Results

![Initial field](aslamsharp/initial.png)

![Extrapolation toward gas](aslamsharp/extrapolated1.png)

![Extrapolation toward liquid](aslamsharp/extrapolated2.png)

![Heaviside function without interface](aslamsharp/heaviside.png)

## References

~~~bib
@article{aslam2004partial,
  title={A partial differential equation approach to multidimensional extrapolation},
  author={Aslam, Tariq D},
  journal={Journal of Computational Physics},
  volume={193},
  number={1},
  pages={349--355},
  year={2004},
  publisher={Elsevier}
}
~~~
*/

