/** 
# Testing `initial_conditions_dimonte_fft2.h`

In this example we initialize a scalar field `f` using an annular spectrum 
with amplitude 1 and wavenumbers between $k_{min}=25$ and $k_{max}=75$. In this 
case, we use `isvof=1`. The result should look like this:

![Initialized field](test_init_fft3/test_init_fft3.png)
*/

#include "grid/octree.h"
#include "view.h"
#include "initial_conditions_dimonte_fft2.h"

#define r2 (sq(x) + sq(y))
int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 5;
  init_grid(N);

#if TREE
  refine(((z < L0/16) && (z > -L0/16)) && level < 9);
#endif

  scalar f[];
  {
    vertex scalar phi[];
    initial_condition_dimonte_fft2(phi, amplitude=0.01, NX=512, NY=512, kmin = 10, kmax = 50, isvof=1);
    fractions(phi, f);
  }

  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g %g %g\n", t, s.sum, s.min, s.max, s.stddev, s.volume);

  view(camera="iso");
  draw_vof ("f", color="z");
  save("test_init_fft3.png");
}