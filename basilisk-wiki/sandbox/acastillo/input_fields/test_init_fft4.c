/** 
# Testing `initial_conditions_dimonte_fft1.h`

In this example we initialize a scalar field `f` using a 1D spectrum 
with amplitude 0.1 and wavenumbers between $k_{min}=25$ and $k_{max}=75$. In this 
case, we use `isvof=1`. The result should look like this:

![Initialized field](test_init_fft4/test_init_fft4.png)
*/

#include "view.h"
#include "initial_conditions_dimonte_fft1.h"

#define r2 (sq(x) + sq(y))
int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 9;
  init_grid(N);

#if TREE
  refine(((r2 < sq(0.25)) && (r2 > sq(0.1))) && level < 9);
#endif

  scalar f[];
  {
    vertex scalar phi[];
    initial_condition_dimonte_fft(phi, 0.01, kmin = 25, kmax = 75);
    fractions(phi, f);
  }

  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g %g %g\n", t, s.sum, s.min, s.max, s.stddev, s.volume);

  box();
  squares("f", linear = false);
  save("test_init_fft4.png");
}