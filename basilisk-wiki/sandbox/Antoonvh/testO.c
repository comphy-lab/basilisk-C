/**
We test if optimzation is actually turned on by the [Makefile](Makefile) by checking for output:

![It is here. the sandbox server is just slower than I expected](testO/s.png)

*/
#include "utils.h"
scalar s[];
int main() {
  system ("rm s.png s2.png");
#ifdef __OPTIMIZE__
  init_grid (256);
  foreach()
    s[] = x + sq(y);
  output_ppm (s, file = "s.png");
#else
  init_grid (128);
  foreach()
    s[] = sq(x) + y;
  output_ppm (s, file = "s2.png");
#endif
}
/**
The non `__OPTIMIZE__`d output does not exist:

![see?](testO/s2.png)
*/
