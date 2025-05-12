/**
# Another assertion by `refine_embed_linear`

Apart from [this example](embed_and_refine.c), an assertion can be
triggered when running this script in parallel with 4 threads... 

It does *not* trigger when 

* Running with a quadtree 
* Running with 1,2,3 or 5 threads 
* doing `init_grid(N)` and `refine()` in the `main()` function (no call to `run()`)
* When using `adapt_wavelet()` iteratively to refine.

It seems the `refine()` function allows non-proper trees? 
*/

#include "grid/octree.h" 
#include "embed.h"
#include "run.h"

scalar s[];
int maxlevel = 7;

double M_top = 5;
double base_height = 10;
#define MOUNTAIN (base_height + M_top*exp(-((sq(x + L0/4.)/sq(10.)) + sq(z)/sq(30.))) - y)

int main() {
  s.refine = refine_embed_linear;
  L0 = 100;
  X0 = Z0 = -L0/2.;
  run();
}

event init (t = 0) {
  refine (fabs(MOUNTAIN) < 2.5 && level < maxlevel);
}
