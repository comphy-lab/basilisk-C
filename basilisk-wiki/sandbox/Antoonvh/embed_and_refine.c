/**
# Embed + refine

`Refine_embed_linear()` may trigger an assertion
*/
#include "embed.h"
#include "run.h"

#define GAUSS (exp(-(sq(x))) - y)

scalar s[];
int main() {
  L0 = 100.;
  X0 = -L0/2.;
  init_grid (1 << 8);
  run();
}

event init (t = 0) {
  s.refine = refine_embed_linear;
  refine (fabs(GAUSS) < 0.1 && level < 9); //This goes well
  foreach()
    printf ("%g %g\n", x, y);
  fflush (stdout);
  /** 
      Refine and log the progress untill it crashes
*/
  FILE * fp2 = fopen ("iterated", "w");
  refine (fabs(GAUSS) < 0.1 && level < 10 &&
	  fprintf(fp2, "%g %g\n", x, y) && !fflush(fp2));//This goes wrong
}

/**


   
~~~gnuplot A wierd pattern?
set xr [-1.5:3.5]
set yr [-2 : 3]
set size square
plot 'out' t 'all', 'iterated'
~~~
 */
