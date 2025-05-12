/**
# Check compatibility of embedded boundaries, axi and VOF */

#include "embed.h"
#include "axi.h"
#define ro 1.2
#define veloc(r) (1./(r))
#define fpos(t) sqrt(sq(ro) + 2.*(t))

#include "advection.h"
#include "vof.h"
#include "curvature.h"

scalar f[], * tracers = NULL, * interfaces = {f};

int main()
{
  Y0 = 1.0;
  L0 = 1.0 [0];
  DT = HUGE [0];
  N = 64;
  run();
}

event init (i = 0)
{
#if EMBED  
  foreach()
    cs[] = 1;
  foreach_face()
    fs.x[] = 1;
#endif // EMBED
  
  fraction (f, ro - y);
  foreach_face(y)
    u.y[] = veloc(y)*fm.y[];
}

/**
## Results

We just compare the analytical and numerical results... */

event prof_pos (i += 5; t <= 0.8) {
  scalar pos[];
  position (f, pos, {0,1.});
  fprintf (stderr, "%g %g %g\n", t, statsf(pos).max, fpos(t));
}

/**
~~~gnuplot Interface position
set xlabel 't'
set ylabel 'interface position'
set key top left
plot 'log' u 1:2  t 'Numerical', 'log' u 1:3 w line t 'Theoretical'
~~~
*/
