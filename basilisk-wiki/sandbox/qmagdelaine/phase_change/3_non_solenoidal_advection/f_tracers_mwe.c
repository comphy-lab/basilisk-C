/**
# f.tracers errors - minimum working example

In order to handle binary mixture evaporation, we want to use tracers associated
to $f$. This allows the transport of a tracer in a phase without any exchange
with the other one, as it is done in [momentum.h](/src/momentum.h). */

#include "navier-stokes/centered.h"
#include "vof.h"

/**
We allocate several scalar and vector fields to describe both the
interface and the concentration field. */

scalar f[], * interfaces = {f};
face vector uf_save[];
scalar concentration[];

const double cs = 1.0;
const double tend = 100.;
const double L = 2.; // size of the box
int LEVEL = 7;

/**
The main function of the program: */

int main()
{
  size (2.*L);
  origin (- L, - L);
  N = 1 << LEVEL;
  init_grid (N);
  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x,y) (sq(1.) - sq(x) -sq(y))

event init (t = 0) {
  fraction (f, circle(x,y));
  foreach() {
    concentration[] = f[]*cs;
  }
  boundary({concentration});
}

/**
We overload the *vof()* event to transport consistently the two
concentrations associated to each phase. */

static scalar * interfaces1 = NULL;

event vof (i++) {

  /**
  f.tracers: *concentration* is associated to $f$. */

  concentration.inverse = false;
  concentration.gradient = minmod2;
  f.tracers = {concentration};
  
  /**
  Begin of added lines */
#if 1
  vof_advection({f},i);
  interfaces1 = interfaces, interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces1;
#endif
/**
End of added lines */
}
  
#if 0
event gfsview (i += 10) {
  foreach()
    f[] = clamp(f[],0,1);
  boundary ({f});
  static FILE * fpview = popen ("gfsview2D", "w");
  output_gfs (fpview);
}
#endif

/**
# Comments

We had to add the *added lines*. Without them, the program aborts immediatly
with this error message : "/basilisk/src/vof.h:112:error: Program received
signal SIGSEGV, Segmentation fault." Line 112 in [vof.h](/src/vof.h) is just:
"if (t.inverse)", where $t$ points sucessively to the the different members of
f.tracers.

I do not understand why we have to add these lines, since the only line in the
*vof()* event of [vof.h](/src/vof.h) is *vof_advection({f},i)*.
*/
