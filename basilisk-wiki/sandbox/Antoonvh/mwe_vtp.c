/**
# A minimal example for a crashing call to the `reposition()` function.

Thanks to Omri Argov for [discovering this
flaw](https://groups.google.com/forum/#!topic/basilisk-fr/ANaemjZ-BCw).

We take some setup inspiration from this [test](reversedp.c).

Plot the last known locations:

~~~gnuplot This is the intended result
set size ratio -1
set key right outside
plot 'out' w l t 'Cells', 'last_known' , #\
   #   'last_facets' w l lw 3 
~~~

The crash can be reproduced by reverting this patch:

[http://basilisk.fr/sandbox/Antoonvh/vof-tracer-particles.h?changes=20200421105835]()

You will find that it corresponds to particles moving
out the domain. 

## Setup

Omri correctly expected that the issue arises when particles move 
out the domain. We set that up, and output some data for visualization.
*/
#include "vof-tracer-particles.h" //We need the reposition function
#include "vof.h"                  //Which expects a vof field f.

Particles Pin;
double R = 0.2, yo = 0.5, xo = 0.5;
#define circle(x,y) (sq(R) - (sq(x - xo) + sq(y - yo)))

face vector uf[];

scalar f[], * interfaces = {f};

int main() {
  N = 32;
  run();
}

event init (t = 0) {
  fraction (f, circle (x, y));
  Pin = new_vof_tracer_particles (10,1);
  struct Init_P I = {10, xo, yo, R};
  place_in_circle (Pin, I);
  foreach_face(x)
    uf.x[] = 1;
  output_cells();
}

event set_dtmax (i++)
  dt = dtnext (0.005);

event trace_pos (i++) {
  FILE * fp = fopen ("last_known", "w"); 
  foreach_particle()
    fprintf (fp, "%g %g\n", x, y);
  fclose (fp);
  
  FILE * fpf = fopen ("last_fac", "w");
  output_facets (f, fp = fpf);
  fclose (fpf);
  }

event stop (t = 0.75);

