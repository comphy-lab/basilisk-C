/**
# Minimum working example with scalar attribute for a scalar

I add to the scalar fields a attribute which is a scalar field itself. Although
the attribute is not used, it is modified after its initialization. This minimum
working example show the different factors I found to have a influence on this
strange and unexpected (at least by me) behaviour.

I define some flags. If MAIN is activated, the initialization of the attribute
is done in the main() function, where we set usually f.sigma for example. with
DEFAULTS, it will be done in the defaults() event, and with MY_INIT, it will be 
done in the init() event.

With BEFORE, f[] (field to which we add a new attribute) is declared before 
#include "navier-stokes/centered.h", if not, it is done after.

With DO_BOUNDARY, we add a free flow boundary on top.

The results are at the end. */

#define MAIN 0
#define DEFAULTS 0
#define MY_INIT 1

#define BEFORE 0
#define DO_BOUNDARY 1

#include "fractions.h"

#if BEFORE
attribute {
  scalar associated_field;
}

scalar f[];
#endif

#include "navier-stokes/centered.h"

#if !BEFORE
attribute {
  scalar associated_field;
}

scalar f[];
#endif

#define L 1. // size of the box
#define LEVEL 3 // 7
#define T_END 1. //4
#define DELTA_T (T_END/100.)

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);
  
  FILE * fpc = fopen ("cells", "w");
  output_cells (fpc);
  fclose (fpc);
  
  /**
  If f.associated_field is set in main(), it is erased. */ 

#if MAIN
  scalar field = f.associated_field;
  foreach()
    field[] = 1.;
  boundary({field});
#endif
  FILE * fp = fopen ("0_field_main", "w");
  scalar field1 = f.associated_field;
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, field1[]);
  fclose (fp);

  run();
}

#if DO_BOUNDARY
p[top] = dirichlet(0.);
u.n[top] = neumann(0.);
#endif

event defaults (i = 0) {
#if DEFAULTS
  scalar field = f.associated_field;
  foreach()
    field[] = 1.;
  boundary({field});
#endif
  FILE * fp = fopen ("1_field_defaults", "w");
  scalar field1 = f.associated_field;
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, field1[]);
  fclose (fp);
}

event init (i = 0) {
  fraction (f, x - y);
  boundary({f});

#if MY_INIT
  scalar field = f.associated_field;
  foreach()
    field[] = 1.;
  boundary({field});
#endif
  FILE * fp = fopen ("2_field_init", "w");
  scalar field1 = f.associated_field;
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, field1[]);
  fclose (fp);
}


#if 1
event current_time (i++; i <= 10) {
  fprintf(stderr, "%d\t%g\n", i, t);
  fflush(stderr);
}
#endif

event projection (i = 0) {
  FILE * fp = fopen ("3_field_projection", "w");
  scalar field1 = f.associated_field;
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, field1[]);
  fclose (fp);
}

event adapt (i = 0) {
  FILE * fp = fopen ("4_field_adapt", "w");
  scalar field1 = f.associated_field;
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, field1[]);
  fclose (fp);
}

event this_is_the_end (t = end) {
  FILE * fp = fopen ("5_field_end", "w");
  scalar field1 = f.associated_field;
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, field1[]);
  fclose (fp);
}

/**
# Results

f.associated_field is initialized to 1. everywhere, and it should keep this value.
It is not the case.

For each situation (different initializations, with BEFORE or not, with free
flow condition or not), I look where the the field f.associated_field is modified.
The two tables below sum up the results, but I comment them after.

          main                          defaults        init
--------  ----------------------------  --------------  -----------
before    0 in defaults 0 and 1 after   random 0 and 1  ok
after     0 everywhere                  ok              ok

Table:  Without free flow boundary condition

          main                          defaults                          init
--------  ----------------------------  --------------------------------  -----------
before    0 in defaults 0 and 1 after   random 0 and 1                    ok
after     0 everywhere                  1 untill projection ~1e-24 after  1 untill projection ~1e-24 after

Table:  With free flow boundary condition

I really do not understand this bug, nevertheless:

* if f.associated_field is initialized in main(), it is set to 0. before
the defaults() event,
* if $f$ is declared before the inclusion of "centered.h", initializing
f.associated_field in the init() event seems ok, but we cannot do so if we
include "two-phase.h" (which has to be included after "centered.h" and which
declares $f$),
* I tried to see where in the projection() event f.field_attribute was modified
but to do so I needed to access to f.associated_field in projection(), which 
required to declare $f$ before "centered.h", and in this case, f.associated_field
is not modified anymore during the projection() event...

Here is an example of f.associated_field at the beginning of the adapt event,
i.e. just after the projection() event.

~~~gnuplot Associated field at the beginning of the adapt() event
set terminal @PNG enhanced size 640,640 font ",8"
set size ratio -1
unset key 
unset border
unset tics
plot 'cells' w l, '4_field_adapt' u 1:2:(sprintf("%.1e", $3)) with labels
~~~
*/
