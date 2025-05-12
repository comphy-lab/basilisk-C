/**
# Refine + axi.h + tracer.h collapse 

This is the axiadvection test peeled-off and with a refine. 
I wonder if this bug it is related to the previous one involving mask().
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tracer.h"

scalar f1[];
scalar * tracers = {f1}, * interfaces = NULL;

int main()
{
  X0 = -0.5;
  N = 64;
  TOLERANCE = 1e-12;
  f1.gradient = minmod2;
  run();
}

u.n[left] = dirichlet(1);
u.t[left] = dirichlet(0);
p[left]   = neumann(0);

u.n[top] = neumann(0);
p[top]   = dirichlet(0);
pf[top]  = dirichlet(0);

#define ellipse(xc, yc, a, b) (sq((x - xc)/(a)) + sq((y - yc)/(b)) - 1.)

event init (i = 0) {
  refine (sq(x) + sq ( y-0.3) - sq(0.12) < 0 && level < 7);
  foreach()
    u.x[] = 1.;
  fraction (f1, - ellipse (0, 0.3, 0.1, 0.1));
}

event output (i = 1)
  printf("Refine done\n");
