/**
# A Gaussian bump

As a hommage to the [tutorial](/Tutorial), we study the evolution of a
Gaussian bump in a heavy fluid.

![The case](diag_displace/f.mp4)
 */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "profile6.h"

int main() {
  const face vector grav[] = {0, -1};
  a = grav;
  rho1 = 2;
  rho2 = 1;
  L0 = 10;
  X0 = Y0 = -L0/2.;
  run(); 
}

event init (t = 0) {
  fraction(f, 2*exp(-sq(x)) - y);
}

event mov (t += 0.1) {
  output_ppm (f, file = "f.mp4" , n = 300);
}

event stop (t = 10);

/** 
## Diagnosis

We want to diagnose some field average over a translated
interface. Unfortunately, `vof.h` does not include an advection
function with the following syntax:

`vof_advection (scalar f, double dt, (const) face vector uf)`

Infact, we must apply some tricks to achieve something similar.
*/

void vof_advection2 (scalar f2, double dt2, face vector uf2, int i) {
  face vector uf_stored = uf;
  double dt_stored = dt;
  uf = uf2;
  dt = dt2;
  vof_advection ({f2}, i);
  dt = dt_stored;
  uf = uf_stored; 
}

/**
Now we setup the diagnosis itself. Say we wish to make a profile of
the $u_y^2$ field, along vertically displaced versions of the vof
interface.
 */

event diag (t += 1) {
  scalar uysq[], f2[];
  scalar_clone (f2, f);
  foreach() {
    f2[] = f[];
    uysq[] = sq(u.y[]);
  }
  boundary ({uysq, f2});
  
  char fname[99];
  sprintf (fname, "prof_above%d", (int)(t + 0.5));
  FILE * fp = fopen (fname, "w");
  
  face vector uuf[];  coord ufv = {0, 1};
  foreach_face()
    uuf.x[] = ufv.x;
  boundary ((scalar*){uuf});

  double dt2 = CFL*L0/(1 << grid->maxdepth);
  double yp = 0, displacement = 3;
  while (yp < displacement) {
    double val[1] = {0};
    output_ppm (f2, file = "f2.mp4", n = 300);
    interface_average ({uysq}, val, f2);
    if (pid() == 0)
      fprintf (fp,"%g %g\n", yp, val[0]);
    vof_advection2 (f2, dt2, uuf, i);
    yp += dt2*ufv.y;
  }
  fclose (fp);
}

/**

## Results

![Displaced interfaces](diag_displace/f2.mp4)

~~~gnuplot
set size square
set grid
set xlabel 'u_y^2 [m^2/s^2]
set ylabel 'y_d [m]'
plot 'prof_above1' u 2:1 w l t 't = 1',\
'prof_above6' u 2:1 w l t 't = 6',\
 'prof_above9' u 2:1 w l t 't = 9'  
~~~

*/
