/**
# "Slave" solver for the coupling example

This should be read in combination with [master.c](). Aside from the
`slave.h` include, this file is almost identical to [karman.c](). */

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"

#include "slave.h"

double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

double D = 0.5, U0 = 1.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{
  solid (cs, fs, sqrt(sq(x) + sq(y)) - D/2.);
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}

event logfile (i++)
  fprintf (stderr, "slave %d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);

/**
We start the slave movie at $t = 15$ so that the master and slave
movies are synchronized (see how the `slave_start` variable is used in
[master.c]()). */

event movies (t = 15; t += 0.05)
{
  scalar omega[], m[];
  foreach() {
    omega[] = (u.y[1] - u.y[-1] + u.x[0,-1] - u.x[0,1])/(2.*Delta);
    m[] = cs[] - 0.5;
  }
  output_ppm (omega, file = "vort-slave.mp4",
	      min = -5, max = 5, linear = true, mask = m, n = 256);
}

int main()
{
  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = 128;
  mu = muv;
  run();
}
