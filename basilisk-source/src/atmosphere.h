#include "utils.h"

face vector u[], un[];
scalar h[], hn[], zb[];

// Default parameters
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;
// Viscosity
double NU = 0.;

trace
void advection_centered (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((f[] + f[-1,0])*u.x[] - 
	    (f[] + f[1,0])*u.x[1,0] +
	    (f[] + f[0,-1])*u.y[] - 
	    (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
}

trace
void advection_upwind (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] - 
	    (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
	    (u.y[] < 0. ? f[] : f[0,-1])*u.y[] - 
	    (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
}

trace
double timestep (void)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach(reduction(min:dtmax)) {
    Delta *= Delta;
    if (h[] > 0.) {
      double dt = Delta/(G*h[]);
      if (dt < dtmax) dtmax = dt;
    }
    foreach_dimension()
      if (u.x[] != 0.) {
	double dt = Delta/sq(u.x[]);
	if (dt < dtmax) dtmax = dt;
      }
  }
  return sqrt (dtmax)*CFL;
}

trace
void momentum (vector u, scalar h, vector du)
{
  scalar ke[];
  vertex scalar psi[];
  scalar dux[], dvy[];
  vector d;
  d.x = dux; d.y = dvy;

  foreach() {
#if 1
    ke[] = (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.;
#else
    double uc = u.x[]*u.x[] + u.x[1,0]*u.x[1,0];
    double vc = u.y[]*u.y[] + u.y[0,1]*u.y[0,1];
    ke[] = (uc + vc)/4.;
#endif
    foreach_dimension()
      d.x[] = (u.x[1,0] - u.x[])/Delta;
  }
  foreach_vertex()
    psi[] = (u.y[] - u.y[-1,0] + u.x[0,-1] - u.x[])/Delta;  

  coord f = {1.,-1.};
  foreach_face()
    du.x[] = 
      - (G*(h[] + zb[]) + ke[] - G*(h[-1,0] + zb[-1,0]) - ke[-1,0])/Delta
      + f.x*(((psi[] + psi[0,1])/2. + F0)*
	     (u.y[] + u.y[0,1] + u.y[-1,0] + u.y[-1,1])/4.)
      + NU*(u.x[0,1] + u.x[0,-1] - 2.*u.x[])/sq(Delta)
      + NU*(d.x[] - d.x[-1,0])/Delta;
}

trace
void advance (double t, scalar * f, scalar * df)
{
  vector u = {f[0], f[1]}, du = {df[0], df[1]};
  scalar h = f[2], dh = df[2];

  advection_centered (h, u, dh);
  momentum (u, h, du);
}

void update (double t, scalar * f)
{
}

event defaults (i = 0)
{
  foreach()
    h[] = 1.;
}

event init (i = 0)
{
}

void run (void)
{
  init_grid (N);

  timer start = timer_start();
  iter = 0, t = 0;
  while (events (true)) {
    double dt = dtnext (timestep ());
#if 1
    advection_centered (h, u, hn);
    foreach()
      h[] += hn[]*dt;
    momentum (u, h, un);
    foreach_face()
      u.x[] += un.x[]*dt;
#else /* unstable! */
    scalar f[3] = { u, v, h };
    scalar df[2][3] = {{ un,  vn,  hn },
		       { un1, vn1, hn1 }};
    runge_kutta (2, t, dt, 3, f, df, advance, update);
#endif
    iter = inext, t = tnext;
  }
  timer_print (start, iter, 0);

  free_grid ();
}
