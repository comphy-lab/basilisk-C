/**
# Oldroyd-B lid-driven cavity

The viscoelastic fluid is confined in a unit-square cavity. A
time-dependent tangential velocity is imposed on the top boundary.
$$
u_{top} = 8 \left[ 1 + \tanh \left(t - \frac{1}{2} \right) \right] x^2 (1-x)^2 
$$
This test case has been proposed by [Fattal \& Kupferman (2005)](#fattal2005).

We assume that the solvent and polymer viscosities are equal ($\beta =
0.5$) and the Weissenberg number $Wi$ is equal to 1. */

#define DT_MAX 5e-4
// Total viscoelastic viscosity 
#define MU0 1.
// Ratio of the solvent viscosity to the total viscosity
#define BETA 0.5
// Polymer viscosity
#define MUP ((1. - BETA)*MU0)
// Solvent viscosity 
#define MUS (BETA*MU0)
// Weissberger number
#define WI 1.
#define uwall(x,t) (8.*(1. + tanh(8.*(t - 0.5)))*sq(x)*sq(1. - x))

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "log-conform.h"

int main()
{
  L0 = 1.[0];
  DT = DT_MAX;
  N = 64;
  init_grid (N);
  lambda[] = WI;
  mup[] = MUP;
  mu[] = {MUS,MUS};
  stokes = true ;
  run();
}

u.t[top]    = dirichlet(uwall(x,t));
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

event init (i = 0) 
{
  scalar s = tau_p.x.x;
  s[left] = dirichlet (0.);
  s[right] = dirichlet (0.);
  s = tau_p.y.y;
  s[bottom] = dirichlet (0.);	
  foreach()
    u.x[] = 0.;
}

static double energy()
{
  double se = 0.;
  if (u.x.face)
    foreach(reduction(+:se))
      se += (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.*sq(Delta);
  else // centered
    foreach(reduction(+:se))
      se += (sq(u.x[]) + sq(u.y[]))/2.*sq(Delta);
  return se;
}

event kinetic_energy (i += 50)
{
  fprintf (stderr, "%g %g\n", t, energy());
}

event profile (t = 8.)
{
  FILE * fp = fopen("yprof", "w");
  for (double y = 0; y <= 1.; y += 0.01)
    fprintf (fp, "%g %g \n", y, interpolate (u.x, 0.5, y));
  fclose (fp);
  
  fp = fopen("xprof", "w");
  for (double x = 0; x <= 1; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0.75));
  fclose (fp);
}

/**
~~~gnuplot Time evolution of the total kinetic energy
set xlabel 't'
set ylabel 'kinetic energy'
plot 'lid-oldroydb.kinetic' w l t 'Fattal and Kupferman (2005)', \
     'log' w l t 'Basilisk'
~~~
~~~gnuplot Velocity profile u_x for x=0.5
set ylabel 'y'
set xlabel 'u_x'
plot 'lid-oldroydb.ux' w l t 'Fattal and Kupferman (2005)', \
     'yprof' u 2:1 t 'Basilisk'
~~~
~~~gnuplot Velocity profile u_y for y=0.75
set xlabel 'x'
set ylabel 'u_y'
plot 'lid-oldroydb.uy' w l t 'Fattal and Kupferman (2005)', 'xprof' t 'Basilisk'
~~~

## References

~~~bib
@article{fattal2005,
  title={Time-dependent simulation of viscoelastic flows at 
         high {W}eissenberg number using the log-conformation representation},
  author={Fattal, Raanan and Kupferman, Raz},
  journal={Journal of Non-Newtonian Fluid Mechanics},
  volume={126},
  number={1},
  pages={23--37},
  year={2005},
  publisher={Elsevier}
}
~~~

## See also

* [Lid-driven cavity for a Newtonian fluid](lid.c)
*/
