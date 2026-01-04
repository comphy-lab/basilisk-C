/**
# Soluble gas diffusing from a static bubble

This is the example discussed in section 3.3.1 of [Farsoiya et al.,
2021](#farsoiya2021).

The concentration at a point inside and outside the bubble is compared
with the analytical solution provided in [Farsoiya et al.,
2021](#farsoiya2021).

~~~pythonplot Concentration field inside the bubble
import numpy as np
import matplotlib.pyplot as plt
	
plt.figure()
t,Gb,Gl = np.loadtxt('static_bubble.sol',delimiter=' ',unpack=True);
plt.plot(t/40, Gb,'k',label='Analytical')

ts,cb9,cl9 = np.loadtxt('log',delimiter=' ',unpack=True)
plt.plot(ts/40,cb9,'k--',label=r'$d_0/\Delta x \approx 102$');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t\; \mathscr{D}_l/d_0^2$')
plt.ylabel(r'$c_b/c_{b0}$')
plt.tight_layout()

plt.savefig('p001cbt.svg')
~~~

~~~pythonplot Concentration field outside the bubble
plt.figure()

plt.plot(t/40, Gl,'k',label='Analytical')

plt.plot(ts/40,cl9,'k--',label=r'$d_0/\Delta x \approx 102$');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t\; \mathscr{D}_l/d_0^2$')
plt.ylabel(r'$c_l/c_{b0}$')
plt.tight_layout()

plt.savefig('p001clt.svg')
~~~

## References

~~~bib
@hal{farsoiya2021, hal-03227997}
~~~
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "henry.h"

scalar c[], * stracers = {c};
double bubble_radius = 1.;
double box_size = 10.;
double conc_liq1 = 0, conc_gas1 = 1.;

int MAXLEVEL = 9;

int main (int argc, char **argv)
{
  size (box_size);
	
  N = 1 << MAXLEVEL;

  rho1 = 1.;
  rho2 = 0.01;
  c.alpha = 0.001;
  TOLERANCE = 1e-4;
 
  c.D1 = 0.1;
  c.D2 = 1.;
  DT = 1.;
  run();
}

event init (t = 0)
{
  fraction (f, - (sq(bubble_radius) - sq(x - box_size*0.5) - sq(y)));
  foreach()
    c[] = conc_liq1*f[] + conc_gas1*(1. - f[]);
}
	
event extract (i++; t <= 48)
{		
  if (i == 0)
    fprintf (stderr, "# t ci co\n");
  fprintf (stderr, "%g %g %g\n",
	   t, interpolate(c,5.0,0.5), interpolate(c,5.0,1.5));
}
