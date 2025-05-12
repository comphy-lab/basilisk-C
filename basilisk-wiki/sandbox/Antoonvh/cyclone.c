/**
# A Cyclone

A circular Gaussian vortex ($\omega = \omega_0 e^{-r^2 R^{-2}}$) exist
in a fluid placed in a rotating tank of depth $H$, kept inside the
basin walls by the acceleration of gravity $g$. In the co-moving
frame, the Coriolis parameter is $f$ and the fluid has viscosity
$\nu$.

On this page we consider Ekman pumping in a cyclonic vortex:

$$Ro = \frac{\omega_0}{f} = 20 $$
$$Fr = \frac{\omega_0 R}{\sqrt{gH}} = 0.1, $$
$$Re = \frac{\omega_0 R^2}{\nu} = 1000,$$
$$\Pi = \frac{R}{H} = 1.$$

![Evolution of the vorticity in the midplane](cyclone/omg.mp4)

The vortex decays relatively quick due to the secondary circulation.

~~~gnuplot x-z Slice: Upwelling in the vortex centre
set size ratio -1
set key outside
set xlabel 'x'
set ylabel 'z'
plot 'out' u 1:2:3:5 with vectors filled head lw 3 t 'Flow vectors'
~~~

~~~gnuplot Convergent flow in the bottom layer
 set ylabel 'y'
 plot "log" using 1:2:3:4 with vectors filled head lw 3 t 'Flow vectors'
~~~
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

double Ro = 20, Fr = 0.1, Re = 1000, Pi = 1;

double R = 1, omega_0 = 1.;
double H ,f , muz;
int main () {
  H = R/Pi;
  G = sq((omega_0*R/Fr))/H;
  f = omega_0/Ro;
  muz = omega_0*sq(R)/Re;

  for (vector u in ul) {
    u.n[left]   = -radiation (H);
    u.n[right]  =  radiation (H);
    u.n[bottom] = -radiation (H);
    u.n[top]    =  radiation (H);
  }
  L0 =  12.*R;
  X0 = Y0 = -L0/2.;
  nl = 10;
  N = 64;
  run();
}

#define RAD (sqrt(sq(x) + sq(y)))

event init (t = 0) {
  foreach() {
    vector u;
    scalar h;
    for (h, u, in hl, ul) {
      h[] = H/nl;
      u.x[] = -omega_0/2.*y*exp(-sq(RAD/R));
      u.y[] =  omega_0/2.*x*exp(-sq(RAD/R)); 
    }
  }
}
/**
The folowing physical processes are added to the solver

* Bottom drag
* Horizontal diffusion
* The Coriolis force

*/
event bottom_drag (i++) {
  vector u = ul[0];
  scalar h = hl[0];
  foreach() {
    foreach_dimension()
      u.x[] *= exp(-dt*muz/sq(h[]));
  }
}

#include "diffusion.h"
event horizontal_diff (i++) {
  const face vector nu[] = {muz, muz};
  vector u;
  scalar w; 
  for (u, w in ul, wl) {
    diffusion (u.x, dt, nu);
    diffusion (u.y, dt, nu);
    diffusion (w, dt, nu);
  }
}

event pressure (i++) {
  coord cor = {-f, f};
  foreach_face() {
    face vector av, uf;
    vector u;
    for (av, uf, u in al, ufl, ul) {
      foreach_dimension() {
	double ax = cor.x*(u.y[] + u.y[-1])/2;
	uf.x[] += ax*dt;
	av.x[] += ax;
      }
    }
  }
}


event outputs (t += 0.1) {
  scalar omg[], omgb[];
  double U = omega_0 * R;
  vorticity (ul[nl/2], omg);
  vorticity (ul[0], omgb);
  output_ppm (omg,  file = "omg.mp4", n = 300,
	      min = -1.5*omega_0, max = 1.5*omega_0);
  output_ppm (omgb, file = "omgb.mp4", n = 300,
	      min = -1.5*omega_0, max = 1.5*omega_0);
#ifdef NH
  output_ppm (wl[nl/2], file = "w.mp4", n = 300,
	      min = -0.05*U, max = 0.05*U);
#endif
  output_ppm (eta,  file = "eta.mp4", min = H*(1 - Fr), max = H*(1 + Fr));
}

event slice (t = 20) {
  double width = 3*R;
  vector u;
  scalar w;
  int np = 20;
  coord c[np];
  double xp = -width/2. * (1 - 1./(np));
  for (int j = 0; j < np; j++) {
    c[j] = (coord){xp, 0};
    xp += width/np;
  }
  double V[np*3];
  int lay = 0;
  for (u, w in ul, wl) {
    interpolate_array ({u, w}, c, np, V, true); 
    xp = -width/2. * (1 - 1./(np));
    for (int n = 0; n < np; n++) {
      printf ("%g %g %g %g %g\n", xp, (lay + 0.5)*H/nl, V[n*3],
	      V[n*3 + 1], V[n*3 + 2]);
      xp += width/np;
    }
    lay++;
  }

  foreach() {
    vector u = ul[0];
    if (RAD < width/2) 
      fprintf (stderr," %g %g %g %g\n", x, y, u.x[], u.y[]);
  }
}

