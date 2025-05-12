/**
# A Simple surface energy budget boundary condition

For an atmospheric medium with a viscosity $\nu$, the buoyancy budget
at the surface may read,

$$G + B = Q_n,$$

With $G$, the soil flux, expressed using a lumped-parameter (i.e. a
linear feedback) model:

$$G = \Lambda (b_{\mathrm{surf}} - b_0).$$

with $b_0$ a reference buoyancy value and $\Lambda$ the lumped
 soil conduction parameter. For the sensible buoyancy flux $B$, per
definition,

$$B = \kappa \frac{\partial b}{\partial \mathbf{n}},$$

with $\kappa$ the heat conductivity of the atmospheric
medium. Finally, for a simple diurnal cycle with period $T$, $Q_n$
is written as,

$$Q_n = \mathrm{max} \left[ B_0\ 
\mathrm{sin}\left( \frac{2\pi t}{T} \right),\ B_1 \right],$$

with $B_0$ and $B_1$ buoyancy flux scales. Giving a [Robin boundary
condtion](https://en.wikipedia.org/wiki/Robin_boundary_condition) for
$b$ at the bottom surface,

$$\Lambda b_{\mathrm{surf}} + \kappa \frac{\partial b}{\partial \mathbf{n}} = 
Q_n + \Lambda b_0.$$

Initially, the atmospheric fluid is at rest with a constant (vertical
$y$) stable stratification:

$$b(t=0) = b_0 + N^2y,$$

From the 7 system parameters (in order of appearance: $\{
\nu, \Lambda, \kappa, B_0, T, B_1, N\}$), we can construct 5 dimensionless
parameters:

$$\Pi_1 = \frac{B_0}{B_1} \approx -6,$$
$$\Pi_2 = TN \approx 2000 \rightarrow 500,$$
$$\Pi_4 = \frac{\sqrt{B_0T}}{\Lambda} \approx 5000 \rightarrow 2500,$$
$$Pr = \frac{\nu}{\kappa} \approx 1,$$
$$Re = \frac{B_0}{\nu} \left( \frac{2T}{\pi N^2} \right) ^{2/3} \approx 10^8 \rightarrow 5000.$$

See [Van Hooft et al (2019)](https://doi.org/10.1175/JAS-D-19-0023.1)
for details.

## Results

![The evolution of the buoyancy field and its horizontally-averaged
 profile](lumpedrob/output.mp4)

![The buoyancy field structure via $\mathrm{log}\left(\| \nabla
 b\| + 1\right)$](lumpedrob/db.mp4)
 */
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta)))	\
		      + ((neumann (0))*					\
			 ((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))

#define QN (max(B0*sin(2*pi*t/T), B1))

double T = 24*60, Lc = 1;
double Pi1 = -6, Pi2 = 500, Pi4 = 2500, Re = 5000;
double LAMBDA, B0, B1, sqN, NU;

int maxlevel = 8;

scalar b[], * tracers = {b};
b[bottom] = robin (LAMBDA, NU, QN);
b[top] = dirichlet (sqN*y);

FILE * gp;
face vector av[];

int main() {
  sqN    = sq(Pi2/T);
  B0 = sq(Lc)*sqN*pi/(2.*T);
  B1     = B0/Pi1;
  LAMBDA = sqrt(B0*T)/Pi4;
  NU = B0*pow(2*T/(pi*sqN), 2./3.)/Re;
  L0     = 3*Lc;
  X0     = -L0/2;
  periodic (left);
  u.t[bottom] = dirichlet (0.);
#if (dimension == 3)
  u.r[bottom] = dirichlet (0.);
  periodic (back);
#endif
  const face vector muc[] = {NU, NU, NU};
  mu = muc;
  a = av;
  run();
}

event init (t = 0) {
  DT = 1./(sqrt(sqN)*10.);
  foreach() {
    b[] = sqN*y;
    foreach_dimension()
      u.x[] = noise()/1000.;
  }
  boundary ({b, u});
  gp = popen ("gnuplot", "w");
  fprintf(gp,
	  "set term pngcairo size 500, 500\n"
	  "set xr [%g: %g]\n"
	  "set yr [0 : %g]\n"
	  "set grid\n"
	  "set size square"
	  "set key topleft"
	  "set xlabel 'b'\n"
	  "set ylabel 'y/L_c'\n",
	  B1/LAMBDA, sqN*Lc*2., Lc*2.);
}

event acceleration (i++) {
  coord grav = {0, 1, 0};
  boundary({b});
  foreach_face()
    av.x[] = grav.x*face_value(b, 0);
}

event tracer_diffusion(i++) 
  diffusion (b, dt, mu);

event adapt (i++) {
  double be = sqN*Lc/40., ue = 0.05;
  adapt_wavelet ({b, u}, (double[]){be, ue, ue, ue}, maxlevel, 4);
}

event damp (i++) {
  foreach() {
    if (y > 2.*L0/3.) {
      foreach_dimension()
	u.x[] *= exp (-dt*(y - 2*L0/3));
      b[] -= (b[] - sqN*y)*(1. - exp (-dt*(y - 2*L0/3)));
    }
  }
  boundary ({b, u});
}

#include "profile5c.h"
int frame = 0;
event movies (t += T/1000.) {
  scalar db[];
  foreach() {
    db[] = 0;
    foreach_dimension() 
      db[] += sq((b[1] - b[-1])/(2*Delta));
    db[] = db[] > 0 ? log(sqrt(db[]) + 1.) : 0;
  }
  boundary ({db});
  output_ppm (b, file = "b.mp4", linear = true,
	      n = 500, min = 0, max = sqN*Lc);
  output_ppm (db, file = "db.mp4", linear = true,
	      n = 500, min = 0, max = 3*log(sqN + 1));
  double by, yp = 1e-6;
  fprintf (gp,
	   "set output 'plot%d.png'\n"
	   "set title 'Time since sunrise:  %02d:%02d (HH:MM)'\n"	
	   "plot '-' w l lw 3 t 'buoyancy'\n",
	   frame, (int)(t/60.), (int)floor(fmod(t, 60.)));
  while (yp < L0) {
    by = 0;
    double dy = average_over_yp ({b}, &by, yp);
    fprintf (gp,"%g %g\n", by, yp);
    yp += dy;
  }
  fprintf (gp,"e\n");
  frame++;
}

event stop (t = 2.*T) {
  pclose (gp);
  system ("rm mov.mp4");
  system ("ffmpeg -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  sleep (10);
  system ("ffmpeg -y -i b.mp4 -i mov.mp4 -filter_complex hstack output.mp4");
  return 1;
}
