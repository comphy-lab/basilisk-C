/**
# Heat transfer with an ellipse

For complex flow phenomena it is desirable to have a description with
a reduced complexity for the most prominent flow features. Here we
focus on the heat exchange between a ellipse and the surrounding flow.

The ellipse ($r_{min} = 1, r_{max} = 5r_{min}$ is in a flow with $Re =
\frac{Ur_{min}}{\nu} = 50$ and has a normalized scalar value at its
surface. We are interested in the effect of its orientation on the
flow and scalar surface flux.

## Results

The evolution of the "heat" field

![The two flows are rather different](ellipse/s.mp4)

Next, we consider the evolution of the scalar surface flux in the
latest quartile of each simulation.

~~~gnuplot The differ by 20%
set yr [0:1.3]
set xlabel 'tUR^{-1}_{min} [-]'
set ylabel 'flux [[s]R^3T^{-1}]'
set grid
set size square
plot 'diag1' u 1:(-$2) w l lw 3 t 'Exp. 1' ,	\
'diag2' u 1:(-$2) w l lw 3 t 'Exp. 2'
~~~
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar s[], * tracers = {s};

u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
s[embed]   = dirichlet (1.);

double r1 = 1, r2 = 5;
#define ELLIPSE (sq(x/r1) + sq(y/r2) - 1.)

int maxlevel = 10;
double ue = 0.05, be = 0.05;

face vector muc[];
double Re = 50;

bool turned = false;

int main() {
  periodic (left);
  L0 = 120;
  X0 = Y0 = -L0/2;
  mu = muc;
  run();
  r1 = r2, r2 = 1.;
  turned = true;
  run();
}

event properties (i++) {
  foreach_face() 
    muc.x[] = fs.x[]/Re;
  boundary ((scalar*){muc});
}

event init (t = 0) {
  refine (ELLIPSE <  2.5 && level  <  maxlevel - 1);
  refine (ELLIPSE > -0.5 && ELLIPSE <  0.5 &&
  	  level  <  maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = ELLIPSE;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach()
    u.x[] = cs[] > 0;
  boundary ({u.x});
}

event tracer_diffusion (i++) {
  diffusion (s, dt, muc);
}

event force (i++) {
  double FB = L0/5., tau = 1;
  foreach() {
    if (x < X0 + FB) {
      s[] -= s[]*dt/tau;
      u.y[] -= u.y[]*dt/tau;
      u.x[] -= (u.x[] - 1.)*dt/tau;
    }
  }
  boundary ({s, u});      
}

event adapt (i++) {
  adapt_wavelet ({cs, s, u},
		 (double[]){1e-5, be, ue,ue}, maxlevel, 5);
}

event mov (t += 1) {
  scalar m[];
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
  output_ppm (s, file = "s.mp4", n = 512, mask = m,
	      linear = true, max = .5, min = -.5, map = cool_warm,
	      box = {{X0 + 15, -15},{X0 + L0, 15}});
}

event diag_flux (t = 300; t += 1) {
  double flx = 0;
  foreach(reduction(+:flx)) {
    double val = 0, e = embed_flux (point, s, mu, &val);
    if (val)
      flx += (val - e*s[])*sq(Delta);
  }
  if (!turned) {
    static FILE * fp = fopen ("diag1", "w");
    fprintf (fp, "%g %g\n", t, flx);
  } else {
    static FILE * fp = fopen ("diag2", "w");
    fprintf (fp, "%g %g\n", t, flx);
  }
}

event stop (t = 400);
