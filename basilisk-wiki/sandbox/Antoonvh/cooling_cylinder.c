/**
# A cooling cylinder

A round (radius $R%) cylinder cools by long wave radiation ($Q =
-40Wm^{-2}$) which in turn cools the air via its duffusivity $\kappa =
1.5 \times 10^{-5}$, the resulting surface buoyancy flux induces a
flow of the fluid with viscousity $\nu = \kappa$. Furthermore, a
windspeed prevails with speed $U$.

Because $Q, \nu$ and $\kappa$ are fixed, the dimensional analysis will
not reduce the degrees of freedom, and we stay in the dimensinal
framework. We are curious about the cooling at the interface as a
function of $U$ and $R$.

![Temperature](cooling_cylinder/b.mp4)

~~~gnuplot Cooling
set xlabel 'Radius [m]'
set ylabel 'wind [m/s]'
set zlabel 'Cooling [K]'
set dgrid3d splines
set grid z
splot 'temp.RUb' with lines

~~~
 */
#define RAD (sqrt(sq(x) + sq(y)))
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "profile6.h"
scalar b[], * tracers = {b};

double R, U, B = -1.8e-3;
double  nuv = 1.5e-5, kappav = 1.5e-5;
b[embed] = neumann (B/kappav); // azimuthal Flux condition

int maxlevel = 8;
double be = 5e-3, ue = 1.5e-3;
face vector nu[], av[];
int row, row2;
int main() {
  periodic (left);
  mu = nu;
  a = av;
  DT = 0.01;
  for (U = 0; U <= 0.2; U += 0.04) { 
    for (R = 0.002; R <= 0.01; R += 0.004)  { //2 mm to 1 cm
      maxlevel = R < 0.005 ? 7 : 8;
      L0 = 20*R;
      X0 = Y0 = -L0/2;
      run();
    }
    row++;
  }
}

event properties (i++) {
  foreach_face()
    nu.x[] = fs.x[]*nuv;
  boundary((scalar*){nu});
}

event init (t = 0) {
  refine (RAD < 2*R && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = RAD - R;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach()
    u.x[] = U*(cs[] > 0);
  boundary ({u.x});
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[0,-1] + b[])/2.;
}
  
event heating_forcing (i++) {
  double tau = dt*20;
  foreach() {
    if (y < Y0 + 3*R || x < X0 + 3*R || x > X0 + L0 - 3*R) {
      b[] -= b[]*dt/tau;
      u.x[] -= (u.x[] - U)*dt/tau;
    }
  }
  boundary ({b, u.x});
}

event adapt (i++) {
  adapt_wavelet ({cs, b, u}, (double[]){1e-9, be, ue, ue, ue}, maxlevel);
}

event tracer_diffusion (i++)
  diffusion (b, dt, nu);

event mov (t += 0.05) {
  double bs;
  scalar bd[];
  foreach()
    bd[] = x < 0 ? b[] : nodata;
  interface_average ({b}, &bs, cs, fs);
  squares ("bd", min = bs - 0.01, max = -bs + 0.01, map = cool_warm );
  translate (z = 1e-5)
    draw_vof ("cs", "fs", filled = -1, fc = {0.3, 0.3, 0.3});
  translate (z = -1e-3)
    cells();
  char str[99];
  sprintf (str, "U = %g m/s, R = %g m ", U, R);
  draw_string (str, 15);
  save ("b.mp4");
}
 
int row2;
event stop (t = 5) {
  static FILE * fp = fopen ("temp_mat", "w");
  static FILE * fp2 = fopen ("temp.RUb", "w");
  double bs = 0;
  if (row2 != row) {
    fputc ('\n', fp);
    row2 = row;
  }
  interface_average({b}, &bs, cs, fs);
  fprintf (fp, "%g ", 27*bs);
  fprintf (fp2, "%g %g %g\n", R, U, 27*bs);
  fflush (fp);
  fflush (fp2);
}
