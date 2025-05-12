/**
# Diffusion of a Gaussian pulse in a uniform domain

We consider the diffusion of a Gaussian pulse:

![The embedded region is marked with a black line. Notice that noting
 special seems to happen there](embed_gaussian/mov.mp4)

~~~gnuplot Convergence is OK
set logscale xy
set grid
set xr [16:512]
set yr [1e-5:0.1]
set xlabel 'N'
set ylabel 'Error'
plot 'out' u 1:2 t 'L_1', '' u 1:3 t 'Max', 5*x**(-2) t 'Second order'
~~~
 */
#include "flux_embed.h"
#include "utils.h"
#include "view.h"
double total_flux_embed (Point point, coord p, coord n) {
  return 0.;
}

double t0 = 0.25;
#define sqR (sq(x - 0.1234) + sq(y + 0.1357))
#define Gaussian (exp(-sqR/(4*(t + t0)))/(4*pi*(t+t0)))

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2;
  
  for (N = 32; N <= 256; N *= 2) {
    DT = 0.00025;
    run();
  }
}

event init (t = 0) {
  foreach() 
    s1[] = s2[] = Gaussian;
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x - 1) + sq(y - .5)) - 1;
  fractions (phi, cs, fs);
  output_ppm (s1, file = "s1.png", n = 300);
  output_ppm (cs, file = "cs.png", n = 300);
}

event dt_setter (i++) {
  dt = dtnext(DT);
}

event movie (t += 0.1) {
  scalar s[];
  foreach() {
    s[] = cs[]*s1[] + (1 - cs[])*s2[];
  }
  squares ("s", min = -0.1, max = 0.1);
  draw_vof ("cs", "fs");
  save ("mov.mp4");
    
}
     

event stop (t = 1) {
  double err = 0, emax = -1;
  foreach (reduction (+:err) reduction (max:emax)){
    double el = nodata;
    if (cs[] == 1) 
      el = fabs(s1[] - Gaussian);
    else if (cs[] == 0)
       el = fabs(s2[] - Gaussian);
    if (el != nodata) {
      if (el > emax)
	emax = el;
      err += sq(Delta)*el;
    }
  }
  printf ("%d %g %g\n", N, err, emax);
}  
