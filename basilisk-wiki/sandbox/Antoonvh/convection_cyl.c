/**
![Some planets are smaller than others (e.g. King Kai's planet). Image via [the super official Dragon Ball website.](https://aminoapps.com/c/dragonball-super-7815629/page/item/king-kai-planet/BB1j_VRcmIJre75Wk5EpG5Lr84W64jnN3J)](http://pm1.narvii.com/6607/9e3791bf3cbb854033b00450e3e0c2089c15a48e_00.jpg)

# Atmospheric convection on a sphere

Rather than using an exotic mapping to allow a spherical metric, we
may *embed* a sphere in a 3D domain. Here we test the
concept by embedding a cylindrical earth into a quadtree. 

## The case

We take inspiration from Van Heerwaarden and Mellado (2016) and
consider a stably stratified atmosphere (BV-freq $N$) over a warm
earth surface with radius $R$. 

![The resulting movie](convection_cyl/movie.mp4)
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "profile6.h"
#include "view.h"

scalar s[];
scalar * tracers = {s};

double R = .5;
#define RADIUS (sqrt(sq(x) + sq(y) + sq(z)))
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
#if dimension == 3
u.r[embed] = dirichlet (0.);
#endif
s[embed]  = dirichlet (RADIUS + 1.);
s[left]   = dirichlet (RADIUS);
s[right]  = dirichlet (RADIUS);
s[bottom] = dirichlet (RADIUS);
s[top]    = dirichlet (RADIUS);

face vector muc[], av[];
double nu = 0.0001;

int main() {
  L0 = 4;
  X0 = Y0 = -L0/2.;
  N = 256;
  mu = muc;
  a = av;
  run(); 
}

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = RADIUS - R;
  fractions (phi, cs, fs);
  foreach()
    s[] = (RADIUS + 0.005*noise())*cs[];
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu;
  boundary ((scalar*){muc});
}

event tracer_diffusion(i++)
  diffusion (s, dt, mu);

/**
## Gravity

Gravity is oriented in the radial direction. */

event acceleration (i++) {
  foreach_face() {
    coord r = {x, y, z};
    av.x[] = (s[] + s[-1])*r.x/RADIUS;
  }
}

event movie (t += 0.1) {
  draw_vof ("cs", "fs", filled = -1, fc = {0.9, 0.9, 0.9});
  squares ("s", min = R, max = R + 1.5, linear = true);
  save ("movie.mp4");
}

event prof (t = {0, 10, 50, 100}) {
  char fname[99];
  sprintf(fname, "prof%g", t);
  profile ({s}, RADIUS, fname); 
}

/**
Here are some radial profiles of the $s$ field.

~~~gnuplot
set yr [0.5:2.2]
set xr [0.4:2.4]
set xlabel 'Buoyancy'
set ylabel 'Radius'
set key top left
set size ratio -1
plot 'prof0' u 2:1 w l lw 2 t 't = 0' ,\
     'prof10' u 2:1 w l lw 2 t 't = 10' ,\
     'prof50' u 2:1 w l lw 2 t 't = 50',\
     'prof100' u 2:1 w l lw 2 t 't = 100'
~~~
*/

event adapt (i++)
  adapt_wavelet ((scalar*){cs,s, u}, (double[]){0.01, 0.05, 0.02, 0.02}, 9, 5);
