/**
![Antique cars often use a so-called
 [thermosiphon](https://en.wikipedia.org/wiki/Thermosiphon) instead of
 a water pump. Hence the large radiators. Image of a Ford Model T via
 [Wikipedia](https://en.wikipedia.org/wiki/Antique_car)](https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/1916_Ford_Model_T_touring_car.JPG/250px-1916_Ford_Model_T_touring_car.JPG)

# Thermal circulation

We model a circulation in a torus with major radius $R_t$ and minor
radius $r_t$. The tube is diferentially heated so that a
buoyancy-driven flow starts.
*/
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#define BVIEW 1
#include "particles.h"

scalar b[], *tracers = {b};
face vector av[], muc[];

/**
   The torus is `#define`d using these variables and macros:
 */
double minor_rt = 1., major_rt = 3.;

#define R      (sqrt(sq(x) + sq(y)))
#define SINPSI (y/R)
#define COSPSI (x/R)
#define TUBE   (R <= 0.01 ? major_rt				\
		: (sqrt(sq(x - major_rt*COSPSI) +		\
			sq(y - major_rt*SINPSI) + sq(z))))
/**
The donut is a no-slip wall.
*/

u.n[embed] = dirichlet (0);
u.t[embed] = dirichlet (0);
#if (dimension == 3)
u.r[embed] = dirichlet (0);
#endif
b[embed]   = dirichlet (x);

int main() {
  L0 = 10;
  X0 = Y0 = Z0 = -L0/2.;
  N = 64;
  TOLERANCE = 1.e-4;
  mu = muc;
  a = av;
  run();
}
/**
   The geometry is initialized and tracer particles are seeded.
*/
event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex() 
    phi[] = minor_rt - TUBE;
  fractions (phi, cs, fs);
  
  n_part = 400;
  int j = 0;
  loc = malloc (n_part*sizeof(coord));
  do {
    double a = noise();
    double b = noise();
    if (sq(a) + sq(b) < sq(0.8*minor_rt)) {
      loc[j].x = a + major_rt;
      loc[j].y = 0;
      loc[j].z = b;
    } 
  } while (j++ < n_part);
}

event properties (i++) {   //diffusivity 
  foreach_face()
    muc.x[] = fm.x[]/200.;
  boundary ((scalar*){muc});
}

event acceleration (i++) { //Gravity
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event tracer_diffusion (i++) 
  diffusion (b, dt, mu);	 

/**
The maximum resolution corresponds to a $64^3$ cell equidistant grid.
 */

event adapt (i++) 
  adapt_wavelet ((scalar*) {cs, b, u},
		 (double[]){0.03, 0.4, 0.3, 0.3, 0.3}, 6);

event stop (t = 50.);

/**
## Results

We output a movie, displaying the buoyancy field, the velocity field,
particles and a poor reconstruction of a part of the embedded
boundary.

![Thermal convection in a tube](tube/mov.mp4)

The particles reveal an initial, short-lived recurring
(i.e. compensating) flow in the centre of the tube. This is when a
shallow ring of warm fluid rises arround the neutral fluid in the
center. Furthermore, an attractive secondary circulation is visible!

Naturally, we wish to study the evolution of the azimutal velocity
($u_{\theta}$), averaged over the minor-radial coordinate ($r$). The
results are plotted below.

~~~gnuplot Radial profiles of the azimuthal velocity.
set xr [0:1.1]
set yr [-0.1:1.5]
set xlabel 'minor radius'
set ylabel 'u_{theta}'
set size ratio 1
set grid
plot 'prof1' u 1:2 w l lw 2 t 't = 1',		\
     'prof3' u 1:2 w l lw 2 t 't = 3' ,		\
     'prof5' u 1:2 w l lw 2 t 't = 5' ,		\
     'prof10' u 1:2 w l lw 2 t 't = 10' ,	\
     'prof20' u 1:2 w l lw 2 t 't = 20' ,	\
     'prof35' u 1:2 w l lw 2 t 't = 35' ,	\
     'prof50' u 1:2 w l lw 2 t 't = 50'
~~~
*/


event movie (t += 0.1) {
#if (dimension == 3)
  scalar U[], t[], bd[];
  face vector tf[];
  double slice_a = -1;
  foreach_cell() {
    bd[] = b[];
    t[] = cs[];
    if (x > slice_a)
      t[] = 0;
    if (cm[] < 1e-6) {
      U[] = nodata;
      bd[] = nodata;
    }
    else {
      U[] = 0;
      foreach_dimension()
	U[] += sq(u.x[]);
      U[] = sqrt(U[]);
    }
  }
  foreach_face() {
    tf.x[] = fm.x[];
    if (x > slice_a)
      tf.x[] = 0;
  }
  boundary_flux ({tf});
  view (phi = 0.1, psi = 0., theta = 0.8);
  cells();
  squares ("bd");
  squares ("U", n = {1,0,0}, alpha = slice_a + 0.1, min = 0);
  draw_vof ("t", "tf");
  scatter (loc, s = 40);
  box();
#else
  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  cells();
  squares ("b");
#endif
  save ("mov.mp4");
}

#include "profile6.h"
event profs (t += 1.) {
  char fname[99];
  sprintf (fname, "prof%g", t);
  scalar Uth[];
  foreach() 
    Uth[] = (u.y[] * COSPSI - u.x[]*SINPSI)*cm[];
  boundary ({Uth});
  vertex scalar phi[];
  foreach_vertex()
    phi[] = TUBE;
  profiles ({Uth}, phi, fname, 1.1*minor_rt, 0, 0.25);
}
