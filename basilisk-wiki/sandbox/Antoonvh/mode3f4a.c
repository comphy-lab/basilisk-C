/**
# A mode 3 instability 

Unstable vortex following Flierl 

~~~gnuplot Inital profile
set xr [0:2]
set yr [-0.1: 1.1]
set size square
set ylabel 'v_{theta} [-]'
set xlabel 'r [-]'
set grid
plot 'prof' u 1:2 w l lw 2 t 'Piecewise',\
         '' u 1:3 w l lw 3 t 'Initialized'
~~~

![Vorticity and meterial line](mode3f4a/mov.mp4)

~~~gnuplot Stetching of material line
set xr [0:30]
set yr [0: 5]
set size square
set xlabel 'Time [-]'
set ylabel 'Stetch factor'
set grid
plot 'log' u 2:(1 +($9- 2.45*3.1415)/(2.45*3.1415)) w l lw 2 t 'Stretch'
~~~

The `TOLERANCE` controls the divergence, which can be diagnosed
without any discrete approximation.

~~~gnuplot Maximum divergence can differ slightly from the projection redidual after grid refinement 
set xr [0:30]
set yr [1e-8: 2e-4]
set logscale y
set size square
set xlabel 'Time [T]'
set ylabel 'Divergence [T^{-1}]'
set grid
set key right outside
plot 'log' u 2:10 w l lw 4 t 'Max Divergence', '' u 2:11 w l lw 2 t 'mgp2.resa', 1e-5 w l lw 2 t 'TOLERANCE'
~~~
 */
#define RKORDER (3)
#include "nsf4t.h"
#include "profile6.h"
vector uv[];
#define u uv
#define interpolate_linear interpolate_qubic_vertex
#include "tracer-particles.h"
#undef u
#include "view.h"
#include "scatter2.h"

scalar * tracers = NULL;
Particles parts;

double ue = 5e-5;
long unsigned int np = 1e5;
double P = 1e-5, m = 3; // Perturbation and mode
double b = 1.45; //Vortex parameter

int maxlevel = 8;

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 8.;
  origin (-pi*4./3., -2*exp(1) + 1.2); // Not centered
  N = 1 << maxlevel;
  //for (ue = 1e-4; ue >= 1e-6; ue /= pow(10, 1./3.))
    run();
}

#define RAD (sqrt((sq(x) + sq(y))))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M))))

event init (t = 0) {
  TOLERANCE = 1e-5;
  parts = new_tracer_particles (np);
  double Theta = 0;
  foreach_particle_in(parts) {
    p().x = (b + 1.)/2.*cos(Theta);
    p().y = (b + 1.)/2.*sin(Theta);
    p().z = 0.01;
    Theta += 2*pi/(double)(np + 1);
  }
  refine (sq(x) + sq(y) < sq(b*1.2) && level < maxlevel + 2 );
  printf ("# Initializing using %ld cells...\n", grid->tn);
  double betav = 1./(1 - sq(b)), alphav = -betav*sq(b);
  vector uc[], ucf[];
  foreach() {
    double rp = RAD*RADP(P,m), vr = 0;
    if (rp <= 1.)
      vr = rp;
    else if (rp > 1 && rp <= b) 
      vr = alphav/rp + betav*rp;
    uc.x[] = -y/rp*vr;
    uc.y[] =  x/rp*vr;
    ucf.x[] = ucf.y[] = 0;
  }
  boundary ((scalar*){uc, ucf});
  /**
     The Helmholtzfilter is applied to smoothen the piecewise flow
     profile.
  */
  const face vector alphaf[] = {-sq(0.4/(2.*pi)), -sq(0.4/(2.*pi))};
  foreach_dimension()
    poisson (ucf.x, uc.x, alphaf, unity);  
  foreach_face()
    u.x[] = ((-ucf.x[-2] + 7*(ucf.x[-1] + ucf.x[]) - ucf.x[1])/12.);
  scalar vtb[], vta[], * ab = {vtb, vta};
  foreach() {
    vtb[] = x/RAD*uc.y[] - y/RAD*uc.x[];
    vta[] = x/RAD*ucf.y[] - y/RAD*ucf.x[];
  }
  profile (ab, sqrt(sq(x) + sq(y)), "prof");
  project (u, p2);
}

event vertex_field (i++) {
  foreach_dimension() {
    uv.x.restriction = restriction_vert;
    uv.x.prolongation = refine_vert5;
  }
  foreach() {
    foreach_dimension()
      uv.x[] = FACE_TO_VERTEX_4(u.x);
  }
  boundary ((scalar*){uv});
}


event adapt (i++) {
  adapt_flow (ue, 99, 1);
}

event logger (i += 5) {
  scalar div[];
  foreach() {
    div[] = 0;
    foreach_dimension()
      div[] += (u.x[1] - u.x[0]);
    div[] = fabs(div[]/Delta);
  }
  stats ds = statsf (div);
  fprintf (stderr, "%d %g %d %d %d %d %ld %d %g %g %g %g\n", i, t, mgp.i, 
           mgp.nrelax, mgp2.i, mgp2.nrelax, grid->tn,
	   grid->maxdepth, plength (parts), ds.max, mgp2.resa, ue);
}

event mov (t += 0.5) {
  scalar omg[];
  vorticityf (u, omg);
  view (fov = 19, width = 900, height = 450);
  translate (x = -L0*.505) {
    cells();
  }
  translate (x = L0*.505) {
    squares ("omg", map = cool_warm);
    draw_curve (parts, lw = 10);
  }
  save ("mov.mp4");
  output_ppm (omg, min = 1.8, max = 2.2, file = "o.mp4");
}

event stop (t = 30);

/**
# Reference

G. R. Flierl, "On the instability of geostrophic vortices", *J. Fluid
Mech. 197, 349*
*/
