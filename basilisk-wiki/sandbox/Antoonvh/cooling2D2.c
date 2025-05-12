/**
# Thermal inertia

Outputs include:

![Visualization showing the evolution of 1) the temperature field (top
 left), 2) the radially averaged temperature profile (top right), 3)
 The vorticity field (bottom left) and 4) the grid cell distribution
 (bottom right)](cooling2D2/output2.mp4)

~~~gnuplot
set size square
set ylabel 'Temp [^oC]'
set xlabel 'time [s]'
set grid
plot 'out' w l lw 2 t 'Cylinder centre',\
'' u 1:3 w l lw 2 t 'Air'
~~~

## Setup

The Navier-Stokes equations are solved with a tracer field that is
subject to diffusion. For the flow field around the cylinder, we use
embedded boundaries. Meaning that the cylinder' surface cuts the
cells.
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "diffusion.h"
#include "tracer.h"
#include "view.h"
/**
The system is defined by a Gaussian-in-time inflow (i.e. the jet). 
 */
#define U_in (U_jet*exp(-sq((t - t_jet)/tau_jet)) + U_base)
#define T_in (T_jet*exp(-sq((t - t_jet)/tau_jet))) 
/**
A bunch of physical parameters are set.
 */
double T_jet = 5., U_jet = 1, t_jet = 20, tau_jet = 10;
double tend = 100, U_base = 0.1;
double R = 0.005;

double muv = 1.7e-5;
double lambda_a = 0.026, lambda_b = 0.6;
double rho_Cp_a = 1.*1006., rho_Cp_b = 1000.*4200.;
/**
Including numerical parameters for the maximum grid refinement level,
and the refinement criteria.
 */
int maxlevel = 7;
double Te = 0.05, ue = 0.05;
/**
Now we let the solver(s) know about our setup and parameters.  The
temperature field `T` is declared and we tell the solver we wish to
advect it with the flow. Further, the variable viscosity field (which
is defined on faces) is declared.
 */
scalar T[], * tracers = {T};
face vector muf[];

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

u.n[left] = dirichlet (U_in);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
T[left] = dirichlet (T_in);

p[right] = dirichlet (0.);
u.n[right] = neumann (0.);

int main() {
  L0 = 25*R;
  X0 = -L0/3.5;
  Y0 = -L0/2.;
  mu = muf;   
  N = 128;
  run();
}
/**
During initialization we overload the tracer attributes because it is
now also defined inside the embedded boundary.
 */
event init (t = 0) {
  for (scalar s in tracers) {
    s.refine = s.prolongation = refine_bilinear;
    s.coarsen = s.restriction = restriction_average;
  }
  // Embed a cylinder using a distance field
  refine (sqrt(sq(x) + sq(y)) < 2*R && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x) + sq(y));
  fractions (phi, cs, fs, R);
  //Initialize `u.x` and `T` 
  foreach() {
    u.x[] = U_in * (cs[] > 0);
    T[] = T_in;
  }
  boundary ({T, u.x});
} 
/**
The momentum does not diffuse into the solid.
 */
event properties (i++) {
  foreach_face() 
    muf.x[] = muv*fm.x[];
  boundary ((scalar*){muf});
}
/**
## Temperature diffusion

The temperature diffuses over the entire field. We set a field for
$\rho C_p$ and $\lambda$, which are weighted averages of the
properties of the corresponding material fractions.
 */
event tracer_diffusion(i++) {
  face vector kappaf[];
  scalar Cp[];
  foreach_face() 
    kappaf.x[] = fs.x[]*lambda_a + lambda_b*(1. - fs.x[]);
  foreach()
    Cp[] = cs[]*rho_Cp_a + rho_Cp_b*(1. - cs[]);
  boundary ({Cp, kappaf});
  diffusion (T, dt, kappaf, theta = Cp);
}

/**
We automatically focus the cells near the (evolving) boundary layer, and
wake.
 */

event adapt (i++) {
  ue = U_in/5.;
  adapt_wavelet ({cs, T, u}, (double[]){1e-3, Te, ue, ue, ue},
		 maxlevel, maxlevel - 3);
  unrefine (x > X0 + 4*L0/5); //Graceful outflow
}

event stop (t = tend);

/**
## Outputs and movie

Movies displaying the temperature field (`T.mp4`, using output_ppm),
the vorticity field and cells (`vor.mp4`, using bview), and
temperature profile (`mov.mp4`, via gnuplot and `ffmpeg`) are
generated. They are merged at the end of the simulation with `ffmpeg`.

But first, we generate a simple data file tracking the core
temperature.
 */
event track (t += 0.1) 
  printf ("%g %g %g\n", t, interpolate (T, 0, 0), T_in);


FILE * gp, * fpm, * fpv;

event init (t = 0) {
  gp = popen ("gnuplot", "w");
  fpm = popen("ppm2mp4 T.mp4", "w");
  fpv = popen("ppm2mp4 vor.mp4", "w");
  fprintf(gp,
	  "set term pngcairo size 400, 400\n"
	  "set xr [%g: %g]\n"
	  "set yr [%g : %g]\n"
	  "set grid\n"
	  "set size square\n"
	  "set xlabel 'Distance'\n"
	  "set key off\n"
	  "set ylabel 'T [K]'\n",
	  -R, 2.5*R,
	  -1. , T_jet + 1);
}

#undef EMBED //Do not use interpolate_embed
#include "profile6.h"
#define EMBED 1

int frame = 0;
event mov (t += 0.1) {
  output_ppm (T, fpm, n = 400, max = T_jet, min = 0.,
	      linear = true);

  view (width = 800, height = 400, fov = 19, tx = -1.5/7.);
  translate (x = -L0/2) {
    draw_vof ("cs", "fs", filled = -1, fc = {0.5, 0.5, 0.5});
    scalar omg[];
    vorticity (u, omg);
    squares ("omg", map = cool_warm);
  }
  translate (x = L0/2) {
    cells();
  }
  save (fp = fpv);
  
  fprintf (gp,
	   "set output 'plot%d.png'\n"
	   "set title 't = %d sec\n"
	   "plot '-' w l lw 3\n"
	   "%g %g\n", frame, (int)t,
	   -R, interpolate (T, 0, 0));
  double Tr, rp = -R*0.8;
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x) + sq(y)) - R; //Distance function
  scalar css[];
  face vector fss[];
  while (rp < 3*R) {
    fractions (phi, css, fss, rp);
    double dy = interface_average ({T}, &Tr, css, fss);
    fprintf (gp,"%g %g\n", rp, Tr);
    rp += dy;
  }
  fprintf (gp,"e\n");
  frame++;
}

event movie_maker (t = end) {
  pclose (gp);
  pclose (fpm);
  pclose (fpv);
  system ("rm mov.mp4");
  system ("ffmpeg -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  system ("ffmpeg -y -i T.mp4 -i mov.mp4 -filter_complex hstack output.mp4");
  system ("ffmpeg -y -i output.mp4 -i vor.mp4 -filter_complex vstack output2.mp4");
  return 1;
}
