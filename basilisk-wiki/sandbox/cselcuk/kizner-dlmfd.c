/**
# Vortex Ejection from a mode 3 instability with dlmfd (distributed Lagrange multipliers/Fictitious domain method)

This page is the fictitious-domain version of [Antoon's simulation with embeded geometry](http://basilisk.fr/sandbox/Antoonvh/kizner2.c). It is copy/pasted here and adapted for the fictitious domain method.
According to [Kizner et
al. (2013)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/instabilities-of-the-flow-around-a-cylinder-and-emission-of-vortex-dipoles/2327C2CFA76059D27461A2BC6A09F146),
a flow around a no-slip cylinder with radius $R$ maybe unstable and eject three
dipolar vortex pairs. We study the flow using the fictitious-domain method and
the Navier-Stokes solver (without the double-projection scheme). Furthermore, we will `view` our results.
*/



/* Physical parameters */
# define rhoval (1.) // fluid density
# define radi (1.) // particle radius
# define Re 30000.
# define tval (1./Re)
# define MAXLEVEL 12
# define mydt 0.04
# define maxtime  100.

/* Criteria for adaptivity */
# define adaptive 1
# define UxAdaptCrit (4.E-3)
# define UyAdaptCrit (4.E-3)

# include "dlmfd.h"
# include "view.h"


/**
The maximum resolution is set to $\Delta_{min}=R/100$. This allows to run on the sandbox server.
*/

int main() {
  L0 = 40.;
  init_grid(64);
  /**
Rather than choosing stress-free outer-domain boundaries, periodic
boundaries are used. This markably increased the congergence
properties of the iterative Multigrid strategy applied to the
advection and viscous problems.
   */
  foreach_dimension()
    periodic(left);
    /*
      Since the perturbed initialized solution is slightly inconsistent, the
      Poisson solver is tuned to be extra robust for the first ten
      timesteps.
    */
  CFL = 0.6;
  DT = mydt;
  TOLERANCE = 1E-4;
  NITERMIN = 5;
  run();
}



/**
The cylinder is defined and the flow field is initialized c.f. Kizner
et al. (2013) with an $m=3$ perturbation.
 */
#define RAD (pow(sq(x-0.5*L0) + sq(y-0.5*L0), 0.5))
#define THETA(M) (M*asin((x-0.5*L0)/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M)))/(pow(1 + 0.5*sq(P), 0.5)))


event init(t = 0) {
   origin (0., 0.);
  
  double a1 = 1.5 , b1 = 2.25 ;
  double P = 0.005, m = 3;
  double gamma = (sq(a1) - 1.)/(sq(b1) - sq(a1));
  refine (RAD < b1  && level < (MAXLEVEL-1));
  refine (RAD < 1.05 && RAD > 0.95  && level < MAXLEVEL);
  foreach() {
    double r = RAD;
    double r1 = RADP(P,m)*r;
    double vr;
    if (r1 > 0.9 && r1 < a1)
      vr = r1 - 1./r1;
    else if (r1 >= a1 && r1 <= b1)
      vr = -gamma*r1 + ((1 + gamma)*sq(a1) - 1.)/r1;
    else // (0.9 > r || r > b)
      vr = 0;
    u.x[] = cm[]*0.5*vr*r*-(y-L0/2.)/(sq((x-L0/2.)) + sq((y-L0/2.)));
    u.y[] = cm[]*0.5*vr*r* (x-L0/2.)/(sq((x-L0/2.)) + sq((y-L0/2.)));
  }
  
  
  /* initial condition: particles/fictitious-domains positions */
  particle * p = particles;
 
  for (int k = 0; k < NPARTICLES; k++) {
    GeomParameter gp = {0};
    gp.center.x = L0/2; 
    gp.center.y = L0/2;
    gp.radius = radi;
    p[k].g = gp;

  }
}

event relax_a_little (i = 10) {
  NITERMIN = 1;
}

/**
## Ouput and Results

Movies are generated that display the vorticity dynamics and the
grid-cell structure. 

<video width="600" height="600" controls>
<source src="http://www.basilisk.fr/sandbox/cselcuk/kizner-dlmfd/kizner12.mp4" type="video/mp4">
</video> 
<video width="600" height="600" controls>
<source src="http://www.basilisk.fr/sandbox/cselcuk/kizner-dlmfd/kizner_cells12.mp4" type="video/mp4">
</video> 
 */
event movie (t += 0.4; t <= maxtime) {
  p.nodump = false;
  
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  dump("dump");
  view (fov = 5, quat = {0,0,0,1}, tx = -0.5*radi, ty = -0.5*radi, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  
  clear();
  squares ("omega", min = -0.75, max = 0.75,
	   map = cool_warm, linear = true);
  save ("kizner12.mp4");
  clear();
  cells();
  save("kizner_cells12.mp4");
  clear();
  view (fov = 2, quat = {0,0,0,1}, tx = -0.5*radi, ty = -0.5*radi, bg = {1,1,1}, width = 1080, height = 1080, samples = 1);
  cells();
  squares ("omega", min = -0.75, max = 0.75,
	   map = cool_warm, linear = true);
  save ("kizner_cells_zoom.mp4");
}

/**
Also there is this movie:

![A zoom](kizner-dlmfd/kizner_cells_zoom.mp4)(width="50%")
*/

event lastdump (t = end) {
  dump (file="dump");
  particle * p = particles;
  dump_particle(p);
  
  /* performances display/output */
  output_dlmfd_perf (dlmfd_globaltiming, i, p);
}


/**
~~~gnuplot Number of constrained cells, Lagrange multipliers, total cells.
set yr [0 : 25000]
set xlabel 'time iteration'
set ylabel 'number of Cells or Lagrange multipliers'
plot "dlmfd-cells" u 1:3 w l lw 2 title "Constrained cells","dlmfd-cells" u 1:2 w l lw 2 title "Lagrange multpliers", "dlmfd-cells" u 1:4 w l lw 2 title "Total cells"
~~~

~~~gnuplot Number of iterations for the fictitious domain solver.
reset
set xlabel 'time iteration'
set ylabel 'number of iterations'
plot "converge-uzawa" u 1:2 w l lw 2 title "dlmfd-solver"
~~~

~~~gnuplot error ||u_dlmfd - u_particle||_2 vs time iteration (here u_particle = 0 as the cylinder is immobile)
reset
set yr [0 : 0.000012]
set xlabel 'time iteration'
set ylabel 'dlmfd-solver error'
f(x) = 0.00001
plot "converge-uzawa" u 1:3 w l lw 2 title "||u_dlmfd - u_particle||_2",f(x) w l lw 2 title "convergence-criterion"
~~~

The performances of the dlmfd-solver is stored in the "dlmfd-perf" file. 

##Reference

Kizner, Z., Makarov, V., Kamp, L., & Van Heijst, G. (2013). *Instabilities of the flow around a cylinder and emission of vortex dipoles*. Journal of Fluid Mechanics, 730, 419-441. doi:10.1017/jfm.2013.345
*/
