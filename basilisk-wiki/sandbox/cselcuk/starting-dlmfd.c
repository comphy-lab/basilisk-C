/**
# Starting flow around a cylinder with DLMFD

This is a canonical case of complex boundary layer separation,
inspired by the experiments of [Bouard & Coutanceau,
1980](#bouard1980). Notable early numerical simulations include the
results of [Koumoutsakos and Leonard, 1995](#koumoutsakos1995),
hereafter K & L, which will be used in the comparisons below.

We will solve the Navier--Stokes equations and we will be using the [fictititous-domain](DLMFD_reverse_Uzawah.h) method.

Note that we could also use the "double projection" method that is apparently necessary for the
cut-cell method but we dont. */

#include "DLMFD_reverse_Uzawa.h"
#include "navier-stokes/perfs.h"
#include "view.h"

/**
High-resolution is needed to resolve the boundary layers
properly. [Mohaghegh et al., 2017](#mohaghegh2017) propose to use a
maximum resolution of order $D/10/\sqrt{Re}$, with $Re$ the Reynolds
number and $D$ the cylinder diameter. The cylinder diameter will be
set to unity, and the domain size to 18, so that the corresponding
levels of refinement are approximately 12 and 16 for $Re=1000$ and
$Re=9500$, respectively. */

int maxlevel = 12;  // 15/16 for Re = 9500, 12 for Re = 1000
double Re = 1000;   // or 9500
double cmax = 3e-3; // 1e-3 for Re = 9500, 3e-3 for Re = 1000
double deltau;
scalar un[];


/**
We set a unit velocity inflow on the left and an outflow on the
right. */

u.n[left] = dirichlet(1);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

/**
Command line arguments can be used to change the default
parameters. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    Re = atof (argv[2]);
  if (argc > 3)
    cmax = atof (argv[3]);
  
  /**
  The domain is $18\times 18$ and we model the full cylinder. */
  
  size (18);

  /**
  We tune the Poisson solver. */
   
  TOLERANCE = 1e-4;
  NITERMIN = 2;
  
  run();
}


/**
## Initial conditions
*/

event init (t = 0)
{

  /* Set dynamic viscosity */
  const face vector muc[] = {1/Re, 1/Re};
  mu = muc;

  /* set density of the flow */ 
  const scalar rhoc[] = 1.;
  rho = rhoc;

  /* initial condition: particles position */
  particle * pp = particles;

  GeomParameter gp = {0};
    
  pp[0].iscube = 0;
  pp[0].iswall = 0;
  /* write the position with respect to 0 */
  gp.center.x = L0/2;
  gp.center.y = L0/2;
  gp.radius   = 0.5;
      
  pp[0].g = gp;
   
  /* particle id */
  pp[0].pnum = 0;

  /**
  We can restart from a dump file. */
  if (!restore ("restart")) {

    init_file_pointers(pdata, fdata, 0); // 0 - initial simulation / 1 - restarting simulation
 
    /**
       Otherwise, we first create a mesh initially refined only around
       the cylinder. */
    
    refine (sqrt (sq(x-gp.center.x) + sq(y-gp.center.y)) < 1.2*gp.radius && sqrt (sq(x-gp.center.x) + sq(y-gp.center.y)) > 0.8*gp.radius && level < maxlevel);

  }
  else { // restart
  
    init_file_pointers(pdata, fdata, 1);
 
  }

  foreach()
    un[] = u.x[];
}



/**
## Positions of the separation points

We would like to track the positions with time of the separation points
on the surface of the cylinder, as done by [K & L,
1995](#koumoutsakos1995).

We first define a function which computes (and interpolates) the vorticity at the surface
of the cylinder and returns an array of $(\theta,\omega)$ pairs, with
$\theta$ the angular coordinate and $\omega$ the corresponding value
of vorticity. */

typedef struct {
  double theta, omega;
} ThetaOmega;

int compar_theta (const void * a, const void * b)
{
  const ThetaOmega * p1 = a, * p2 = b;
  return p1->theta > p2->theta ? 1 : -1;
}

void theta_omega_dlmfd (const particle * pp, ThetaOmega ** to)
{
  scalar vortz[];
  vorticity (u, vortz);
   
  Array * a = array_new();
  coord center = (pp->g).center;
  for (int ii = 0; ii < pp->s.m; ii++) {

    ThetaOmega t;
    t.omega = interpolate (vortz, pp->s.x[ii], pp->s.y[ii]);
    t.theta = atan2 (pp->s.y[ii] - center.y, pp->s.x[ii] - center.x);
    array_append (a, &t, sizeof (ThetaOmega));
  }
  
  qsort (a->p, a->len/sizeof(ThetaOmega), sizeof(ThetaOmega), compar_theta);
  ThetaOmega t = {nodata, nodata};
  array_append (a, &t, sizeof (ThetaOmega));
 
  ThetaOmega * p = a->p;
  free (a);
  *to =  p;
}

/**
The zeros of the function approximated by the $(\theta,\omega)$ array
are then recorded, together with the corresponding time, in the file
pointed to by *fp*.  */

void omega_zero_dlmfd (FILE * fp, const particle * pp)
{
   /* fixme: this function will not work with MPI because of the
      interpolations which are thread-dependent */
  ThetaOmega * a = NULL;
  theta_omega_dlmfd (pp, &a);

  for (ThetaOmega * o = a; (o + 1)->theta != nodata; o++) {
    ThetaOmega * o1 = o + 1;
    if (o1->omega*o->omega < 0.)
      fprintf (fp, "%g %g\n", t,
	       o->theta + o->omega*(o1->theta - o->theta)/
	       (o->omega - o1->omega));
  }
  free (a);
  fflush (fp);
}


/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.x, un);
  fprintf (stderr, "log output %d %g %d %d %g %g %g %ld\n", i, t, mgp.i, mgu.i, mgp.resa, mgu.resa, deltau, grid->tn);

}

/**
## Images and animations

We display the vorticity field and the corresponding adaptive mesh. */

void display_omega (int width, int height)
{
 view (fov = 2.16967, quat = {-0.00418115,0.0542593,-0.000125816,0.998518}, tx = -0.53112, ty = -0.499909, bg = {1,1,1}, width = width, height = height, samples = 1);

  squares ("omega", min = -12, max = 12, linear = 1, map = cool_warm);
}

/**
Still images are created at times matching those in various papers. */

event snapshots (t = 1.; t <= 3.; t += 1.)
{
  scalar omega[];
  vorticity (u, omega);

  if (t == 1. || t == 3.) {
    // Fig. 4.
    view (fov = 2.16967, quat = {-0.00418115,0.0542593,-0.000125816,0.998518}, tx = -0.53112, ty = -0.499909, bg = {1,1,1}, width = 966, height = 775, samples = 1);
    
    double max = Re == 9500 ? 40 : 12;
    squares ("omega", min = -max, max = +max, linear = 1, map = cool_warm);
    char name[80];
    sprintf (name, "zoom-%g.png", t);
    save (name);

    squares ("level");
    sprintf (name, "cells-%g.png", t);
    save (name);
  }
  
  if (t == 3.) {
    // Fig. 3.
    display_omega (640, 480);
    save ("omega-3.png");
  }
  
  p.nodump = false;
  char name [80];
  sprintf (name, "dump-%g", t);
  dump (name);
}

/**
An animation is created. The iteration interval is adjusted depending
on the maximum spatial resolution. */

event movie (i += 10*(1 << (maxlevel - 12)))
{
  scalar omega[];
  vorticity (u, omega);
  display_omega (1280, 960);
  save ("omega.mp4");
}


event movie1 (i += 10*(1 << (maxlevel - 12)))
{
  scalar omega[];
  vorticity (u, omega);
  display_omega (1280, 960);
  cells();
  save ("movie1.mp4");
}

/**
## Surface vorticity profiles

We are also interested in the details on the surface of the cylinder,
in particular surface vorticity. */

void cpout_dlmfd (FILE * fp, particle * pp)
{
  scalar omega[];
  vorticity(u, omega);
  
  for (int ii = 0; ii < pp->s.m; ii++) {
    GeomParameter gc = pp->g;
    coord c = gc.center;
    double w = atan2 (pp->s.y[ii] - c.y, pp->s.x[ii] - c.x);
    if (w < 0) w += 2*pi; 
    fprintf (fp, "%g %g %g %g \n",
	     pp->s.x[ii], // 1
	     pp->s.y[ii], // 2
	     w, // 3
	     interpolate (omega, pp->s.x[ii], pp->s.y[ii])
	     );
  }
}





/**
## Adaptive mesh refinement

The mesh is adapted according to the embedded boundary and velocity
field. */

event adapt (i++) {

  particle * pp = particles;
  double rhoval = 1.;
  
  DLMFD_subproblem (pp, i, rhoval);

  /* Save forces acting on particles before the adapting the mesh */
  sumLambda (pp, fdata, t, dt, flagfield, DLM_lambda, index_lambda, rhoval);

  static FILE * fp = fopen ("omega", "w");
  omega_zero_dlmfd (fp, pp);

  if (t == 0.5|| t == 0.9|| t ==1.5 || t == 2.5) {
    printf("inside cp block\n");
    char name[80];
    sprintf (name, "cp-%g", t);
    FILE * fp2 = fopen (name, "w");
    cpout_dlmfd (fp2, pp);
    fclose (fp2);
  }
  
  /* Free particle structures (we dont need them anymore) */
  free_particles (pp, NPARTICLES);

  /* Save particles trajectories */
  particle_data(pp, t, i, pdata);

  /* char name[80] = "test-init-starting"; */
  /* scalar * list = {p, index_lambda.x, flagfield}; */
  /* vector * vlist = {u}; */
  /* save_data(list, vlist, i, t, name); */
  /* return 1; */
  
  adapt_wavelet ({flagfield,u}, (double[]){1e-4,cmax,cmax}, maxlevel, 5);
}

event surface_profiles (t = {0.5,0.9,1.5,2.5})
{
  printf("forcing output data at t= {0.5,0.9,1.5,2.5}\n")
}


/**
## Results

### Re = 1000

![Animation of the vorticity field and adaptive mesh, Re = 1000.](starting-dlmfd/omega.mp4)(width="640" height="480")

![Animation of the adaptive mesh, Re = 1000.](starting-dlmfd/movie1.mp4)(width="640" height="480")

The final state at $tU/D = 3$ can be compared with figure 3 (top row) of
[Mohahegh et al. 2017](#mohahegh2017).

~~~gnuplot Time history of drag coefficient. Re = 1000. See Fig. 1a of [Mohahegh et al. 2017](#mohahegh2017).
set xlabel 'tU/D'
set ylabel 'C_D'
set grid
set pointsize 0.5
plot [][0:2]\
     "`echo $BASILISK`/test/starting/level-13/log" u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk-cut-cells (13 levels)',\
     "`echo $BASILISK`/test/starting/fig1a.SIM" u ($1/2.):2 w l lw 2 t 'SIM (Mohahegh et al 2017)', \
     "`echo $BASILISK`/test/starting/fig1a.KL" u ($1/2.):2 pt 7 t 'K and L. 1995', \
     "`echo $BASILISK`/test/starting/fig19.f" u ($1/2.):2 pt 9 t 'friction, K and L. 1995', \
     "`echo $BASILISK`/test/starting/fig19.p" u ($1/2.):2 pt 11 t 'pressure, K and L. 1995', \
     'sum_lambda-0' every 10 u 1:($2*2) w l lw 2 t 'Basilisk-dlmfd (12 levels)'
~~~

Note that the points of Figure 4 of [K. & L. 1995](#koumoutsakos1995)
do not seem to match the data in Fig. 5a and 5b of the same paper,
which explains the disagreement in the figure below. This agreement
should be much better as can be seen on the more detailed surface
vorticity figures below.

~~~gnuplot Location of the points of zero surface vorticity. Re = 1000. See Fig. 4 of [K. & L. 1995](#koumoutsakos1995).
reset
set xlabel 'tU/D'
set ylabel 'Angle/pi'
set grid
set key bottom right
set pointsize 0.5
plot "`echo $BASILISK`/test/starting/level-13/omega" u 1:($2/pi) pt 5 t 'Basilisk-cut-cell (13 levels)', \
     "`echo $BASILISK`/test/starting/fig4.1000.KL" u ($1/2.):2 pt 7 ps 0.9 t 'K and L. 1995', \
     'omega' every 10 u 1:($2/pi) pt 5 t 'Basilisk-dlmfd (12 levels)'
~~~

Note that the vorticity field is not constrained within the ficititous
domain problem and an additional interpolation is needed to get the
vorticity values at the surface of the cylinder on the location of the
Lagrange multipliers. This results in a loss of accuracy in the post-process.

~~~gnuplot Surface vorticity at $tU/D=0.5$. Re = 1000. See Fig. 5a of [Mohahegh et al. 2017](#mohahegh2017).
reset
set xlabel 'theta/pi'
set ylabel 'omega_sD/U'
set grid

plot [0:2] '< sort -k3,4 "$BASILISK/test/starting/level-13/cp-0.5" |awk -f "$BASILISK/test/starting/surface.awk"' w l lw 2 t 'Basilisk-cut-cells (13 levels)',\
     "`echo $BASILISK`/test/starting/fig5a.SIM" every 20 pt 7 t 'SIM (Mohahegh et al 2017)',\
     "`echo $BASILISK`/test/starting/fig5a.KL" w l lw 2 t 'K and L 1995', \
     'cp-0.5' every 10 u ($3/pi):($4*2) w l lw 2 \
     t 'Basilisk-dlmfd (12 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=1.5$. Re = 1000. See Fig. 5b of [Mohahegh et al. 2017](#mohahegh2017).

plot [0:2] '< sort -k3,4 "$BASILISK/test/starting/level-13/cp-1.5" |awk -f "$BASILISK/test/starting/surface.awk"' w l lw 2 t 'Basilisk-cut-cells (13 levels)',\
     "`echo $BASILISK`/test/starting/fig5b.SIM" every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     "`echo $BASILISK`/test/starting/fig5b.KL" w l lw 2 t 'K and L 1995', \
     'cp-1.5'  every 10 u ($3/pi):($4*2) w l lw 2 \
     t 'Basilisk-dlmfd (12 levels)'
~~~

### Re = 9500

![Animation of the vorticity field and adaptive mesh, Re = 9500, 15 levels.](starting-dlmfd/level-15/omega.mp4)(width="640" height="480")

![Animation of the adaptive mesh, Re = 9500.](starting-dlmfd/level-15/movie1.mp4)(width="640" height="480")
The final state at $tU/D = 3$ can be compared with figure 3 (bottom row) of
[Mohahegh et al. 2017](#mohahegh2017).

~~~gnuplot Time history of drag coefficient. Re = 9500. See Fig. 1b of [Mohahegh et al. 2017](#mohahegh2017).
set xlabel 'tU/D'
set ylabel 'C_D'
set grid
plot [][0:2.5] \
     "`echo $BASILISK`/test/starting/fig1b.SIM" u ($1/2.):2 w l lw 2 t 'SIM (Mohahegh et al 2017)', \
     "`echo $BASILISK`/test/starting/fig1b.KL" u ($1/2.):2 pt 7 t 'K and L. 1995', \
     "`echo $BASILISK`/test/starting/level-15/log" u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk (15 levels)', \
     "`echo $BASILISK`/test/starting/level-16/log" u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk (16 levels)',\
     'level-15/sum_lambda-0' every 10 u 1:($2*2) w l lw 2 t 'Basilisk-dlmfd (15 levels)'
~~~

~~~gnuplot Location of the points of zero surface vorticity. Re = 9500. See also Fig. 4 of [K. & L. 1995](#koumoutsakos1995).
reset
set term pngcairo enhanced font ",10"
set output 'loc9500.png'
set xlabel 'tU/D'
set ylabel 'Angle/pi'
set grid
set key bottom right
plot "`echo $BASILISK`/test/starting/level-15/omega" u 1:($2/pi) pt 7 ps 0.25 t 'Basilisk (15 levels)', \
     "`echo $BASILISK`/test/starting/level-16/omega" u 1:($2/pi) pt 5 ps 0.25 t 'Basilisk (16 levels)',\
     'level-15/omega' every 10 u 1:($2/pi) pt 5 t 'Basilisk-dlmfd (15 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=0.9$. Re = 9500. See Fig. 6a of [Mohahegh et al. 2017](#mohahegh2017).
reset
set term @SVG
set xlabel 'theta/pi'
set ylabel 'omega_sD/U'
set grid

plot [0:2] "`echo $BASILISK`/test/starting/fig6a.SIM" every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     "`echo $BASILISK`/test/starting/fig6a.KL" w l lw 2 t 'K and L 1995', \
     '< sort -k3,4 "$BASILISK/test/starting/level-15/cp-0.9" | awk -f "$BASILISK/test/starting/surface.awk"' w l lw 2 \
     t 'Basilisk (15 levels)', \
     '< sort -k3,4 "$BASILISK/test/starting/level-16/cp-0.9" | awk -f "$BASILISK/test/starting/surface.awk"' w l lw 2 \
     t 'Basilisk (16 levels)',\
     'level-15/cp-0.9'  every 10 u ($3/pi):($4*2) w l lw 2 \
     t 'Basilisk-dlmfd (15 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=2.5$. Re = 9500. See Fig. 6b of [Mohahegh et al. 2017](#mohahegh2017).
plot [0:2] "`echo $BASILISK`/test/starting/fig6b.SIM" every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     "`echo $BASILISK`/test/starting/fig6b.KL" w l lw 2 t 'K and L 1995', \
     '< sort -k3,4 "$BASILISK/test/starting/level-15/cp-2.5" | awk -f "$BASILISK/test/starting/surface.awk"' w l lw 2 \
     t 'Basilisk (15 levels)', \
     '< sort -k3,4 "$BASILISK/test/starting/level-16/cp-2.5" | awk -f "$BASILISK/test/starting/surface.awk"' w l lw 2 \
     t 'Basilisk (16 levels)',\
     'level-15/cp-2.5'  every 10 u ($3/pi):($4*2) w l lw 2 \
     t 'Basilisk-dlmfd (15 levels)'
~~~

## References

~~~bib
@article{bouard1980,
  title={The early stage of development of the wake behind an 
         impulsively started cylinder for 40 < {Re} < 10^4^},
  author={Bouard, Roger and Coutanceau, Madeleine},
  journal={Journal of Fluid Mechanics},
  volume={101},
  number={3},
  pages={583--607},
  year={1980},
  publisher={Cambridge University Press}
}

@article{koumoutsakos1995,
  title={High-resolution simulations of the flow around an 
         impulsively started cylinder using vortex methods},
  author={Koumoutsakos, Petros and Leonard, A},
  journal={Journal of Fluid Mechanics},
  volume={296},
  pages={1--38},
  year={1995},
  publisher={Cambridge University Press}
}

@article{mohaghegh2017,
  title={Comparison of sharp and smoothed interface methods for simulation
         of particulate flows II: Inertial and added mass effects},
  author={Mohaghegh, Fazlolah and Udaykumar, HS},
  journal={Computers \& Fluids},
  volume={143},
  pages={103--119},
  year={2017},
  publisher={Elsevier}
}
~~~
*/
