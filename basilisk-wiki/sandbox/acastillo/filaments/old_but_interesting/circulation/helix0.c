#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "lambda2.h"
#include "view.h"
#include "PointTriangle.h"

#define sortme(x) a = x[i]; x[i] = x[j]; x[j] = a;
#define MINLEVEL 5

/**
# Simulated Navier-Stokes helical vortex

The evolution of a helical vortex is simulated as in
[Antoon's sandbox](http://basilisk.fr/sandbox/Antoonvh/helical.c). The goal
of the physical space initialisation is to map an analytically defined trefoil
vortex onto an Eulerian (static) numerical mesh. The helical trajectory
discussed in here is defined by:
$$
x = R \cos{(t)}
$$
$$
y = R \sin{(t)}
$$
$$
z = \frac{H}{2\pi} t
$$
where $t$ ranges between 0 and $2\pi$, $R$ is the radius and $H$ the pitch.

![Space-curve of a helical vortex relative to the domain size](helix/helix0.png)
*/

int n_seg;
double Gamma= 1.0, as = 0.1, Hs = 1.0, Rs = 0.5, n_turns;

#include "../filaments.h"
#include "../filaments.c"
#include "../../output_fields/output_matrix_normal.h"

void helical_filament(vector omega, double Gamma){
  // Define the helix as the spacecurve c
  double delta_t0 = (2 + n_turns)*2*pi/((double)n_seg-1);
  double t0[n_seg];
  coord c[n_seg], dc[n_seg], d2c[n_seg];

  for (int i = 0; i < n_seg; i++){
    t0[i] = delta_t0 * (double)i - 2*pi;
    c[i].x = Rs * cos(t0[i]);
    c[i].y = Rs * sin(t0[i]);
    c[i].z = (Hs/(2*pi)) * t0[i] - L0/2 ;

    dc[i].x = -Rs * sin(t0[i]);
    dc[i].y =  Rs * cos(t0[i]);
    dc[i].z = (Hs/(2*pi));

    d2c[i].x = -Rs * cos(t0[i]);
    d2c[i].y = -Rs * sin(t0[i]);
    d2c[i].z = 0.;
  }

  view (camera="iso", fov=40);
  draw_space_curve(n_seg, c);
  box();
  save ("helix0.png");

  // Use the first and second derivatives to compute a Frenet-Serret Frame
  coord tvec[n_seg], nvec[n_seg], bvec[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec[i].x =  dc[i].x/sqrt(vecdot( dc[i],  dc[i]));
      nvec[i].x = d2c[i].x/sqrt(vecdot(d2c[i], d2c[i]));
    }
    bvec[i] = vecdotproduct(tvec[i], nvec[i]);
  }

  // Initialize the vorticity field
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega = get_vorticity_filament2(1.0, 1.0, pcar, n_seg, as, t0, c, tvec, nvec, bvec, 0);
    foreach_dimension()
      omega.x[] = Gamma*val_omega.x;
  }
}

/**
For this example, we use the (compressible) Navier-Stokes equations inside a
triply periodic box. Periodicity implies that the net circulation contained
within the domain must be precisely zero, which can trigger an artificial
centrifugal instability of an isolated vortex (see,
[Pradeep and Hussain](DOI: 10.1017/S002211200400076X)).
We use the approach by [Otheguy, Chomaz, and Billant](doi:10.1017/S0022112006008901)
and work in a rotating-frame with angular velocity $\Omega^F$, such that
the circulation in the rotating frame is zero.
$$
\Omega^F = \Gamma_0/2L_xL_y
$$
*/

coord omega_mean={0,0,0} ;
#include "moments3d.h"



int main(){
  L0 = 2;
  X0 = Y0 = Z0 = -L0/2;
  omega_mean = (coord) {0,0,(-0.5*1.0/sq(L0))};

  init_grid (1 << MINLEVEL);

  //TOLERANCE = 1e-6;

  n_turns = L0/Hs;
  n_seg = 64 * n_turns;

  periodic(top);
  periodic(left);
  periodic(front);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-4,1e-4,1e-4}, MAXLEVEL, MINLEVEL);
}

event output (t = 5){
  printf( "Done!... \n" );
}

event init (t = 0){

  FILE * fp;
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Atot\t [4]Gamma\t [5]mu_x\t [6]mu_y\t [7]mu_z\t [8]M20\t [9]M02\t [10]M11\t [11]a\t [12]b\t [13]c\t [14]e\t [15]maxvor \n");
  fp = fopen("vortex_x0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_y0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_x0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_y0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("globals.asc", "w"); fclose (fp);
  fp = fopen("closed_paths.asc", "w"); fclose (fp);

#ifdef _INIT
  printf( "Creating initial conditions... \n" );

  refine  ((sqrt(sq(x) + sq(y))) < 16*Rs && level < MAXLEVEL-3);
  refine  ((sqrt(sq(x) + sq(y))) < 8*Rs && level < MAXLEVEL-2);
  refine  ((sqrt(sq(x) + sq(y))) < 4*Rs && level < MAXLEVEL-1);
  refine  ((sqrt(sq(x) + sq(y))) < 2*Rs && level < MAXLEVEL);

  vector c[], b[];
  helical_filament(b, Gamma);

  double bsum = 0.;
  foreach(reduction(+:bsum)) {
    bsum += b.z[]*dv();
  }
  fprintf (stderr, "bsum (before): %g\n", bsum);
  fprintf (stderr, "omega_mean (before): %g\n", omega_mean.z);

  omega_mean = (coord){0,0,-0.5*bsum/pow(L0,3)};
  fprintf (stderr, "omega_mean: %.15g\n", omega_mean.z);

  foreach()
    foreach_dimension()
      c.x[] = 0.;

  double bsum_x = 0., bsum_y = 0., bsum_z = 0.;
  foreach(reduction(+:bsum_x), reduction(+:bsum_y), reduction(+:bsum_z)) {
    double r = sqrt(sq(x) + sq(y));
    c.z[] -= sq(r)*omega_mean.z/2.0;
    b.z[] += 2.0*omega_mean.z;
    foreach_dimension()
      bsum_x += b.x[]*dv();
  }
  boundary ((scalar *){c,b});

  foreach_dimension()
    fprintf (stderr, "bsum: %g\n", bsum_x);

  /** The Poisson equation is solved. */
  foreach_dimension(){
    timer t = timer_start();
    mgstats s = poisson (c.x, b.x, tolerance=1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
  }
  boundary ((scalar *){c});

  foreach(){
    u.x[] = -((c.z[0,1,0] - c.z[0,-1,0]) - (c.y[0,0,1] - c.y[0,0,-1]))/(2.*Delta);
    u.y[] = -((c.x[0,0,1] - c.x[0,0,-1]) - (c.z[1,0,0] - c.z[-1,0,0]))/(2.*Delta);
    u.z[] = -((c.y[1,0,0] - c.y[-1,0,0]) - (c.x[0,1,0] - c.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  vector omega[];
  vorticity3d (u, omega);

  vector psi[];
  foreach()
    foreach_dimension()
      psi.x[] = c.x[];
  boundary ((scalar *){psi, omega});

  foreach_dimension(){
    timer t = timer_start();
    mgstats s = poisson (psi.x, omega.x, tolerance = 1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
  }

  scalar e1[], e2[];
  foreach(){
    e1[] = c.z[]-psi.z[];
    e2[] = b.z[]-omega.z[];
  }

  dump ("dump0", list = (scalar *){u});

  /**
  The solution is displayed using bview.
  */


  squares ("c.z", spread = -1);
  save ("a.png");
  clear();

  squares ("b.z", spread = -1);
  save ("b.png");
  clear();

  squares ("psi.z", spread = -1);
  save ("psi.png");
  clear();

  squares ("omega.z", spread = -1);
  save ("omega.png");
  clear();

  squares ("u.x", spread = -1);
  save ("u.png");
  clear();

  squares ("u.y", spread = -1);
  save ("v.png");
  clear();

  squares ("u.z", spread = -1);
  save ("w.png");
  clear();

  squares ("e1", spread = -1);
  save ("e1.png");
  clear();

  squares ("e2", spread = -1);
  save ("e2.png");
  clear();

  cells();
  save ("cells.png");
  clear();

  scalar l2[];
  lambda2 (u, l2);
  squares ("l2", linear = false);
  isosurface ("l2", -1e-3); cells();
  save ("lambda2.png");

  printf( "Done! \n" );

  stats s0 ;
  s0 = statsf (c.z); fprintf (stderr, "a: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (b.z); fprintf (stderr, "b: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (psi.z); fprintf (stderr, "psi: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (omega.z); fprintf (stderr, "omega: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (e1); fprintf (stderr, "e(psi): %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (e2); fprintf (stderr, "e(omega): %g %g %g \n", s0.min, s0.max, s0.sum);
  fprintf (stderr, "\n");

#else
  printf( "Loading initial conditions... \n" );
  restore("dump0");
  printf( "Done! \n" );

#endif
}


