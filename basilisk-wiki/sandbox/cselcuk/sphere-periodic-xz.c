/**
# Sphere advected in the (x+z)-directions in a tri-periodic domain with DLMFD
*/
# define LEVEL 6 
# include "grid/octree.h"
# define DLM_Moving_particle 1
# define NPARTICLES 1

# define adaptive 1
# define MAXLEVEL (LEVEL + 1)

/** 
# Physical parameters 
*/
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (1.) // particle diameter
# define ReD  200. // Reynolds number based on the particle's diameter
# define fs_density_ratio  2. // fluid solid density ratio
# define Ld_ratio 5. // box size-particle diameter ratio
# define Ldomain (Ld_ratio*diam)
# define rhosolid (fs_density_ratio*rhoval) //particle density
# define tval (rhoval*Uc*diam/ReD) // fluid dynamical viscosity
 

/**
# output and numerical parameters */
# define Tc (diam/Uc)  // caracteristic time scale
# define mydt (Tc/200) // maximum time-step 
# define maxtime (1.)
# define tsave (Tc/200.)

/** 
We include the fictitious-domain implementation with a toy model granular solver
*/
# include "dlmfd-toygs.h"
# include "view.h"

double deltau;
scalar un[];

int main() {
  L0 = Ldomain;

  /* set time step */
  DT = mydt;
     
  /* initialize grid */
  init_grid(1 << (LEVEL));
  
  /* boundary conditions: periodicity everywhere */
  foreach_dimension() 
    periodic(top);

  /* Convergence criteria */
  TOLERANCE = 1e-3;
  
  run();

  
}

/**
We initialize the fluid and particle variables. */

event init (i = 0) {
  /* set origin */
  origin (0., 0., 0.);

  /* Initialize acceleration (face) vectors for pressure gradient to drive the flow */
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
  
  /* fluid initial condition: */
  foreach() {
    /* foreach_dimension() */
      u.x[] = 0;
      un[] = u.x[];
  }

  /* initial condition: particles position */
  particle * p = particles;
   
  for (int k = 0; k < NPARTICLES; k++) {
    GeomParameter gp;
    gp.center.x = L0- (diam/2);
    gp.center.y = L0/2;
    gp.center.z = L0- (diam/2);
    gp.radius = diam/2;
    p[k].g = gp;
   
    /* initial condition: particle's velocity */
    coord c = {0., 0., 0.};
    p[k].U = c;
    p[k].w = c;
  }
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.x, un);
  fprintf (stderr, "log output %d %g %d %d %g %g %g\n", i, t, mgp.i, mgu.i, mgp.resa, mgu.resa, deltau);
  
}


event acceleration(i++) {
  face vector av = a;
  coord dp = {20./L0/rhoval, 0., 20./L0/rhoval};
  foreach_face(){
    av.x[] = dp.x;
  }
}


event output_data (t += tsave; t < maxtime) {
/* event output_data (i++; i < 20) {  */

  /* vue sur le plan z */
  /* view (fov = 29.0823, quat = {-0.5,0.5,-0.5,0.5}, tx = 0.5, ty = 0.5, bg = {0.3,0.4,0.6}, width = 804, height = 748, samples = 1); */

  /* vue 3D */
  /* view (fov = 26.8568, quat = {-0.733356,-0.326189,0.262887,0.535427}, tx = -0.66604, ty = 0.560201, bg = {0.3,0.4,0.6}, width = 886, height = 810, samples = 1); */
  /* vue 3D, xy */

  view (fov = 32.0833, quat = {-0.567402,-0.204857,-0.260493,-0.753812}, tx = -0.601431, ty = -0.431815, bg = {0.3,0.4,0.6}, width = 886, height = 810, samples = 1);
  stats statsvelox;

  scalar unorm[];
  foreach(){
    unorm[] = sqrt( sq(u.x[]) + sq(u.z[]));
  }
  statsvelox = statsf (unorm); 

  clear();
  box();
  squares ("unorm", n = {1,0,0},alpha = 0, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  squares ("unorm", n = {0,0,1},alpha = 0, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  squares ("unorm", n = {0,0,1},alpha = L0, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  squares ("unorm", n = {0,1,0},alpha = L0/2, map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  cells(n = {0,1,0},alpha = L0/2);
  save ("movie.mp4");
  dump(file = "dump");
}


/** 
# Results
~~~gnuplot particle's z-position
plot "particle-data-0" u 1:4 w l
~~~

~~~gnuplot particle's z-velocity
plot "particle-data-0" u 1:7 w l
~~~

~~~gnuplot particle's z-angular velocity
plot "particle-data-0" u 1:10 w l
~~~

~~~gnuplot particle's x-position
plot "particle-data-0" u 1:2 w l
~~~

~~~gnuplot particle's x-velocity
plot "particle-data-0" u 1:5 w l
~~~


~~~gnuplot particle's x-angular velocity
plot "particle-data-0" u 1:8 w l
~~~

~~~gnuplot particle's y-position
plot "particle-data-0" u 1:3 w l
~~~

~~~gnuplot particle's y-velocity
plot "particle-data-0" u 1:6 w l
~~~

~~~gnuplot particle's y-angular velocity
plot "particle-data-0" u 1:9 w l
~~~


*/

/**
# Animation

## Scalar field unorm = sq(u_x) + sq(u_z)
<video width="890" height="862" controls>
<source src="sphere-periodic-xz/movie.mp4" type="video/mp4">
</video>
*/
