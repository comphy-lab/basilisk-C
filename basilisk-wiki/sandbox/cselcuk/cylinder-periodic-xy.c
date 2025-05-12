/**
# Cylinder moving in the (x+y)-directions in a fully periodic domain with DLMFD
*/

/** A cylinder twice as heavy as the surrounding fluid is advected in
  the xy-directions by a pressure gradient */

# define LEVEL 7 
# include "grid/quadtree.h"
# define DLM_Moving_particle 1
# define adaptive 1
# define NPARTICLES 1

# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif

/**
# Physical parameters 
*/
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (1.) // particle diameter
# define ReD  200.0 // Reynolds number based on the particle's diameter
# define fs_density_ratio  2. // fluid solid density ratio
# define Ld_ratio 5. // box size-particle diameter ratio
# define Ldomain (Ld_ratio*diam)
# define rhosolid (fs_density_ratio*rhoval) //particle density
# define tval (rhoval*Uc*diam/ReD) // fluid dynamical viscosity
 
/** 
# Output and numerical parameters 
*/
# define Tc (diam/Uc) // caracteristic time scale
# define mydt (Tc/400) // maximum time-step (the time step is also adaptive in time but it won't exceed this value)
# define maxtime (4.)
# define tsave (Tc/200.)

/** 
We include the ficitious-domain implementation 
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
  
  /* boundary conditions */
  periodic(right); 
  periodic(top);

  /* Convergence criteria */
  TOLERANCE = 1e-3;
  
  run(); 
}

/**
We initialize the fluid and particle variables. */

event init (i = 0) {
  /* set origin */
  origin (0, 0);

  /* Initialize acceleration (face) vectors for pressure gradient to drive the flow */
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
  
  /* If new simulation: set fluid initial condition and initialise
     data file pointers */
  if (!restore (file = "dump")) {
    /* fluid initial condition: */
    foreach() {
      foreach_dimension()
	u.x[] = 0;
      un[] = u.x[];
    }

    /* initial condition: particles position */
    particle * p = particles;
    for (int k = 0; k < NPARTICLES; k++) {
      GeomParameter gp;
      gp.center.x = L0-diam/2;
      gp.center.y = gp.center.x;
      gp.radius = diam/2;
      p[k].g = gp;
	
      /* initial condition: particle's velocity */
      coord c = {0., 0., 0.};
      /* Translational velocity U */
      p[k].U = c;
      /* Rotational velocity w */
      p[k].w = c;
    }
  }
  else // Restart of a simulation
    {
      /* the default init=0 event will take of it */
    }
}

/* We log the number of iterations of the 
    multigrid solver for pressure and viscosity */
event logfile (i++) {
  deltau = change(u.x,un);
  fprintf (stderr, "%d %g %d %d %g\n", i, t, mgp.i, mgu.i, deltau);
}

event acceleration (i++) {
  face vector av = a;
  coord dp = {20./L0/rhoval, 20./L0/rhoval, 0};
  foreach_face(){
    av.x[] = dp.x;
  }
}

event output_data (t += tsave; t < maxtime) {
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  stats statsvelox;
  scalar normu[];
  foreach() {
    normu[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }
  statsvelox = statsf (normu); 
  clear();
  squares ("normu", map = cool_warm, min = statsvelox.min, max = statsvelox.max);
  cells();
  save ("movie.mp4");
}


event lastdump (t = maxtime) {
  dump(file = "dump");
}


/** 
# Results
~~~gnuplot particle's x-position
plot "particle-data-0" u 1:2 w l
~~~

~~~gnuplot particle's streamwise velocity
plot "particle-data-0" u 1:4 w l
~~~

~~~gnuplot particle's y-position
plot "particle-data-0" u 1:3 w l
~~~

~~~gnuplot particle's normal-streamwise velocity
plot "particle-data-0" u 1:5 w l
~~~

~~~gnuplot particle's angular velocity
plot "particle-data-0" u 1:6 w l
~~~
*/

/**
# Animation
norm of the velocity field
<video width="890" height="862" controls>
<source src="cylinder-periodic-xy/movie.mp4" type="video/mp4">
</video>
*/
