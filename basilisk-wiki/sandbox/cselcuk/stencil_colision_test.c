/**
# Two cylinders freely rotating in a sheared-flow with DLMFD
*/

/** Two cylinders, next to each other and twice as heavy as the
surrounding fluid, are freely rotating in a shear-flow.

We place the cylinders this way to test how the stencils are
constructed when they close to colision. */

# define LEVEL 7 
# include "grid/quadtree.h"
# define DLM_Moving_particle 1
# define TRANSLATION 0
# define ROTATION 1
# define DLM_alpha_coupling 1
# define adaptive 1
# define NPARTICLES 2

# if adaptive
# define MAXLEVEL (LEVEL)
# endif

/** # Physical parameters */
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (0.5) // particle diameter
# define ReD  (10) // Reynolds number based on the particle's diameter
# define Ldomain 1.0
# define rhosolid 2.0 //particle density
# define tval (rhoval*Uc*Ldomain/ReD) //fluid dynamical viscosity
 

/** # Output and numerical parameters */
# define Tc (diam/Uc) // caracteristic time scale
# define mydt (Tc/200) // maximum time-step
# define maxtime (1.0)
# define tsave (Tc/1.)

/** 
We include the ficitious-domain implementation 
*/
# include "dlmfd.h"
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
  u.n[left] = dirichlet(-Uc + 2*y*Uc/L0);
  u.t[left] = dirichlet(0);

  u.n[right] = dirichlet(-Uc + 2*y*Uc/L0);
  u.t[right] = dirichlet(0);
	
  u.n[top]   = dirichlet(0);
  u.t[top]   = dirichlet(Uc);

  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(-Uc);

   
  /* Convergence criteria */
  TOLERANCE = 1e-4;

  run();
}

/**
   We initialize the fluid and particle variables. */

event init (i = 0) {

  /* set origin */
  origin (0, 0);
  
  if (!restore (file = "dump")){
    /* fluid initial condition: */
    foreach() {
      u.x[] = -Uc + 2*y*Uc/L0;
      un[] = u.x[];
    }

    /* initial condition: particles position */
    particle * p = particles;

    for (int k = 0; k < NPARTICLES; k++) {
      GeomParameter gp = {0};
      /** We place the cylinders close to each others to test how the
	  stencils are constructed when they close to colision. */
      if (k == 0) {
	gp.center.x = L0/4;
	gp.center.y = L0/2;
	gp.radius   = diam/2;
      }
      else {
	gp.center.x = 3*L0/4;
	gp.center.y = L0/2;
	gp.radius   = diam/2;

      }
      p[k].g = gp;
      /* initial condition: particle's velocity */
      coord c = {0., 0., 0.};
      p[k].w = c;
    }
  } else {
    // restart run, the default init event will take care of it
  }

}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.x, un);
  fprintf (stderr, "log output %d %g %d %d %g %g %g %ld\n", i, t, mgp.i, mgu.i, mgp.resa, mgu.resa, deltau, grid->tn);
}

event output_data (t += 0.01; t < maxtime) {
  stats statsvelox;
  /* scalar omega[]; */
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  /* vorticity(u, omega); */
  statsvelox =  statsf (flagfield);
  clear();
  squares ("flagfield", map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells();
  save ("movie.mp4");
}


event output_data_2 (t += 0.01; t < maxtime) {
  stats statsvelox;
  /* scalar omega[]; */
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  /* vorticity(u, omega); */
  statsvelox =  statsf (u.x);
  clear();
  squares ("u.x", map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells();
  save ("movie1.mp4");
}

/** 
# Results
~~~gnuplot particles angular velocities

set grid
show grid
plot "particle-data-0" u 1:6 w l, "particle-data-1" u 1:6 w l 
~~~

~~~gnuplot Torque of the particles

set grid
show grid
plot "sum_lambda-0" u 1:3 w l, "sum_lambda-1" u 1:3 w l

~~~
*/

/**
# Annimation
## Constrained cells by the Lagrange multipliers
<video width="890" height="862" controls>
<source src="stencil_colision_test/movie.mp4" type="video/mp4">
</video>

## Streamwise velocity u.x
<video width="890" height="862" controls>
<source src="stencil_colision_test/movie1.mp4" type="video/mp4">
</video>
*/
