/**
# Cylinder freely rotating in a sheared-flow with DLMFD
*/

/** A cylinder twice as heavy as the surrounding fluid is freely rotating in a shear-flow. */

# define LEVEL 7 
# include "grid/quadtree.h"
# define DLM_Moving_particle 1
# define TRANSLATION 0
# define ROTATION 1
# define DLM_alpha_coupling 1
# define NPARTICLES 1

# define adaptive 1
# define MAXLEVEL (LEVEL)


/** # Physical parameters */
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (0.4) // particle diameter
# define ReD  (10) // Reynolds number based on the particle's diameter
# define Ldomain 1.0
# define rhosolid 2.0 //particle density
# define tval (rhoval*Uc*Ldomain/ReD)
 

/** # Output and numerical parameters */
# define Tc (diam/Uc) // caracteristic time scale
# define mydt (Tc/200) // maximum time-step
# define maxtime (2.5)
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
 
  if (!restore (file = "dump")) {
    /* fluid initial condition: */
    foreach() {
      u.x[] = -Uc + 2*y*Uc/L0;
      un[] = u.x[];
    }
    /* initial condition: particles position */
    particle * p = particles;
    for (int k = 0; k < NPARTICLES; k++) {
      GeomParameter gp = {0};
       
      gp.center.x = L0/2;
      gp.center.y = L0/2;
      gp.radius   = diam/2;
   
      p[k].g = gp;
   
      /* initial condition: particle's velocity */
      coord c = {0., 0., 0.};
      p[k].w = c;

    }
  }
  else { // restart of a run, the default init event will take care of
	 // it

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
  scalar omega[];
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  vorticity(u, omega);
  statsvelox =  statsf (u.x);
  clear();
  squares ("u.x", map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells();
  save ("movie.mp4");
}

/** 
# Results
~~~gnuplot particle's angular velocity
set grid
show grid
plot "particle-data-0" u 1:6 w l
~~~

~~~gnuplot Torque
set grid
show grid
plot "sum_lambda-0" u 1:4 every::10 w l
~~~
*/

/**
# Annimation

## Streamwise velocity component: u.x
<video width="890" height="862" controls>
<source src="shear2D_rot/movie.mp4" type="video/mp4">
</video>
*/
