/**
# Three overlapping cylinders freely rotating in a sheared-flow with DLMFD
*/

/** Three cylinders, overlapping each other and twice as heavy as the
surrounding fluid, are freely rotating in a linear shear-flow.

 */

# define LEVEL 6 
# include "grid/quadtree.h"
# define DLM_Moving_particle 1
# define TRANSLATION 0
# define ROTATION 1
# define DLM_alpha_coupling 0
# define adaptive 1
# define NPARTICLES 3

# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif


/** # Physical parameters */
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (0.5) // particle diameter
# define ReD  (10) // Reynolds number based on the particle's diameter
# define Ldomain 1.0
# define rhosolid 2.0 //particle density
# define tval (rhoval*Uc*Ldomain/ReD) //fluid's dynamic viscosity
 

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
  TOLERANCE = 1e-3;

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
    /* particles initial position */
    particle * p = particles;

    for (int k = 0; k < NPARTICLES; k++) {
      GeomParameter gp = {0};
    
      /** We place the cylinders close to each others to test how the
	  stencils are constructed when they are close to colision. */
      if (k == 0) {
	gp.center.x = L0/4;
	gp.center.y = L0/4;
	gp.radius   = diam/2;
      }

      if (k == 1) { 
	gp.center.x = 3*L0/4;
	gp.center.y = L0/4;
	gp.radius   = diam/2;

      }

      if (k == 2) { 
	gp.center.x = L0/2;
	gp.center.y = L0/4;
	gp.radius   = diam/2;

      }
    
      p[k].g = gp;
       
      /* particle's initial angular velocity */
      coord c = {0., 0., 0.};
      p[k].w = c;
    }
  }
  else // Restart of a simulation
    {
     /* the default init=0 event will take of it */
  }
}



event output_data (t += maxtime/50; t < maxtime) {
  stats statsvelox;
  view (fov = 22.3366, quat = {0, 0, 0, 1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);

  scalar constrained_cells[];
  int aa;
  foreach() {
    aa = 0;
    if ((int)index_lambda.y[] > -1 && flagfield[] < 1)
      aa = 1;
    constrained_cells[] = aa;
    constrained_cells[] += 3*flagfield[];
  }
  statsvelox =  statsf (constrained_cells);
  clear();
  squares ("constrained_cells", map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells();
  save ("movie.mp4");
}


event output_data_2 (t += maxtime/50; t < maxtime) {
  stats statsvelox;
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  statsvelox =  statsf (u.x);
  clear();
  squares ("u.x", map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells();
  save ("movie1.mp4");
}

event output_data_3 (t += maxtime/50; t < maxtime) {
  stats statsvelox;
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  statsvelox =  statsf (index_lambda.y);
  clear();
  squares ("index_lambda.y", map = cool_warm, min = statsvelox.min, max = statsvelox.max); 
  cells();
  save ("movie2.mp4");
}

/** 
# Results
~~~gnuplot particles angular velocities

set grid
show grid
set xlabel "time"
set ylabel "Omegaz"
plot "particle-data-0" u 1:6 w l, "particle-data-1" u 1:6 w l, "particle-data-2" u 1:6 w l 
~~~

~~~gnuplot Torque of the particles
set grid
show grid
set xlabel "time"
set ylabel "Tz"
plot "sum_lambda-0" u 1:3 w l, "sum_lambda-1" u 1:3 w l, "sum_lambda-2" u 1:3 w l

~~~

~~~gnuplot Number of iterations to solve the saddle problem with the Uzawa method

set grid
show grid
set xlabel "Time iteration"
set ylabel "Iterations to solve the saddle problem"
plot "converge-uzawa" u 1:2 w l
~~~
*/

/**
# Annimation
## Constrained cells by the Lagrange multipliers

Red cells: constrained by the set of Lagrange multipliers which lie on
the surface of the particle.

Light blue cells: constrained by Lagrange multipliers that are inside
the particle.

<video width="890" height="862" controls>
<source src="2D-3multidomain_overlapping_test-2/movie.mp4" type="video/mp4">
</video>

## Streamwise velocity u.x
<video width="890" height="862" controls>
<source src="2D-3multidomain_overlapping_test-2/movie1.mp4" type="video/mp4">
</video>

## Fictitious domains of the particles
<video width="890" height="862" controls>
<source src="2D-3multidomain_overlapping_test-2/movie2.mp4" type="video/mp4">
</video>
*/
