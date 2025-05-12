/**
Test of 1 light sphere rising in a box with Grains3D
*/

# include "grid/octree.h"
# define DLM_Moving_particle 1
# define NPARTICLES 1
# define DLM_alpha_coupling 1

/* Parameters related to the spatial resolution*/
# define LEVEL 4
# define adaptive 1

/* if adaptivity is enabled define the maximum level of refinement MAXLEVEL */
# if adaptive
# define MAXLEVEL (LEVEL + 4)
# endif
# define FlagAdaptCrit (1.E-5)
# define UxAdaptCrit (1.E-2)
# define UyAdaptCrit (1.E-2)
# define UzAdaptCrit (1.E-2)

/* Physical parameters */
# define rhoval 960. // fluid density
# define tval 0.058 // fluid dynamic viscosity
# define Uc 0.3431 // characteristic velocity
# define b_explicit_added_mass true // for buoyant and lighter particles
# define gravity_x 0. // x component of gravity acceleration
# define gravity_y 0. // x component of gravity acceleration
# define gravity_z -9.81 // x component of gravity acceleration

/* output and numerical parameters */
# define mydt (2.E-4) //time-step
# define maxtime (20000*mydt) // simulation time

/* 1st order Marshuk-Yanenko implementation of the DLM-FD method */
# include "dlmfd-grains.h"

int main () {

  origin (0., 0., 0.);

  L0 = 0.16; // lenght of the fluid-domain

  /* set time step */
  DT = mydt;

  /* initialise a uniform grid */
  init_grid (1 << LEVEL);

  /* boundary conditions */
  u.t[left]   = dirichlet(0.); //v
  u.r[left]   = dirichlet(0.); //w
  u.n[left]   = dirichlet(0.); //u

  u.t[right]  = dirichlet(0.); //v
  u.r[right]  = dirichlet(0.); //w
  u.n[right]  = dirichlet(0.); //u

  u.n[bottom] = dirichlet(0.); //v
  u.t[bottom] = dirichlet(0.); //w
  u.r[bottom] = dirichlet(0.); //u

  u.n[top]    = dirichlet(0.); //v
  u.t[top]    = dirichlet(0.); //w
  u.r[top]    = dirichlet(0.); //u

  u.n[front]  = dirichlet(0.); //w 
  u.t[front]  = dirichlet(0.); //u 
  u.r[front]  = dirichlet(0.); //v 

  u.n[back]   = dirichlet(0.); //w 
  u.t[back]   = dirichlet(0.); //u 
  u.r[back]   = dirichlet(0.); //v 

  //periodic(top);
  
  // Convergence criteria for the N-S solver
  TOLERANCE = 1e-3;

  run();
}

event init (i = 0) {
  fprintf (stderr, "Init: Nothing special to do \n");  
}

event adapt (i++) {
  fprintf (stderr, "Adapt: Nothing special to do \n"); 
}

/* output in terms of iteration */
event output_data (t += maxtime/25; t <= maxtime) {
  char name[80] = "allfields";

  scalar * list = {p};
  vector * vlist = {u};
  save_data (list, vlist, i, t, name);

  dump (file = "dump");
  particle * p = particles;
  dump_particle (p);

  if (pid() == 0)
    SaveResults_Grains (i, &restarted_simu) ;
}

/* event lastdump (t = end) { */
/*   dump(file = "dump"); */
/*   particle * pp = particles; */
/*   dump_particle(pp); */

/*   char name[80] = "allfields"; */
/*   scalar * list = {p}; */
/*   vector * vlist = {u}; */
/*   save_data(list, vlist, i, t, name); */

/*   if (pid() == 0) */
/*     SaveResults_Grains(i, &restarted_simu); */
/* } */

/* We log the number of iterations of the 
    multigrid solver for pressure and viscosity */
event logfile (i++) {
  fprintf (stderr, "%d %g %d %d \n", i, t, mgp.i, mgu.i);
}
