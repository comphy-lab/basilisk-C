/** 
# Coupling Grains3D (a c++ granular solver) with basilisk
## Flow past a cube at Re  = 20 (based on cube's length)
*/

# include "grid/octree.h"
# define DLM_Moving_particle 0
# define NPARTICLES 1

/* Coupling Interface for Grains3D */
# include "InterfaceGrains.h"

/* Fictitious domain implementation */
# include "DLMFD_reverse_Uzawa.h"

/* Additional helper functions for the coupling Grains3D */
# include "BasiliskGrainsCouplingFunctions.h"

/* Parameters related to the patial resolution*/
# define LEVEL 6
# define adaptive 1

/* if adaptivity is enabled define the maximum level of refinement MAXLEVEL */
# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif

/* Physical parameters */
# define rhoval 1. // fluid density
# define tval 0.01 // fluid dynamic viscosity 
# define Uc 1 // caracteristic velocity 

/* output and numerical parameters */
# define mydt  (0.002) //time-step
# define maxtime (1500*mydt) // simulation time

double deltau;
scalar un[];
int restarted_simu = 0;

struct BasiliskDataStructure BasiliskData[NPARTICLES];

int main () {
  
  origin (0., 0., 0.);
 
  L0 = 1.;
   
  /* set time step */
  DT = mydt;

  /* initialise a uniform grid */
  init_grid (1 << LEVEL);

  /* boundary conditions */
  u.t[left]   = dirichlet(0.); //v
  u.r[left]   = dirichlet(0.); //w
  u.n[left]   = dirichlet(Uc); //u
  
  u.t[right]  = dirichlet(0.); //v
  u.r[right]  = dirichlet(0.); //w
  u.n[right]  = dirichlet(Uc); //u
 
  u.n[bottom] = dirichlet(0.); //v
  u.t[bottom] = dirichlet(0.); //w
  u.r[bottom] = dirichlet(Uc); //u
  
  u.n[top]    = dirichlet(0.); //v
  u.t[top]    = dirichlet(0.); //w
  u.r[top]    = dirichlet(Uc); //u
  
  u.n[front]  = dirichlet(0.); //w
  u.t[front]  = dirichlet(Uc); //u
  u.r[front]  = dirichlet(0.); //v
  
  u.n[back]   = dirichlet(0.); //w
  u.t[back]   = dirichlet(Uc); //u
  u.r[back]   = dirichlet(0.); //v
  
  // Convergence criteria for the N-S solver
  TOLERANCE = 1e-3; 
  
  run();
}


event init (i = 0) {
  
  particle * p = particles;
  
  const scalar rhoc[] = rhoval;
  rho = rhoc;

  const face vector muc[] = {tval, tval, tval};
  mu = muc;

  /* If new simulation: set fluid's initial
     condition and initialise data file pointers */
  if (!restore (file = "dump")) {

    /* Initial condition for the fluid */
    foreach() {
      foreach_dimension() {
	u.x[] = Uc;
      }
      un[] = u.x[];
    }

    /* Initialise data file pointers */
    init_file_pointers (pdata, fdata, 0);
  }
  else { 
    restore_particle (p);

    restarted_simu = 1;
    
    /* Re-initialize file poiters */
    init_file_pointers(pdata, fdata, 1);
    fprintf (ferr, "simulation restored: it has to go for t_end = %f\n", maxtime);
  }

  /* Initialize Grains with its parameters */
  bool b_restart = false;
  bool b_intitializeClonePer = false;
  double grid_size = 0.;
  bool is_solidsolver_parallel = false;
  int my_rank = pid();
  int nb_procs = npe();

  /* Grains runs in sequential */
  if (pid() == 0) {
    Init_Grains (rhoval, b_restart,
		 b_intitializeClonePer,
		 grid_size,
		 is_solidsolver_parallel, my_rank, nb_procs);
    
    /* Transfer the data to the common C structure */
    Data_GrainsToCstruct (&BasiliskData[0], NPARTICLES);
      
    /* Check that Paraview writer is activated (is this needed ?) */
    checkParaviewPostProcessing_Grains ("Grains/Res");
  }
      
  /* Update Basilisk particle structure */ 
  UpdateParticlesBasilisk (&BasiliskData[0], p, NPARTICLES);

  /* Unallocate the BasiliskDataStructure used for Grains-Basilisk
     communication.  At this point Basilisk has all the particle data
     in the structure particles */
  unallocateBasiliskDataStructure(&BasiliskData[0], NPARTICLES);

  /* Set the particle types (spheres and cubes only for now) */
  for (int k = 0; k < NPARTICLES; k++) {   
    p[k].iswall = 0;
    GeomParameter * gg;
    gg = &(p[k].g);
    if (gg->ncorners == 8)
      p[k].iscube = 1;
    
    p[k].pnum = k;
#if DLM_Moving_particle
    p[k].gravity.x = 0.;
    p[k].gravity.y = 0.;
    p[k].gravity.z = 0.;
#endif
  }

}

/** We log the number of iterations of the multigrid solver for
pressure and viscosity. */

event logfile (i++) {
  deltau = change (u.y, un);
  fprintf (stderr, "%d %g %d %d %g \n", i, t, mgp.i, mgu.i, deltau);
  if (t > maxtime){
    return 1;
  }
}

/** Overloading the viscount event to add an explicit coupling term
 that improves the resolution of two (uncoupled) sub-problems */
#if DLM_alpha_coupling
event viscous_term (i++) {
  foreach() {
    foreach_dimension() {
      u.x[] += -dt*DLM_explicit.x[]/(rho[]*dv());
    }
  }
  boundary((scalar*){u});
}
#endif

/** Overloading the adapt event to plug-in the granular and fictitious
 domain problem */
event adapt (i++) {
  
  particle * pp = particles;

  /* Predictor step: pure granular problem with Grains (Grains works
     in sequential only) */
  if (pid() == 0) {
    printf("calling Grains at iteration %d\n", i);
    Setdt_Grains(dt);
    Simu_Grains (true, false, false);
    
    Data_GrainsToCstruct (&BasiliskData[0], NPARTICLES);
  }
  UpdateParticlesBasilisk (&BasiliskData[0], pp, NPARTICLES);

  /* Fictitious domain problem */
  DLMFD_subproblem (pp, i, rhoval);

  UpdateBasiliskStructure (&BasiliskData[0], pp, NPARTICLES);
  
  if (pid() == 0)
    Update_Velocity_Grains (&BasiliskData[0]);

  unallocateBasiliskDataStructure(&BasiliskData[0], NPARTICLES);
  
  /* Save the forces acting on particles before adapting the mesh */
  sumLambda (pp, fdata, t, dt, flagfield, DLM_lambda, index_lambda, rhoval);

  /* We free particle's dlmfd points (we dont need them anymore) */
  free_particles (pp, NPARTICLES);
  
  /* We save all particles trajectories */
  particle_data (pp, t, i, pdata);

  /**
     Call the adapt_wavelet function to adapt the mesh
  */
#if adaptive
#if DLM_Moving_particle
  astats s = adapt_wavelet ((scalar *){flagfield_mailleur, u}, (double[]){1e-5, 1e-2, 1e-2, 1e-2}, maxlevel = MAXLEVEL);
#else
  astats s = adapt_wavelet ((scalar *){flagfield, u}, (double[]){1e-5, 1e-2, 1e-2, 1e-2}, maxlevel = MAXLEVEL);
#endif
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
#endif
}


event output_data (t += 0.0125; t <= maxtime)
{
  char name[80] = "allfields";
  scalar * list = {p, index_lambda.x};
  vector * vlist = {u};
  save_data(list, vlist, i, t, name);

  dump(file = "dump");
  particle * p = particles;
  dump_particle(p);

  if (pid() == 0)
    SaveResults_Grains(i, &restarted_simu);
}

event output_perf (i+=100; t <= maxtime) {
  particle * p = particles;
  if (i>10) output_dlmfd_perf (dlmfd_globaltiming, i, p);
 

}
event lastdump (t = end) {
  dump(file = "dump");
  particle * p = particles;
  dump_particle(p);

  if (pid() == 0)
    SaveResults_Grains(i, &restarted_simu);
}

/**
## Results
### Annimation
<video width="1280" height="640" controls>
<source src="http://www.basilisk.fr/sandbox/cselcuk/pics/flow_past_a_cube_Re_20.mp4" type="video/mp4">
</video>
*/
