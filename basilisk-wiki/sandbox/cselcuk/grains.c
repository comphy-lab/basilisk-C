/** 
# Coupling Grains3D (a c++ granular solver) with basilisk
*/

# include "grid/octree.h"
# define DLM_Moving_particle 1
# define NPARTICLES 4

/* Coupling Interface for Grains3D */
# include "InterfaceGrains.h"

/* Fictitious domain implementation */
# include "DLMFD_reverse_Uzawa.h"

/* Additional helper functions for the coupling Grains3D */
# include "BasiliskGrainsCouplingFunctions.h"

/* Parameters related to the patial resolution*/
# define LEVEL 6
# define adaptive 1

# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif

/* Physical parameters */
# define rhoval 1.
# define tval 0.05
# define Uc 1

/* output and numerical parameters */
# define mydt  (0.002)
# define maxtime (1500*mydt)

double deltau;
scalar un[];
int restarted_simu = 0;

struct BasiliskDataStructure BasiliskData[NPARTICLES];

int main () {
  
  origin (0., 0., 0.);
 
  L0 = 1.;
   
  // set time step
  DT = mydt;
  
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
  
  // Convergence criteria
  TOLERANCE = 1e-3; 
  
  run();
}


event init (i = 0) {
  
  particle * p = particles;
  
  const scalar rhoc[] = rhoval;
  rho = rhoc;

  const face vector muc[] = {tval, tval, tval};
  mu = muc;

  /* if new simulation fluid initial
     condtion and data file pointers */
  if (!restore (file = "dump")) {

    /* Initial condition for the fluid */
    foreach() {
      foreach_dimension() {
	u.x[] = Uc;
      }
      un[] = u.x[];
    }

    /* Initial condition: particles/fictitious-domains positions */
    init_file_pointers (pdata, fdata, 0);
  }
  else { 
    restore_particle (p);

    restarted_simu = 1;
    
    /* Initialize file poiters */
    init_file_pointers(pdata, fdata, 1);
    fprintf (ferr, "simulation restored: it has to go for t_end = %f\n", maxtime);
  }

  /* Initialize Grains with its needed parameters */
  bool b_restart = false;
  bool b_intitializeClonePer = false;
  double grid_size = 0.;
  bool is_solidsolver_parallel = false;
  int my_rank = pid();
  int nb_procs = npe();

  if (pid() == 0) {
    Init_Grains (rhoval, b_restart,
		 b_intitializeClonePer,
		 grid_size,
		 is_solidsolver_parallel, my_rank, nb_procs);
      
    /*Set viscosity for grains3D */
    Setviscosity_Grains(tval);
     
    /* Transfer the data to the common C structure */
    Data_GrainsToCstruct (&BasiliskData[0], NPARTICLES);
      
    /* Activate the Paraview writer for grains3D */
    checkParaviewPostProcessing_Grains ("Grains/Res");
  }
      
  /* Update Basilisk particle structure (mpi ok) */ 
  UpdateParticlesBasilisk (&BasiliskData[0], p, NPARTICLES);

  /* Unallocate the BasiliskDataStructure used for Grains-Basilisk
     communication.  At this point Basilisk has all the particle data
     in the structure particles */
  unallocateBasiliskDataStructure(&BasiliskData[0], NPARTICLES);

  /* Set the particle types (for now only spheres and cubes are implemented) */
  for (int k = 0; k < NPARTICLES; k++) {   
    p[k].iswall = 0;
    GeomParameter * gg;
    gg = &(p[k].g);
    if (gg->ncorners == 8)
      p[k].iscube = 1;
    p[k].pnum = k;
  }

}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++) {
  deltau = change (u.y, un);
  fprintf (stderr, "%d %g %d %d %g \n", i, t, mgp.i, mgu.i, deltau);
  if (t > maxtime){
    return 1;
  }
}


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

event adapt (i++) {
  
  particle * pp = particles;

  /* Predictor step: pure granular problem with Grains (Grains turn in
     sequential only) */
  if (pid() == 0) {
    printf("calling Grains at iteration %d\n", i);
    Setdt_Grains(dt);
    Simu_Grains (true, false, false);
    
    Data_GrainsToCstruct (&BasiliskData[0], NPARTICLES);
  }

  UpdateParticlesBasilisk (&BasiliskData[0], pp, NPARTICLES);
   
  unallocateBasiliskDataStructure(&BasiliskData[0], NPARTICLES);

  /* Fictitious domain problem */
  DLMFD_subproblem (pp, i, rhoval);

  /* Save the forces acting on particles before adapting the mesh */
  sumLambda (pp, fdata, t, dt, flagfield, DLM_lambda, index_lambda, rhoval);

  /* We free particle's dlmfd points (we dont need them anymore) */
  free_particles (pp, NPARTICLES);
  
  /* We save all particles trajectories */
  particle_data (pp, t, i, pdata);

  /* Call the adapt_wavelet function to adapt the mesh */
#if adaptive
  astats s = adapt_wavelet ((scalar *){flagfield_mailleur, u}, (double[]){1e-3, 1e-2, 1e-2, 1e-2}, maxlevel = MAXLEVEL);
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


/**
## Results

### Annimation
<video width="1280" height="640" controls>
<source src="http://www.basilisk.fr/sandbox/cselcuk/pics/three_cubes_sphere_rho_ratio_8_falling.mp4" type="video/mp4">
</video>

*/
