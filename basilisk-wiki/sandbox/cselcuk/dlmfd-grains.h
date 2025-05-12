/* The general DLMFD plugin */ 
# include "dlmfd-plugin.h"

/* Coupling Interface for Grains3D */
# include "InterfaceGrains.h"

/* Additional helper functions for the coupling with Grains3D */
# include "BasiliskGrainsCouplingFunctions.h"

struct BasiliskDataStructure BasiliskData[NPARTICLES];

/* Here we overload the generic events defined in the general DLMFD plugin
   DLMFD.h such that it uses Grains3D as a granular solver */


/** Overloading of the granular solver init event */
// -------------------------------------------------
event GranularSolver_init (t < -1.)
{
    // Initialize Grains with its parameters 
  bool b_restart = false;
  bool b_intitializeClonePer = false;
  double grid_size = 0.;
  bool is_solidsolver_parallel = false;
  int my_rank = pid();
  int nb_procs = npe();

  // Grains runs in sequential 
  if (pid() == 0)
    {
      // Output the call to Grains3D
      printf ("Grains3D\n");
    
      // Initialize Grains
      Init_Grains (rhoval, b_restart, b_intitializeClonePer, grid_size,
		    is_solidsolver_parallel, my_rank, nb_procs);

      // Activate explicit added mass if solid density < fluid density
      if (b_explicit_added_mass)
	{
	  ActivateExplicitAddedMass (b_restart); 
	  InitializeExplicitAddedMass (b_restart, "Grains/Res");
	}
    
      // Transfer the data to the common C structure
      Data_GrainsToCstruct (&BasiliskData[0], NPARTICLES);

      // Check that Paraview writer is activated (is this needed ?)
      checkParaviewPostProcessing_Grains ("Grains/Res");
    }


  // Update Basilisk particle structure
  UpdateParticlesBasilisk (&BasiliskData[0], particles, NPARTICLES);

  // Unallocate the BasiliskDataStructure used for Grains-Basilisk
  // communication.  At this point Basilisk has all the particle data
  // in the structure particles 
  unallocateBasiliskDataStructure (&BasiliskData[0], NPARTICLES);

  
} 

/** Overloading of the granular solver predictor event */
// ------------------------------------------------------
event GranularSolver_predictor (t < -1.)
{
  // Predictor step: pure granular problem solved by Grains (Grains works
  // in serial only )
  if (pid() == 0) 
    {
      // Output the call to Grains3D
      printf ("# granular solver predictor step at iter %d: calling Grains3D\n", i);
       
      // Set the fluid time step magnitude in Grains3D
      Setdt_Grains (dt);
    
      // Run the granular simulation
      Simu_Grains (true, false, b_explicit_added_mass);

      // Transfer the data to the common C structure
      Data_GrainsToCstruct (&BasiliskData[0], NPARTICLES);
    
      // Set dt for explicit mass calculation in Grains3D at next time step
      if (b_explicit_added_mass) Setdt_AddedMassGrains (dt) ;
    }

  /* printBasiliskDataStructure(&BasiliskData[0]); */
    
  // Update Basilisk particle structure
  UpdateParticlesBasilisk (&BasiliskData[0], particles, NPARTICLES);
}

/** Overloading of the granular solver velocity update event */
// ------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.)
{
  // Output the call to Grains3D
  if (pid() == 0) printf ("# granular solver velocity update at iter %d: sending velocity to Grains3D\n", i);

  // Update Basilisk particle structure  
  UpdateBasiliskStructure (&BasiliskData[0], particles, NPARTICLES);

  // Update particles velocity on the granular side
  if (pid() == 0)
    Update_Velocity_Grains (&BasiliskData[0], b_explicit_added_mass);

  // Unallocate the BasiliskDataStructure used for Grains-Basilisk communication
  unallocateBasiliskDataStructure (&BasiliskData[0], NPARTICLES);
}
