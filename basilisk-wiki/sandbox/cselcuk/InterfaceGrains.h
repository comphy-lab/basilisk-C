/** 
# Interface functions for Grains/Basilisk
*/

#ifndef INTERFACEGRAINS_H 
#define INTERFACEGRAINS_H 

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __cplusplus
#include "BasiliskInterfaceDataStructure.h"
#endif

  void Init_Grains (double fluid_density, const bool b_restart,
        const bool b_initializeClonePer,
        const double grid_size,
	const bool is_solidsolver_parallel,
	const int my_rank, const int nb_procs);
  
  void Simu_Grains (bool predictor, const bool isPredictorCorrector, const bool explicit_added_mass);
  
  void Data_GrainsToCstruct (struct BasiliskDataStructure * b, const int m);

  void Data_CstructToGrains (struct BasiliskDataStructure * b);

  void Setdt_Grains (double const dtfluid);

  void Setdt_AddedMassGrains (double dtfluid);
  
  void Setviscosity_Grains (double const viscosity);

  void SaveResults_Grains (const int counter, int * restarted_simu);

  void checkParaviewPostProcessing_Grains (char * solid_resDir );

  void Update_Velocity_Grains (struct BasiliskDataStructure * b, bool explicit_added_mass);

  void ActivateExplicitAddedMass (bool restart);
  
  void InitializeExplicitAddedMass (bool b_restart, char * solid_resDir);

#ifdef __cplusplus
}
#endif


#endif
