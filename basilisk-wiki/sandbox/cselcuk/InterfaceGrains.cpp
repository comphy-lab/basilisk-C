/**
# Coupling interface for Grains3D and Basilisk 
*/

#include "Grains_BuilderFactory.H"
#include "InterfaceGrains.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  static GrainsCoupledWithFluid * grains = NULL;
  
  void Init_Grains (double fluid_density, const bool b_restart,
        const bool b_initializeClonePer,
        const double grid_size,
	const bool is_solidsolver_parallel,
	const int my_rank, const int nb_procs) {

    string simulation_file = "Grains/simul.xml";
    
    ReaderXML::initialize();
        
    string simulation_file_exe = Grains_BuilderFactory::init (simulation_file, my_rank, is_solidsolver_parallel ? nb_procs : 1);
  
    DOMElement* rootNode = ReaderXML::getRoot (simulation_file_exe);
    
    // cout << "simulation_file_exe=:"<< simulation_file_exe  <<endl;
    // cout << "rootNode: " << rootNode << " on pid()= " << my_rank  << endl;
    
    grains = Grains_BuilderFactory::createCoupledWithFluid (rootNode, fluid_density, grid_size);

    // cout << "grain = Grains_BuilderFactory:" << grains << "pid()=" << my_rank  << endl;

    if ( b_restart ) grains->setReloadSame();
    grains->Construction (rootNode);
    // cout << "grain->Construction ok" << "on pid()= " << my_rank  << endl;
    grains->Forces (rootNode);
    // cout << "grain->Forces ok" << "on pid()= " << my_rank  << endl;
    grains->Chargement (rootNode);
    // cout << "grain->Chargement ok " << "on pid()= " << my_rank  << endl;
    if ( b_initializeClonePer ) grains->initializeClonesPeriodiques();

    ReaderXML::terminate();
     
    cout << "Construction of Grains completed" << endl;
  }

  void Simu_Grains (const bool predictor, const bool isPredictorCorrector, const bool explicit_added_mass) {
    grains->Simulation (predictor, isPredictorCorrector, explicit_added_mass);
	
  }

  void Data_GrainsToCstruct (struct BasiliskDataStructure * b, const int m) {
    grains->WriteParticulesInFluid (b);
  }

  void Data_CstructToGrains (struct BasiliskDataStructure * b) {};

  void Setdt_Grains (double const dtfluid) {
    grains->set_timeStep (dtfluid);
  }
  
  void Setdt_AddedMassGrains (double dtfluid) {
    grains->set_ExplicitAddedMasstimeStep (dtfluid);
  } 

  void Setviscosity_Grains (double const viscosity) {
    grains->setFluidViscosity (viscosity);
  }

  void SaveResults_Grains (const int counter, int * restarted_simu) {
    int cc = counter ;
    if ((cc == 0) || (*restarted_simu)) {
      if (*restarted_simu)
	cc--;
      grains->setInitialCycleNumber (cc);
      grains->InitialPostProcessing();

      *restarted_simu = 0;
    }
    else
      grains->doPostProcessing();

  }

  /* Check that Paraview writer is activated 
     ------------------------------------------*/
  void checkParaviewPostProcessing_Grains (char * solid_resDir ) {
    
    grains->checkParaviewPostProcessing ("grains", solid_resDir, true);
  }

  void Update_Velocity_Grains (struct BasiliskDataStructure * b, bool explicit_added_mass) {
    grains->UpdateParticulesVelocities (b, explicit_added_mass);
  }

   void ActivateExplicitAddedMass (bool restart) {
    grains->AddExplicitAddedMass( restart, "" );
  }
  
  void InitializeExplicitAddedMass (bool restart, char * solid_resDir) {
    grains->InitializeExplicitAddedMassRestart (restart, solid_resDir);
  }
  
#ifdef __cplusplus
}
#endif
