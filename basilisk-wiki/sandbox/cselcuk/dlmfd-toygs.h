/* The general DLMFD plugin */ 
# include "dlmfd-plugin.h"

/* Here we overload the generic events defined in the global DLMFD plugin
   dlmfd-plugin.h such that it uses my toy granular solver */


/** Overloading of the granular solver init event */
// -------------------------------------------------
event GranularSolver_init (t < -1.) {
  particle * pp = particles;
  for (int k = 0; k < NPARTICLES; k++) {
    /* Contact model parameters needed to setup the granular time-step */
    pp[k].wished_ratio = 0.1;
    pp[k].en = 1;
    pp[k].vzero = 1;
    compute_wo (&pp[k]);

    /* add this term to make sure that the gravity is added only once in the granular subproblem (it is already present in the granular solver) */
    foreach_dimension()
      pp[k].adforce.x = -pp[k].gravity.x;
  }
} 

/** Overloading of the granular solver predictor event */
// ------------------------------------------------------
event GranularSolver_predictor (t < -1.) {
  if ( pid() == 0 ) printf("Prediction: my rk-4 toy granular solver\n");
  particle * pp = particles;
  granular_subproblem (pp, 1, dt, rhoval);
}

/** Overloading of the granular solver velocity update event */
// ------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.) {
  if ( pid() == 0 ) printf("Correction: my rk-4 toy granular solver\n");
  particle * pp = particles;
  granular_subproblem (pp, 0, dt, rhoval);
}
