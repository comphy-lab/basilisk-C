/* The general DLMFD plugin */ 
# include "dlmfd-plugin.h"

/* Here we can overload the generic events defined in the general
   DLMFD plugin dlmfd-plugin.h such that it uses a granular solver */


/** Overloading of the granular solver init event */
// -------------------------------------------------
event GranularSolver_init (t < -1.)
{

} 

/** Overloading of the granular solver predictor event */
// ------------------------------------------------------
event GranularSolver_predictor (t < -1.)
{
}

/** Overloading of the granular solver velocity update event */
// ------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.)
{
}
