/**
# Rain term in saint venant

*/

// Intensity of the rainfall
scalar rain[];
// Its horizontal velocity
vector urain[]; // Not yet implemented

void updaterain(scalar * evolving, scalar * sources, double dtmax, int numbersource ){
  // Updates for evolving quantities
  scalar dsh = sources[0];
 
  // Computing the source term
  foreach(){
    dsh[] += rain[];
  }
  
  // Calling the next source term
  numbersource++;
  updatesource[numbersource](evolving,sources,dtmax,numbersource);
}

// Initialisation
event initrain(i = 0){
  updatesource[numbersource] = updaterain;
  numbersource++;
  updatesource[numbersource] = fnull ;
}
