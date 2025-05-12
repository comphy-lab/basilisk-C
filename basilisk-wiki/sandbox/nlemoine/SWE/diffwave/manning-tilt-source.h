// Manning coefficient
//double n = 0.025;  
scalar nmanning[]; // n variable spatialement

// Tilt
coord tilt = {.x = 0., .y = 0.};

void updatemanningtilt(scalar * evolving, scalar * sources, double dtmax, int numbersource ){
  // We first recover the evolving fields
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };
//  double s, h2inv;
//  double eps = 

  // Updates for evolving quantities
  vector dshu = { sources[1], sources[2] }; 

  foreach(){
    if(h[] > dry){
//      double s = dtmax*n*n*G*norm(u)/pow(h[],4./3.);	// sans dimension
//    FIXME : désingulariser h pour obtenir un ressuyage plus rapide des zones faiblement submergées ?
      double s = dtmax*nmanning[]*nmanning[]*G*norm(u)/pow(h[],4./3.);	// Equation (5) du mémo
//      h2inv = 
      
      foreach_dimension()
	// Translate it in an explicit form
//	dshu.x[] += G*h[]*tilt/(1+s) - h[]*u.x[]*s/(1+s)/dtmax;
	dshu.x[] += G*h[]*tilt.x/(1+s) - h[]*u.x[]*s/(1+s)/dtmax; // Equation (6) du mémo
    }
  }
  boundary((scalar *){dshu});

  // Calling of the next source term
  numbersource++;
  updatesource[numbersource](evolving,sources,dtmax,numbersource);
}

// Overloading 
event initmanningtilt(i = 0){
  updatesource[numbersource]=updatemanningtilt;
  numbersource++;
  updatesource[numbersource] = fnull;
}