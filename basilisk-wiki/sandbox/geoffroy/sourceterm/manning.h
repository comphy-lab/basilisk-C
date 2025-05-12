/**
# Friction term : Manning in saint venant

When the Reynolds number is high (>2000), the stream becomes turbulent. In this case, the friction term can't be easily solved  analytically. Manning proposed an empirical law describing this term, it can be written in its full form as : 
$$ 
   Cf = - n^2 g \frac{q|q|}{h^{7/3}} 
$$
where n is a free-parameter wich depends on the nature of the soil.

The overloading process is fully explained in  [poiseuille.h](poiseuille.h)*/

// Manning coefficient
double n = 0.025;

/**
We define the function which will replace the update function in the
predictor-corrector
*/

void updatemanning(scalar * evolving, scalar * sources, double dtmax, int numbersource ){
  // We first recover the evolving fields
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };
  
  // Updates for evolving quantities
  vector dshu = { sources[1], sources[2] }; 

  foreach(){
    if(h[] > dry){
      
      /** 
      We Compute the new field u with an implicit scheme. The
      $u^2$ term is linearised */

      double s = dtmax*n*n*G*norm(u)/pow(h[],4./3.);
      foreach_dimension()
	// Translate it in an explicit form
	dshu.x[] -= h[]*u.x[]*s/(s+1)/dtmax;
    }
  }
  // Calling of the next source term
  numbersource++;
  updatesource[numbersource](evolving,sources,dtmax,numbersource);
}

// Overloading 
event initmann(i = 0){
  updatesource[numbersource]=updatemanning;
  numbersource++;
  updatesource[numbersource] = fnull;
}
