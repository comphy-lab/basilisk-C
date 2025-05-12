/**
# Friction term : Darcy in saint venant

When the Reynolds number is high (>2000), the stream becomes turbulent. In this case, the friction term can't be easily solved  analytically. Darcy and Weisbach found an empirical law describing this term, it can be written in its full form as : 
$$ 
   Cf = \frac{f}{8} \frac{q|q|}{h^2}, 
$$
where f is a free-parameter wich depends on the nature of the soil.

The overloading process is fully explained in  [poiseuille.h](poiseuille.h)

*/

// Darcy coefficient
double f = 0.5; 

/**
We define the function which will replace the update function in the
predictor-corrector
*/

void updatedarcy(scalar * evolving, scalar * sources, double dtmax, int numbersource ){

  // We first recover the evolving fields
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };
  
  // Updates for evolving quantities
  vector dshu = { sources[1], sources[2] }; 


  foreach(){
    if(h[] > dry){
      /** Computing the new field u with an implicit scheme. The term $u^2$ is linearised.
       */
      double s = dtmax*norm(u)*f/(8*h[]);
      foreach_dimension()
	dshu.x[] -= h[]*u.x[]*s/(1+s)/dtmax;
      }
    }
  }
  // Calling of the next source term
  numbersource++;
  updatesource[numbersource](evolving,sources,dtmax,numbersource);

}

// Overloading
event initdarcy (i = 0){
  updatesource[numbersource]=updatedarcy;
  numbersource++;
  updatesource[numbersource] = fnull;
}
