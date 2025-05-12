/**
# Friction term : Poiseuille in saint venant

## Poiseuille flow

For a low Reynolds number, the stream is laminar. In this case, we can
compute an exact theoretical friction term in the momentum equation of
saint-venant, supposing a Poiseuille stream. We call this term the
"Poiseuille friction term". It can be written as : 

$$ 
Cf = 3 \nu \frac{q}{h^2} 
$$ 

where $\nu$ is the dynamic viscosity of the fluid.

*/
// Dynamic viscosity of the water
double nu = 1e-6; 

/**
We define the function which will replace the update function in the
predictor-corrector

 */
void updatepoiseuille(scalar * evolving, scalar * sources, double dtmax, int numbersource ){

  // We first recover the evolving fields
  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };
  
  // Updates for evolving quantities
  vector dshu = { sources[1], sources[2] }; 
/**
## Implicit treatment with explicit scheme

We chose to compute all our source term with an explicit term for
conveniance. If $\Psi (U)$ is the function used to solved the $x$
differential term and $S_i(U)$ the $i^{th}$ source term in the
right-hand side. Then, our scheme can be written as :

$$
\begin{array}{c} U^* =U^n+ \frac{\Delta t}{\Delta x} \Psi (U^n)\\
U^{n+1}= U^* +\Delta t\Sigma_{i=1}^N S_i(U^*) \end{array} 
$$

But, even with an explicit scheme, we can compute the friction term in
an implicit way. To do this, we first compute implicitly an
intermediary field . We then compute the $S_i$ explicit term using the
intermediary field as follows :
*/

  foreach(){
    if(h[] > dry){
      // Compute the new field u with an implicit scheme
      double s = dtmax*3.*nu/(h[]*h[]);
      foreach_dimension()
	// Translate it in an explicit form
	dshu.x[] -= h[]*u.x[]*s/(s+1)/dtmax;
    }
  }
  // We finish by calling the next source term
  numbersource++;
  updatesource[numbersource](evolving,sources,dtmax,numbersource);
}

/**
## Overloading 

In an initial event, we simply overload the $i^{th}$ updatesource[]
fonction with the function defined before. We also overload the
$i+1^{th}$ function with the fnull function.*/

event initpoiseuille (i = 0){
  updatesource[numbersource] = updatepoiseuille;
  numbersource++;
  updatesource[numbersource] = fnull;
}
