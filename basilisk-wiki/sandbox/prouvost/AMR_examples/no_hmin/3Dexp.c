/**
This case is 3D.

We search an adapted mesh containing $N_{obj}$ elements which minimizes the total error (error between the numerical solution of the Basilisk Poisson-Helmholtz solver) and the (known) analytical solution.

The [solution](./exp_3D.h) is a 3D an exponential function.

 */

#include "grid/octree.h"
#include "poisson.h" // solver

#include "./exp_3D.h"  // ''problem''
#include "./../../AMR_tools/amr.h"  // AMR
#include "./../../utils/gauss_quadrature.h"  // compute total error

#include "utils.h"

double TOL = 1.e-6;

const face vector alp[] = {D,D,D};
const scalar lam[] =s;

scalar psi[]; // the numerical solution
psi[left] = dirichlet(exact(x,y,z));
psi[right] = dirichlet(exact(x,y,z));
psi[top] = dirichlet(exact(x,y,z));
psi[bottom] = dirichlet(exact(x,y,z));
psi[front] = dirichlet(exact(x,y,z));
psi[back] = dirichlet(exact(x,y,z));

scalar psi_exact[];  // the analytical solution
psi_exact[left] = dirichlet(exact(x,y,z));
psi_exact[right] = dirichlet(exact(x,y,z));
psi_exact[top] = dirichlet(exact(x,y,z));
psi_exact[bottom] = dirichlet(exact(x,y,z));
psi_exact[front] = dirichlet(exact(x,y,z));
psi_exact[back] = dirichlet(exact(x,y,z));


int main() {

  FILE * fpglobal = fopen("error","w");

  L0=1.;
  origin(-0.5,-0.5);
  
  int mylev=5;
  int Nobj = pow(2,mylev*3);  // initial objective number of element

  init_grid (1 << (mylev-1));

  scalar rhs[]; // RHS term for Poisson-Helmholtz equation

/**
We do a loop to obtain several adapted meshes having $N_{obj}$ elements.

For each $N_{obj}$, we do a loop to obtain the objective number of elements.

*/

  for (int j=0; Nobj <= pow(2,3*6)/1.6; j++){

    for (int ki=0; ( fabs((double)(grid->tn - Nobj)/Nobj)) > 0.03 ; ki++){

      foreach ()
        rhs[] = src(x,y,z);

      poisson (psi, rhs, alp, lam, tolerance = TOL);   // poisson solver 

/**
We update an epsilon criterion to obtain the objective number of elements
*/
      
      // AMR criterion 
      if (ki < 5 && j == 0) {  // estimate AMReps the first run with uniform refinement
        struct PreFactorData cd = compute_prefactors (2, {psi} );
        AMReps =  cd.cuniform*pow(Nobj,-2./3.)/pow(Nobj,1./2.); 
      } else {
        AMReps = update_epsilon_control(Nobj);
      }

/**
We adapt the mesh.
*/
      
      adapt_metric( {psi} );
    }


/**
We compute total error and interpolation error.

*/

    foreach ()
      rhs[] = src(x,y,z);
    poisson (psi, rhs, alp, lam, tolerance = TOL);

    // TOTAL ERROR
    double errtot = norm_gauss_5p (psi, user_norm, exact);  // compute total error :   ||u_num - u_exact||

    // INTERPOLATION ERROR estimate with metric-based formulae
    double Interr = Interpolation_error(user_norm);

    // INTERPOLATION ERROR exact
    foreach ()
      psi_exact[] = exact(x,y,z);

    double Int_err_exact = norm_gauss_5p (psi_exact, user_norm, exact);   // compute exact interpolation error :   ||u_interp - u_exact||

    // theoretical optimal and uniform errors
    struct PreFactorData cd;
    cd = compute_prefactors (2, {psi}, NULL );

    fprintf(fpglobal,"%ld %g %.10g %.10g %.10g %g\n", grid->tn, errtot, Int_err_exact, Interr, cd.copt*pow(grid->tn,-2./3.), cd.cuniform*pow(grid->tn,-2./3.));


/**

We plot the total error and the interpolation error.
We see that the measured total error is far from the expected optimal, whereas the interpolation error is close to the optimal.
Thus, we show that, in this case, the obtained meshes minimizes the interpolation error, but not the total error.

~~~gnuplot Error
reset
set term pngcairo enhanced size 500,500 
set output 'error.png'

set logscale
set xtics (16,32,64,128,256)
set format y "10^{%T}"

p "error" u (($1)**(1./3.)):6 w l t "uniform",\
"error" u (($1)**(1./3.)):5 w l t "optimal",\
"error" u (($1)**(1./3.)):2 w p t "total error",\
"error" u (($1)**(1./3.)):3 w p t "interp error"
~~~


*/


    
/**
We verify that we obtained the expected result
*/

    assert ( fabs((double)(grid->tn - Nobj)/Nobj) <= 0.03 && errtot>5.*cd.copt*pow(grid->tn,-(2./3.)) );

    
/**
loop update

*/
    
    Nobj *= 1.5;

  }

  free_grid();
  fclose(fpglobal);
 
}



/**

 */










