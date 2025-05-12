/**
This case is 3D.

We search an adapted mesh containing $N_{obj}$ elements which minimizes the total error (error between the numerical solution of the Basilisk Poisson-Helmholtz solver) and the (known) analytical solution.

The [solution](./../metric_amr_no_constraint/exp_3D.h) is a 3D an exponential function.

 */

#include "grid/octree.h"
#include "poisson.h" // solver

#include "./../no_hmin/exp_3D.h"  // ''problem''
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
We adapt the mesh with an additionnal constraint on the minimal size computed from an estimation of the optimal conmpression ratio.
*/
      /* restriction on the minimum grid size */
      double etaopt = estimate_eta_opt(2, {psi});
      maxlevel = floor(1./dimension*log(Nobj/etaopt)/log(2.));
      
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

We plot the total error and the interpolation error and we compare with the case without constraint.
We plot the error in function of the cubic root of the number of element (equivalent to mean cell size in 3D).
We see that the total error is lower with the constraint than without constraint, but the interpolation error is slightly higher with the constraint.

~~~gnuplot Error
reset
set term pngcairo enhanced size 700,400 
set output 'error.png'

set logscale
set xtics (16,32,64,128,256)
set format y "10^{%T}"

set key below

set multiplot layout 1,2

set title "total error"
p "error" u (($1)**(1./3.)):6 w l t "uniform",\
"error" u (($1)**(1./3.)):5 w l t "optimal",\
"error" u (($1)**(1./3.)):2 w p t "total error with constraint",\
"./../../no_hmin/3Dexp/error" u (($1)**(1./3.)):2 w p t "total error no constraint"

set title  "interpolation error"
p "error" u (($1)**(1./3.)):6 w l t "uniform",\
"error" u (($1)**(1./3.)):5 w l t "optimal",\
"error" u (($1)**(1./3.)):3 w p t "interp error with constraint",\
"./../../no_hmin/3Dexp/error" u (($1)**(1./3.)):3 w p t "interp error no constraint"
unset multiplot
~~~


*/


    
/**
We verify that we obtained the expected result
*/

    assert ( fabs((double)(grid->tn - Nobj)/Nobj) <= 0.03);
    if (Nobj>100000)
      assert( errtot<2.*cd.copt*pow(grid->tn,-(2./3.)) );

    
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










