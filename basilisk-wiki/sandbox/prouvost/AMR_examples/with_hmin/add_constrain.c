#include "poisson.h"   // solver

#include "./../no_hmin/exponentiel.h"  // ''problem''
#include "./../../AMR_tools/amr.h"   // AMR
#include "./../../utils/gauss_quadrature.h"  // compute total error

/**
We do a small mesh adaptation experiment in which we start from a mesh minimizing the interpolation error.
Then, we impose a constrain on the minimal element size and we iterate until the final mesh has the same number of element than the first mesh.

We will see how the total error (difference between the continuous solution and the discrete solution) evolves.
*/


double TOL = 1.e-7;  // poisson solver tolerance

const face vector alp[] = {D,D};
const scalar lam[] =s;

scalar psi[];   // numerical solution
psi[left] = dirichlet(exact(x,y,0.));
psi[right] = dirichlet(exact(x,y,0.));
psi[top] = dirichlet(exact(x,y,0.));
psi[bottom] = dirichlet(exact(x,y,0.));

scalar psi_exact[];    // point-wise exact solution
psi_exact[left] = dirichlet(exact(x,y,0.));
psi_exact[right] = dirichlet(exact(x,y,0.));
psi_exact[top] = dirichlet(exact(x,y,0.));
psi_exact[bottom] = dirichlet(exact(x,y,0.));


int main() {

  int mylev=10;

  init_grid (1 << mylev);

  scalar rhs[];


/**
We obtain a mesh minimizing the interpolation error in $L^2$-norm.
*/
  int norm=2;
  
  for (int j=0; j<=10; j++){

    foreach ()
      rhs[] = src(x,y,0);
    poisson (psi, rhs, alp, lam, tolerance = TOL);
    
    AMReps = 0.6e-4;
    adapt_metric( {psi} );
  }


/**
We register the total and interpolation error on this mesh
*/

  int it=0;  // compteur

  foreach ()
    rhs[] = src(x,y,0);
  poisson (psi, rhs, alp, lam, tolerance = TOL); // compute numerical solution

  foreach ()
    psi_exact[] = exact(x,y,0);   // compute point-wise theorical solution

  double errtot = norm_gauss_5p (psi, norm, exact);  // compute total error :   ||u_num - u_exact||

  double Int_err_exact = norm_gauss_5p (psi_exact, norm, exact);   // compute exact interpolation error :   ||u_interp - u_exact||

  fprintf(stderr,"%ld %.10g %.10g %i\n", grid->tn, errtot, Int_err_exact, it);

/**
we also register the local total error and the level of refinement
*/  
  
  scalar err_total_local[];  

  struct RealError rer;
  rer.Lnorm=2;
  rer.f=psi;
  norm_gauss_5p_local (rer,err_total_local);

  scalar lev[];
  foreach()
    lev[]=level;

  FILE * fpvtk = fopen("fields_0","w");
  foreach()
    fprintf(fpvtk, "%g %g %g %g\n", x, y, lev[], err_total_local[]);
  fclose(fpvtk);


/**
Then, we add a constrain on the minimal element size based on the optimal compression ratio, and we do several adaptations to obtain a constrained mesh containing almost the same number of element.
*/

  int Nobj = 0;   // number of element to reach
  foreach()
    Nobj++;

  int testi=1;  // we want to be sure to do the first iteration
  for (it=1; (testi || (fabs((double)(grid->tn - Nobj)/Nobj)) > 0.03) ; it++){ // we loop until the objective number of element is reached
    testi=0;

    // AMR criterion 
    double etaopt = estimate_eta_opt(2, {psi});
    maxlevel = 0.5*log(Nobj/etaopt)/log(2.);   // we impose the maxlevel contrain
    
    AMReps = update_epsilon_control(Nobj);   // we update epsilon to obtain Nobj
    adapt_metric( {psi} );

/**
We register the total/interpolation error at each iteration    
*/
    
    foreach ()
      rhs[] = src(x,y,0);
    poisson (psi, rhs, alp, lam, tolerance = TOL); // compute numerical solution

    foreach ()
      psi_exact[] = exact(x,y,0);   // compute point-wise theorical solution

    double errtot = norm_gauss_5p (psi, norm, exact);  // compute total error :   ||u_num - u_exact||

    double Int_err_exact = norm_gauss_5p (psi_exact, norm, exact);   // compute exact interpolation error :   ||u_interp - u_exact||

    fprintf(stderr,"%ld %.10g %.10g %i\n", grid->tn, errtot, Int_err_exact, it);
    
  }

/**
At the end, we register the local fields
*/  

  norm_gauss_5p_local (rer,err_total_local);

  foreach()
    lev[]=level;

  fpvtk = fopen("fields_1","w");
  foreach()
    fprintf(fpvtk, "%g %g %g %g\n", x, y, lev[], err_total_local[]);
  fclose(fpvtk);
    
  free_grid();
}






/**

We see that the final constrain mesh has a total error reduced by more than a decade in comparison with the unconsrtain (initial) mesh.
This is done at the cost of a slight increase of the interpolation error.

~~~gnuplot Evolution of the total and interpolation errors when adding a correct minimal cell size constrain
set term pngcairo enhanced size 400,400 
set output 'error.png'

set logscale y
set xlabel "iteration"
set ylabel "L^2 error"
set format y "10^{%T}"

p 'log' u 4:3 w lp t 'interpolation error', 'log' u 4:2 w lp t 'total error'
~~~


We also plot the mesh levels and the total local error to see that, on the non-constrained mesh, the max of this error is localised in region where the element size changes, whereas is localized in the region of maximal resolution in the constrained mesh.

~~~gnuplot Left: mesh level.  Right: local total error. Top: non-constrained mesh. Bottom: constrained mesh.
reset
set term pngcairo enhanced size 800,800
set output 'fields.png'

unset xtics
unset ytics

set multiplot layout 2,2
set view map

sp 'fields_0' u 1:2:3 w p pt 5 ps 3 palette not
#set logscale cb
sp 'fields_0' u 1:2:4 w p pt 5 ps 3 palette not
unset logscale
sp 'fields_1' u 1:2:3 w p pt 5 ps 3 palette not
#set logscale cb
sp 'fields_1' u 1:2:4 w p pt 5 ps 3 palette not
~~~




 */









