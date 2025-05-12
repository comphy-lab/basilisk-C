/**

# AMR algorithm

The objective is to perform AMR simulations using $L^p$-norm metric-based error estimates.
When dealing with numerical solutions obtained using a solver, at least norm 2 is recommanded.

This file provide some global variables which can be modified:

* user_norm: $L^p$ norm. By default equal to 2.

* maxlevel: the maximum level of refinement. By default equal to 100, meaning unconstrained mesh adaptation.

* AMReps: the epsilon criterion used to refine/coarsen the elements.


First, we include the file which computes metric-based errors.
*/

#include "error_metric.h"

scalar AMRerror[];     // note: only one error field, but it may contain the
double AMReps = 1.e-5; // error from several sources
                       // fixme: to detail better, and cf 
                       // compute_metric_error_isotropic

int maxlevel = 100;

int user_norm = 2;
double Cadapt=1.;  // if necessary, set Cadapt<=1

void normed_restriction(Point point, scalar s) {

   double sum=0.;
   foreach_child()
      sum += pow(s[],user_norm);
   s[] = pow(sum,1./user_norm);

}

void normed_prolongation(Point point, scalar s) {    //  parent_cell_prolongation would be a better name

   double sum=0.;
   foreach_child()
      sum += pow(s[],user_norm);
   foreach_child()
      s[] = pow(sum,1./user_norm);   
   
}


/** The goal of the adaptation is to obtain a mesh containing $N$ elements. 
To do so, an epsilon criterion is used to find the cells which should be refined/coarsen.
This criterion must be iteratively adjusted in order to obtain the desired number of elements.

We want to find a good $eps$ criterion to use in the adaptation method.

The number of cells we an refine is the current number of
element minus the maximum number of element obtainable 
with a given max level: $N_{cellref} = N_{cells} - N_{cellmax}$.

The total error is supposed to be at order 2: $E = C_0\ h^2$. 
We want an equirepartition of the error on the elements, that means we ideally want that $eps = E/N$.

<span style="color:red"> TO DO: clearly explain all this section, and link to the article when it is published</span>

  $eps * N_{cellref} = C_0 dx^2 = C (Refvol/N_{cellref})^{2/dim}$
  
  $eps = C (Refvol)^{2/dim} / N_{cellref}^{2/dim + 1}$
  
  $log(eps) = log (C (Refvol)^{2/dim}) - {2/dim + 1} log(N_{cells} - N_{cellsmax})$
  
  $dlog(eps)/dN_{cells} = - {2/dim + 1}/(N_{cells} - N_{cellsmax})$
  
  $log(eps_{new}) = log(eps_{old}) + d(log(eps))/dN (N_{obj} - N_{cells})$
  
  $eps_{new} = eps_{old}* exp( - {2/dim + 1}/(N_{cells} - N_{cellsmax}) (N_{obj} - N_{cells})$
*/

double update_epsilon_control (int Nobj) {

  int Nmaxlevel = 0;
  if (maxlevel < 100) {
    foreach (reduction(+:Nmaxlevel)) { 
      if (level == maxlevel) {
        Nmaxlevel++;
      }  
    }
    Nmaxlevel = min(Nmaxlevel, 0.9*grid->tn);
  }

  if (AMReps == 0) {
    foreach (reduction(+:AMReps)) {
      AMReps+=AMRerror[];
    }
    AMReps /= Nobj;
  }


  return AMReps*exp(-Cadapt/user_norm*(Nobj - grid->tn)/(grid->tn- Nmaxlevel));

}


/**
The adaptation loop presents slight modification of the adapt_wavelet() Basilisk function and has the same objective: it refine the element for which the error is greater than the epsilon criterion, coarsen the elements which must be coarsen, and respects a maximum level presecription (which is not necessarily user-imposed).
It also enforce the 2:1 (1-irregularity) element size constrain.
*/

trace
astats adapt_metric_generic ( int maxlevel ) // fixme: either remove int maxlevel or add , double AMReps, it is weird to call for one of the two variables but not the other...
{

   restriction({AMRerror});
   
   astats st = {0, 0};
   scalar * listc = NULL;
   for (scalar s in all)
      if (!is_constant(s) && s.restriction != no_restriction)
	 listc = list_add (listc, s);

   // refinement
   tree->refined.n = 0;
   static const int refined = 1 << user, too_fine = 1 << (user + 1);
   foreach_cell() {
      if (is_active(cell)) {
	 static const int too_coarse = 1 << (user + 2);
	 if (is_leaf (cell)) {
	    if (cell.flags & too_coarse) {
	       cell.flags &= ~too_coarse;
	       refine_cell (point, listc, refined, &tree->refined);
	       st.nf++;
	    }
	    continue;
	 }
	 else { // !is_leaf (cell)
	    if (cell.flags & refined) {
	       // cell has already been refined, skip its children
	       cell.flags &= ~too_coarse;
	       continue;
	    }
	    // check whether the cell or any of its children is local
	    bool local = is_local(cell);
	    if (!local) {
	       foreach_child() {
		  if (is_local(cell)) {
		     local = true;
                     break;
                  }
               }
            }
	    if (local) {
	       static const int just_fine = 1 << (user + 3);
		  foreach_child() {
		     if (AMRerror[] > AMReps && level < maxlevel) {
			cell.flags &= ~too_fine;
			cell.flags |= too_coarse;
		     }
		     else if ((AMRerror[] <= AMReps/4. || level > maxlevel) &&   
			      !(cell.flags & (too_coarse|just_fine))) {
			   cell.flags |= too_fine;
		     }
		     else if (!(cell.flags & too_coarse)) {
			cell.flags &= ~too_fine;
			cell.flags |= just_fine;
		     }
		  }
	       foreach_child() {
		  cell.flags &= ~just_fine;
		  if (!is_leaf(cell)) {
		     cell.flags &= ~too_coarse;
		     if (level >= maxlevel)
			cell.flags |= too_fine;
		  }
		  else if (!is_active(cell))
		     cell.flags &= ~too_coarse;
	       }
	    }
	 }
      }
      else // inactive cell
	 continue;
   }
   mpi_boundary_refine (listc);
  
   // coarsening
   // the loop below is only necessary to ensure symmetry of 2:1 constraint
   for (int l = depth(); l >= 0; l--) {
      foreach_cell()
	 if (!is_boundary(cell)) {
	    if (level == l) {
	       if (!is_leaf(cell)) {
		  if (cell.flags & refined)
		     // cell was refined previously, unset the flag
		     cell.flags &= ~(refined|too_fine);
		  else if (cell.flags & too_fine) {
		     if (is_local(cell) && coarsen_cell (point, listc))
			st.nc++;
		     cell.flags &= ~too_fine; // do not coarsen parent
		  }
	       }
	       if (cell.flags & too_fine)
		  cell.flags &= ~too_fine;
	       else if (level > 0 && (aparent(0).flags & too_fine))
		  aparent(0).flags &= ~too_fine;
	       continue;
	    }
	    else if (is_leaf(cell))
	       continue;
	 }
      mpi_boundary_coarsen (l, too_fine);
   }
   free (listc);

   mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
   mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
   if (st.nc || st.nf)
      mpi_boundary_update (all);


   return st;
}



/**
A user interface performs the adaptation. It is similar (but not identical to) the adapt wavelet interface:
the user simply need to replace 
*/
/*
int maxl = 42;
double eps=0.01;
adapt_wavelet({psi}, (double[]){eps}, maxlevel=maxl);
*/
/**
by
*/
/*
      // restriction on the minimum grid size based on compression ratio (optional) 
      double etaopt = estimate_eta_opt(2, {psi});
      maxlevel = 0.5*log(Nobj/etaopt)/log(2.);   // maxlevel: global variable. Optional.

      // epsilon criteria for cell refinement/coarsening (mandatory) 
      AMReps = 0.01;                             // AMReps: global variable

      // AMR 
      adapt_metric( {psi} );                     // user interface similar to adapt_wavelet()
*/

trace
astats adapt_metric ( struct Adapt p )
{

  AMRerror.prolongation = normed_prolongation;
  AMRerror.restriction = normed_restriction;

  compute_metric_error_isotropic (user_norm, p, AMRerror);

  astats st = adapt_metric_generic (maxlevel);

  return st;

}


/** 
Global interpolation error estimation:
compute the global interpolation error assuming the local interpolation error is known and stored in the scalar field AMRerror[].
This global interpolation error is not used in the mesh adaptation procedure, it is mainly a post-treatment tool.

FOR LATER: I think this is not the right place for this function, and I'm even not sure it is very useful. For now, it is defined here for practical purposes (the scalar field AMRerror[] is defined in this file).

 */
double Interpolation_error (int norm) {

  double err = 0.;
  foreach (reduction(+:err))
    err += pow(AMRerror[],norm);

  return pow(err, 1./norm);

}



