#include <fftw3.h>
#include <complex.h>
 
#define M_PI acos(-1.0)

int get_LHS_blocks( int nsink, int nbound, int ncoef, int ntheta, double ** _THETA, 
                    double ** _ReZstart, double ** _ImZstart, double ** _ReZend, double ** _ImZend,
                    double ** _UL, double ** _UC, double ** _UR, double ** _LL, double ** _LC, double ** _LR,
                    double ** _ReZ, double ** _ImZ )
{
        int nslit = nsink + nbound;
        int nrows_J = ntheta*nslit; 

	// DYNAMIC ALLOCATIONS FOR FORTRAN SUBROUTINE

	double * ReCn = (double *) malloc(nslit*ncoef*sizeof(double));
	double * ImCn = (double *) malloc(nslit*ncoef*sizeof(double));
	double * Q    = (double *) malloc(nslit*sizeof(double));
	(*_ReZ) = (double *) malloc(nrows_J*sizeof(double));		
	(*_ImZ) = (double *) malloc(nrows_J*sizeof(double));		

	double * PHI  = (double *) malloc(nrows_J*sizeof(double));		
	double * PSI  = (double *) malloc(nrows_J*sizeof(double));
 	double * Vx   = NULL, * Vy   = NULL;

        int ndof = (2*ncoef+1)*nslit;
	// Enable computation of Jacobian for potential and velocity
        int varout[3] = {0,1,1};
   	// Allocate arrays for Jacobian of potential	
	double * J_PHI = (double *) malloc(nrows_J*ndof*sizeof(double));
	double * J_PSI = (double *) malloc(nrows_J*ndof*sizeof(double));
   	// Allocate arrays for Jacobian of velocity	
	double * J_Vx = (double *) malloc(nrows_J*ndof*sizeof(double));
	double * J_Vy = (double *) malloc(nrows_J*ndof*sizeof(double));
   	// Allocate work arrays	
	double * J_PHI1 = (double *) malloc(ndof*sizeof(double));
	double * J_PSI1 = (double *) malloc(ndof*sizeof(double));
	double * J_Vx1 = (double *) malloc(ndof*sizeof(double));
	double * J_Vy1 = (double *) malloc(ndof*sizeof(double));
	double * Xdof = (double *) malloc(ndof*sizeof(double));

        double v0x = 0., v0y = 0., Phi0 = 0.;

	(void) conformsystem_ ( (*_THETA),&ntheta,
				(*_ReZstart),(*_ImZstart),(*_ReZend),(*_ImZend),&nslit,   
				ReCn,ImCn,&ncoef,Q,&v0x,&v0y,&Phi0,
				varout,(*_ReZ),(*_ImZ),
				PHI,PSI,Vx,Vy,J_PHI,J_PSI,J_Vx,J_Vy,
				J_PHI1,J_PSI1,J_Vx1,J_Vy1,
				Xdof,&ndof);
        free(PHI);free(PSI);
        free(ReCn);free(ImCn);
        free(Q);
        free(J_PHI1);free(J_PSI1);
        free(J_Vx1);free(J_Vy1);
        free(Xdof);
        free(J_PSI);

        int row,col,row_J,col_J,rec;

        // UPPER-LEFT BLOCK [ nsink*ntheta , nsink ]:
        // Jacobian of Phi at stream control points, w.r.t. sink terms
        (*_UL) = (double *) malloc(nsink*ntheta*nsink*sizeof(double));

        for(row=0;row<(nsink*ntheta);row++){
          row_J = row;
          for(int sink=0;sink<nsink;sink++){
            col = sink;
            col_J = (2*ncoef+1)*sink + 2*ncoef;

            (*_UL)[row+(nsink*ntheta)*col] = J_PHI[row_J+nrows_J*col_J]; 
        }}

        fprintf(stderr,"Done computing UPPER-LEFT block.\n");
	fflush(stderr);

        // UPPER-CENTER BLOCK [ nsink*ntheta , nsink*ncoef ]:
        // Jacobian of Phi at stream control points, w.r.t. Re{C_n} coefficients at stream slits
        // Rq: the Im{C_n}'s are zero for ensuring continuity of real potential
        (*_UC) = (double *) malloc(nsink*ntheta*(nsink*ncoef)*sizeof(double));

        for(row=0;row<(nsink*ntheta);row++){
          row_J = row;
          for(int sink=0;sink<nsink;sink++){
          for(int coef=0;coef<ncoef;coef++){
            col = ncoef*sink + coef;
            col_J = (2*ncoef+1)*sink + 2*coef;

            (*_UC)[row+(nsink*ntheta)*col] = J_PHI[row_J+nrows_J*col_J]; 
          }}
	}

        fprintf(stderr,"Done computing UPPER-CENTER block.\n");
	fflush(stderr);

        // UPPER-RIGHT BLOCK [ nsink*ntheta , 2*nbound*ncoef+3 ]:
        // Jacobian of Phi at stream control points, w.r.t. C_n coefficients at boundary slits, and v0x / v0y / Phi0
        (*_UR) = (double *) malloc(nsink*ntheta*(2*nbound*ncoef+3)*sizeof(double));

        for(row=0;row<(nsink*ntheta);row++){
          row_J = row;
          for(int bound=0;bound<nbound;bound++){
          for(int coef=0;coef<(2*ncoef);coef++){

            col_J = (2*ncoef+1)*(nsink+bound) + coef;
            col = (2*ncoef)*bound + coef;

            (*_UR)[row+(nsink*ntheta)*col] = J_PHI[row_J+nrows_J*col_J];

          }}

          col = (2*nbound*ncoef);
          (*_UR)[row+(nsink*ntheta)*col] = +1.0; // w.r.t. Phi0

          col++;
          (*_UR)[row+(nsink*ntheta)*col] = -1.0*((*_ReZ)[row_J]-X0); // w.r.t. v0x

          col++;
          (*_UR)[row+(nsink*ntheta)*col] = -1.0*((*_ImZ)[row_J]-Y0);  // w.r.t. v0y
        }

        fprintf(stderr,"Done computing UPPER-RIGHT block.\n");
	fflush(stderr);

        // LOWER-LEFT BLOCK [ nbound*ntheta , nsink ]:
        // Jacobian of vn (normal velocity component) at boundary control points, w.r.t. sink terms
        (*_LL) = (double *) malloc(nbound*ntheta*nsink*sizeof(double));

        for(row=0;row<(nbound*ntheta);row++){
          row_J = (nsink*ntheta)+row;
          // Compute inward normal to boundary at control point
          int slit = nsink + (row / ntheta);
          double complex vv = I*((*_ReZend)[slit]-(*_ReZstart)[slit] + I*(*_ImZend)[slit]-I*(*_ImZstart)[slit]);
          vv /= cabs(vv);
          double nx = creal(vv);
          double ny = cimag(vv);
          for(int sink=0;sink<nsink;sink++){
            col = sink;
            col_J = (2*ncoef+1)*sink + 2*ncoef;
            (*_LL)[row+(nbound*ntheta)*col] = nx*J_Vx[row_J+nrows_J*col_J]
                                            + ny*J_Vy[row_J+nrows_J*col_J]; 
        }}

        fprintf(stderr,"Done computing LOWER-LEFT block.\n");
	fflush(stderr);

        // LOWER-CENTER BLOCK [ nbound*ntheta , nsink*ncoef ]:
        // Jacobian of vn (normal velocity component) at boundary control points, w.r.t. Re{C_n} coefficients at stream slits
        // Rq: the Im{C_n}'s are zero for ensuring continuity of real potential
        (*_LC) = (double *) malloc(nbound*ntheta*(nsink*ncoef)*sizeof(double));

        for(row=0;row<(nbound*ntheta);row++){
          row_J = (nsink*ntheta)+row;
          // Compute inward normal to boundary at control point
          int slit = nsink + (row / ntheta);
          double complex vv = I*((*_ReZend)[slit]-(*_ReZstart)[slit] + I*(*_ImZend)[slit]-I*(*_ImZstart)[slit]);
          vv /= cabs(vv);
          double nx = creal(vv);
          double ny = cimag(vv);
          for(int sink=0;sink<nsink;sink++){
          for(int coef=0;coef<ncoef;coef++){

            col_J = (2*ncoef+1)*sink + 2*coef;
            col = ncoef*sink + coef;

            (*_LC)[row+(nbound*ntheta)*col] = nx*J_Vx[row_J+nrows_J*col_J]
                                            + ny*J_Vy[row_J+nrows_J*col_J];
          }}
        }

        fprintf(stderr,"Done computing LOWER-CENTER block.\n");
	fflush(stderr);

        // LOWER-RIGHT BLOCK [ nbound*ntheta , 2*nbound*ncoef+3 ]:
        // Jacobian of vn (normal velocity component) at boundary control points, w.r.t. C_n coefficients at boundary slits, and v0x / v0y / Phi0
        (*_LR) = (double *) malloc(nbound*ntheta*(2*nbound*ncoef+3)*sizeof(double));

        for(row=0;row<(nbound*ntheta);row++){

          row_J = (nsink*ntheta)+row;
          // Compute inward normal to boundary at control point
          int slit = nsink + (row / ntheta);
          double complex vv = I*((*_ReZend)[slit]-(*_ReZstart)[slit] + I*(*_ImZend)[slit]-I*(*_ImZstart)[slit]);
          vv /= cabs(vv);
          double nx = creal(vv);
          double ny = cimag(vv);

          for(int bound=0;bound<nbound;bound++){
          for(int coef=0;coef<(2*ncoef);coef++){

            col_J = (2*ncoef+1)*(nsink+bound) + coef;
            col = (2*ncoef)*bound + coef;

            (*_LR)[row+(nbound*ntheta)*col] = nx*J_Vx[row_J+nrows_J*col_J]
                                            + ny*J_Vy[row_J+nrows_J*col_J];

          }}

          col = (2*nbound*ncoef);
          (*_LR)[row+(nbound*ntheta)*col] = 0.; // w.r.t. Phi0

          col++;
          (*_LR)[row+(nbound*ntheta)*col] = nx; // w.r.t. v0x

          col++;
          (*_LR)[row+(nbound*ntheta)*col] = ny; // w.r.t. v0y
        }

        fprintf(stderr,"Done computing LOWER-RIGHT block.\n");
	fflush(stderr);

        free(J_Vx);
        free(J_Vy);
        free(J_PHI);

        return(0);
}

int prune_tree(	 double ** _ReZstart, double ** _ImZstart, double ** _Phistart, 
		 double ** _ReZend, double ** _ImZend, double ** _Phiend,
		 double ** _Aup, double ** _dA, 
                 int ** _DownID, int * _nsink,
                 double Aup_min)
{
  for(int s=0;s<(*_nsink);s++)
  {
    if( (*_Aup)[s] < Aup_min )  // merge with downstream reach subcatchment
    {     
	int down = (*_DownID)[s];
	(*_dA)[down] += (*_Aup)[s] ;
    }
  }

  double * ReZstart = NULL, * ImZstart = NULL, * Phistart = NULL,
         * ReZend = NULL, * ImZend = NULL, * Phiend = NULL, * dA = NULL;
  
// Discard slits draining less than Aup_min

  int nsink = 0;
  double sumdA = 0.;

  for(int s=0;s<(*_nsink);s++)
  {
     if( (*_Aup)[s] >= Aup_min )
     {
        ReZstart = (double *) realloc(ReZstart,(nsink+1)*sizeof(double));
        ReZstart[nsink] = (*_ReZstart)[s];
        ImZstart = (double *) realloc(ImZstart,(nsink+1)*sizeof(double));
        ImZstart[nsink] = (*_ImZstart)[s];
        Phistart = (double *) realloc(Phistart,(nsink+1)*sizeof(double));
        Phistart[nsink] = (*_Phistart)[s];

        ReZend = (double *) realloc(ReZend,(nsink+1)*sizeof(double));
        ReZend[nsink] = (*_ReZend)[s];
        ImZend = (double *) realloc(ImZend,(nsink+1)*sizeof(double));
        ImZend[nsink] = (*_ImZend)[s];
        Phiend = (double *) realloc(Phiend,(nsink+1)*sizeof(double));
        Phiend[nsink] = (*_Phiend)[s];

        dA = (double *) realloc(dA,(nsink+1)*sizeof(double));
        dA[nsink] = (*_dA)[s];
        
        sumdA += dA[nsink];
        nsink++;
     }
  }

  fprintf(stderr," Sum of reach subcatchments after pruning: %.6lfe6\n",sumdA/1.e6);  

  (*_nsink) = nsink;

  free(*_ReZstart);
  *_ReZstart = ReZstart;
  free(*_ImZstart);
  *_ImZstart = ImZstart;
  free(*_Phistart);
  *_Phistart = Phistart;

  free(*_ReZend);
  *_ReZend = ReZend;
  free(*_ImZend);
  *_ImZend = ImZend;
  free(*_Phiend);
  *_Phiend = Phiend;

  free(*_dA);
  *_dA = dA;

  return(0);
}

int solve_correcting_dipoles( 	int nsink, int nbound, int ncoef, int ntheta,
                   		double ** _UC, double ** _UR, double ** _LC, double ** _LR,
               			double ** _ReZ, double ** _ImZ,
// RHS ingredients:
                   		double ** _SlitLength, double ** _PHI_topo,
				double ** _PHI_areal1, double ** _Vn_areal1,
				double ** _UL, double ** _LL, double gamma0, double ** _dA,
                   		double ** _Q, double ** _ReCn, double ** _ImCn, double * _v0x, double * _v0y, double * _Phi0,
                   		double ** _PHI )
{
   double * A, * B; // for DGELS
   char trans = 'N';
   int nrhs = 1;
   int LWORK, INFO;
   double * WORK;

   bool fit_trend =  false;

   /////////////////////////////////////////////////////////////////////////
   /////////////// Fill matrix of the system, and RHS vector ///////////////
   /////////////////////////////////////////////////////////////////////////
   
   int ncnt = ntheta*(nsink+nbound);
   int ndof = fit_trend ? (nsink*ncoef+2*nbound*ncoef+3) : (nsink*ncoef+2*nbound*ncoef+1); // γ0, Re{Cn}'s and Im{Cn}'s for boundary segments, Phi0, (v0x, v0y)

   A = (double *) malloc(ncnt*ndof*sizeof(double));
   B = (double *) malloc(ncnt*sizeof(double));

   fprintf(stderr,"Fill matrix of the system with %d degrees of freedom and %d control points... ",ndof,ncnt);
   fflush(stderr);

   int nrows_U = ntheta*nsink, nrows_L = ntheta*nbound;

   // Fill upper part of A and B (constraints on potential in the thalweg network)

   for(int r=0;r<nrows_U;r++)
   {
      // RHS:
      B[r] = (*_PHI_topo)[r] - fabs(gamma0) * (*_PHI_areal1)[r]; // RHS is residual topographic elevation

      for(int s=0;s<nsink;s++)
         B[r] -= fabs(gamma0) * ((*_UL)[r+nrows_U*s])*((*_dA)[s]);

      // LHS:
      for(int c=0;c<(nsink*ncoef);c++)
         A[r+ncnt*c] = (*_UC)[r+nrows_U*c];

      int coffset = nsink*ncoef , cmax = ndof-coffset;  

      for(int c=0;c<cmax;c++)
         A[r+ncnt*(coffset+c)] = (*_UR)[r+nrows_U*c];
   }

   // Fill lower part of A and B (constraints on normal velocity on the boundary)
   
   for(int r=0;r<nrows_L;r++)
   {
      int slit = nsink + (r / ntheta);
      double Lw = sqrt(((double)nrows_U)/((double)nrows_L)) * (*_SlitLength)[slit];

      // RHS:
      B[nrows_U+r] = 0. - fabs(gamma0) * (*_Vn_areal1)[r]; // RHS is residual normal velocity component
 
      for(int s=0;s<nsink;s++)
        B[nrows_U+r] -= fabs(gamma0) * ((*_LL)[r+nrows_L*s])*((*_dA)[s]); 

      B[nrows_U+r] *= Lw;

      // LHS:
      for(int c=0;c<(nsink*ncoef);c++)
         A[(nrows_U+r)+ncnt*c] = Lw * (*_LC)[r+nrows_L*c];

      int coffset = nsink*ncoef , cmax = ndof-coffset;  

      for(int c=0;c<cmax;c++)
         A[(nrows_U+r) + ncnt*(coffset+c)] = Lw * (*_LR)[r+nrows_L*c];
   }

   fprintf(stderr,"Done.\n");
   fflush(stderr);

   fprintf(stderr,"Solving system... ");
   fflush(stderr);

   ////////////////////////////////////////////////////////////////////
   /////////////// Solve system and retrieve parameters ///////////////
   ////////////////////////////////////////////////////////////////////

   LWORK = 4*ndof; // MN = min(nrows, ncols) and LWORK >= max( 1, MN + max( MN, NRHS )*NB ) = ndof + ndof*OptimBlockSize
   WORK = (double *) malloc(LWORK*sizeof(double));

//   fprintf(stderr,"Solving system with %d degrees of freedom and %d control points...\n",ndof,nrows);
//   fflush(stderr);

   (void) dgels_ ( &trans, &ncnt, &ndof, &nrhs, A, &ncnt,B, &ncnt, WORK, &LWORK, &INFO);
   // SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,INFO )

   fprintf(stderr,"Done.\n");
   fflush(stderr);
   
   free(A);
   free(WORK);
   
   (*_Phi0) = fit_trend ? B[ndof-3] : B[ndof-1];
   (*_v0x)  = fit_trend ? B[ndof-2] : 0.;
   (*_v0y)  = fit_trend ? B[ndof-1] : 0.;

   // Store discharge for sink slits

//   for(int s=0;s<nsink;s++)
//     (*_Q)[s] = 0.;
  
   // Set coefficients of Laurent series to zero for sink slits

   for(int c=0;c<(nsink*ncoef);c++)
   {
      (*_ReCn)[c] = B[c]; 
      (*_ImCn)[c] = 0.;    
   }

   // Set sink term to zero for boundary slits

   for(int b=0;b<nbound;b++)
      (*_Q)[nsink+b] = 0.;

   // Store coefficients of Laurent series at boundary slits

   for(int c=0;c<(nbound*ncoef);c++)
   {
      (*_ReCn)[(nsink*ncoef)+c] = B[(nsink*ncoef)+2*c]; 
      (*_ImCn)[(nsink*ncoef)+c] = B[(nsink*ncoef)+2*c+1];    
   }

   // Update value of potential at all sink control points
   // (even at unactive ones because they could be further reactivated)

   for(int row=0;row<(nsink*ntheta);row++)
   {
      (*_PHI)[row] = fabs(gamma0) * (*_PHI_areal1)[row];

      for(int s=0;s<nsink;s++)
        (*_PHI)[row] += fabs(gamma0) * ((*_UL)[row+nrows_U*s])*((*_dA)[s]);

      for(int colUC=0;colUC<(nsink*ncoef);colUC++)
        (*_PHI)[row] += (*_UC)[colUC*nrows_U+row]*B[colUC];

      int coffset = nsink*ncoef, cmax = ndof-coffset;

      for(int colUR=0;colUR<cmax;colUR++)
        (*_PHI)[row] += (*_UR)[colUR*nrows_U+row]*B[coffset+colUR];
   }

   free(B);

   return(0);
}

int get_LHS_blocks_p (	int np, double ** _Xp, double ** _Yp, 
	                int nsink, int nbound, int ncoef,
                        double ** _ReZstart, double ** _ImZstart, double ** _ReZend, double ** _ImZend,
                        double ** _ReCn, double ** _ImCn, double ** _Q,
                        double ** _PL, double ** _PC, double ** _PR, double ** _PHI_areal1_p 
                      )
{
  int one=1;
  int col,col_J;
  int varout[3] = {0,1,0};
  double PHI1,PSI1,Vx1,Vy1;
  double * J_Vx = NULL, * J_Vy = NULL;
  double v0x = 0., v0y = 0., Phi0 = 0.;

  int nslit = nsink+nbound;
  double * J_PHI = (double *) malloc((2*ncoef+1)*nslit*sizeof(double));
  double * J_PSI = (double *) malloc((2*ncoef+1)*nslit*sizeof(double));

  (*_PL) = malloc( np *       nsink        * sizeof(double)); // Jacobian block w.r.t. Q's at sink slits
  (*_PC) = malloc( np *   (ncoef*nsink)    * sizeof(double)); // Jacobian block w.r.t. Re{Cn}'s at sink slits
  (*_PR) = malloc( np * (2*ncoef*nbound+3) * sizeof(double)); // Jacobian block w.r.t. Re{Cn}'s and Im{Cn}'s at boundary slits and Phi0, v0x, v0y

  (*_PHI_areal1_p) = malloc(np*sizeof(double));

  for(int row=0;row<np;row++)
  {
     (void) slit_ ( (*_Xp)+row,(*_Yp)+row, &one,
	   	    (*_ReZstart),(*_ImZstart),(*_ReZend),(*_ImZend),&nslit,   
		    (*_ReCn),(*_ImCn),&ncoef,(*_Q),&v0x,&v0y,&Phi0,
		    varout,
		    &PHI1,&PSI1,&Vx1,&Vy1,J_PHI,J_PSI,J_Vx,J_Vy);

     for(int sink=0;sink<nsink;sink++){
       col = sink;
       col_J = (2*ncoef+1)*sink + 2*ncoef;
       (*_PL)[row+np*col] = J_PHI[col_J]; 
     }

     for(int sink=0;sink<nsink;sink++){
     for(int coef=0;coef<ncoef;coef++){
       col = ncoef*sink + coef;
       col_J = (2*ncoef+1)*sink + 2*coef;   // Re{Cn} only
       (*_PC)[row+np*col] = J_PHI[col_J]; 
     }}

     for(int bound=0;bound<nbound;bound++){
     for(int coef=0;coef<(2*ncoef);coef++){
       col = (2*ncoef)*bound + coef;
       col_J = (2*ncoef+1)*(nsink+bound) + coef;
       (*_PR)[row+np*col] = J_PHI[col_J];
     }}

     int col = (2*nbound*ncoef);
     (*_PR)[row+np*col] = +1.0; // w.r.t. Phi0

     col++;
     (*_PR)[row+np*col] = -1.0*((*_Xp)[row]-X0); // w.r.t. v0x

     col++;
     (*_PR)[row+np*col] = -1.0*((*_Yp)[row]-Y0);  // w.r.t. v0y

     double gamma0 = -1.0, A;

     (void) arealsink_ ( (*_ReZstart)+nsink,(*_ImZstart)+nsink,&nbound,   
                         (*_Xp)+row,(*_Yp)+row,&one,&gamma0,
                         &PHI1,&Vx1,&Vy1,&A);

     (*_PHI_areal1_p)[row] = PHI1;
  }

  free(J_PHI);free(J_PSI);

  return(0);
}

/** 
# Solve full system with areal sink and borehole data
*/

int solve_recharge_p ( int nsink, int nbound, int ncoef, int ntheta, int np,
                       double ** _UL, double ** _UC, double ** _UR,
                       double ** _LL, double ** _LC, double ** _LR,
                       double ** _PL, double ** _PC, double ** _PR,
                       double ** _PHI_areal1, double ** _Vn_areal1,  
                       double ** _PHI_areal1_p,  
                       double ** _ReZ, double ** _ImZ, double ** _Xp, double ** _Yp,
                       double ** _SlitLength, double ** _PHI_topo, double ** _PHI_p,
                       double ** _dA, double ** _Aup,
                       double * _gamma0, double ** _Q, double ** _ReCn, double ** _ImCn, double * _v0x, double * _v0y, double * _Phi0,
                       double ** _PHI )  
{
   double * A, * B; // for DGELS
   char trans = 'N';
   int nrhs = 1;
   int LWORK, INFO;
   double * WORK;

   bool fit_trend =  false;

   /////////////////////////////////////////////////////////////////////////
   /////////////// Fill matrix of the system, and RHS vector ///////////////
   /////////////////////////////////////////////////////////////////////////
   
   int nslit = nsink+nbound;
   int ncnt = ntheta*(nsink+nbound)+np;
   int ndof = 1+(ncoef*nsink)+(2*ncoef*nbound); // γ0, Re{Cn}'s for sinks, Re{Cn}'s and Im{Cn}'s for boundary segments
   ndof += fit_trend ? 3 : 1; // add Phi0, as well as v0x and v0y if needed

   A = (double *) malloc(ncnt*ndof*sizeof(double));
   B = (double *) malloc(ncnt*sizeof(double));

   fprintf(stderr,"Fill matrix of the system with %d degrees of freedom and %d control points... ",ndof,ncnt);
   fflush(stderr);

   int nrows_U = ntheta*nsink, nrows_L = ntheta*nbound;
   int coffset, cmax;

   // Fill upper part of A and B (constraints on potential in the thalweg network)

   for(int r=0;r<nrows_U;r++)
   {
      double w = sqrt(((double)ncnt)/((double)nrows_U));

      B[r] = w * (*_PHI_topo)[r] ; // RHS is topographic elevation
      double Ar1 = (*_PHI_areal1)[r]; // Jacobian w.r.t γ0 in first column

      for(int s=0;s<nsink;s++)
         Ar1 += ((*_UL)[r+nrows_U*s])*((*_dA)[s]);

      A[r] = Ar1 * w;

      // In the upper center block, just copy/paste UC

      coffset = 1;
      cmax = (nsink*ncoef);  

      for(int c=0;c<cmax;c++)
         A[r+ncnt*(coffset+c)] = w * (*_UC)[r+nrows_U*c];

      // In the upper right block, just copy/paste UR

      coffset = 1+(nsink*ncoef);
      cmax = fit_trend ? (2*ncoef*nbound)+3 : (2*ncoef*nbound)+1;  

      for(int c=0;c<cmax;c++)
         A[r+ncnt*(coffset+c)] = w * (*_UR)[r+nrows_U*c];
   }

   // Fill lower part of A and B (constraints on normal velocity on the boundary)
   
   for(int r=0;r<nrows_L;r++)
   {
      int slit = nsink + (r / ntheta);
      double Lw = sqrt(((double)ncnt)/((double)nrows_L)) * (*_SlitLength)[slit];
//      double Lw = 1000.;

      B[nrows_U+r] = 0. ; // RHS is zero normal velocity component
      double Ar1 = (*_Vn_areal1)[r]; // Jacobian w.r.t γ0 in first column
      
      for(int s=0;s<nsink;s++)
        Ar1 += ((*_LL)[r+nrows_L*s])*((*_dA)[s]); 

      // Store it in position A(r+nrows_U,1) with row-major order, weighted by Lw
      A[nrows_U+r] = Lw * Ar1;

      // In the lower center block, just copy/paste LC 

      coffset = 1;
      cmax = (nsink*ncoef);  

      for(int c=0;c<cmax;c++)
         A[(nrows_U+r) + ncnt*(coffset+c)] = Lw * (*_LC)[r+nrows_L*c];

      // In the lower right block, just copy/paste LR 

      coffset = 1+(nsink*ncoef);
      cmax = fit_trend ? (2*ncoef*nbound)+3 : (2*ncoef*nbound)+1;

      for(int c=0;c<cmax;c++)
         A[(nrows_U+r) + ncnt*(coffset+c)] = Lw * (*_LR)[r+nrows_L*c];         
   }

   // Fill lowest part of A and B (constraints on potential at well data points)

   for(int r=0;r<np;r++)
   {
      double w = sqrt(0.5*((double)ncnt)/((double)np));

      B[(nslit*ntheta)+r] = w * (*_PHI_p)[r] ; // RHS is groundwater level data
      double Ar1 = (*_PHI_areal1_p)[r]; // Jacobian w.r.t γ0 in first column

      for(int s=0;s<nsink;s++)
         Ar1 += ((*_PL)[r+np*s])*((*_dA)[s]);

      A[(nslit*ntheta)+r] = Ar1 * w;

      // In the lowest center block, just copy/paste PC

      coffset = 1;
      cmax = (nsink*ncoef);  

      for(int c=0;c<cmax;c++)
         A[(nslit*ntheta)+r+ncnt*(coffset+c)] = w * (*_PC)[r+np*c];

      // In the lowest right block, just copy/paste PR

      coffset = 1+(nsink*ncoef);
      cmax = fit_trend ? (2*ncoef*nbound)+3 : (2*ncoef*nbound)+1;  

      for(int c=0;c<cmax;c++)
         A[(nslit*ntheta)+r+ncnt*(coffset+c)] = w * (*_PR)[r+np*c];
   }

   fprintf(stderr,"Done.\n");
   fflush(stderr);

   fprintf(stderr,"Solving system... ");
   fflush(stderr);

   ////////////////////////////////////////////////////////////////////
   /////////////// Solve system and retrieve parameters ///////////////
   ////////////////////////////////////////////////////////////////////

   LWORK = 4*ndof; // MN = min(nrows, ncols) and LWORK >= max( 1, MN + max( MN, NRHS )*NB ) = ndof + ndof*OptimBlockSize
   WORK = (double *) malloc(LWORK*sizeof(double));

//   fprintf(stderr,"Solving system with %d degrees of freedom and %d control points...\n",ndof,nrows);
//   fflush(stderr);

   (void) dgels_ ( &trans, &ncnt, &ndof, &nrhs, A, &ncnt,B, &ncnt, WORK, &LWORK, &INFO);
   // SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,INFO )

   fprintf(stderr,"Done.\n");
   fflush(stderr);
   
   free(A);
   free(WORK);
   
   (*_gamma0) = -1.*B[0];

   (*_Phi0) = fit_trend ? B[ndof-3] : B[ndof-1];
   (*_v0x)  = fit_trend ? B[ndof-2] : 0.;
   (*_v0y)  = fit_trend ? B[ndof-1] : 0.;

   // Store discharge for sink slits

   double sumQ = 0.;

   for(int s=0;s<nsink;s++)
   {
     (*_Q)[s] = B[0] * (*_dA)[s];
     sumQ += (*_Q)[s];
   }
  
   fprintf(stderr,"sumQ = %.4e, gamma0 = %.4e\n",sumQ,(*_gamma0));
   fflush(stderr);

   // Set coefficients of Laurent series to zero for sink slits

   for(int c=0;c<(nsink*ncoef);c++)
   {
      (*_ReCn)[c] = B[1+c]; 
      (*_ImCn)[c] = 0.;    
   }

   // Set sink term to zero for boundary slits

   for(int b=0;b<nbound;b++)
      (*_Q)[nsink+b] = 0.;

   // Store coefficients of Laurent series at boundary slits

   for(int c=0;c<(nbound*ncoef);c++)
   {
      (*_ReCn)[(nsink*ncoef)+c] = B[1+(nsink*ncoef)+2*c]; 
      (*_ImCn)[(nsink*ncoef)+c] = B[1+(nsink*ncoef)+2*c+1];    
   }

   // Update value of potential at all sink control points
   // (even at unactive ones because they could be further reactivated)

   for(int row=0;row<(nsink*ntheta);row++)
   {
      (*_PHI)[row] = B[0]*(*_PHI_areal1)[row];

      for(int s=0;s<nsink;s++)
        (*_PHI)[row] += (*_UL)[s*nrows_U+row]*((*_Q)[s]);

      for(int colUC=0;colUC<(nsink*ncoef);colUC++)
        (*_PHI)[row] += (*_UC)[colUC*nrows_U+row]*B[1+colUC];

      cmax = fit_trend ? (2*ncoef*nbound)+3 : (2*ncoef*nbound)+1;

      for(int colUR=0;colUR<cmax;colUR++)
        (*_PHI)[row] += (*_UR)[colUR*nrows_U+row]*B[1+(nsink*ncoef)+colUR];
   }

   FILE * fgw = fopen("gw_rmse.csv","w");
   fprintf(fgw,"X;Y;Hobs;PHIp\n");

   double rmse = 0.;

   for(int row=0;row<np;row++)
   {
      double PHIp = B[0]*(*_PHI_areal1_p)[row];

      for(int s=0;s<nsink;s++)
        PHIp += (*_PL)[s*np+row]*((*_Q)[s]);

      for(int colPC=0;colPC<(nsink*ncoef);colPC++)
        PHIp += (*_PC)[colPC*np+row]*B[1+colPC];

      cmax = fit_trend ? (2*ncoef*nbound)+3 : (2*ncoef*nbound)+1;

      for(int colPR=0;colPR<cmax;colPR++)
        PHIp += (*_PR)[colPR*np+row]*B[1+(nsink*ncoef)+colPR];

      fprintf(fgw,"%.2lf;%.2lf;%.2lf;%.2lf\n",(*_Xp)[row],(*_Yp)[row],(*_PHI_p)[row],PHIp);
      fflush(fgw);

      rmse += sq(PHIp-(*_PHI_p)[row]);
   }
   fclose(fgw);
   fprintf(stderr,"RMSE piezo = %.2lf m\n",sqrt(rmse/np));

   free(B);

   return(0);
}