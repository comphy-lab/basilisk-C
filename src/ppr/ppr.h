/*
   NPOS  no. edges in old grid.
   NNEW  no. edges in new grid.
   NVAR  no. discrete variables to remap.
   NDOF  no. degrees-of-freedom per cell.
   XPOS  old grid edge positions. XPOS is a length NPOS array.
   XNEW  new grid edge positions. XNEW is a length NNEW array.
   FDAT  grid-cell moments on old grid. FNEW has SIZE = NDOF-by-NVAR-by-NPOS - 1.
   FNEW  grid-cell moments on new grid. FNEW has SIZE = NDOF-by-NVAR-by-NNEW - 1.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#pragma autolink -L$BASILISK/ppr -lppr -lgfortran

#define p1e_method 100
#define p3e_method 101
#define p5e_method 102

#define pcm_method 200
#define plm_method 201
#define ppm_method 202
#define pqm_method 203

#define null_limit 300
#define mono_limit 301
#define weno_limit 302

#define bcon_loose 400
#define bcon_value 401
#define bcon_slope 402

void my_remap (int * npos, int * nnew, int * nvar, int * ndof,
	       double * xpos, double * xnew,
	       double * fdat, double * fnew,
	       int * edge_meth, int * cell_meth, int * cell_lim);
