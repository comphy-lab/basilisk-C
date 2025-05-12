
/**
Here you will find some tools to manipulate tensors, and in particular some linear algebra tools (eigen values and vectors of a symmetric matrix).

*/



/**
Pseudo tensor and vector
========================

These structures are similar to basilisk tensors and vectors, but they do not store all the fields, they are only used locally in foreach() loops. Useful to use less memory.

*N.B.* : these structures are originally defined in [log-conform](http://basilisk.fr/src/log-conform.h), but with events which occurs every iteration. So include the log-conform library would result in creating and calculating a lot of variables for nothing.

*/

#if dimension==2
typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;
#endif
#if dimension==3
typedef struct { double x, y, z;}   pseudo_v;
typedef struct { pseudo_v x, y, z;} pseudo_t;
#endif


/**
Function declaration

*/
void symmetrize_tensor (tensor t);
void multiply_tensor (Point point, tensor a, tensor b, tensor result);
void transpose_tensor (Point point, tensor a, tensor t_a);

void matrix_from_tensor (double mat[dimension][dimension], pseudo_t t);
void tensor_from_matrix (pseudo_t t, double mat[dimension][dimension]);

void matrix_from_pseudo_t (double mat[dimension][dimension], pseudo_t t);
pseudo_t pseudo_t_from_matrix (double mat[dimension][dimension]);


void diagonal_tensor_from_vector (Point point, vector v, tensor t);
void trace_tensor (Point point, tensor t, scalar trace);



/**
Eigeenvalues and vectors of symmetric matrix
============================================

   They are computed in basilisk in 
   [src/lambda2.h](http://basilisk.fr/src/lambda2.h), in order to compute the
   [lambda2](https://en.wikipedia.org/wiki/Lambda2_method) criterion on vortex.

   However, write #include "lambda2.h" does not work. This is due to the function 
   lambda2 defined in the file, which uses the .z component of a vector : that 
   makes that the code can't compile in 2D due to a function we don't use... **So, 
   for now, the functions eigsrt and eigenvalues are just copied from 
   [src/lambda2.h](http://basilisk.fr/src/lambda2.h)**, and put here.
   (copy 20/01/2020).

   **begin of copy**

   */




static void eigsrt (double d[dimension],
		    double v[dimension][dimension])
{
   int k, j, i;
   double p;

   for (i = 0; i < dimension - 1; i++) {
      p = d[k = i];

      for (j = i + 1; j < dimension; j++)
	 if (d[j] >= p) 
	    p = d[k = j];
      if (k != i) {
	 d[k] = d[i];
	 d[i] = p;
	 for (j = 0; j < dimension; j++) {
	    p = v[j][i];
	    v[j][i] = v[j][k];
	    v[j][k] = p;
	 }
      }
   }
}

#define ROTATE(a,i,j,k,l) {						\
      g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}


/**
   eigenvalues:
   @a: a symmetric matrix.
   @d: a vector.
   @v: another matrix.
   Fills @d (resp. @v) with the eigenvalues (resp. eigenvectors) of
   matrix @a.
*/


void eigenvalues (double a[dimension][dimension],
		  double d[dimension],
		  double v[dimension][dimension])
{
   int j, iq, ip, i;
   double tresh, theta, tau, t, sm, s, h, g, c, b[dimension], z[dimension];

   for (ip = 0; ip < dimension; ip++) {
      for (iq = 0; iq < dimension; iq++)
	 v[ip][iq] = 0.0;
      v[ip][ip] = 1.0;
   }

   for (ip = 0; ip < dimension; ip++) {
      b[ip] = d[ip] = a[ip][ip];
      z[ip] = 0.0;
   }

   for (i = 1; i <= 50; i++) {
      sm = 0.0;
      for (ip = 0; ip < dimension - 1; ip++) {
	 for (iq = ip + 1; iq < dimension; iq++)
	    sm += fabs (a[ip][iq]);
      }
      if (sm == 0.0) {
	 eigsrt (d, v);
	 return;
      }
      if (i < 4)
	 tresh = 0.2*sm/(dimension*dimension);
      else
	 tresh = 0.0;
      for (ip = 0; ip < dimension - 1; ip++) {
	 for (iq = ip + 1; iq < dimension; iq++) {
	    g = 100.0*fabs (a[ip][iq]);
	    if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&
		fabs(d[iq]) + g == fabs(d[iq]))
	       a[ip][iq] = 0.0;
	    else if (fabs (a[ip][iq]) > tresh) {
	       h = d[iq] - d[ip];
	       if (fabs(h) + g == fabs(h))
		  t = a[ip][iq]/h;
	       else {
		  theta = 0.5*h/a[ip][iq];
		  t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
		  if (theta < 0.0) t = -t;
	       }
	       c = 1.0/sqrt (1 + t*t);
	       s = t*c;
	       tau = s/(1.0 + c);
	       h = t*a[ip][iq];
	       z[ip] -= h;
	       z[iq] += h;
	       d[ip] -= h;
	       d[iq] += h;
	       a[ip][iq] = 0.0;
	       for (j = 0; j <= ip - 1; j++)
		  ROTATE (a, j, ip, j, iq);
	       for (j = ip + 1; j <= iq - 1; j++)
		  ROTATE (a, ip, j, j, iq);
	       for (j = iq + 1; j < dimension; j++)
		  ROTATE(a, ip, j, iq, j);
	       for (j = 0; j < dimension; j++)
		  ROTATE(v, j, ip, j, iq);
	    }
	 }
      }
      for (ip = 0; ip < dimension; ip++) {
	 b[ip] += z[ip];
	 d[ip] = b[ip];
	 z[ip] = 0.0;
      }
   }
   /* Too many iterations */
   for (i = 0; i < dimension; i++) {
      for (j = 0; j < dimension; j++)
	 fprintf (stderr, "%10.3g ", a[i][j]);
      fprintf (stderr, "\n");
   }
   assert (false);
}


/**
 **End of the copy**

 */

/**
Tensor eigen values and vectors
-------------------------------

   The above function for eigenvalues and vectors is written for matrix and not 
   basilisk tensors. 

*/

void tensor_eigen(pseudo_t t, pseudo_v * eig_val, pseudo_t * eig_vec) {

   double matrix[dimension][dimension];
   double eig_vec_matrix[dimension][dimension];
   double eig_val_vector[dimension];

   // store tensor in matrix
   matrix_from_pseudo_t(matrix, t);

   // compute eigenvalues
   eigenvalues (matrix, eig_val_vector, eig_vec_matrix);

   // store eigen values in basilisk vector
   eig_val->x = eig_val_vector[0];
#if dimension > 1
   eig_val->y = eig_val_vector[1];
#endif
#if dimension > 2
   eig_val->z = eig_val_vector[2];
#endif

   // store rotation tensor (i.e. eigen vectors matrix as they are orthonormal)
   *eig_vec = pseudo_t_from_matrix(eig_vec_matrix);

}

/**
Other tools
===========

*/




/**
Symmetric local tensor
----------------------

   $H_R(u) = \frac{H^*(u) + H^*(u)^t}{2}$

   *N.B.* : we use the fact that the diagonal elements are not changing. 

   It modify the tensor in place.

*/


pseudo_t local_symmetrize_pseudo_t (pseudo_t t) {

   pseudo_t tmp;

   // same diagonal elements
   foreach_dimension()
      tmp.x.x = t.x.x;

   // mean of non-diagonal elements
#if dimension > 1
   tmp.x.y = ( t.x.y + t.y.x ) / 2.;
   tmp.y.x = t.x.y;
#if dimension > 2
   tmp.x.z = ( t.x.z + t.z.x ) / 2.;
   tmp.z.x = t.x.z;
   tmp.y.z = ( t.y.z + t.z.y ) / 2.;
   tmp.z.y = t.y.z;
#endif
#endif

   return tmp;
}





/**
Multiplication of local tensor
------------------------------

   Multiplication of a pseudo-tensor by a pseudo-tensor in a pseudo-tensor.

*/

pseudo_t multiply_local_tensor (pseudo_t a, pseudo_t b) {
	
   pseudo_t result;

   foreach_dimension() {
#if dimension == 1
      result.x.x = a.x.x*b.x.x;
#endif
#if dimension == 2
      result.x.x = a.x.x*b.x.x + a.x.y*b.y.x;
      result.x.y = a.x.x*b.x.y + a.x.y*b.y.y;
#endif
#if dimension == 3
      result.x.x = a.x.x*b.x.x + a.x.y*b.y.x + a.x.z*b.z.x;
      result.x.y = a.x.x*b.x.y + a.x.y*b.y.y + a.x.z*b.z.y;
      result.x.z = a.x.x*b.x.z + a.x.y*b.y.z + a.x.z*b.z.z;
#endif
   }

   return result;
}





/**
Transposition of local tensor
-----------------------------

   Transposition of a pseudo-tensor.

*/


pseudo_t transpose_pseudo_t (pseudo_t t) {

   pseudo_t tmp;

   foreach_dimension() {
      tmp.x.x = t.x.x;
#if dimension > 1
      tmp.x.y = t.y.x;
#endif
#if dimension > 2
      tmp.x.z = t.z.x;
#endif
   }

   return tmp;
}






/**
Tensor - matrix - local tensor transformations
----------------------------------------------

   Transform basilisk tensor in a pseudo-tensor.

*/

pseudo_t create_local_tensor (Point point, tensor t) {

   pseudo_t tmp;

   tmp.x.x = t.x.x[];
#if dimension > 1
   tmp.y.y = t.y.y[];
   tmp.x.y = t.x.y[];
   tmp.y.x = t.y.x[]; 
#endif
#if dimension > 2
   tmp.z.z = t.z.z[];
   tmp.x.z = t.x.z[];
   tmp.y.z = t.y.z[];
   tmp.z.x = t.z.x[];
   tmp.z.y = t.z.y[];
#endif

   return tmp;
}

/**
   Transform basilisk tensor in a matrix.

*/

void matrix_from_pseudo_t (double mat[dimension][dimension], pseudo_t t) {
	
   mat[0][0] = t.x.x;
#if dimension > 1
   mat[1][1] = t.y.y;
   mat[0][1] = t.x.y;
   mat[1][0] = t.y.x; 
#endif
#if dimension > 2
   mat[2][2] = t.z.z;
   mat[0][2] = t.x.z;
   mat[1][2] = t.y.z;
   mat[2][0] = t.z.x;
   mat[2][1] = t.z.y;
#endif
}


/**
   Transform a matrix in a basilisk tensor

*/

pseudo_t pseudo_t_from_matrix (double mat[dimension][dimension]) {
	
   pseudo_t t;

   t.x.x = mat[0][0];
#if dimension > 1
   t.y.y = mat[1][1];
   t.x.y = mat[0][1];
   t.y.x = mat[1][0];
#endif
#if dimension > 2
   t.z.z = mat[2][2];
   t.x.z = mat[0][2];
   t.y.z = mat[1][2];
   t.z.x = mat[2][0];
   t.z.y = mat[2][1];
#endif

   return t;
}




/**
Create diagonal tensor from a vector
------------------------------------

   Useful to create the diagonal matrix containing the eigen values of a tensor, to 
   use tensor multiplication : $M = R\Lambda R^{-1}$.

*/



pseudo_t diagonal_pseudo_t_from_pseudo_v (pseudo_v v) {

   pseudo_t tmp;

   foreach_dimension() {
      tmp.x.x = v.x;
#if dimension > 1
      tmp.x.y = 0.;
#endif
#if dimension > 2
      tmp.x.z = 0.;
#endif
   }

   return tmp;
}




/**
Trace of a tensor
-----------------

   Compute the trace of a pseudo-tensor.

*/

double trace_pseudo_t (pseudo_t t) {

  double vtrace = 0.;
   foreach_dimension() 
      vtrace += t.x.x;
   return vtrace;
}










