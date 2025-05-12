/**
# Higher-order-accurate functions and definitions

## Additional attributes for prolongation and refinement.
   
An implementation of 1st upto 5th-order accurate prolongation methods
for centered scalar fields in upto three dimensions, called
`refine_1st;` upto `refine_5th;`, respectively, via `_2nd, _3rd and
_4th`. The uneven-order methods are conservative.

The 1D interplation methods for averaged values:
*/
static inline double int_1pt (double * a) { //mimics Basilisk's `injection`
  return a[0];
}

static inline double int_2pt (double * a) { //linear (mimics Basilisk's `bilinear`)
  return (3.*a[1] + a[0])/4.;
}

static inline double int_3pt (double * a) { //quadratic (i.e. Basilisk's `linear` in 1D)
  return (8.*a[1] + a[0] - a[2])/8.;
}

static inline double int_4pt (double * a) { //cubic
  return (-3*a[0] + 17*a[1] + 55*a[2] - 5*a[3])/64.;
}

static inline double int_5pt (double * a) { //quartic
  return   (-3*a[0] + 22.*a[1] + 128.*a[2] - 22*a[3] + 3.*a[4])/128.;
}

/**
A function that iteratively reduces the dimensionality of the neighbor
data, at a given `order`.
 */
static inline double mul_xpt (Point point, scalar s, int order) {
  int start = -order/2, end = (int)(((double)order/2.) + 0.6);
  double a[order], (*int_xpt)(double *);
  if (order < 2)
    int_xpt = int_1pt;
  else if (order < 3)
    int_xpt = int_2pt;
  else if (order < 4)
    int_xpt = int_3pt;
  else if (order < 5)
    int_xpt = int_4pt;
  else 
    int_xpt = int_5pt;
  
#if (dimension == 1)
  for (int j = start; j < end; j++)
    a[j - start] = coarse(s,-j*child.x);
#elif (dimension == 2)
  for (int j = start; j < end; j++) {
    double b[order];
    for (int k = start; k < end; k++)
      b[k - start] = coarse(s,-j*child.x, -k*child.y);
    a[j - start] = int_xpt (b);
  }
#else // dimension == 3
  for (int j = start; j < end; j++) {
    double b[order];
    for (int k = start; k < end; k++) {
      double c[order];
      for (int m = start; m < end; m++)
	c[m - start] = coarse(s,-j*child.x, -k*child.y, -m*child.z);
      b[k - start] = int_xpt (c) ;
    }
    a[j - start] = int_xpt (b);
  }
#endif
  return int_xpt (a);
}
/**
The refinement attributes:
 */
trace
static inline void refine_1st (Point point, scalar s) {
  foreach_child()
    s[] = mul_xpt (point, s, 1);
}
 trace
static inline void refine_2nd (Point point, scalar s) {
  foreach_child()
    s[] = mul_xpt (point, s, 2);
}
trace
static inline void refine_3rd (Point point, scalar s) {
  foreach_child()
    s[] = mul_xpt (point, s, 3);
}
trace
 static inline void refine_4th (Point point, scalar s) {
  foreach_child()
    s[] = mul_xpt (point, s, 4);
}
trace
static inline void refine_5th (Point point, scalar s) {
  foreach_child()
    s[] = mul_xpt (point, s, 5);
}
/**
## 4th-order face refinement in 2D. 

   Following the perspective's of [Rajarshi Roy Chowdhury's
   thesis](https://hal.archives-ouvertes.fr/tel-02056238). We implement a
   high-order face-refinement method.

   See also the [test page](rf4t.c).

   we destinguish between
   
   * the semi-colocated fine faces prolongation  
   * the inner-cell face values 
   
   The first catogory uses "simple" interpolation at a lower
   dimension.

The second catogory requires a staggered grid interpolation before we
may do the actual prolongation.

The work horse is the 5th order 5-point 1D refinement, which is
conservative.
*/
static inline void refine_order_5 (double val[5], double * lr) {
  lr[0] =  (-3*val[0] + 22.*val[1] + 128.*val[2] +
	    -22*val[3] + 3.*val[4])/128.;
  lr[1] = 2*val[2] - lr[0];
}

static inline void colocation_values (double val[4][5], double * intrp) {
  for (int j = 0; j < 5; j++)
      intrp[j] = (9.*(val[2][j] + val[1][j]) - (val[0][j] + val[3][j]))/16.;
  //!=  (7.*(val[2][j] + val[1][j]) - (val[0][j] + val[3][j]))/12.;
}
/**
Alternatively these 4/5th-order accurate 4/5 point methods could be
replaced by 4th-order accurate 5-point ols methods.
 */
static inline void refine_order_4_ols (double val[5], double * lr) {
  double a = //rare coefficients
    -489./4480. * val[0] +
    1153./2240. * val[1] +
      17./35.   * val[2]   +
     383./2240. * val[3]  +
    -279./4480. * val[4];
  lr[0] = a ;
  lr[1] = 2*val[0] - a;
}

static inline void colocation_values_ols (double val[5][5], double * intrp) {
  for (int j = 0; j < 5; j++)
      intrp[j] = (-3.*val[0][j] + 2*val[1][j] + 72.*val[2][j]
		  + 102*val[3][j] - 13.*val[4][j])/160.;
}
/**
   These functions can be used to construct a face-refinement
   method
 */
#if (TREE)
#if (dimension == 2)
foreach_dimension()
void refine_face_4_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    double val[5], lr[1 << (dimension - 1)];
    for (int j = -2; j < 3; j++)
      val[j + 2] = v.x[0,j];
    refine_order_5 (val, &lr[0]);
    fine(v.x,0,0,0) = lr[0];
    fine(v.x,0,1,0) = lr[1];
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {//rhs 
    double val[5], lr[2];
    for (int j = -2; j < 3; j++)
      val[j + 2] = v.x[1,j];
    refine_order_5 (val, &lr[0]);
    fine(v.x,2,0,0) = lr[0];
    fine(v.x,2,1,0) = lr[1];
  }
  if (is_local(cell)) {//cross of faces
    double vals[4][5], val[5], lr[2];
    for (int k = -1; k < 3; k++) 
      for (int j = -2; j < 3; j++)
	vals[k + 1][j + 2] = v.x[k,j];
    colocation_values (vals, &val[0]);
    refine_order_5 (val, &lr[0]);
    fine(v.x,1,0,0) = lr[0];
    fine(v.x,1,1,0) = lr[1];
  }
}
#else // (dimension == 3)
foreach_dimension()
void refine_face_4_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    for (int cx = -1; cx < 2; cx += 2) {
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = - 2; j < 3; j++) {
	  double b[5];
	  for (int k = - 2; k < 3; k++)
	    b[k + 2] = v.x[0, -j*cx, -k*cy];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 0, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
    }
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {//rhs 
    for (int cx = -1; cx < 2; cx += 2) 
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = -2; j < 3; j++) {
	  double b[5];
	  for (int k = -2; k < 3; k++)
	    b[k + 2] = v.x[1, -j*cx, -k*cy];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 2, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
  }
  if (is_local(cell)) { //cross of faces
    double vals[5][5]; //4th-order accurate co-location values
    for (int j = -2; j < 3; j++) 
      for (int k = -2; k < 3; k++) 
	vals[j + 2][k + 2] = (9.*(v.x[0,j,k] + v.x[1,j,k]) - (v.x[-1,j,k] + v.x[2,j,k]))/16.;
    for (int cx = -1; cx < 2; cx += 2) 
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = -2; j < 3; j++) {
	  double b[5];
	  for (int k = -2; k < 3; k++)
	    b[k + 2] = vals[-j*cx + 2][-k*cy + 2];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 1, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
  }
}
#endif
#endif
/**
and a prolongation methods for error-estimation
 */

#if (TREE)
#if (dimension == 2)
foreach_dimension()
void interpolant_face_4_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    double val[5], lr[1 << (dimension - 1)];
  for (int j = -2; j < 3; j++)
    val[j + 2] = v.x[0,j];
  refine_order_5 (val, &lr[0]);
  fine(v.x,0,0,0) = lr[0];
  fine(v.x,0,1,0) = lr[1];
  }
  if (is_local(cell)) {
    double vals[4][5], val[5], lr[2];
    for (int k = -1; k < 3; k++) 
      for (int j = -2; j < 3; j++)
	vals[k + 1][j + 2] = v.x[k,j];
    colocation_values (vals, &val[0]);
    refine_order_5 (val, &lr[0]);
    fine(v.x,1,0,0) = lr[0];
    fine(v.x,1,1,0) = lr[1];
  }
}

#else // (dimension == 3)
foreach_dimension()
void interpolant_face_4_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    for (int cx = -1; cx < 2; cx += 2) {
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = - 2; j < 3; j++) {
	  double b[5];
	  for (int k = - 2; k < 3; k++)
	    b[k + 2] = v.x[0, -j*cx, -k*cy];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 0, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
    }
  }
  if (is_local(cell)) {
    double vals[5][5]; //4th-order accurate co-location values
    for (int j = -2; j < 3; j++) 
      for (int k = -2; k < 3; k++) 
	vals[j + 2][k + 2] = (9.*(v.x[0,j,k] + v.x[1,j,k]) - (v.x[-1,j,k] + v.x[2,j,k]))/16.;
    for (int cx = -1; cx < 2; cx += 2) 
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = -2; j < 3; j++) {
	  double b[5];
	  for (int k = -2; k < 3; k++)
	    b[k + 2] = vals[-j*cx + 2][-k*cy + 2];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 1, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
  }
}
#endif
#endif



/**
And their OLS counterparts:
 */

#if (TREE)
#if (dimension == 2)
foreach_dimension()
void refine_face_4_ols_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    double val[5], lr[1 << (dimension - 1)];
    for (int j = -2; j < 3; j++)
      val[j + 2] = v.x[0,j];
    refine_order_4_ols (val, &lr[0]);
    fine(v.x,0,0,0) = lr[0];
    fine(v.x,0,1,0) = lr[1];
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {//rhs 
    double val[5], lr[2];
    for (int j = -2; j < 3; j++)
      val[j + 2] = v.x[1,j];
    refine_order_4_ols (val, &lr[0]);
    fine(v.x,2,0,0) = lr[0];
    fine(v.x,2,1,0) = lr[1];
  }
  if (is_local(cell)) {//cross of faces
    double vals[5][5], val[5], lr[2];
    for (int k = -2; k < 3; k++) 
      for (int j = -2; j < 3; j++)
	vals[k + 2][j + 2] = v.x[k,j];
    colocation_values_ols (vals, &val[0]);
    refine_order_4_ols (val, &lr[0]);
    fine(v.x,1,0,0) = lr[0];
    fine(v.x,1,1,0) = lr[1];
  }
}
#else // (dimension == 3)
foreach_dimension()
void refine_face_4_ols_x (Point point, scalar s) {
  assert (0); //not implemented
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    for (int cx = -1; cx < 2; cx += 2) {
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = - 2; j < 3; j++) {
	  double b[5];
	  for (int k = - 2; k < 3; k++)
	    b[k + 2] = v.x[0, -j*cx, -k*cy];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 0, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
    }
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {//rhs 
    for (int cx = -1; cx < 2; cx += 2) 
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = -2; j < 3; j++) {
	  double b[5];
	  for (int k = -2; k < 3; k++)
	    b[k + 2] = v.x[1, -j*cx, -k*cy];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 2, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
  }
  if (is_local(cell)) { //cross of faces
    double vals[5][5]; //4th-order accurate co-location values
    for (int j = -2; j < 3; j++) 
      for (int k = -2; k < 3; k++) 
	vals[j + 2][k + 2] = (9.*(v.x[0,j,k] + v.x[1,j,k]) - (v.x[-1,j,k] + v.x[2,j,k]))/16.;
    for (int cx = -1; cx < 2; cx += 2) 
      for (int cy = -1; cy < 2; cy += 2) {
	double a[5];
	for (int j = -2; j < 3; j++) {
	  double b[5];
	  for (int k = -2; k < 3; k++)
	    b[k + 2] = vals[-j*cx + 2][-k*cy + 2];
	  a[j + 2] = int_5pt (b);
	}
	fine(v.x, 1, (cx + 1)/2, (cy + 1)/2) = int_5pt (a);
      }
  }
}
#endif
#endif

/**
### Solenodial face refinement

Popinet's (2009) local projection method still applies. The fact that
`p` is determined with second-order accuracy only, does not mean that
the divergece is not *exactly* zero (provided the solvebility
condition), when refining a solenodial field. It can be viewed as the
circulation (from the inital guess) perserving solenoidal field. In
fact, the addional contraint results in an extra order of accuracy!

Vraag niet hoe het kan, maar profiteer ervan.

The following code should set 5th-order solenoidal refinement (in 2D):

~~~literatec
face vector uf[];
int main() {
   uf.x.refine = refine_face_solenoidal;
   foreach_dimension()
      uf.x.prolongation = refine_face_4_x;
 ...
~~~

## Point-value interpolation   

3rd and 5th-Order-accurate point-value interpolation from
cell-averaged (i.e. centered) scalar quantities. `interpolate_3` and
`interpolate_5` are alternatives to `interpolate_linear`. 

Note: These maybe optimized by not using the coefficient matrices
`BBB`.
*/
static inline double interp_3 (double * s, double xp) {
  static double BBB[3][3] = {{ 2,  5, -1},
			     {-3,  3,  0},
			     { 1, -2,  1}};
  double c[3];
  for (int j = 0; j < 3; j++) {
    c[j] = 0;
    for (int i = 0; i < 3; i++)
      c[j] += BBB[j][i]*s[i];
    c[j] = c[j]/6.;
  }
  return (c[0] + 2.*c[1]*xp + 3.*c[2]*sq(xp)); 
}

static inline double interpolate_quadratic (Point point, scalar v,
					    double xp = 0, double yp = 0, double zp = 0)
{
  coord pos = {x, y, z}, p = {xp, yp, zp};
  foreach_dimension()
    pos.x = (p.x - pos.x)/Delta + 0.5;
  double s[3];
  for (int j = -1; j <= 1; j++) {
#if (dimension < 2)
    s[j + 1] = v[j]; 
#else 
    double s1[3];
    for (int k = -1; k <= 1; k++) {
#if (dimension < 3)
      s1[k + 1] = v[j, k];
#else //(dimension == 3)
      double s2[3];
      for (int l = -1; l <= 1; l++) 
	s2[l + 1] = v[j, k, l];
      s1[k + 1] = interp_3 (s2, pos.z);
#endif //dimension < 3 
    }
    s[j + 1] = interp_3 (s1, pos.y);
#endif //dimension > 1
    }
  return interp_3 (s, pos.x);
}

trace 
double interpolate_3 (scalar v, double xp = 0, double yp = 0, double zp = 0)
{
  Point P = locate (xp, yp, zp);
  if (P.level < 0)
    return nodata;
  return interpolate_quadratic (P, v, xp, yp, zp);
}

static double get_val (Point point, scalar s) { return s[]; }

trace
void interpolate_array_3 (scalar * list, coord * a, int n, double * v, bool linear) {
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point P = locate (a[i].x, a[i].y, a[i].z);
    if (P.level >= 0) {
      for (scalar s in list)
	v[j++] = !linear ? get_val (P, s) :
	  interpolate_quadratic (P, s, a[i].x, a[i].y, a[i].z);
    }
    else
      for (scalar s in list)
	v[j++] = nodata;
  }
#if _MPI
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
		MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
		MPI_MIN, 0, MPI_COMM_WORLD);
#endif
}


static inline double interp_5 (double * s, double xp) {
  static double BBB[5][5] = {{-6,  54,  94, -26,  4},
			     { 5, -75,  75, -5,   0},
			     { 5,  10, -40,  30, -5},
			     {-5,  15, -15,  5,   0},
			     { 1, -4,   6,  -4,   1}};
  double c[5];
  for (int j = 0; j < 5; j++) {
    c[j] = 0;
    for (int i = 0; i < 5; i++)
      c[j] += BBB[j][i]*s[i];
    c[j] = c[j]/120.;
  }
  return (c[0] + 2.*c[1]*xp + 3.*c[2]*sq(xp) +
	  4.*c[3]*cube(xp) + 5.*c[4]*sq(xp)*sq(xp)); 
}

static inline double interpolate_quartic (Point point, scalar v,
					  double xp = 0, double yp = 0, double zp = 0)
{
  coord pos = {x, y, z}, p = {xp, yp, zp};
  foreach_dimension()
    pos.x = (p.x - pos.x)/Delta + 0.5;
  double s[5];
  for (int j = -2; j <= 2; j++) {
#if (dimension < 2)
    s[j + 2] = v[j]; 
#else 
    double s1[5];
    for (int k = -2; k <= 2; k++) {
#if (dimension < 3)
      s1[k + 2] = v[j, k];
#else //(dimension == 3)
      double s2[5];
      for (int l = -2; l <= 2; l++) 
	s2[l + 2] = v[j, k, l];
      s1[k + 2] = interp_5 (s2, pos.z);
#endif //dimension < 3 
    }
    s[j + 2] = interp_5 (s1, pos.y);
#endif //dimension > 1
    }
  return interp_5 (s, pos.x);
}

trace 
double interpolate_5 (scalar v, double xp = 0, double yp = 0, double zp = 0)
{
  Point P = locate (xp, yp, zp);
  if (P.level < 0)
    return nodata;
  return interpolate_quartic (P, v, xp, yp, zp);
}

trace
void interpolate_array_5 (scalar * list, coord * a, int n, double * v, bool linear) {
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point P = locate (a[i].x, a[i].y, a[i].z);
    if (P.level >= 0) {
      for (scalar s in list)
	v[j++] = !linear ? get_val (P, s) :
	  interpolate_quartic (P, s, a[i].x, a[i].y, a[i].z);
    }
    else
      for (scalar s in list)
	v[j++] = nodata;
  }
#if _MPI
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
		MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
		MPI_MIN, 0, MPI_COMM_WORLD);
#endif
}


static inline double interp_FD_4 (double * s, double xp) {
  static double BBB[4][4] =   {{ 0. , 6. , 0.,  0.},
			       {-2., -3., 6., -1.},
			       { 3.,  -6., 3.,  0.},
			       {-1.,   3.,-3.,  1.}};
  double c[4];
  for (int j = 0; j < 4; j++) {
    c[j] = 0;
    for (int i = 0; i < 4; i++)
      c[j] += BBB[j][i]*s[i];
    c[j] = c[j]/6.;
  }
  return (c[0] + c[1]*xp + c[2]*sq(xp) + c[3]*cube(xp)); 
}

static inline double interpolate_qubic_vertex (Point point, scalar v,
					       double xp = 0., double yp = 0., double zp = 0.) {
  coord pos = {x, y, z};
  coord p = {xp, yp, zp};
  foreach_dimension()
    pos.x = (p.x - pos.x)/Delta + 0.5;
  double s[4];
  for (int j = -1; j <= 2; j++) {
#if (dimension < 2)
    s[j + 1] = v[j]; 
#else 
    double s1[5];
    for (int k = -1; k <=2; k++) {
#if (dimension < 3)
      s1[k + 1] = v[j, k];
#else //(dimension == 3)
      double s2[5];
      for (int l = -1; l <= 2; l++) 
	s2[l + 1] = v[j, k, l];
      s1[k + 1] = interp_FD_4 (s2, pos.z);
#endif //dimension < 3 
    }
    s[j + 1] = interp_FD_4 (s1, pos.y);
#endif //dimension > 1
    }
  return interp_FD_4 (s, pos.x);
}

double interpolate_vertex_4 (scalar v,
			     double xp = 0, double yp = 0, double zp = 0) {
  Point P = locate (xp, yp, zp);
  if (P.level < 0)
    return nodata;
  return interpolate_qubic_vertex (P, v, xp, yp ,zp);
}

foreach_dimension()
void face_to_vertex_x (scalar u, vertex scalar m) {
  foreach_vertex()
    m[] = (-u[0,-2] + 7*u[0,-1] + 7*u[0] - u[0,1])/12;
}
/**
## Compact schemes

### Face-average to vertex in 2D

Note that vector v should already have suitable prolongation
attributes for vertices.
 */
int compact_iters = 10;
void reduce_average (scalar sa, scalar sp) {
  scalar rhs[];
  foreach() { 
    sp[] = ((-sa[-2] + 7.*(sa[-1] + sa[]) - sa[1])/12.);
    rhs[] =  29./36.*(sa[] + sa[-1]) + 1./36.*(sa[-2] + sa[1]);
  }
  boundary ({sp});
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      sp[] = rhs[] - 1./3.*(sp[-1] + sp[1]);
    }
    boundary ({sp});
  }
}

void cell_to_face (scalar s, face vector f) {
  vector rhs[];
  for (int l = 1; l <= depth(); l++) {
    foreach_level(l) {
      foreach_dimension() {
	f.x[] = ((-s[-2] + 7.*(s[-1] + s[]) - s[1])/12.);
	rhs.x[] =  29./36.*(s[] + s[-1]) + 1./36.*(s[-2] + s[1]);
      }
    }
    boundary_level ((scalar*){f}, l);
    for (int it = 0; it < compact_iters; it++) {
      foreach_level(l) {
	foreach_dimension() {
	    f.x[] = rhs.x[] - 1./3.*(f.x[-1] + f.x[1]);
	}
      }
    }
  }
  boundary ((scalar*){f});
}

void cell_tendency_to_face (vector v, vector f, vector u) {
  vector rhs[];
    for (int l = 1; l <= depth(); l++) {
    foreach_level(l) {
      foreach_dimension() {
	f.x[] = (-v.x[-2] + 7.*(v.x[-1] + v.x[]) - v.x[1])/12.; //4th order guess
	rhs.x[] = u.x[] > 0. ? 251./152.*v.x[-1] + 1./8.*v.x[-2] : 
	                       251./152.*v.x[0]  + 1./8.*v.x[1];
      }
    }
    boundary_level ((scalar*){f}, l);
    for (int it = 0; it < compact_iters; it++) {
      foreach_level(l) {
	foreach_dimension() {
	  double das = u.x[] > 0. ? - 413./456.*f.x[-1] + 23./152.*f.x[1]  - 5./228*f.x[2] :
	    - 413./456.*f.x[1]  + 23./152.*f.x[-1] - 5./228*f.x[-2];
	  f.x[] = (f.x[] + 8.*(rhs.x[] + das))/9.; // 8/9th under relaxation
	}
      }
    }
    }
    boundary ((scalar*){f});
}

void cell_to_face_vector (vector s, face vector f) {
  vector rhs[];
  for (int l = 1; l <= depth(); l++) {
    foreach_level(l) {
      foreach_dimension() {
	f.x[] = ((-s.x[-2] + 7.*(s.x[-1] + s.x[]) - s.x[1])/12.);
	rhs.x[] =  29./36.*(s.x[] + s.x[-1]) + 1./36.*(s.x[-2] + s.x[1]);
      }
    }
    boundary_level ((scalar*){f}, l);
    for (int it = 0; it < 0; it++) {
      foreach_level(l) {
	foreach_dimension() {
	    f.x[] = rhs.x[] - 1./3.*(f.x[-1] + f.x[1]);
	}
      }
    }
  }
  boundary ((scalar*){f});
}

#include "my_vertex.h"
void face_to_vertex (vector f, scalar v) {
  vector vv[], rhs[];
  foreach_dimension() { 
    vv.x.restriction = restriction_vert;
    vv.x.prolongation = refine_vert5;
  }
  boundary ((scalar*){f});
  for (int l = 1; l <= depth(); l++) {
    foreach_level(l) {
      foreach_dimension() {
	vv.x[] = ((-f.x[0,-2] + 7.*(f.x[0,-1] + f.x[]) - f.x[0,1])/12.);
	rhs.x[] =  (29./36.*(f.x[] + f.x[0,-1]) + 1./36.*(f.x[0,-2] + f.x[0,1]));
      }
    }
    boundary_level ((scalar*){vv}, l);
    for (int it = 0; it < compact_iters; it++) {
      foreach_level(l) {
	foreach_dimension()
	  vv.x[] = rhs.x[] - 1./3.*(vv.x[0,-1] + vv.x[0,1]);
      }
    }
  }
  foreach() 
    v[] = ((vv.y[]) + (vv.x[]))/2.;;
  boundary ({v});
}


void compact_face_av_to_vertex (face vector f, vector v) {
  vector rhs[];
  for (int l = 1; l <= depth(); l++) {
    foreach_level(l) { 
      foreach_dimension() {
	v.x[] = ((-f.x[0,-2] + 7*(f.x[0,-1] + f.x[]) - f.x[0,1])/12.);
	rhs.x[] =  29./36*(f.x[] + f.x[0,-1]) + 1./36.*(f.x[0,-2] + f.x[0,1]);
      }
    }
    boundary_level ((scalar*){v}, l);
    for (int it = 0; it < compact_iters; it++) {
      foreach_level(l) {
	foreach_dimension()
	  v.x[] = rhs.x[] - 1./3.*(v.x[0,-1] + v.x[0,1]);
      }
    }
  }
}

void vertex_vector_to_face (vector v, face vector f) {
  face vector rhs[];
  foreach() {
    foreach_dimension() {
      f.x[] =  (-v.x[0,-1] + 13.*(v.x[] + v.x[0,1]) - v.x[0,2])/24.;
      rhs.x[] = 27./38.*(v.x[] + v.x[0,1]) + 3./38.*(v.x[0,-1] + v.x[0,2]) ;
    }
  }
  boundary ((scalar*){f});
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      foreach_dimension()
	f.x[] = rhs.x[] - 11./38.*(f.x[0,-1] + f.x[0,1]);
    }
    boundary ((scalar*){f});
  }
}
/**
### Vertex to face average
 */
void compact_vertex_to_face_av (vector v, face vector f) {
  face vector rhs[];
  foreach_face() { 
    f.x[] =  (-v.x[0,-1] + 13.*(v.x[] + v.x[0,1]) - v.x[0,2])/24.;
    rhs.x[] = 27./38.*(v.x[] + v.x[0,1]) + 3./38.*(v.x[0,-1] + v.x[0,2]) ;
  }
  boundary ((scalar*){f});
  for (int it = 0; it < compact_iters; it++) {
    foreach_face() 
      f.x[] = rhs.x[] - 11./38.*(f.x[0,-1] + f.x[0,1]);
  }
}


/**
### Finite first derivative difference
*/
void compact_first_derivative (scalar * sl, vector * dsl) {
  vector * rhsl;
  rhsl = (vector*)list_clone ((scalar*)dsl);
  double alphaP = 1./4., aP = 1.5;
  foreach() {
    scalar s;
    vector ds, rhs;
    for (s, ds, rhs in sl, dsl, rhsl) {
      foreach_dimension() {
	ds.x[] = ((s[1] - s[-1])/(2*Delta));
	rhs.x[] = aP*(s[1] - s[-1])/(2.*Delta);
      }
    }
  }
  boundary((scalar*)dsl);
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      vector ds, rhs;
      for  (ds, rhs in dsl, rhsl) {
	foreach_dimension()	
	  ds.x[] = rhs.x[] - alphaP*(ds.x[-1] + ds.x[1]);
      }
    }
    boundary ((scalar*)dsl);
  }
  delete((scalar*)rhsl);
  free (rhsl); 
  rhsl = NULL;
}
/**
~~~gnuplot Modified wavenumbers for upwind schemes
set size square
set grid
set xr [0 : 3.14]
set yr [0 : 3.14]
a = 165./197.
b =  -192./197.
c =  27./197.
alpha = -76./197
beta =  17./197.
i = {0.0, 1.0}
f4(x) = (c*exp(-i*2*x) + b*exp(-i*x) + a)/(1 + beta*exp(i*2.*x) + alpha*exp(i*x))

alpha6 = 43./123.
beta6 = 35./123.
gamma6 = -1./123.;
a6 = -11./369.
b6 = -33/41.
c6 = 3./41.
d6 = 281./369;
f6(x) = (a6*exp(-i*2*x) + b6*exp(-i*x) + c6 + d6*exp(i*x))/(1 + alpha6*exp(-i*x) + beta6*exp(i*x) + gamma6*exp(2*i*x))

plot exp(-real(f4(x))), imag(f4(x)), exp(-real(f6(x))), imag(f6(x)), x t 'exact'
~~~
*/
void compact_upwind (scalar * sl, vector * dsl, vector v) {
  double alpha = -76./197, beta = 17/197.;
  double a = 165./197., b = -192/197., c = 27./197.;
  vector * rhsl;
  rhsl = (vector*)list_clone ((scalar*)dsl);
  foreach() {
    scalar s;
    vector ds, rhs;
    for  (s, ds, rhs in sl, dsl, rhsl) {
      foreach_dimension() {
	rhs.x[] = (v.x[] > 0 ? (a*s[] + b*s[-1] + c*s[-2])/Delta
		   : -(a*s[] + b*s[1] + c*s[2])/Delta);
	ds.x[] = (8*(s[1] - s[-1]) + s[-2] - s[2])/(12*Delta);
	//ds.x[] = 0.;
      }
    }
  }
  boundary ((scalar*)dsl);
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      vector ds, rhs;
      for  (ds, rhs in dsl, rhsl) {
	foreach_dimension() {
	  ds.x[] =  (v.x[] > 0 ? rhs.x[] - alpha*ds.x[1] - beta*ds.x[2] :
		     rhs.x[] - alpha*ds.x[-1] - beta*ds.x[-2]);
	}
      }
    }
    boundary ((scalar*)dsl);
  }
  delete((scalar*)rhsl);
  free (rhsl); 
  rhsl = NULL;
}

void compact_upwind5 (scalar * sl, vector * dsl, vector v) {
  double alpha = 62./111., beta = -1./37.;
  double a = 10./333., b = -13./37., c = -34./37., d = 413./333.;
  vector * rhsl;
  rhsl = (vector*)list_clone ((scalar*)dsl);
  foreach() {
    scalar s;
    vector ds, rhs;
    for  (s, ds, rhs in sl, dsl, rhsl) {
      foreach_dimension() {
	rhs.x[] = (v.x[] > 0 ? (a*s[-2] + b*s[-1] + c*s[] + d*s[1])/Delta
		   : -(a*s[2] + b*s[1] + c*s[] + d*s[-1])/Delta);
	ds.x[] = (8*(s[1] - s[-1]) + s[-2] - s[2])/(12*Delta);
	//ds.x[] = 0.;
      }
    }
  }
  boundary ((scalar*)dsl);
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      vector ds, rhs;
      for  (ds, rhs in dsl, rhsl) {
	foreach_dimension() {
	  double n = (v.x[] > 0 ? rhs.x[] - alpha*ds.x[1] - beta*ds.x[2]:
		     rhs.x[] - alpha*ds.x[-1] - beta*ds.x[-2]);
	  ds.x[] = (ds.x[] + 2*n)/3.;
	}
      }
    }
    boundary ((scalar*)dsl);
  }
  delete((scalar*)rhsl); free (rhsl); rhsl = NULL;
}

void compact_upwind6 (scalar * sl, vector * dsl, vector v) {
  double alpha = 43./123., beta = 35./123., gamma = -1./123.;
  double a = -11./369., b = -33/41., c = 3./41., d = 281./369;
  vector * rhsl;
  rhsl = (vector*)list_clone ((scalar*)dsl);
  foreach() {
    scalar s;
    vector ds, rhs;
    for  (s, ds, rhs in sl, dsl, rhsl) {
      foreach_dimension() {
	rhs.x[] = (v.x[] > 0 ? (a*s[-2] + b*s[-1] + c*s[] + d*s[1])/Delta
		   : -(a*s[2] + b*s[1] + c*s[] + d*s[-1])/Delta);
	ds.x[] = (8*(s[1] - s[-1]) + s[-2] - s[2])/(12*Delta);
	//ds.x[] = 0.;
      }
    }
  }
  boundary ((scalar*)dsl);
  for (int it = 0; it < compact_iters; it++) {
    foreach() {
      vector ds, rhs;
      for  (ds, rhs in dsl, rhsl) {
	foreach_dimension() {
	  ds.x[] =  (v.x[] > 0 ? rhs.x[] - alpha*ds.x[-1] - beta*ds.x[1] - gamma*ds.x[2]:
		     rhs.x[] - alpha*ds.x[1] - beta*ds.x[-1] - gamma*ds.x[-2]);
	}
      }
    }
    boundary ((scalar*)dsl);
  }
  delete((scalar*)rhsl);
  free (rhsl); 
  rhsl = NULL;
}

void compact_first_derivative6 (scalar * sl, vector * dsl) {
  vector * rhsl;
  rhsl = (vector*)list_clone ((scalar*)dsl);
  double alphaP = 1./3., aP = 14./9., bP = 1./9.;
  for (int l = 1; l <= depth(); l++) {
    foreach_level(l) {
      scalar s;
      vector ds, rhs;
      for  (s, ds, rhs in sl, dsl, rhsl) {
	foreach_dimension() {
	  ds.x[] = ((8*(s[1] - s[-1]) + s[-2] - s[2])/(12*Delta));
	  rhs.x[] = aP*(s[1] - s[-1])/(2.*Delta) + bP*(s[2] - s[-2])/(4.*Delta);
	}
      }
    }
    boundary_level ((scalar*)dsl, l);
    for (int it = 0; it < compact_iters; it++) {
      foreach_level(l) {
	vector ds, rhs;
	for  (ds, rhs in dsl, rhsl) {
	  foreach_dimension()	
	    ds.x[] = rhs.x[] - alphaP*(ds.x[-1] + ds.x[1]);
	}
      }
    }
  }
  delete((scalar*)rhsl);
  free (rhsl); 
  rhsl = NULL;
}

@define n_x neighbor.i
@define n_y neighbor.j
@define n_z neighbor.k


@define layer_nr_x (n_x < GHOSTS ? (GHOSTS - n_x) : n_x - (1 << level) - GHOSTS + 1)
@define layer_nr_y (n_y < GHOSTS ? (GHOSTS - n_y) : n_y - (1 << level) - GHOSTS + 1)
@define layer_nr_z (n_z < GHOSTS ? (GHOSTS - n_z) : n_z - (1 << level) - GHOSTS + 1)
    
@define dirichlet_x(a) (val(_s,0,0,0) + layer_nr_x*2.*(a - val(_s,0,0,0))) //2nd order
@define dirichlet_y(a) (val(_s,0,0,0) + layer_nr_y*2.*(a - val(_s,0,0,0))) //2nd order
@define dirichlet_z(a) (val(_s,0,0,0) + layer_nr_z*2.*(a - val(_s,0,0,0))) //2nd order

@define neumann_x(a)   (val(_s,0,0,0) + layer_nr_x*Delta*a)
@define neumann_y(a)   (val(_s,0,0,0) + layer_nr_y*Delta*a)
@define neumann_z(a)   (val(_s,0,0,0) + layer_nr_z*Delta*a)
   
#define neumann_x_homogeneous (val(_s,0,0,0))
#define neumann_y_homogeneous (val(_s,0,0,0))
#define neumann_z_homogeneous (val(_s,0,0,0))

@define dirichlet_vert_left(a) (layer_nr_x == 1 ? (3*a - 3*val(_s,1,0,0) + val(_s,2,0,0)) : (6*val(_s,0,0,0) - 8*val(_s,1,0,0) + 3*val(_s,2,0,0)))
@define dirichlet_vert_right(a) (layer_nr_x == 1 ? (a) : 3.*((a) - val(_s,0,0,0)) + val(_s,-1,0,0))
@define dirichlet_vert_bottom(a) (layer_nr_y == 1 ? (3*a - 3*val(_s,0,1,0) + val(_s,0,2,0)) :(6*val(_s,0,0,0) - 8*val(_s,0,1,0) + 3*val(_s,0,2,0)))
@define dirichlet_vert_top(a) (layer_nr_y == 1 ? (a) : 3.*((a) - val(_s,0,0,0)) + val(_s,0,-1,0)) 


/**
## 3rd order accurate boundaries
 */
#define dirichlet_left(a) (layer_nr_y == 1 ? (3.*(a) -  5.*val(_s,0,0,0)/2. +    val(_s,1,0,0)/2.) : \
			   (9.*(a) - 21.*val(_s,0,0,0)/2. + 5.*val(_s,1,0,0)/2.))
#define dirichlet_right(a) (layer_nr_y == 1 ? (3.*(a) -  5.*val(_s,0,0,0)/2. +    val(_s,-1,0,0)/2.) : \
			    (9.*(a) - 21.*val(_s,0,0,0)/2. + 5.*val(_s,-1,0,0)/2.))
#define dirichlet_bottom(a) (layer_nr_y == 1 ? (3.*(a) -  5.*val(_s,0,0,0)/2. + val(_s,0,1,0)/2.) : \
			     (9.*(a) - 21.*val(_s,0,0,0)/2. + 5.*val(_s,0,1,0)/2.))
#define dirichlet_top(a) (layer_nr_y == 1 ? (3.*(a) -  5.*val(_s,0,0,0)/2. +    val(_s,0,-1,0)/2.) : \
			  (9.*(a) - 21.*val(_s,0,0,0)/2. + 5.*val(_s,0,-1,0)/2.))
#define dirichlet_back(a) (layer_nr_y == 1 ? (3.*(a) -  5.*val(_s,0,0,0)/2. +    val(_s,0,0,1)/2.) : \
			   (9.*(a) - 21.*val(_s,0,0,0)/2. + 5.*val(_s,0,0,1)/2.))
#define dirichlet_front(a) (layer_nr_y == 1 ? (3.*(a) -  5.*val(_s,0,0,0)/2. +    val(_s,0,0,-1)/2.) : \
			    (9.*(a) - 21.*val(_s,0,0,0)/2. + 5.*val(_s,0,0,-1)/2.))

#define neumann_left(a) (layer_nr_y == 1 ? val(_s,0,0,0) - (a)*Delta : val(_s,1,0,0) - 3.*(a)*Delta);
#define neumann_right(a) (layer_nr_y == 1 ? val(_s,0,0,0) - (a)*Delta : val(_s,-1,0,0) - 3.*(a)*Delta);
#define neumann_bottom(a) (layer_nr_y == 1 ? val(_s,0,0,0) - (a)*Delta : val(_s,0,1,0) - 3.*(a)*Delta);
#define neumann_top(a) (layer_nr_y == 1 ? val(_s,0,0,0) - (a)*Delta : val(_s,0,-1,0) - 3.*(a)*Delta);
#define neumann_back(a) (layer_nr_y == 1 ? val(_s,0,0,0) - (a)*Delta : val(_s,0,0,1) - 3.*(a)*Delta);
#define neumann_front(a) (layer_nr_y == 1 ? val(_s,0,0,0) - (a)*Delta : val(_s,0,0,-1) - 3.*(a)*Delta);

#define symmetric_top (layer_nr_y == 1 ? val(_s,0,0,0) : val(_s,0,-1,0)) 
#define symmetric_bottom (layer_nr_y == 1 ? val(_s,0,0,0) : val(_s,0,1,0)) 

#define dirichlet_left_homogeneous (dirichlet_left(0))
#define dirichlet_right_homogeneous (dirichlet_right(0))
#define dirichlet_bottom_homogeneous (dirichlet_bottom(0))
#define dirichlet_top_homogeneous (dirichlet_top(0))
#define dirichlet_back_homogeneous (dirichlet_back(0))
#define dirichlet_front_homogeneous (dirichlet_front(0))
#define neumann_left_homogeneous (neumann_left(0))
#define neumann_right_homogeneous (neumann_right(0))
#define neumann_bottom_homogeneous (neumann_bottom(0))
#define neumann_top_homogeneous (neumann_top(0))
#define neumann_back_homogeneous (neumann_back(0))
#define neumann_front_homogeneous (neumann_front(0)) 


/**
## 4th order accurate dirichlet boundaries

using 3 or even 4 points!
 */
#define dirichlet_left4(a) (layer_nr_x == 1 ? (4.*(a) + (-13.*val(_s,0,0,0) +  5*val(_s,1,0,0) - val(_s,2,0,0))/3.) : \
			    (16.*(a) + (-70.*val(_s,0,0,0) + 32.*val(_s,1,0,0) - 7*val(_s,2,0,0))/3.))
#define dirichlet_right4(a) (layer_nr_x == 1 ? (4.*(a) + (-13.*val(_s,0,0,0) +  5*val(_s,-1,0,0) - val(_s,-2,0,0))/3.) : \
			    (16.*(a) + (-70.*val(_s,0,0,0) + 32.*val(_s,-1,0,0) - 7*val(_s,-2,0,0))/3.))
#define dirichlet_bottom4(a) (layer_nr_y == 1 ? (4.*(a) + (-13.*val(_s,0,0,0) +  5*val(_s,0,1,0) - val(_s,0,2,0))/3.) : \
			    (16.*(a) + (-70.*val(_s,0,0,0) + 32.*val(_s,0,1,0) - 7*val(_s,0,2,0))/3.))
#define dirichlet_top4(a) (layer_nr_y == 1 ? (4.*(a) + (-13.*val(_s,0,0,0) +  5.*val(_s,0,-1,0) - val(_s,0,-2,0))/3.) : \
			   (16.*(a) + (-70.*val(_s,0,0,0) + 32.*val(_s,0,-1,0) - 7*val(_s,0,-2,0))/3.))
#define dirichlet_back4(a) (layer_nr_z == 1 ? (4.*(a) + (-13.*val(_s,0,0,0) +  5*val(_s,0,0,1) - val(_s,0,0,2))/3.) : \
			    (16.*(a) + (-70.*val(_s,0,0,0) + 32.*val(_s,0,1,0) - 7*val(_s,0,2,0))/3.))
#define dirichlet_front4(a) (layer_nr_z == 1 ? (4.*(a) + (-13.*val(_s,0,0,0) +  5*val(_s,0,0,-1) - val(_s,0,0,-2))/3.) : \
			     (16.*(a) + (-70.*val(_s,0,0,0) + 32.*val(_s,0,0,-1) - 7*val(_s,0,-,-2))/3.))
@define dirichlet_vert_left4(a)   (level > 0 ? layer_nr_x == 1 ? (4.*(a) - 6.*val(_s,1,0,0) + 4.*val(_s,2,0,0) - val(_s,3,0,0))   : (10.*(a) - 20.*val(_s,1,0,0) + 15.*val(_s,2,0,0) - 4*val(_s,3,0,0)) : 0);
@define dirichlet_vert_right4(a)  (level > 0 ? layer_nr_x == 1 ? (a) : (4.*(a) - 6.*val(_s,0,0,0) + 4.*val(_s,-1,0,0) - val(_s,-2,0,0)) : 0)
@define dirichlet_vert_bottom4(a) (level > 0 ? layer_nr_y == 1 ? (4.*(a) - 6.*val(_s,0,1,0) + 4.*val(_s,0,2,0) - val(_s,0,3,0))   : (10.*(a) - 20.*val(_s,0,1,0) + 15.*val(_s,0,2,0) - 4*val(_s,0,3,0)) : 0)
@define dirichlet_vert_top4(a)    (level > 0 ? layer_nr_y == 1 ? (a) : (4.*(a) - 6.*val(_s,0,0,0) + 4.*val(_s,0,-1,0) - val(_s,0,-2,0)) : 0)

#if dimension == 2
double Gauss6_x (double x, double y, double Delta, double (* myfun)(double x, double y)) {
  double w1 = 4./9., w2 = 5./18.;
  double yw = sqrt(3./5.)/2.*Delta;
  return w1*myfun (x, y) + w2*(myfun (x, y - yw) + myfun (x, y + yw));
}

double Gauss6_y (double x, double y, double Delta, double (* myfun)(double x, double y)) {
  double w1 = 4./9., w2 = 5./18.;
  double xw = sqrt(3./5.)/2.*Delta;
  return w1*myfun (x, y) + w2*(myfun (x - xw, y ) + myfun (x + xw, y));
}

double Gauss6 (double x, double y, double Delta, double (* myfun)(double x, double y)) { //cell integral
   double w1 = 25./324., w2 = 10./81., w3 = 16./81.;
  double d = sqrt(3./5.)/2.*Delta;
  return (w3*myfun (x, y) +
	  w1*(myfun (x - d, y - d) + myfun (x - d, y + d) +
	      myfun (x + d, y + d) + myfun (x + d, y - d)) +
	  w2*(myfun (x - d, y) + myfun (x + d, y) +
	      myfun (x, y - d) + myfun (x, y + d)));
}

#else
double Gauss6_x (double x, double y, double z, double Delta, double (* myfun)(double x, double y, double z)) {
  double w1 = 25./324., w2 = 10./81., w3 = 16./81.;
  double d = sqrt(3./5.)/2.*Delta;
  return (w3*myfun (x, y, z) +
	  w1*(myfun (x, y - d, z - d) + myfun (x, y + d, z - d) +
	      myfun (x, y + d, z + d) + myfun (x, y - d, z + d)) +
	  w2*(myfun (x, y, z - d) + myfun (x, y, z + d) +
	      myfun (x, y - d, z) + myfun (x, y + d, z)));
}

double Gauss6_y (double x, double y, double z, double Delta, double (* myfun)(double x, double y, double z)) {
  double w1 = 25./324., w2 = 10./81., w3 = 16./81.;
  double d = sqrt(3./5.)/2.*Delta;
  return (w3*myfun (x, y, z) +
	  w1*(myfun (x - d, y, z - d) + myfun (x + d, y, z - d) +
	      myfun (x + d, y, z + d) + myfun (x - d, y, z + d)) +
	  w2*(myfun (x, y, z - d) + myfun (x, y, z + d) +
	      myfun (x - d, y, z) + myfun (x + d, y, z)));

}

double Gauss6_z (double x, double y, double z, double Delta, double (* myfun)(double x, double y, double z)) {
  double w1 = 25./324., w2 = 10./81., w3 = 16./81.;
  double d = sqrt(3./5.)/2.*Delta;
  return (w3*myfun (x, y, z) +
	  w1*(myfun (x - d, y - d, z) + myfun (x - d, y + d, z) +
	      myfun (x + d, y + d, z) + myfun (x + d, y - d, z)) +
	  w2*(myfun (x - d, y, z) + myfun (x + d, y, z) +
	      myfun (x, y - d, z) + myfun (x, y + d, z)));
}
#endif

#if (dimension == 2)
#define VERTEX_TO_FACE_4(s) ((-s[0,-1] + 13.*(s[] + s[0,1]) - s[0,2])/24.)
#define FACE_TO_VERTEX_4(s) ((-s[0,-2] + 7.*(s[0,-1] + s[]) - s[0,1])/12.)
#else 
#define VERTEX_TO_FACE_4(s) ((1.*  (s[0,-1,-1] + s[0,-1,2] + s[0,2,2]  + s[0,2,-1]) + \
			      -13.*(s[0,-1,0]  + s[0,-1,1] + s[0,2,0]  + s[0,2,1]   + \
				    s[0,0,2]   + s[0,1,2]  + s[0,1,-1] + s[0,0,-1]) + \
			      169.*(s[0,0,1]   + s[0,1,0]  + s[0,1,1]  + s[]))/576.)
#define FACE_TO_VERTEX_4(s) ((1.*(s[0,-2,-2] + s[0,-2,1]  + s[0,1,1]   + s[0,1,-2]) + \
			      -7.*(s[0,-2,0]  + s[0,-2,-1] + s[0,1,0]   + s[0,1,-1]  + \
				   s[0,-1,1]  + s[0,0,1]   + s[0,-1,-2] + s[0,0,-2]) + \
			      49.*(s[0,-1,-1] + s[0,-1,0]  + s[0,0,0]   + s[0,0,-1]))/144.)
#endif
/**
## Tests

* [2D scalar refinement](refine2D.c)  
* [3D scalar refinement](refine3D.c)  
* [2D face vector refinement](rf4t.c)  
* [Interpolation in 1D](thi.c)  
* [Interpolation in 3D](thi3d.c)  
 */
