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
  return     (-3*a[0] + 17*a[1] + 55*a[2] - 5*a[3])/64.;
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
/**
   We also define a 4-th order staggering interpolator. Doing 5
   interpolations from a 4x5-point stencil.
 */

static inline void colocation_values (double val[4][5], double * intrp) {
  for (int j = 0; j < 5; j++)
      intrp[j] = (9.*(val[2][j] + val[1][j]) - (val[0][j] + val[3][j]))/16.;
  //!=  (7.*(val[2][j] + val[1][j]) - (val[0][j] + val[3][j]))/12.;
}
/**
   These two functions can be used to construct a face-refinement
   method
 */
#if (QUADTREE)
foreach_dimension()
static void refine_face_4_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {//lhs
    double val[5], lr[2];
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
#endif
/**
### Solenodial refinement

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

### To do

Three dimensions.

## `#define`-itions

Some more higher-order definitions: 
*/

foreach_dimension() {
#define face_value_4_x(s,i)      (7.*(s[i] + s[i - 1]) - (s[i - 2] + s[i + 1]))/12.
#define face_der_4_x(s,i)        (15.*(s[i] - s[i - 1]) - (s[i + 1] - s[i - 2]))/(12.*Delta)
#define face_2ndder_x(s,i) =     (s[i - 2] - s[i - 1] - s[i] + s[i + 1])/(2.*sq(Delta))
#define centered_der_4_x(s,i)    ((face_value_4_x(s,i + 1) - face_value_4_x(s,i))/Delta) 
#define centered_2ndder_4_x(s,i) ((face_der_4_x(s,i + 1) - face_der_4_x(s,i))/Delta) 
#define centered_3rdder_x(s,i)   ((face_2ndder_x(s,i + 1) - face_2ndder_x(s,i))/Delta)
}

/**
## Point-value Interpolation   

A 5th-order-accurate point-value interpolation from cell-averaged
(i.e. centered) scalar quantities. 
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


static inline double interpolate_quadratic (Point point, struct _interpolate p) {
  scalar v = p.v;
  coord pos = {x, y, z};
  foreach_dimension()
    pos.x = (p.x - pos.x)/Delta + 0.5;
  double s[3];
  for (int j = -1; j <= 1; j++) {
#if (dimension < 2)
    s[j + 2] = v[j]; 
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
double interpolate_3 (struct _interpolate p)   {
  Point point = locate (p.x, p.y, p.z);
  if (point.level < 0)
    return nodata;
  return interpolate_quadratic (point, p);
}

trace
void interpolate_array_3 (scalar * list, coord * a, int n, double * v, bool linear) {
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate (a[i].x, a[i].y, a[i].z);
    if (point.level >= 0) {
      for (scalar s in list)
	v[j++] = !linear ? s[] :
	  interpolate_quadratic (point,
				 (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
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

static inline double interpolate_quartic (Point point, struct _interpolate p) {
  scalar v = p.v;
  coord pos = {x, y, z};
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
double interpolate_5 (struct _interpolate p)   {
  Point point = locate (p.x, p.y, p.z);
  if (point.level < 0)
    return nodata;
  return interpolate_quartic (point, p);
}

trace
void interpolate_array_5 (scalar * list, coord * a, int n, double * v, bool linear) {
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate (a[i].x, a[i].y, a[i].z);
    if (point.level >= 0) {
      for (scalar s in list)
	v[j++] = !linear ? s[] :
	  interpolate_quartic (point,
			       (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
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

static inline double interpolate_qubic_vertex (Point point, struct _interpolate p) {
  scalar v = p.v;
  coord pos = {x, y, z};
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

double interpolate_vertex_4 (struct _interpolate p) {
  Point point = locate (p.x, p.y, p.z);
  if (point.level < 0)
    return nodata;
  return interpolate_qubic_vertex (point, p);
}

foreach_dimension()
void face_to_vertex_x (scalar u, vertex scalar m) {
  foreach_vertex()
    m[] = (-u[0,-2] + 7*u[0,-1] + 7*u[0] - u[0,1])/12;
}

/**
## Tests

* [2D scalar refinement](refine2D.c)  
* [3D scalar refinement](refine3D.c)  
* [2D face vector refinement](rf4t.c)  
* [Interpolation in 1D](thi.c)  
* [Interpolation in 3D](thi3d.c)  
 */
