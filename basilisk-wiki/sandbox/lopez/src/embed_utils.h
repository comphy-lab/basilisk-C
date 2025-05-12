/**
# Embedded boundary utilities
*/

#if 1
foreach_dimension()
double embed_face_avg_gradient_t1_x (Point point, scalar a, int i) {
  double up = nodata, down = nodata;
  double xp = (1.-fs.y[0,i])/2.*sign(fs.y[1,i]-fs.y[-1,i]);
  if (cs[0,i])
    up = (fs.x[0,i] && fs.x[1,i] ? ((xp+1.)*a[1,i] - 2.*xp*a[0,i] + (xp -1.)*a[-1,i])/(2.*Delta) :
	  fs.x[1,i] && fs.x[2,i] ? ((xp-1.)*a[2,i]  -2.*(xp-2.)*a[1,i] + (xp-3.)*a[0,i])/(2.*Delta) :
	  fs.x[0,i] && fs.x[-1,i] ? ((xp+1.)*a[-2,i] -2.*(xp+2.)*a[-1,i] + (xp+3.)*a[0,i])/(2.*Delta) :
	  fs.x[1,i] ? (a[1,i] - a[0,i])/Delta :
	  fs.x[0,i]  ? (a[0,i] - a[-1,i])/Delta : nodata);
  if (cs[0,i-1])
    down = (fs.x[0,i-1] && fs.x[1,i-1] ? ((xp+1.)*a[1,i-1] - 2.*xp*a[0,i-1] + (xp -1.)*a[-1,i-1])/(2.*Delta) :
	  fs.x[1,i-1] && fs.x[2,i-1] ? ((xp-1.)*a[2,i-1]  -2.*(xp-2.)*a[1,i-1] + (xp-3.)*a[0,i-1])/(2.*Delta) :
	  fs.x[0,i-1] && fs.x[-1,i-1] ? ((xp+1.)*a[-2,i-1] -2.*(xp+2.)*a[-1,i-1] + (xp+3.)*a[0,i-1])/(2.*Delta) :
	    fs.x[1,i-1] ? (a[1,i-1] - a[0,i-1])/Delta :
	    fs.x[0,i-1] ? (a[0,i-1] - a[-1,i-1])/Delta : nodata);
 
  return (up == nodata && down == nodata ? 0. :
	  up == nodata ? down :
	  down == nodata ? up :
	  fs.y[0,i] ? (down + up)/2. : 0.);
}

foreach_dimension()
double embed_face_avg_gradient_t2_x (Point point, scalar a, int i) {
  double up = nodata, down = nodata;
//  double xp = 0.;//(1.-fs.y[0,i])/2.*sign(fs.y[1,i]-fs.y[-1,i]);
  if (cs[0,0,i])
    up = (fs.x[0,0,i] && fs.x[1,0,i] ? (a[1,0,i] - a[-1,0,i])/(2.*Delta) :
	  fs.x[1,0,i] && fs.x[2,0,i] ? (-a[2,0,i] + 4.*a[1,0,i] - 3.*a[0,0,i])/(2.*Delta) :
	  fs.x[0,0,i] && fs.x[-1,0,i] ? (a[-2,0,i] - 4.*a[-1,0,i] + 3.*a[0,0,i])/(2.*Delta) :
	  fs.x[1,0,i] ? (a[1,0,i] - a[0,0,i])/Delta :
	  fs.x[0,0,i]  ? (a[0,0,i] - a[-1,0,i])/Delta : 0.);
  if (cs[0,0,i-1])
    down = (fs.x[0,0,i-1] && fs.x[1,0,i-1] ? (a[1,0,i-1] - a[-1,0,i-1])/(2.*Delta) :
	    fs.x[1,0,i-1] && fs.x[2,0,i-1] ? (-a[2,0,i-1] + 4.*a[1,0,i-1] - 3.*a[0,0,i-1])/(2.*Delta) :
	    fs.x[0,0,i-1] && fs.x[-1,0,i-1] ? (a[-2,0,i-1] - 4.*a[-1,0,i-1] + 3.*a[0,0,i-1])/(2.*Delta) :
	    fs.x[1,0,i-1] ? (a[1,0,i-1] - a[0,0,i-1])/Delta :
	    fs.x[0,0,i-1] ? (a[0,0,i-1] - a[-1,0,i-1])/Delta : 0.);
  return (up == nodata && down == nodata ? 0. :
	  up == nodata ? down :
	  down == nodata ? up :
	  fs.z[0,0,i] ? (down + up)/2. : 0.);
}
#endif

#if 0
foreach_dimension()
double embed_face_avg_gradient_t1_x (Point point, scalar a, int i) {
  double up = (fs.x[0,i] && fs.x[1,i] ? (a[1,i] - a[-1,i])/(2.*Delta) :
	       fs.x[1,i] && fs.x[2,i] ? (-a[2,i] + 4.*a[1,i] - 3.*a[0,i])/(2.*Delta) :
	       fs.x[0,i] && fs.x[-1,i] ? (a[-2,i] - 4.*a[-1,i] + 3.*a[0,i])/(2.*Delta) :
	       fs.x[1,i] ? (a[1,i] - a[0,i])/Delta :
	       fs.x[0,i]  ? (a[0,i] - a[-1,i])/Delta : 0.);
  double down = (fs.x[0,i-1] && fs.x[1,i-1] ? (a[1,i-1] - a[-1,i-1])/(2.*Delta) :
	       fs.x[1,i-1] && fs.x[2,i-1] ? (-a[2,i-1] + 4.*a[1,i-1] - 3.*a[0,i-1])/(2.*Delta) :
	       fs.x[0,i-1] && fs.x[-1,i-1] ? (a[-2,i-1] - 4.*a[-1,i-1] + 3.*a[0,i-1])/(2.*Delta) :
	       fs.x[1,i-1] ? (a[1,i-1] - a[0,i-1])/Delta :
	       fs.x[0,i-1] ? (a[0,i-1] - a[-1,i-1])/Delta : 0.);
  return (!up ? down : (!down ? up : (down + up)/2.));			       
}

foreach_dimension()
double embed_face_avg_gradient_t2_x (Point point, scalar a, int i) {
  double up = (fs.x[0,0,i] && fs.x[1,0,i] ? (a[1,0,i] - a[-1,0,i])/(2.*Delta) :
	       fs.x[1,0,i] && fs.x[2,0,i] ? (-a[2,0,i] + 4.*a[1,0,i] - 3.*a[0,0,i])/(2.*Delta) :
	       fs.x[0,0,i] && fs.x[-1,0,i] ? (a[-2,0,i] - 4.*a[-1,0,i] + 3.*a[0,0,i])/(2.*Delta) :
	       fs.x[1,0,i] ? (a[1,0,i] - a[0,0,i])/Delta :
	       fs.x[0,0,i]  ? (a[0,0,i] - a[-1,0,i])/Delta : 0.);
  double down = (fs.x[0,0,i-1] && fs.x[1,0,i-1] ? (a[1,0,i-1] - a[-1,0,i-1])/(2.*Delta) :
	       fs.x[1,0,i-1] && fs.x[2,0,i-1] ? (-a[2,0,i-1] + 4.*a[1,0,i-1] - 3.*a[0,0,i-1])/(2.*Delta) :
	       fs.x[0,0,i-1] && fs.x[-1,0,i-1] ? (a[-2,0,i-1] - 4.*a[-1,0,i-1] + 3.*a[0,0,i-1])/(2.*Delta) :
		 fs.x[1,0,i-1] ? (a[1,0,i-1] - a[0,0,i-1])/Delta :
		 fs.x[0,0,i-1] ? (a[0,0,i-1] - a[-1,0,i-1])/Delta : 0.);
  return (!up ? down : (!down ? up : (down + up)/2.));			       
}
#endif

#if 0
foreach_dimension()
double embed_face_avg_gradient_t1_x (Point point, scalar a, int i) {
  return (fs.x[0,i] && fs.x[1,i] && fs.x[0,i-1] && fs.x[1,i-1] ?
  	  (a[1,i] + a[1,i-1] - a[-1,i-1] - a[-1,i])/(4.*Delta) : //6 cells availables
  	  fs.x[1,i]  && fs.x[1,i-1] ?
	  (a[1,i] + a[1,i-1] - a[0,i] - a[0,i-1] )/(2.*Delta) : //4 cells availables
  	  fs.x[0,i]  && fs.x[0,i-1] ?
	  (a[0,i] + a[0,i-1] - a[-1,i] - a[-1,i-1] )/(2.*Delta) : //4 cells availables
  	  fs.x[0,i] ?  (a[0,i] - a[-1,i])/Delta :
	  fs.x[0,i-1] ?  (a[0,i-1] - a[-1,i-1])/Delta : // 1 cell available
  	  fs.x[1,i] ?  (a[1,i] - a[0,i])/Delta :
	  fs.x[1,i-1] ?  (a[1,i-1] - a[0,i-1])/Delta : 0.); // 1 cell available
}

foreach_dimension()
double embed_face_avg_gradient_t2_x (Point point, scalar a, int i) {
  return (fs.x[0,0,i] && fs.x[1,0,i] && fs.x[0,0,i-1] && fs.x[1,0,i-1] ?
  	  (a[1,0,i] + a[1,0,i-1] - a[-1,0,i-1] - a[-1,0,i])/(4.*Delta) : //6 cells availables
  	  fs.x[1,0,i]  && fs.x[1,0,i-1] ?
	  (a[1,0,i] + a[1,0,i-1] - a[0,0,i] - a[0,0,i-1] )/(2.*Delta) : //4 cells availables
  	  fs.x[0,0,i]  && fs.x[0,0,i-1] ?
	  (a[0,0,i] + a[0,0,i-1] - a[-1,0,i] - a[-1,0,i-1] )/(2.*Delta) : //4 cells availables
  	  fs.x[0,0,i] ?  (a[0,0,i] - a[-1,0,i])/Delta :
	  fs.x[0,0,i-1] ?  (a[0,0,i-1] - a[-1,0,i-1])/Delta : // 1 cell available
  	  fs.x[1,0,i] ?  (a[1,0,i] - a[0,0,i])/Delta :
	  fs.x[1,0,i-1] ?  (a[1,0,i-1] - a[0,0,i-1])/Delta : 0.); // 1 cell available
}
#endif

#undef face_avg_gradient_t1_x
#define face_avg_gradient_t1_x(a,i) embed_face_avg_gradient_t1_x (point, a, i)
#undef face_avg_gradient_t1_y
#define face_avg_gradient_t1_y(a,i) embed_face_avg_gradient_t1_y (point, a, i)
#undef face_avg_gradient_t1_z
#define face_avg_gradient_t1_z(a,i) embed_face_avg_gradient_t1_z (point, a, i)

#undef face_avg_gradient_t2_x
#define face_avg_gradient_t2_x(a,i) embed_face_avg_gradient_t2_x (point, a, i)
#undef face_avg_gradient_t2_y
#define face_avg_gradient_t2_y(a,i) embed_face_avg_gradient_t2_y (point, a, i)
#undef face_avg_gradient_t2_z
#define face_avg_gradient_t2_z(a,i) embed_face_avg_gradient_t2_z (point, a, i)

/**
The function below computes the gradient of an scalar *s* on the
position *o* of the barycenter of a embed fragment. The gradient is
calculated from a bilinear expression 
$$
s(x,y) = a_o + a_1 x + a_2 y + a_3 xy
$$
or a trilinear one in 3D,
$$
s(x,y,z) = a_o + a_1 x + a_2 y + a_3 z + a_4 xy + a_5 xz + a_6 yz
$$
whose coefficients are calculated from neighboring values 
(if available).
 */

#if dimension == 2 
#define NCELLS 4 //A macro like this is already defined?
#else
#define NCELLS 8
#endif

double bilinear_embed_gradient (Point point, scalar s, coord * grad, coord * n)
{
  assert (cs[] > 0. && cs[] < 1.); //only in mixed cells
  double area = 0.;
  coord o;
  area = embed_geometry (point, &o, n);
  if (metric_embed_factor)
    area *= metric_embed_factor (point, o);
  bool dirichlet;
  double vb = s.boundary[embed] (point, point, s, &dirichlet);
  
  /**
  By definition the normal *n* points out outward, i.e. from fluid to
  solid, we change the sign because we are interested into locate the
  largest adjacent fluid cells. */
  
  int i = sign(-n->x), j = sign(-n->y);
  
  i = !i ? (cs[-1] >= cs[1] ? -1 : 1)  : i; 
  j = !j ? (cs[0,-1] >= cs[0,1] ? -1 : 1)  : j;

#if dimension > 2  
  int k = sign(-n->z);
  k = !k ? (cs[0,0,-1] >= cs[0,0,1] ? -1 : 1)  : k;
  
  if(cs[i] && cs[0,j] && cs[i,j] &&
     cs[0,0,k] && cs[0,j,k] && cs[i,0,k] && cs[i,j,k]) {
#else
  if(cs[i] && cs[0,j] && cs[i,j]) {
#endif
    
    double **m, b[NCELLS];
    m = (double **) matrix_new (NCELLS, NCELLS, sizeof(double));  //Allocate!
    
    /** 
    In case of neighbouring mixed cells, I should ask Stephane if
    below we assume *s* is defined in the cell center or in the
    centroid of the fluid. If *s* value is defined at the centroids a
    correction with the centroid coordinates should/could be applied
    below. 

    The system of reference is moved to the center of the embeded
    interfacial segment. Therefore, coefficients results of the resolution of
    algebraic equations of the form (in 2D).

    $$
    \left(
    \begin{array}{c}
    s_{10} \\ s_{01} \\ s_{11} \\ v_b
    \end{array}
    \right) = \left(
    \begin{array}{cccc}
    1 & x_{10} & y_{10} & x_{10} y_{10} \\
    1 & x_{01} & y_{01} & x_{01} y_{01} \\
    1 & x_{11} & y_{11} & x_{11} y_{11} \\
    1 &  0 & 0 & 0 \\
    \end{array}
    \right) \left(
    \begin{array}{c}
    a_o \\ a_1 \\ a_2 \\ a_3
    \end{array}
    \right)
    $$

    where we have assumed that the BC is of Dirichlet type. If the BC
    is of Neumann type, the last equation would write...
    $$
     v_b
     = \left(
    \begin{array}{cccc}
    0 &  n_x & n_y & 0 \\
    \end{array}
    \right) \left(
    \begin{array}{c}
    a_o \\ a_1 \\ a_2 \\ a_3
    \end{array}
    \right)
    $$
    */
    
    m[0][0] = 1.; m[0][1] = i-o.x; m[0][2] =  -o.y; b[0] = s[i];
    m[1][0] = 1.; m[1][1] =  -o.x; m[1][2] = j-o.y; b[1] = s[0,j];
    m[2][0] = 1.; m[2][1] = i-o.x; m[2][2] = j-o.y; b[2] = s[i,j];

#if dimension > 2
    m[0][3] = m[1][3] = m[2][3] = -o.z;
    m[3][0] = 1.; m[3][1] =  -o.x; m[3][2] =  -o.y; m[3][3] = k-o.z; b[3] = s[0,0,k];
    m[4][0] = 1.; m[4][1] = i-o.x; m[4][2] =  -o.y; m[4][3] = k-o.z; b[4] = s[i,0,k];
    m[5][0] = 1.; m[5][1] =  -o.x; m[5][2] = j-o.y; m[5][3] = k-o.z; b[5] = s[0,j,k];
    m[6][0] = 1.; m[6][1] = i-o.x; m[6][2] = j-o.y; m[6][3] = k-o.z; b[6] = s[i,j,k];
#endif 
    
    for(int ii = 0; ii < NCELLS; ii++) {
      if(ii < NCELLS-1) {
	m[ii][dimension+1] = m[ii][1]*m[ii][2];
#if dimension > 2
	m[ii][5] = m[ii][1]*m[ii][3];
	m[ii][6] = m[ii][2]*m[ii][3];
#endif	
      }
      m[NCELLS-1][ii] = 0.;
    }
    
    b[NCELLS-1] = vb;
    if(dirichlet)
      m[NCELLS-1][0] = 1.;
    else {
      m[NCELLS-1][1] = n->x; m[NCELLS-1][2] = n->y;
#if dimension > 2
      m[NCELLS-1][3] = n->z;
#endif	
    }
    double pivmin = matrix_inverse (m, NCELLS, 1e-10);
    grad->x = grad->y = 0.;
#if dimension > 2
     grad->z = 0.
#endif       
    if(pivmin) 
       for(int ii = 0; ii < NCELLS; ii++) {
	 grad->x += m[1][ii]*b[ii];
	 grad->y += m[2][ii]*b[ii];
#if dimension > 2
	 grad->z += m[3][ii]*b[ii];
#endif       
    }
    matrix_free (m);
  }
  
  /** 
  In case we do not have available enough contiguous fluid cells, we
  assume that the tangential component of the derivative is zero. The
  normal one is calculated with *dirichlet_gradient* or taken directly
  from the boundary condition if the BC is Neumann. */

  else {
    double coef = 0.;
    if(dirichlet)
      vb = dirichlet_gradient (point, s, cs, *n, o, vb, &coef);
    foreach_dimension()
      grad->x = vb*n->x;
  }
  
  /**
  Let's return the real derivative by dividing with the size of the
  cell. */
  
  foreach_dimension()
    grad->x /= Delta;
  return area;
}

/**
The function below applies a cuadratic correction of the form  
$$
s(x,y) = a_o + a_1 x + a_2 y + a_3 xy 
{\color{blue} + a_4 x^2 + a_5 y^2}
$$
for 2D or of the form,
$$
s(x,y,z) = a_o + a_1 x + a_2 y + a_3 z + a_4 xy + a_5 xz + a_6 yz 
{\color{blue} + a_7 x^2 + a_8 y^2 + + a_9 z^2}
$$
for 3D, can be obtained by including some additional cells
of the 5x5 stencil. */

double bilinear_corrected_embed_gradient (Point point, scalar s,
					  coord * grad, coord * n)
{
  assert (cs[] > 0. && cs[] < 1.); //only in mixed cells
  double area = 0.;
  coord o;
  area = embed_geometry (point, &o, n);
  if (metric_embed_factor)
    area *= metric_embed_factor (point, o);

  bool dirichlet;
  double vb = s.boundary[embed] (point, point, s, &dirichlet);
  
  /**
  By definition the normal *n* points out outward, i.e. from fluid to
  solid, we change the sign because we are interested into locate the
  largest adjacent fluid cells. */
  
  int i = sign(-n->x), j = sign(-n->y), i2 = 0, j2 = 0;
  
  i = !i ? (cs[-1] >= cs[1] ? -1 : 1)  : i; 
  j = !j ? (cs[0,-1] >= cs[0,1] ? -1 : 1)  : j;

  /**
  The extra cells for the correction (if exists) is also the largest
  one. */

  i2 = cs[-i] > cs[2*i] ? -i : 2*i;
  j2 = cs[0,-j] > cs[0,2*j] ? -j : 2*j;
  //  fprintf(stderr, "%d %d %d %d %g %g %g \n", i, i2, j, j2, o.x, o.y, o.z);
#if dimension > 2  
  int k = sign(-n->z), k2 = 0;
  k = !k ? (cs[0,0,-1] >= cs[0,0,1] ? -1 : 1)  : k;
  k2 = cs[0,0,-k] > cs[0,0,2*k] ? -k : 2*k;
  
  if(cs[i] && cs[0,j] && cs[i,j] &&
     cs[0,0,k] && cs[0,j,k] && cs[i,0,k] && cs[i,j,k]) {
#else
  if(cs[i] && cs[0,j] && cs[i,j]) {
#endif

    /**
    The size of system enlarges by *dimension* with the correction. */
    
    double **m, b[NCELLS + dimension];
    m = (double **) matrix_new (NCELLS + dimension, NCELLS + dimension, sizeof(double));  //Allocate!
    
    /** 
    In case of neighbouring mixed cells, I should ask Stephane if
    below we assume *s* is defined in the cell center or in the
    centroid of the fluid. If *s* value is defined at the centroids a
    correction with the centroid coordinates should/could be applied
    below. */
    
    m[0][0] = 1.; m[0][1] = i-o.x; m[0][2] =  -o.y; b[0] = s[i];
    m[1][0] = 1.; m[1][1] =  -o.x; m[1][2] = j-o.y; b[1] = s[0,j];
    m[2][0] = 1.; m[2][1] = i-o.x; m[2][2] = j-o.y; b[2] = s[i,j];

    /**
    Add information from extra cells at rows below. */
   
    m[NCELLS  ][0] = 1.; m[NCELLS  ][1] = i2-o.x; m[NCELLS  ][2] =   -o.y; b[NCELLS  ] = s[i2];
    m[NCELLS+1][0] = 1.; m[NCELLS+1][1] =   -o.x; m[NCELLS+1][2] = j2-o.y; b[NCELLS+1] = s[0,j2];

#if dimension > 2
    m[0][3] = m[1][3] = m[2][3] = m[NCELLS][3] = m[NCELLS+1][3] = -o.z;
    m[3][0] = 1.; m[3][1] =  -o.x; m[3][2] =  -o.y; m[3][3] = k-o.z; b[3] = s[0,0,k];
    m[4][0] = 1.; m[4][1] = i-o.x; m[4][2] =  -o.y; m[4][3] = k-o.z; b[4] = s[i,0,k];
    m[5][0] = 1.; m[5][1] =  -o.x; m[5][2] = j-o.y; m[5][3] = k-o.z; b[5] = s[0,j,k];
    m[6][0] = 1.; m[6][1] = i-o.x; m[6][2] = j-o.y; m[6][3] = k-o.z; b[6] = s[i,j,k];
    
    m[NCELLS+2][0] = 1.;
    m[NCELLS+2][1] =-o.x;
    m[NCELLS+2][2] =-o.y;
    m[NCELLS+2][3] = k2-o.z;
    b[NCELLS+2] = s[0,0,k2];
#endif

    /**
    For all rows except that of NCELLS-1 where the boundary condition
    is applied we complete the crossed and the correction terms. The
    row NCELLS-1 is initialized to 0.  */
    
    for(int ii = 0; ii < NCELLS + dimension; ii++) { 
      if(ii < NCELLS-1 || ii > NCELLS -1) { 
	m[ii][dimension+1] = m[ii][1]*m[ii][2]; // x*y terms
	m[ii][NCELLS] = sq(m[ii][1]); //x^2 terms
	m[ii][NCELLS+1] = sq(m[ii][2]); //y^2 terms
#if dimension > 2
	m[ii][5] = m[ii][1]*m[ii][3]; // x*z terms
	m[ii][6] = m[ii][2]*m[ii][3]; // y*z terms
	m[ii][NCELLS+2] = sq(m[ii][3]); //z^2 terms
#endif
      }
      m[NCELLS-1][ii] = 0.;
    }

    b[NCELLS-1] = vb;
    if(dirichlet)
      m[NCELLS-1][0] = 1.;
    else {
      m[NCELLS-1][1] = n->x; m[NCELLS-1][2] = n->y;
#if dimension > 2
      m[NCELLS-1][3] = n->z;
#endif
    }
    /**
    In case the cuadratic correction could not be performed in some
    direction by the lack of a proper cell, the corresponding term is
    set to zero.*/

    if(!i2) {
      for (int ii = 0; ii < NCELLS + dimension; ii++)
	m[NCELLS][ii] = 0.;
      m[NCELLS][NCELLS] = 1.; b[NCELLS] = 0.;
    }
    if(!j2) {
      for (int ii = 0; ii < NCELLS + dimension; ii++)
	m[NCELLS+1][ii] = 0.;
      m[NCELLS+1][NCELLS+1] = 1.; b[NCELLS+1] = 0.;
    }
#if dimension > 2
    if(!k2) {
      for (int ii = 0; ii < NCELLS + dimension; ii++)
	m[NCELLS+2][ii] = 0.;
      m[NCELLS+2][NCELLS+2] = 1.; b[NCELLS+2] = 0.;
    }
#endif

    double pivmin = matrix_inverse (m, NCELLS + dimension, 1e-10);
    grad->x = grad->y = 0.;
#if dimension > 2
     grad->z = 0.
#endif
    if(pivmin)
       for(int ii = 0; ii < NCELLS + dimension; ii++) {
	 grad->x += m[1][ii]*b[ii];
	 grad->y += m[2][ii]*b[ii];
#if dimension > 2
	 grad->z += m[3][ii]*b[ii];
#endif
    }
    matrix_free (m);
  }

  /**
  In case we do not have available enough contiguous fluid cells, we
  assume that the tangential component of the derivative is zero. The
  normal one is calculated with *dirichlet_gradient* or taken directly
  from the boundary condition if the BC is Neumann. */

  else {
    double coef = 0.;
    if(dirichlet)
      vb = dirichlet_gradient (point, s, cs, *n, o, vb, &coef);
    foreach_dimension()
      grad->x = vb*n->x;
  }

  /**
  Let's return the real derivative by dividing with the size of the
  cell. */

  foreach_dimension()
    grad->x /= Delta;
  return area;
}

/**
This function is totally equivalent to the above function except it
provides the array of coefficients of the interpolation $a_o,a_1,
a_2...$. The bilinear terms are stored in *A1* while those
corresponding to the cuadratic correction are located at *A2*.
The gradients must be then calculated outside the
function. It can be used to calculate, for example, not only the
gradient at the embed segment but also at the contiguous faces.  */

#if dimension > 2
typedef struct { double x, y, z;}   pseudo_v;
typedef struct { pseudo_v x, y, z;} pseudo_t;
#else
typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;
#endif

double coef_bilinear_corrected_embed_gradient (Point point, scalar s,
					       pseudo_t * A1, pseudo_v * A2,
					       coord * o, coord * n)
{
  assert (cs[] > 0. && cs[] < 1.); //only in mixed cells
  double area = 0.;
  area = embed_geometry (point, o, n);
  if (metric_embed_factor)
    area *= metric_embed_factor (point, *o);
  bool dirichlet;
  double vb = s.boundary[embed] (point, point, s, &dirichlet);
  
  /**
  By definition the normal *n* points out outward, i.e. from fluid to
  solid, we change the sign because we are interested into locate the
  largest adjacent fluid cells. */
  
  int i = sign(-n->x), j = sign(-n->y), i2 = 0, j2 = 0;
  
  i = !i ? (cs[-1] >= cs[1] ? -1 : 1)  : i; 
  j = !j ? (cs[0,-1] >= cs[0,1] ? -1 : 1)  : j;

  /**
  The extra cells for the correction (if exists) is also the largest
  one. */

  i2 = cs[-i] > cs[2*i] ? -i : 2*i;
  j2 = cs[0,-j] > cs[0,2*j] ? -j : 2*j;

  //  fprintf(stderr, "%d %d %d %d %g %g %g \n", i, i2, j, j2, o.x, o.y, o.z);
  
#if dimension > 2  
  int k = sign(-n->z), k2 = 0;
  k = !k ? (cs[0,0,-1] >= cs[0,0,1] ? -1 : 1)  : k;
  k2 = cs[0,0,-k] > cs[0,0,2*k] ? -k : 2*k;
  
  if(cs[i] && cs[0,j] && cs[i,j] &&
     cs[0,0,k] && cs[0,j,k] && cs[i,0,k] && cs[i,j,k]) {
#else
  if(cs[i] && cs[0,j] && cs[i,j]) {
#endif

    /**
    The size of system enlarges by *dimension* with the correction. */
    
    double **m, b[NCELLS + dimension], a[NCELLS + dimension];
    m = (double **) matrix_new (NCELLS + dimension, NCELLS + dimension, sizeof(double));  //Allocate!
    
    /** 
    In case of neighbouring mixed cells, I should ask Stephane if
    below we assume *s* is defined in the cell center or in the
    centroid of the fluid. If *s* value is defined at the centroids a
    correction with the centroid coordinates should/could be applied
    below. */
    
    m[0][0] = 1.; m[0][1] = i-o->x; m[0][2] =  -o->y; b[0] = s[i];
    m[1][0] = 1.; m[1][1] =  -o->x; m[1][2] = j-o->y; b[1] = s[0,j];
    m[2][0] = 1.; m[2][1] = i-o->x; m[2][2] = j-o->y; b[2] = s[i,j];

    /**
    Add information from extra cells at rows below. */
   
    m[NCELLS  ][0] = 1.; m[NCELLS  ][1] = i2-o->x; m[NCELLS  ][2] =   -o->y; b[NCELLS  ] = s[i2];
    m[NCELLS+1][0] = 1.; m[NCELLS+1][1] =   -o->x; m[NCELLS+1][2] = j2-o->y; b[NCELLS+1] = s[0,j2];

#if dimension > 2
    m[0][3] = m[1][3] = m[2][3] = m[NCELLS][3] = m[NCELLS+1][3] = -o->z;
    m[3][0] = 1.; m[3][1] =  -o->x; m[3][2] =  -o->y; m[3][3] = k-o->z; b[3] = s[0,0,k];
    m[4][0] = 1.; m[4][1] = i-o->x; m[4][2] =  -o->y; m[4][3] = k-o->z; b[4] = s[i,0,k];
    m[5][0] = 1.; m[5][1] =  -o->x; m[5][2] = j-o->y; m[5][3] = k-o->z; b[5] = s[0,j,k];
    m[6][0] = 1.; m[6][1] = i-o->x; m[6][2] = j-o->y; m[6][3] = k-o->z; b[6] = s[i,j,k];
    
    m[NCELLS+2][0] = 1.;
    m[NCELLS+2][1] =-o->x;
    m[NCELLS+2][2] =-o->y;
    m[NCELLS+2][3] = k2-o->z;
    b[NCELLS+2] = s[0,0,k2];
#endif

    /**
    For all rows except that of NCELLS-1 where the boundary condition
    is applied we complete the crossed and the correction terms. The
    row NCELLS-1 is initialized to 0.  */
    
    for(int ii = 0; ii < NCELLS + dimension; ii++) { 
      if(ii < NCELLS-1 || ii > NCELLS -1) { 
	m[ii][dimension+1] = m[ii][1]*m[ii][2]; // x*y terms
	m[ii][NCELLS] = sq(m[ii][1]); //x^2 terms
	m[ii][NCELLS+1] = sq(m[ii][2]); //y^2 terms
#if dimension > 2
	m[ii][5] = m[ii][1]*m[ii][3]; // x*z terms
	m[ii][6] = m[ii][2]*m[ii][3]; // y*z terms
	m[ii][NCELLS+2] = sq(m[ii][3]); //z^2 terms
#endif
      }
      m[NCELLS-1][ii] = 0.;
    }

    b[NCELLS-1] = vb;
    if(dirichlet)
      m[NCELLS-1][0] = 1.;
    else {
      m[NCELLS-1][1] = n->x; m[NCELLS-1][2] = n->y;
#if dimension > 2
      m[NCELLS-1][3] = n->z;
#endif
    }
    /**
    In case the cuadratic correction could not be performed in some
    direction by the lack of a proper cell, the corresponding term is
    set to zero.*/

    if(!i2) {
      for (int ii = 0; ii < NCELLS + dimension; ii++)
	m[NCELLS][ii] = 0.;
      m[NCELLS][NCELLS] = 1.; b[NCELLS] = 0.;
    }
    if(!j2) {
      for (int ii = 0; ii < NCELLS + dimension; ii++)
	m[NCELLS+1][ii] = 0.;
      m[NCELLS+1][NCELLS+1] = 1.; b[NCELLS+1] = 0.;
    }
#if dimension > 2
    if(!k2) {
      for (int ii = 0; ii < NCELLS + dimension; ii++)
	m[NCELLS+2][ii] = 0.;
      m[NCELLS+2][NCELLS+2] = 1.; b[NCELLS+2] = 0.;
    }
#endif
    
    double pivmin = matrix_inverse (m, NCELLS + dimension, 1e-10);
    if(pivmin)
      for(int  jj = 0; jj < NCELLS + dimension; jj++) {
	a[jj] = 0;
	for(int ii = 0; ii < NCELLS + dimension; ii++) 
	  a[jj] += m[jj][ii]*b[ii];
      }

    matrix_free (m);
    A1->x.x = a[1];
    A1->y.y = a[2];
    A1->x.y = a[3];
    A1->y.x = A1->x.y;
#if dimension > 2
    A1->x.z = a[4];
    A1->y.z = a[5];
    A1->z.x = A1->x.z;
    A1->z.y = A1->y.z;
    A2->x = a[6];
    A2->y = a[7];
    A2->z = a[8];
#else
    A2->x = a[4];
    A2->y = a[5];
#endif
  }

  /**
  In case we do not have available enough contiguous fluid cells, we
  assume that the tangential component of the derivative is zero. The
  normal one is calculated with *dirichlet_gradient* or taken directly
  from the boundary condition if the BC is Neumann. */

  else {
    double coef = 0.;
    if(dirichlet)
      vb = dirichlet_gradient (point, s, cs, *n, *o, vb, &coef);
    //    for(int ii = 0; ii < NCELLS + dimension; ii++) 
    //      a[ii] = 0.;
    A1->x.x = vb*n->x; A1->y.y = vb*n->y;
#if dimension > 2 
    A1->z.z = vb*n->z;
#endif
  }
  return area;
}
